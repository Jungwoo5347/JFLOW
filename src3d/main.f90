!> JFLOW3D : 3D 압축성 오일러 솔버 (곡선 정렬 격자, 다중 블록, 셀 중심 FVM)
!>   2D 솔버 (src/) 와 같은 구조, 같은 기법 (MUSCL + Rusanov/Roe + SSP-RK3), 같은 입력 형식
!>
!> 빌드:  build.bat 3d release omp       ->  main3d.exe
!> 실행:  .\main3d.exe cases3d\sod       (케이스 폴더, 생략하면 현재 폴더)
!>   케이스 폴더의 grid.in (격자, 경계조건) 과 flow.in (계산 제어, 유동 조건) 을 읽고
!>   결과를 같은 폴더에 저장 (입력 형식은 case_setup.f90 맨 위 설명)
!>     <output>_result.plt / .dat  : Tecplot (블록 하나 = zone 하나)
!>     <output>_line.dat           : 셀 줄 하나 (1D 그래프용)
!>     <output>_residual.dat       : 잔차 이력
!>     <output>_restart.bin        : 재시작 파일 (flow.in 의 restart = 1 로 이어서 계산)
!>   계산 중 케이스 폴더에 STOP 파일을 만들면 그 반복을 마치고 저장한 뒤 끝냄
program main
    use, intrinsic :: iso_fortran_env, only: output_unit
    use grid
    use euler
    use case_setup
    use flow
    use solver
    use implicit, only: lusgs
    use monitor
    use output
    !$ use omp_lib, only: omp_get_max_threads           ! -fopenmp 로 빌드할 때만 사용

    implicit none

    character(len=256) :: arg
    real(dp) :: t, dt, dtb, t0
    integer  :: step, step0, ib
    logical  :: converged, stop_req

    ! ---- 입력 (케이스 폴더) ----
    call get_command_argument(1, arg)
    if (len_trim(arg) == 0) arg = '.'
    t0 = wall_seconds()
    call read_case(trim(arg))
    call print_case()
    !$ write(*,'(a,i0)') ' OpenMP   : threads = ', omp_get_max_threads()

    ! ---- 격자 ----
    call build_blocks(specs, nspecs, nblock_in)
    call check_grid()

    ! ---- 경계조건, 초기조건 ----
    call setup_bc()
    call init_flow()

    if (line_block < 1 .or. line_block > nblk) then
        write(*,'(a,i0,a,i0)') ' Error: line_block = ', line_block, ' must be between 1 and ', nblk
        error stop
    end if
    if (line_row < 1 .or. line_row > blk(line_block)%ny .or. line_k < 1 .or. line_k > blk(line_block)%nz) then
        write(*,'(a,i0,a,i0,a,i0,a,i0)') ' Error: line_row = ', line_row, ', line_k = ', line_k, &
            ' must be within 1..', blk(line_block)%ny, ', 1..', blk(line_block)%nz
        error stop
    end if

    do ib = 1,nblk
        associate (bk => blk(ib))
            allocate(bk%R(nvar,bk%nx,bk%ny,bk%nz), bk%Rt(nvar,bk%nx,bk%ny,bk%nz))
            allocate(bk%Q0(nvar,bk%nx,bk%ny,bk%nz))
        end associate
    end do

    ! ---- 재시작 ----
    step = 0
    t = 0.0_dp
    if (restart) then
        call read_restart(trim(case_path(trim(output_name)//'_restart.bin')), step, t, res0)
        write(*,'(a,i0,a,es12.4)') ' restart  : continuing from step ', step, ',  t = ', t
    end if
    step0 = step

    ! ---- 시간 루프 ----
    call monitor_open(step, t)
    write(*,'(a,a)') ' setup    : ', trim(hms(wall_seconds() - t0))
    flush(output_unit)
    do
        if (steady) then
            if (step >= max_steps) exit
        else
            if (t >= t_end*(1.0_dp - 1.0e-12_dp)) exit
        end if

        call bc()                       ! 고스트 셀 채우기 (모든 블록)
        dt = huge(1.0_dp)
        do ib = 1,nblk
            call compute_dt(ib, dtb)
            dt = min(dt, dtb)
        end do
        if (.not. steady) dt = min(dt, t_end - t)    ! 마지막 스텝은 t_end에 정확히 맞춤
        ! 국소 시간 간격은 정상 문제에서만. 비정상 문제나 local_dt = 0 이면 모든 셀이 같은 dt (전체 최솟값)
        if (.not. (steady .and. local_dt)) then
            do ib = 1,nblk
                blk(ib)%dtl = dt
            end do
        end if

        if (implicit_ts) then
            ! 암시적 LU-SGS: 잔차 한 번 (원시 변수 W 도 여기서 계산) 뒤 블록마다 두 번 스윕해서 갱신
            do ib = 1,nblk
                call residual(ib, blk(ib)%R)
            end do
            do ib = 1,nblk
                call lusgs(ib)
            end do
        else
        ! 시간 전진: SSP-RK3  (단계마다 모든 블록의 잔차를 구한 뒤 갱신)
        do ib = 1,nblk
            call residual(ib, blk(ib)%R)
            blk(ib)%Q0 = blk(ib)%Q(:,1:blk(ib)%nx,1:blk(ib)%ny,1:blk(ib)%nz)
        end do
        do ib = 1,nblk
            call stage(ib, 1.0_dp, 0.0_dp, blk(ib)%R)             ! Q1 = Q0 + dt L(Q0)
        end do
        call bc()
        do ib = 1,nblk
            call residual(ib, blk(ib)%Rt)
        end do
        do ib = 1,nblk
            call stage(ib, 0.25_dp, 0.75_dp, blk(ib)%Rt)          ! Q2 = 3/4 Q0 + 1/4 (Q1 + dt L(Q1))
        end do
        call bc()
        do ib = 1,nblk
            call residual(ib, blk(ib)%Rt)
        end do
        do ib = 1,nblk
            call stage(ib, 2.0_dp/3.0_dp, 1.0_dp/3.0_dp, blk(ib)%Rt)   ! Q = 1/3 Q0 + 2/3 (Q2 + dt L(Q2))
        end do
        end if

        t = t + dt
        step = step + 1

        call check_solution(step, t)                 ! NaN, 음수 밀도면 위치 알리고 멈춤
        call monitor_update(step, t, converged)      ! 잔차 기록, 화면 출력, 수렴 판정
        if (converged) exit

        ! 중간 저장
        if (save_interval > 0) then
            if (mod(step, save_interval) == 0) call save_all(.false.)
        end if

        ! 케이스 폴더에 STOP 파일이 있으면 지우고 끝냄 (아래에서 결과 저장)
        inquire(file=trim(case_path('STOP')), exist=stop_req)
        if (stop_req) then
            call remove_file(trim(case_path('STOP')))
            write(*,'(a,i0)') ' STOP file found: stopping after step ', step
            exit
        end if
    end do
    call monitor_close()

    print '(a,i8,a,es12.4)', ' done: step', step, '   t =', t
    call save_all(.true.)
    write(*,'(a,a,a,i0,a)') ' wall time: ', trim(hms(wall_seconds() - t0)), '  (', step - step0, ' iterations this run)'

contains

    !> RK 한 단계 (블록 b):  Q = a0*Q0 + a1*(Q + dtl*Rs)   (내부 셀만)
    subroutine stage(ib, a1, a0, Rs)
        integer,  intent(in) :: ib
        real(dp), intent(in) :: a1, a0
        real(dp), intent(in) :: Rs(:,:,:,:)
        integer :: ic, jc, kc
        real(dp), pointer, contiguous :: Q(:,:,:,:), Q0(:,:,:,:), dtl(:,:,:)
        Q => blk(ib)%Q;  Q0 => blk(ib)%Q0;  dtl => blk(ib)%dtl
        !$omp parallel do collapse(2) private(ic)
        do kc = 1,blk(ib)%nz
            do jc = 1,blk(ib)%ny
                do ic = 1,blk(ib)%nx
                    Q(:,ic,jc,kc) = a0*Q0(:,ic,jc,kc) + a1*(Q(:,ic,jc,kc) + dtl(ic,jc,kc)*Rs(:,ic,jc,kc))
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine stage

    !> 결과 (plt / ascii) 와 재시작 파일 저장. final 이면 line 파일도
    subroutine save_all(final)
        logical, intent(in) :: final
        real(dp) :: ts
        ts = wall_seconds()
        if (out_plt)   call write_plt(trim(case_path(trim(output_name)//'_result.plt')), t)
        if (out_ascii) call write_tecplot(trim(case_path(trim(output_name)//'_result.dat')), t)
        call write_restart(trim(case_path(trim(output_name)//'_restart.bin')), step, t, res0)
        if (final) call write_line(trim(case_path(trim(output_name)//'_line.dat')), line_block, line_row, line_k)
        write(*,'(a,i0,a,a,a)') ' saved    : step ', step, ' (result + restart, ', trim(hms(wall_seconds() - ts)), ')'
        flush(output_unit)
    end subroutine save_all

    subroutine remove_file(fname)
        character(len=*), intent(in) :: fname
        integer :: u, ios
        open(newunit=u, file=fname, status='old', iostat=ios)
        if (ios == 0) close(u, status='delete')
    end subroutine remove_file

end program main
