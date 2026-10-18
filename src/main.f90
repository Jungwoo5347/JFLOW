!> JFLOW : 2D 압축성 오일러 솔버 (곡선 정렬 격자, 셀 중심 FVM)
!>
!> 실행:  .\main.exe cases\wedge          (케이스 폴더, 생략하면 현재 폴더)
!>   케이스 폴더의 grid.in (격자, 경계조건) 과 flow.in (계산 제어, 유동 조건) 을 읽고
!>   결과를 같은 폴더에 <output>_result.plt (Tecplot 바이너리, output_format 으로 .dat 도), _line.dat, _residual.dat 로 저장
!>   <output>_restart.bin : 재시작 파일 (flow.in 의 restart = 1 로 이어서 계산, save_interval 마다 중간 저장)
!>   (입력 형식은 case_setup.f90 맨 위 설명)
program main
    use grid
    use euler
    use case_setup
    use flow
    use solver
    use monitor
    use output
    use implicit, only: lusgs
    use, intrinsic :: iso_fortran_env, only: output_unit
    !$ use omp_lib, only: omp_get_max_threads           ! -fopenmp 로 빌드할 때만 사용

    implicit none

    character(len=256) :: arg
    real(dp), allocatable :: R(:,:,:)
    real(dp), allocatable :: R_temp(:,:,:)
    real(dp), allocatable :: Q_temp(:,:,:)
    real(dp), allocatable :: dtl(:,:)                  ! 셀마다의 dt (국소 시간 간격)
    real(dp) :: t, dt, t0
    integer  :: step, step0, iter0, ic, jc
    logical  :: converged

    ! dual time (time_scheme = implicit, steady = 0): BDF2 에 필요한 지난 두 시간의 해 (내부 셀)
    real(dp), allocatable :: Qn(:,:,:), Qold(:,:,:)     ! Q^n (이번 스텝 시작), Q^(n-1)
    real(dp) :: dt_old = 0.0_dp                         ! 지난 스텝의 dt (BDF2 가변 간격 계수용)
    logical  :: have_old = .false.                      ! Q^(n-1) 있음 (첫 스텝은 BDF1)
    integer  :: iter_total = 0                          ! 안쪽 반복 누적 수 (잔차 기록의 Iteration, 재시작해도 이어짐)

    ! ---- 입력 (케이스 폴더) ----
    t0 = wall_seconds()
    call get_command_argument(1, arg)
    if (len_trim(arg) == 0) arg = '.'
    call read_case(trim(arg))
    call print_case()
    !$ write(*,'(a,i0)') ' OpenMP   : threads = ', omp_get_max_threads()

    ! ---- 격자 ----
    select case (grid_kind)
    case ('plot3d')
        call read_plot3d(trim(grid_file), grid_form)
    case ('cartesian')
        call init_grid(grid_ni-1, grid_nj-1, grid_x0, grid_x1, grid_y0, grid_y1)
    end select
    call check_grid()

    ! ---- 경계조건, 초기조건 ----
    call setup_bc()
    call init_flow()

    if (line_row < 1 .or. line_row > ny) then
        write(*,'(a,i0,a,i0)') ' Error: line_row = ', line_row, ' must be between 1 and ', ny
        error stop
    end if

    allocate(R(nvar,0:nx+1,0:ny+1))
    allocate(R_temp(nvar,0:nx+1,0:ny+1))
    allocate(Q_temp, mold=Q)
    allocate(dtl(nx,ny))

    if (implicit_ts .and. .not. steady) then
        allocate(Qn(nvar,nx,ny), Qold(nvar,nx,ny))
        Qold = 0.0_dp
    end if

    ! ---- 재시작 ----
    step = 0
    t = 0.0_dp
    if (restart) then
        if (implicit_ts .and. .not. steady) then
            call read_restart(trim(case_path(trim(output_name)//'_restart.bin')), Q, step, t, res0, &
                              Qold, dt_old, have_old, iter_total)
        else
            call read_restart(trim(case_path(trim(output_name)//'_restart.bin')), Q, step, t, res0)
        end if
        write(*,'(a,i0,3a,es12.4)') ' restart  : continuing from step ', step, ',  ', trim(t_label), ' = ', t*t_unit
    end if
    step0 = step
    iter0 = iter_total

    ! ---- 시간 루프 ----
    if (save_series) call make_dir(trim(case_path('res')))   ! 시간 이력 결과 폴더 <case>/res

    call monitor_open(step, t, implicit_ts .and. .not. steady)
    flush(output_unit)
    do
        if (steady) then
            if (step >= max_steps) exit
        else
            if (t >= t_end*(1.0_dp - 1.0e-12_dp)) exit
        end if

        if (implicit_ts .and. steady) then
            ! ---- 정상, 암시적: LU-SGS 반복 한 번 ----
            call bc()
            call compute_dt(Q, dtl, dt)
            if (.not. local_dt) dtl = dt
            call residual(Q, R)
            call implicit_step(R, 0.0_dp)

        else if (implicit_ts) then
            ! ---- 비정상, 암시적: dual time 한 스텝 (BDF2 + 안쪽 LU-SGS 반복) ----
            call dual_time_step()

        else
        call bc()                       ! 고스트 셀 채우기
        call compute_dt(Q, dtl, dt)
        if (.not. steady) dt = min(dt, t_end - t)    ! 마지막 스텝은 t_end에 정확히 맞춤
        ! 국소 시간 간격은 정상 문제에서만 (셀마다 자기 dt 로 전진 -> 시간은 의미 없음, 수렴만 빨라짐)
        !   비정상 문제나 local_dt = 0 이면 모든 셀이 같은 dt (최솟값)
        if (.not. (steady .and. local_dt)) dtl = dt
        call residual(Q, R)

        ! 시간 전진: SSP-Rk3
        Q_temp = Q
        !$omp parallel do private(ic)
        do jc = 1,ny
            do ic = 1,nx
                Q(:,ic,jc) = Q_temp(:,ic,jc) + dtl(ic,jc)*R(:,ic,jc)
            end do
        end do
        !$omp end parallel do
        call bc()                     ! 바뀐 Q 에 맞게 고스트 셀
        call residual(Q, R_temp)      ! R_temp = L(Q1) 또는 L(Q2)

        !$omp parallel do private(ic)
        do jc = 1,ny
            do ic = 1,nx
                Q(:,ic,jc) = 0.75_dp*Q_temp(:,ic,jc) + 0.25_dp*(Q(:,ic,jc)+dtl(ic,jc)*R_temp(:,ic,jc))
            end do
        end do
        !$omp end parallel do
        call bc()                     ! 바뀐 Q 에 맞게 고스트 셀
        call residual(Q, R_temp)      ! R_temp = L(Q1) 또는 L(Q2)

        !$omp parallel do private(ic)
        do jc = 1,ny
            do ic = 1,nx
                Q(:,ic,jc) = 1.0_dp/3.0_dp*Q_temp(:,ic,jc) + 2.0_dp/3.0_dp*(Q(:,ic,jc)+dtl(ic,jc)*R_temp(:,ic,jc))
            end do
        end do
        !$omp end parallel do
        end if
        t = t + dt
        step = step + 1

        call check_solution(Q, step, t)              ! NaN, 음수 밀도면 멈춤
        if (implicit_ts .and. .not. steady) then     ! dual time: 잔차는 안쪽 반복마다 dual_time_step 에서 기록
            converged = .false.
        else
            call monitor_update(step, t, R, converged)   ! 잔차 기록, 화면 출력, 수렴 판정
        end if
        if (converged) exit

        ! 중간 저장
        if (save_interval > 0) then
            if (mod(step, save_interval) == 0) call save_all(.false.)
        end if
    end do
    call monitor_close()

    print '(a,i8,3a,es12.4)', ' done: step', step, '   ', trim(t_label), ' =', t*t_unit
    call save_all(.true.)
    if (implicit_ts .and. .not. steady) then
        write(*,'(a,a,a,i0,a,i0,a)') ' wall time: ', trim(hms(wall_seconds() - t0)), '  (', step - step0, &
            ' time steps, ', iter_total - iter0, ' inner iterations this run)'
    else
        write(*,'(a,a,a,i0,a)') ' wall time: ', trim(hms(wall_seconds() - t0)), '  (', step - step0, ' iterations this run)'
    end if

contains

    !> 암시적 반복 한 번 (implicit_solver 에 따라): dQ/dtau = Rs 를 풀어 Q 갱신
    !>   diag_add : 대각항에 더할 값 / 면적 (dual time: a0/dt,  정상: 0)
    !>   (DADI 등 다른 방법은 같은 인자의 서브루틴을 만들어 여기에 case 추가)
    subroutine implicit_step(Rs, diag_add)
        real(dp), intent(in) :: Rs(:,0:,0:)
        real(dp), intent(in) :: diag_add
        select case (implicit_solver)
        case ('lusgs');  call lusgs(Q, Rs, dtl, diag_add)
        end select
    end subroutine implicit_step

    !> dual time stepping 한 스텝 (실제 시간 dt, BDF2)
    !>   (a0 Q^(n+1) + a1 Q^n + a2 Q^(n-1)) / dt = R(Q^(n+1))  를 가상 시간 반복으로 풂
    !>     dQ/dtau = R* = R(Q) - (a0 Q + a1 Q^n + a2 Q^(n-1))/dt,  R* -> 0 까지 (inner_iter 번 또는 inner_tol 만큼 줄 때까지)
    !>   BDF2 가변 간격 (w = dt/dt_old): a0 = (1+2w)/(1+w), a1 = -(1+w), a2 = w^2/(1+w)   (w = 1: 3/2, -2, 1/2)
    !>   첫 스텝 (Q^(n-1) 없음) 은 BDF1: a0 = 1, a1 = -1, a2 = 0
    !>   안쪽 반복마다 R* 를 monitor_inner 로 기록 (Iteration = iter_total, Step = 이번 시간 스텝)
    !>   끝나면 dt = 이번 스텝 dt
    subroutine dual_time_step()
        real(dp) :: dtp, w, a0, a1, a2, rn, r1
        integer  :: m
        logical  :: done

        dtp = min(dt_phys, t_end - t)                   ! 마지막 스텝은 t_end 에 맞춤 (가변 간격 BDF2 로)
        if (have_old) then
            w  = dtp/dt_old
            a0 = (1.0_dp + 2.0_dp*w)/(1.0_dp + w)
            a1 = -(1.0_dp + w)
            a2 = w*w/(1.0_dp + w)
        else
            a0 = 1.0_dp;  a1 = -1.0_dp;  a2 = 0.0_dp
        end if
        Qn = Q(:,1:nx,1:ny)

        r1 = 0.0_dp;  rn = 0.0_dp
        do m = 1,inner_iter
            call bc()
            call compute_dt(Q, dtl, dt)                  ! 가상 시간 간격: 셀마다 (cfl 은 가상 시간용)
            call residual(Q, R)
            !$omp parallel do private(ic)
            do jc = 1,ny
                do ic = 1,nx
                    R_temp(:,ic,jc) = R(:,ic,jc) - (a0*Q(:,ic,jc) + a1*Qn(:,ic,jc) + a2*Qold(:,ic,jc))/dtp
                end do
            end do
            !$omp end parallel do
            rn = sqrt(sum(R_temp(:,1:nx,1:ny)**2))
            if (m == 1) r1 = rn
            done = m > 1 .and. rn <= inner_tol*r1        ! 충분히 줄었으면 갱신 없이 끝
            iter_total = iter_total + 1
            call monitor_inner(iter_total, step+1, t+dtp, R_temp, done .or. m == inner_iter, m)
            if (done) exit
            call implicit_step(R_temp, a0/dtp)
        end do

        Qold = Qn;  dt_old = dtp;  have_old = .true.
        dt = dtp
    end subroutine dual_time_step

    !> 결과와 재시작 파일 저장. final 이면 line 파일도
    !>   결과: <output>_result.plt (덮어씀).  save_series 이면 중간 저장은 <output>_result_00012000.plt 처럼
    !>   반복 수를 붙인 새 파일 (끝날 때는 둘 다: 번호 파일로 시간 이력의 마지막 장, 번호 없는 파일로 최종 결과)
    subroutine save_all(final)
        logical, intent(in) :: final
        character(len=16) :: tag

        if (.not. save_series .or. final) then
            if (out_plt)   call write_plt(trim(case_path(trim(output_name)//'_result.plt')), Q, t)
            if (out_ascii) call write_tecplot(trim(case_path(trim(output_name)//'_result.dat')), Q, t)
        end if
        if (save_series) then
            write(tag,'(a,i8.8)') '_result_', step
            if (out_plt)   call write_plt(trim(case_path('res/'//trim(output_name)//trim(tag)//'.plt')), Q, t)
            if (out_ascii) call write_tecplot(trim(case_path('res/'//trim(output_name)//trim(tag)//'.dat')), Q, t)
        end if
        if (have_old) then          ! dual time: Q^(n-1) 도 저장 (이어서 계산해도 BDF2 그대로)
            call write_restart(trim(case_path(trim(output_name)//'_restart.bin')), Q, step, t, res0, &
                               Qold, dt_old, iter_total)
        else
            call write_restart(trim(case_path(trim(output_name)//'_restart.bin')), Q, step, t, res0)
        end if
        if (final) call write_line(trim(case_path(trim(output_name)//'_line.dat')), Q, line_row)
        if (save_series) then
            write(*,'(a,i0,3a,es12.4,a,a,a)') ' saved    : step ', step, ',  ', trim(t_label), ' = ', t*t_unit, &
                ' (res/', trim(output_name)//trim(tag), ' + restart)'
        else
            write(*,'(a,i0,a)') ' saved    : step ', step, ' (result + restart)'
        end if
        flush(output_unit)
    end subroutine save_all

    !> 폴더 만들기 (이미 있으면 그대로). Fortran 표준에 mkdir 이 없어서 Windows 명령으로
    subroutine make_dir(path)
        character(len=*), intent(in) :: path
        character(len=len(path)) :: wpath
        integer :: k, stat, cstat

        wpath = path
        do k = 1,len(wpath)                       ! cmd 의 mkdir 은 / 를 못 읽음 -> \
            if (wpath(k:k) == '/') wpath(k:k) = '\'
        end do
        call execute_command_line('if not exist "'//trim(wpath)//'" mkdir "'//trim(wpath)//'"', &
                                  exitstat=stat, cmdstat=cstat)
        if (cstat /= 0 .or. stat /= 0) then
            write(*,'(a,a)') ' Error: cannot create folder ', trim(path)
            error stop
        end if
    end subroutine make_dir

end program main
