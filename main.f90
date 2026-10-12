!> JFLOW : 2D 압축성 오일러 솔버 (곡선 정렬 격자, 셀 중심 FVM)
!>
!> 실행:  .\main.exe cases\wedge          (케이스 폴더, 생략하면 현재 폴더)
!>   케이스 폴더의 grid.in (격자, 경계조건) 과 flow.in (계산 제어, 유동 조건) 을 읽고
!>   결과를 같은 폴더에 <output>_result.dat, _line.dat, _residual.dat 로 저장
!>   (입력 형식은 case_setup.f90 맨 위 설명)
program main
    use grid
    use euler
    use case_setup
    use flow
    use solver
    use monitor
    use output
    !$ use omp_lib, only: omp_get_max_threads           ! -fopenmp 로 빌드할 때만 사용

    implicit none

    character(len=256) :: arg
    real(dp), allocatable :: R(:,:,:)
    real(dp), allocatable :: R_temp(:,:,:)
    real(dp), allocatable :: Q_temp(:,:,:)
    real(dp), allocatable :: dtl(:,:)                  ! 셀마다의 dt (국소 시간 간격)
    real(dp) :: t, dt
    integer  :: step, ic, jc
    logical  :: converged

    ! ---- 입력 (케이스 폴더) ----
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

    ! ---- 시간 루프 ----
    call monitor_open()
    t = 0.0_dp
    step = 0
    do
        if (steady) then
            if (step >= max_steps) exit
        else
            if (t >= t_end*(1.0_dp - 1.0e-12_dp)) exit
        end if

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
        t = t + dt
        step = step + 1

        call check_solution(Q, step, t)              ! NaN, 음수 밀도면 멈춤
        call monitor_update(step, t, R, converged)   ! 잔차 기록, 화면 출력, 수렴 판정
        if (converged) exit
    end do
    call monitor_close()

    print '(a,i8,a,es12.4)', ' done: step', step, '   t =', t

    call write_tecplot(trim(case_path(trim(output_name)//'_result.dat')), Q, t)
    call write_line(trim(case_path(trim(output_name)//'_line.dat')), Q, line_row)

end program main
