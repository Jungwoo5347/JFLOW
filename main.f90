program main
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use param, only: print_interval
    use grid
    use euler
    use init_cond, only: case_name
    use flow
    use solver
    use output

    implicit none

    real(dp), allocatable :: R(:,:,:)
    real(dp) :: t, dt, t_end, res, res0, res_tol
    integer  :: step, max_steps, line_j
    logical  :: steady

    ! ---- 문제 설정 (case_name 은 init_cond.f90 에서 고름) ----
    select case (case_name)
    case ('sod')
        call init_grid()
        bc_type = [BC_outlet, BC_outlet, BC_outlet, BC_outlet]
        steady  = .false.
        t_end   = 0.2_dp
        line_j  = ny/2                  ! 가운데 줄 출력

    case ('wedge')
        call read_plot3d('grids/wedge.x.fmt')
        ! 왼쪽 입구, 오른쪽 출구, 아래 벽(평판 + 쐐기), 위 자유류
        bc_type   = [BC_inlet, BC_outlet, BC_wall, BC_inlet]
        steady    = .true.
        max_steps = 50000
        res_tol   = 1.0e-8_dp           ! 밀도 잔차가 처음의 1e-8 배가 되면 수렴
        line_j    = 1                   ! 벽 바로 위 셀 줄 출력

    case default
        write(*,*) 'Unknown case_name: ', case_name
        error stop
    end select

    call check_grid()
    call init_flow()

    allocate(R(nvar,0:nx+1,0:ny+1))

    ! ---- 시간 루프 (1차 오일러 전진) ----
    t = 0.0_dp
    step = 0
    res0 = 0.0_dp
    do
        if (steady) then
            if (step >= max_steps) exit
        else
            if (t >= t_end*(1.0_dp - 1.0e-12_dp)) exit
        end if

        call bc()                       ! 고스트 셀 채우기
        call compute_dt(Q, dt)
        if (.not. steady) dt = min(dt, t_end - t)    ! 마지막 스텝은 t_end에 정확히 맞춤
        call residual(Q, R)

        Q(:,1:nx,1:ny) = Q(:,1:nx,1:ny) + dt*R(:,1:nx,1:ny)
        t = t + dt
        step = step + 1

        if (any(ieee_is_nan(Q(:,1:nx,1:ny))) .or. any(Q(1,1:nx,1:ny) <= 0.0_dp)) then
            print '(a,i0,a,es12.4)', ' Solution blew up at step ', step, ', t = ', t
            error stop
        end if

        ! 밀도 잔차 (RMS) : 정상 상태 수렴 판단용
        res = sqrt(sum(R(1,1:nx,1:ny)**2)/(nx*ny))
        if (step == 1) res0 = res

        if (mod(step, print_interval) == 0) then
            print '(a,i8,a,es12.4,a,es12.4,a,es10.3)', ' step', step, '   t =', t, '   dt =', dt, &
                  '   res/res0 =', res/res0
        end if

        if (steady .and. res < res_tol*res0) then
            print '(a,es10.3)', ' converged: res/res0 =', res/res0
            exit
        end if
    end do

    print '(a,i8,a,es12.4)', ' done: step', step, '   t =', t

    call write_tecplot('result.dat', Q, t)
    call write_line('line.dat', Q, line_j)

end program main
