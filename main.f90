!> JFLOW : 2D 압축성 오일러 솔버 (곡선 정렬 격자, 셀 중심 FVM)
!>
!> 실행:  .\main.exe cases\wedge          (케이스 폴더, 생략하면 현재 폴더)
!>   케이스 폴더의 grid.in (격자, 경계조건) 과 flow.in (계산 제어, 유동 조건) 을 읽고
!>   결과를 같은 폴더에 <output>_result.dat, _line.dat, _residual.dat 로 저장
!>   (입력 형식은 case_setup.f90 맨 위 설명)
program main
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use grid
    use euler
    use case_setup
    use flow
    use solver
    use output
    !$ use omp_lib, only: omp_get_max_threads           ! -fopenmp 로 빌드할 때만 사용

    implicit none

    character(len=256) :: arg
    real(dp), allocatable :: R(:,:,:)
    real(dp) :: t, dt
    real(dp) :: res(nvar), res0(nvar), res_n(nvar)     ! 잔차 L2 norm, 정규화 기준, 정규화된 값
    integer  :: step, k, ures, nprint, ic, jc
    logical  :: bad

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

    ! ---- 잔차 기록 파일 <output>_residual.dat : 매 반복 기록 (tools/plot_residual.m 으로 그림) ----
    open(newunit=ures, file=trim(case_path(trim(output_name)//'_residual.dat')), status='replace', action='write')
    write(ures,'(a)') '# iter  time  continuity  x-momentum  y-momentum  energy   (L2 norm / max of first 5 iters)'

    ! ---- 시간 루프 (1차 오일러 전진) ----
    t = 0.0_dp
    step = 0
    nprint = 0
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

        ! 갱신 + 잔차 제곱합 + 발산 검사를 한 루프에서 (OpenMP: 합은 reduction 으로 모음)
        res  = 0.0_dp
        bad  = .false.
        !$omp parallel do private(ic) reduction(+:res) reduction(.or.:bad)
        do jc = 1,ny
            do ic = 1,nx
                Q(:,ic,jc) = Q(:,ic,jc) + dt*R(:,ic,jc)
                res = res + R(:,ic,jc)**2
                bad = bad .or. any(ieee_is_nan(Q(:,ic,jc))) .or. Q(1,ic,jc) <= 0.0_dp
            end do
        end do
        !$omp end parallel do
        t = t + dt
        step = step + 1

        if (bad) then
            print '(a,i0,a,es12.4)', ' Solution blew up at step ', step, ', t = ', t
            close(ures)
            error stop
        end if

        ! 잔차 L2 norm (식마다):  sqrt( sum(R^2) / 셀 수 )
        !   정규화: 처음 5번 반복 중 최댓값으로 나눔 (STAR-CCM+ 기본 방식)
        !   처음부터 사실상 0인 식(예: 1D 문제의 y 운동량, 반올림 오차 수준)은 0으로 표시
        !   (반올림 오차로 나누면 의미 없는 큰 값이 나오고 수렴 판정도 막힘)
        res = sqrt(res/real(nx*ny, dp))
        if (step <= 5) res0 = max(res0, res)
        do k = 1,nvar
            if (res0(k) > 1.0e-12_dp*maxval(res0)) then
                res_n(k) = res(k)/res0(k)
            else
                res_n(k) = 0.0_dp
            end if
        end do

        write(ures,'(i8,5es14.6)') step, t, res_n

        if (step == 1 .or. mod(step, print_interval) == 0) then
            if (mod(nprint, 20) == 0) then
                print '(a)', '   Iteration        Time    Continuity    X-momentum    Y-momentum        Energy'
            end if
            print '(i12,es12.4,4es14.4)', step, t, res_n
            nprint = nprint + 1
        end if

        if (steady .and. res_tol > 0.0_dp .and. step > 5 .and. maxval(res_n) < res_tol) then
            print '(i12,es12.4,4es14.4)', step, t, res_n
            print '(a,es9.2)', ' converged: all residuals <', res_tol
            exit
        end if
    end do
    close(ures)

    print '(a,i8,a,es12.4)', ' done: step', step, '   t =', t

    call write_tecplot(trim(case_path(trim(output_name)//'_result.dat')), Q, t)
    call write_line(trim(case_path(trim(output_name)//'_line.dat')), Q, line_row)

end program main
