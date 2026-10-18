!> 계산 감시: 잔차 L2 norm 기록 / 화면 출력 / 수렴 판정, 발산 검사, 계산 시간
!>   main 시간 루프에서
!>     call monitor_open(step0, t0)                ! 시작 전 한 번 (재시작이면 step0, t0 = 이어서 시작하는 반복, 시간)
!>     call check_solution(Q, step, t)             ! 갱신 뒤: NaN 이나 음수 밀도면 멈춤
!>     call monitor_update(step, t, R, converged)  ! 갱신 뒤: 잔차 기록, 수렴하면 converged = .true.
!>     call monitor_close()                        ! 끝난 뒤
module monitor
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: output_unit, int64
    use euler
    use grid
    use case_setup, only: steady, res_tol, print_interval, print_inner, output_name, case_path, max_steps, t_end, restart, &
                          t_unit, t_label, dim_out

    implicit none
    private
    public :: monitor_open, monitor_update, monitor_inner, monitor_close, check_solution, wall_seconds, hms
    public :: res0

    integer  :: ures   = -1          ! <output>_residual.dat
    integer  :: nprint = 0           ! 화면에 찍은 줄 수 (20 줄마다 머리말)
    integer  :: hdr_step = -1        ! dual time: Step 줄을 찍은 마지막 시간 스텝
    integer  :: nstep_hdr = 0        ! dual time: 찍은 Step 줄 수 (열 머리말 10 스텝마다)
    real(dp) :: res0(nvar) = 0.0_dp  ! 정규화 기준 (처음 5번 반복의 최댓값, 재시작 파일에 저장)

    ! 계산 시간: 마지막으로 화면에 찍은 때의 반복, 시간, 벽시계
    integer  :: step_last
    real(dp) :: t_last, clock_last

contains

    !> 프로그램 시작부터의 벽시계 초
    real(dp) function wall_seconds()
        integer(int64), save :: c0 = -1
        integer(int64) :: c, rate
        call system_clock(c, rate)
        if (c0 < 0) c0 = c
        wall_seconds = real(c - c0, dp)/real(rate, dp)
    end function wall_seconds

    !> 초 -> 'h:mm:ss'
    function hms(sec) result(s)
        real(dp), intent(in) :: sec
        character(len=16) :: s
        integer :: t
        t = nint(max(sec, 0.0_dp))
        write(s,'(i0,a,i2.2,a,i2.2)') t/3600, ':', mod(t,3600)/60, ':', mod(t,60)
    end function hms

    !> 잔차 기록 파일 <output>_residual.dat 열기 (매 반복 기록, tools/plot_residual.m 으로 그림, 재시작이면 이어 씀)
    !>   dual (선택): dual time 이면 안쪽 반복마다 기록하는 머리말 (iter  step  time ...)
    subroutine monitor_open(step0, t0, dual)
        integer,  intent(in) :: step0
        real(dp), intent(in) :: t0
        logical,  intent(in), optional :: dual
        logical :: ex, dt_mode
        character(len=512) :: fname

        dt_mode = .false.
        if (present(dual)) dt_mode = dual
        fname = case_path(trim(output_name)//'_residual.dat')
        inquire(file=trim(fname), exist=ex)
        if (restart .and. ex) then
            open(newunit=ures, file=trim(fname), status='old', position='append', action='write')
        else if (dt_mode) then
            open(newunit=ures, file=trim(fname), status='replace', action='write')
            write(ures,'(a)') '# dual time: one line per inner iteration' // &
                              '   iter (all inner iterations)  step (time step)  '//trim(merge('time [s]', 'time    ', dim_out)) // &
                              '  continuity  x-momentum  y-momentum  energy   (L2 norm / max of first 5 iters)' // &
                              '   L2-continuity  L2-x-mom  L2-y-mom  L2-energy   (L2 norm, not normalized)'
        else
            open(newunit=ures, file=trim(fname), status='replace', action='write')
            write(ures,'(a)') '# iter  '//trim(merge('time [s]', 'time    ', dim_out))//'  continuity  x-momentum  y-momentum  energy   (L2 norm / max of first 5 iters)' // &
                              '   L2-continuity  L2-x-mom  L2-y-mom  L2-energy   (L2 norm, not normalized)'
        end if
        nprint = 0
        if (.not. restart) res0 = 0.0_dp
        step_last = step0;  t_last = t0;  clock_last = wall_seconds()
    end subroutine monitor_open

    subroutine monitor_close()
        if (ures > 0) close(ures)
        ures = -1
    end subroutine monitor_close

    !> 잔차 L2 norm (식마다):  sqrt( sum(R^2) / 셀 수 )
    !>   정규화: 처음 5번 반복 중 최댓값으로 나눔 (STAR-CCM+ 기본 방식, 수렴 판정은 이 값으로)
    !>   화면과 파일에 정규화 값과 정규화 안 한 L2 norm 을 둘 다 씀 (파일은 L2 norm 이 맨 뒤 4 열)
    !>   처음부터 사실상 0인 식(예: 1D 문제의 y 운동량, 반올림 오차 수준)은 0으로 표시
    !>   (반올림 오차로 나누면 의미 없는 큰 값이 나오고 수렴 판정도 막힘)
    !>   파일에는 매 반복, 화면에는 print_interval 마다 (+ 반복당 시간, 남은 시간). 정상 문제에서 모두 res_tol 아래면 converged
    subroutine monitor_update(step, t, R, converged)
        integer,  intent(in)  :: step
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: R(:,0:,0:)
        logical,  intent(out) :: converged

        real(dp) :: res(nvar), res_n(nvar), now, rate, eta
        character(len=*), parameter :: fmt = '(i10,es11.3,2x,4es11.3,3x,4es10.2,a)'
        character(len=40) :: timing

        call norms(R, step <= 5, res, res_n)

        write(ures,'(i8,9es14.6)') step, t*t_unit, res_n, res

        converged = steady .and. res_tol > 0.0_dp .and. step > 5 .and. maxval(res_n) < res_tol

        if (step == 1 .or. mod(step, print_interval) == 0 .or. converged) then
            ! 반복당 시간 (지난 출력 이후 평균), 남은 시간 (정상: max_iter 까지, 비정상: t_end 까지)
            now  = wall_seconds()
            rate = (now - clock_last)/max(step - step_last, 1)
            if (steady) then
                eta = rate*(max_steps - step)
            else if (t > t_last) then
                eta = (t_end - t)*(now - clock_last)/(t - t_last)
            else
                eta = 0.0_dp
            end if
            write(timing,'(2x,es8.2,a,a)') rate, ' s/it  left ', trim(hms(eta))
            step_last = step;  t_last = t;  clock_last = now

            if (mod(nprint, 20) == 0) then
                print '(23x,a44,3x,a40)', '<-- normalized (L2/max of first 5 iters) -->', &
                                          '<------ L2 norm sqrt(sum R^2/N) ------>'
                print '(a10,a11,2x,4a11,3x,4a10,a)', 'Iteration', trim(t_label), &
                    'cont', 'x-mom', 'y-mom', 'energy', 'cont', 'x-mom', 'y-mom', 'energy', &
                    '   (time per iteration, time left)'
            end if
            print fmt, step, t*t_unit, res_n, res, trim(timing)
            nprint = nprint + 1
            flush(output_unit)                    ! 화면 출력을 파일로 넘겨도 바로 보이게
            flush(ures)
        end if

        if (converged) then
            print '(a,es9.2)', ' converged: all residuals <', res_tol
        end if
    end subroutine monitor_update

    !> 잔차 L2 norm (식마다) sqrt(sum R^2 / 셀 수) 와 정규화 값
    !>   first : 처음 5번 반복이면 정규화 기준 res0 을 갱신 (최댓값)
    !>   처음부터 사실상 0인 식은 정규화 값 0
    subroutine norms(R, first, res, res_n)
        real(dp), intent(in)  :: R(:,0:,0:)
        logical,  intent(in)  :: first
        real(dp), intent(out) :: res(nvar), res_n(nvar)
        integer :: ic, jc, k

        res = 0.0_dp
        !$omp parallel do private(ic) reduction(+:res)
        do jc = 1,ny
            do ic = 1,nx
                res = res + R(:,ic,jc)**2
            end do
        end do
        !$omp end parallel do
        res = sqrt(res/real(nx*ny, dp))

        if (first) res0 = max(res0, res)
        do k = 1,nvar
            if (res0(k) > 1.0e-12_dp*maxval(res0)) then
                res_n(k) = res(k)/res0(k)
            else
                res_n(k) = 0.0_dp
            end if
        end do
    end subroutine norms

    !> dual time (비정상 + implicit): 안쪽 반복 하나의 잔차 R* 기록 (STAR-CCM+ 처럼 Iteration 과 Time Step 을 나눠서)
    !>   iter   : 안쪽 반복 누적 번호 (모든 시간 스텝을 통틀어, 재시작해도 이어짐)
    !>   step   : 시간 스텝 번호,  t : 그 스텝의 실제 시간 (풀고 있는 시간 n+1)
    !>   last   : 이 시간 스텝의 마지막 안쪽 반복이면 .true. (스텝당 시간, 남은 시간을 이 줄에)
    !>   ninner : 이 시간 스텝 안에서 몇 번째 안쪽 반복인지 (화면 inner 열)
    !>   화면: print_interval 스텝마다 그 스텝의 안쪽 반복을 print_inner 간격으로 한 줄씩 (STAR-CCM+ 처럼, 마지막은 항상)
    !>   파일: 안쪽 반복마다 한 줄  (iter  step  time  정규화 4  L2 4).  정규화 기준은 처음 5번 안쪽 반복의 최댓값
    subroutine monitor_inner(iter, step, t, R, last, ninner)
        integer,  intent(in) :: iter, step, ninner
        real(dp), intent(in) :: t
        real(dp), intent(in) :: R(:,0:,0:)
        logical,  intent(in) :: last

        real(dp) :: res(nvar), res_n(nvar), now, rate, eta
        character(len=40) :: timing
        logical :: show

        call norms(R, iter <= 5, res, res_n)
        write(ures,'(i10,i9,9es14.6)') iter, step, t*t_unit, res_n, res

        ! print_interval 스텝이면 그 스텝의 안쪽 반복을 print_inner 간격으로 화면에 (마지막은 항상, 스텝마다 머리말, 시간은 마지막 줄에)
        show = last
        if (print_inner > 0) show = show .or. mod(ninner, print_inner) == 0
        if ((step == 1 .or. mod(step, print_interval) == 0) .and. show) then
            timing = ''
            if (last) then
                ! 시간 스텝당 시간 (지난 출력 이후 평균), 남은 시간 (t_end 까지)
                now  = wall_seconds()
                rate = (now - clock_last)/max(step - step_last, 1)
                eta  = 0.0_dp
                if (t > t_last) eta = (t_end - t)*(now - clock_last)/(t - t_last)
                write(timing,'(2x,es8.2,a,a)') rate, ' s/step  left ', trim(hms(eta))
                step_last = step;  t_last = t;  clock_last = now
            end if

            if (print_inner > 0) then
                ! 안쪽 반복 여러 줄: 시간 스텝 줄 아래에 안쪽 반복을 들여 써서 (열 머리말은 10 스텝마다)
                if (step /= hdr_step) then
                    if (mod(nstep_hdr, 10) == 0) then
                        print '(25x,a44,3x,a40)', '<-- normalized (L2/max of first 5 iters) -->', &
                                                  '<------ L2 norm sqrt(sum R^2/N) ------>'
                        print '(6x,a10,a7,2x,4a11,3x,4a10,a)', 'Iteration', 'inner', &
                            'cont', 'x-mom', 'y-mom', 'energy', 'cont', 'x-mom', 'y-mom', 'energy', &
                            '   (time per step, time left)'
                    end if
                    print '(a,i8,3a,es12.4)', ' Step ', step, '   ', trim(t_label), ' ', t*t_unit
                    hdr_step  = step
                    nstep_hdr = nstep_hdr + 1
                end if
                print '(6x,i10,i7,2x,4es11.3,3x,4es10.2,a)', iter, ninner, res_n, res, trim(timing)
            else
                ! 마지막 줄만: 시간 스텝마다 한 줄 (머리말 20 줄마다)
                if (mod(nprint, 20) == 0) then
                    print '(32x,a44,3x,a40)', '<-- normalized (L2/max of first 5 iters) -->', &
                                              '<------ L2 norm sqrt(sum R^2/N) ------>'
                    print '(a10,a8,a11,2x,4a11,3x,4a10,a7,a)', 'Iteration', 'Step', trim(t_label), &
                        'cont', 'x-mom', 'y-mom', 'energy', 'cont', 'x-mom', 'y-mom', 'energy', 'inner', &
                        '   (time per step, time left)'
                end if
                print '(i10,i8,es11.3,2x,4es11.3,3x,4es10.2,i7,a)', iter, step, t*t_unit, res_n, res, ninner, trim(timing)
            end if
            nprint = nprint + 1
            flush(output_unit)
            flush(ures)
        end if
    end subroutine monitor_inner

    !> 발산 검사: 내부 셀에 NaN 이나 0 이하 밀도가 있으면 알리고 멈춤
    subroutine check_solution(Q, step, t)
        real(dp), intent(in) :: Q(:,1-ng:,1-ng:)
        integer,  intent(in) :: step
        real(dp), intent(in) :: t

        logical :: bad
        integer :: ic, jc

        bad = .false.
        !$omp parallel do private(ic) reduction(.or.:bad)
        do jc = 1,ny
            do ic = 1,nx
                bad = bad .or. any(ieee_is_nan(Q(:,ic,jc))) .or. Q(1,ic,jc) <= 0.0_dp
            end do
        end do
        !$omp end parallel do

        if (bad) then
            print '(a,i0,3a,es12.4)', ' Solution blew up at step ', step, ', ', trim(t_label), ' = ', t*t_unit
            call monitor_close()
            flush(output_unit)
            error stop
        end if
    end subroutine check_solution

end module monitor
