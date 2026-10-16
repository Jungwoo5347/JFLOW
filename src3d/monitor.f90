!> 계산 감시: 잔차 L2 norm 기록 / 화면 출력 / 수렴 판정, 발산 검사, 계산 시간  (3D, 모든 블록)
!>   main 시간 루프에서
!>     call monitor_open(step0)                    ! 시작 전 한 번 (재시작이면 step0 = 이어서 시작하는 반복)
!>     call check_solution(step, t)                ! 갱신 뒤: NaN 이나 음수 밀도면 위치를 알리고 멈춤
!>     call monitor_update(step, t, converged)     ! 갱신 뒤: 잔차 (blk%R) 기록, 수렴하면 converged = .true.
!>     call monitor_close()                        ! 끝난 뒤
module monitor
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: output_unit, int64
    use euler
    use grid
    use case_setup, only: steady, res_tol, print_interval, output_name, case_path, max_steps, t_end, restart

    implicit none
    private
    public :: monitor_open, monitor_update, monitor_close, check_solution, wall_seconds, hms
    public :: res0

    integer  :: ures   = -1          ! <output>_residual.dat
    integer  :: nprint = 0           ! 화면에 찍은 줄 수 (20 줄마다 머리말)
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

    !> 잔차 기록 파일 <output>_residual.dat 열기 (매 반복 기록, 재시작이면 이어 씀)
    subroutine monitor_open(step0, t0)
        integer,  intent(in) :: step0
        real(dp), intent(in) :: t0
        logical :: ex
        character(len=512) :: fname

        fname = case_path(trim(output_name)//'_residual.dat')
        inquire(file=trim(fname), exist=ex)
        if (restart .and. ex) then
            open(newunit=ures, file=trim(fname), status='old', position='append', action='write')
        else
            open(newunit=ures, file=trim(fname), status='replace', action='write')
            write(ures,'(a)') '# iter  time  continuity  x-momentum  y-momentum  z-momentum  energy' // &
                              '   (L2 norm / max of first 5 iters)' // &
                              '   L2-continuity  L2-x-mom  L2-y-mom  L2-z-mom  L2-energy   (L2 norm, not normalized)'
        end if
        nprint = 0
        if (.not. restart) res0 = 0.0_dp
        step_last = step0;  t_last = t0;  clock_last = wall_seconds()
    end subroutine monitor_open

    subroutine monitor_close()
        if (ures > 0) close(ures)
        ures = -1
    end subroutine monitor_close

    !> 잔차 L2 norm (식마다, 모든 블록):  sqrt( sum(R^2) / 셀 수 ),  R = blk(b)%R (RK 첫 단계)
    !>   정규화: 처음 5번 반복 중 최댓값으로 나눔 (수렴 판정은 이 값으로)
    !>   화면과 파일에 정규화 값과 정규화 안 한 L2 norm 을 둘 다 씀 (파일은 L2 norm 이 맨 뒤 5 열)
    !>   처음부터 사실상 0인 식(예: 2D 문제의 z 운동량)은 0으로 표시
    !>   파일에는 매 반복, 화면에는 print_interval 마다 (+ 반복당 시간, 남은 시간). 정상 문제에서 모두 res_tol 아래면 converged
    subroutine monitor_update(step, t, converged)
        integer,  intent(in)  :: step
        real(dp), intent(in)  :: t
        logical,  intent(out) :: converged

        real(dp) :: res(nvar), res_n(nvar), now, rate, eta
        integer  :: ic, jc, kc, k, b, mx, my, mz
        integer(int64) :: ncell
        real(dp), pointer, contiguous :: R(:,:,:,:)
        character(len=*), parameter :: fmt = '(i10,es11.3,2x,5es11.3,3x,5es10.2,a)'
        character(len=40) :: timing

        res = 0.0_dp;  ncell = 0
        do b = 1,nblk
            R => blk(b)%R;  mx = blk(b)%nx;  my = blk(b)%ny;  mz = blk(b)%nz
                !$omp parallel do collapse(2) private(ic) reduction(+:res)
                do kc = 1,mz
                    do jc = 1,my
                        do ic = 1,mx
                            res = res + R(:,ic,jc,kc)**2
                        end do
                    end do
                end do
                !$omp end parallel do
                ncell = ncell + int(mx,int64)*my*mz
        end do
        res = sqrt(res/real(ncell, dp))

        if (step <= 5) res0 = max(res0, res)
        do k = 1,nvar
            if (res0(k) > 1.0e-12_dp*maxval(res0)) then
                res_n(k) = res(k)/res0(k)
            else
                res_n(k) = 0.0_dp
            end if
        end do

        write(ures,'(i8,11es14.6)') step, t, res_n, res

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
                print '(23x,a55,3x,a50)', '<------ normalized (L2 / max of first 5 iters) ------>', &
                                          '<----------- L2 norm sqrt(sum R^2/N) ----------->'
                print '(a10,a11,2x,5a11,3x,5a10,a)', 'Iteration', 'Time', &
                    'cont', 'x-mom', 'y-mom', 'z-mom', 'energy', 'cont', 'x-mom', 'y-mom', 'z-mom', 'energy', &
                    '   (time per iteration, time left)'
            end if
            print fmt, step, t, res_n, res, trim(timing)
            nprint = nprint + 1
            flush(output_unit)                    ! 화면 출력을 파일로 넘겨도 바로 보이게
            flush(ures)
        end if

        if (converged) then
            print '(a,es9.2)', ' converged: all residuals <', res_tol
        end if
    end subroutine monitor_update

    !> 발산 검사: 내부 셀에 NaN 이나 0 이하 밀도가 있으면 첫 셀의 위치를 알리고 멈춤
    subroutine check_solution(step, t)
        integer,  intent(in) :: step
        real(dp), intent(in) :: t

        logical :: bad
        integer :: ic, jc, kc, b, mx, my, mz
        real(dp), pointer, contiguous :: Q(:,:,:,:)

        do b = 1,nblk
            bad = .false.
            Q => blk(b)%Q;  mx = blk(b)%nx;  my = blk(b)%ny;  mz = blk(b)%nz
                !$omp parallel do collapse(2) private(ic) reduction(.or.:bad)
                do kc = 1,mz
                    do jc = 1,my
                        do ic = 1,mx
                            bad = bad .or. any(ieee_is_nan(Q(:,ic,jc,kc))) .or. Q(1,ic,jc,kc) <= 0.0_dp
                        end do
                    end do
                end do
                !$omp end parallel do
                if (.not. bad) cycle

                print '(a,i0,a,es12.4)', ' Solution blew up at step ', step, ', t = ', t
                do kc = 1,mz
                    do jc = 1,my
                        do ic = 1,mx
                            if (any(ieee_is_nan(Q(:,ic,jc,kc))) .or. Q(1,ic,jc,kc) <= 0.0_dp) then
                                print '(a,i0,a,i0,a,i0,a,i0,a,3es12.4,a)', '   first bad cell: block ', b, &
                                    '  (i, j, k) = (', ic, ', ', jc, ', ', kc, ')   x y z =', &
                                    blk(b)%xc(ic,jc,kc), blk(b)%yc(ic,jc,kc), blk(b)%zc(ic,jc,kc)
                                call monitor_close()
                                flush(output_unit)
                                error stop
                            end if
                        end do
                    end do
                end do
        end do
    end subroutine check_solution

end module monitor
