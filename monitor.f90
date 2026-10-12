!> 계산 감시: 잔차 L2 norm 기록 / 화면 출력 / 수렴 판정, 발산 검사
!>   main 시간 루프에서
!>     call monitor_open()                         ! 시작 전 한 번
!>     call check_solution(Q, step, t)             ! 갱신 뒤: NaN 이나 음수 밀도면 멈춤
!>     call monitor_update(step, t, R, converged)  ! 갱신 뒤: 잔차 기록, 수렴하면 converged = .true.
!>     call monitor_close()                        ! 끝난 뒤
module monitor
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use euler
    use grid
    use case_setup, only: steady, res_tol, print_interval, output_name, case_path

    implicit none
    private
    public :: monitor_open, monitor_update, monitor_close, check_solution

    integer  :: ures   = -1          ! <output>_residual.dat
    integer  :: nprint = 0           ! 화면에 찍은 줄 수 (20 줄마다 머리말)
    real(dp) :: res0(nvar) = 0.0_dp  ! 정규화 기준 (처음 5번 반복의 최댓값)

contains

    !> 잔차 기록 파일 <output>_residual.dat 열기 (매 반복 기록, tools/plot_residual.m 으로 그림)
    subroutine monitor_open()
        open(newunit=ures, file=trim(case_path(trim(output_name)//'_residual.dat')), status='replace', action='write')
        write(ures,'(a)') '# iter  time  continuity  x-momentum  y-momentum  energy   (L2 norm / max of first 5 iters)'
        nprint = 0
        res0   = 0.0_dp
    end subroutine monitor_open

    subroutine monitor_close()
        if (ures > 0) close(ures)
        ures = -1
    end subroutine monitor_close

    !> 잔차 L2 norm (식마다):  sqrt( sum(R^2) / 셀 수 )
    !>   정규화: 처음 5번 반복 중 최댓값으로 나눔 (STAR-CCM+ 기본 방식)
    !>   처음부터 사실상 0인 식(예: 1D 문제의 y 운동량, 반올림 오차 수준)은 0으로 표시
    !>   (반올림 오차로 나누면 의미 없는 큰 값이 나오고 수렴 판정도 막힘)
    !>   파일에는 매 반복, 화면에는 print_interval 마다. 정상 문제에서 모두 res_tol 아래면 converged
    subroutine monitor_update(step, t, R, converged)
        integer,  intent(in)  :: step
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: R(:,0:,0:)
        logical,  intent(out) :: converged

        real(dp) :: res(nvar), res_n(nvar)
        integer  :: ic, jc, k

        res = 0.0_dp
        !$omp parallel do private(ic) reduction(+:res)
        do jc = 1,ny
            do ic = 1,nx
                res = res + R(:,ic,jc)**2
            end do
        end do
        !$omp end parallel do
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

        converged = steady .and. res_tol > 0.0_dp .and. step > 5 .and. maxval(res_n) < res_tol
        if (converged) then
            print '(i12,es12.4,4es14.4)', step, t, res_n
            print '(a,es9.2)', ' converged: all residuals <', res_tol
        end if
    end subroutine monitor_update

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
            print '(a,i0,a,es12.4)', ' Solution blew up at step ', step, ', t = ', t
            call monitor_close()
            error stop
        end if
    end subroutine check_solution

end module monitor
