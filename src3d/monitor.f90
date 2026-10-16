!> 계산 감시: 잔차 L2 norm 기록 / 화면 출력 / 수렴 판정, 발산 검사  (3D)
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

    !> 잔차 기록 파일 <output>_residual.dat 열기 (매 반복 기록)
    subroutine monitor_open()
        open(newunit=ures, file=trim(case_path(trim(output_name)//'_residual.dat')), status='replace', action='write')
        write(ures,'(a)') '# iter  time  continuity  x-momentum  y-momentum  z-momentum  energy' // &
                          '   (L2 norm / max of first 5 iters)' // &
                          '   L2-continuity  L2-x-mom  L2-y-mom  L2-z-mom  L2-energy   (L2 norm, not normalized)'
        nprint = 0
        res0   = 0.0_dp
    end subroutine monitor_open

    subroutine monitor_close()
        if (ures > 0) close(ures)
        ures = -1
    end subroutine monitor_close

    !> 잔차 L2 norm (식마다):  sqrt( sum(R^2) / 셀 수 )
    !>   정규화: 처음 5번 반복 중 최댓값으로 나눔 (수렴 판정은 이 값으로)
    !>   화면과 파일에 정규화 값과 정규화 안 한 L2 norm 을 둘 다 씀 (파일은 L2 norm 이 맨 뒤 5 열)
    !>   처음부터 사실상 0인 식(예: 2D 문제의 z 운동량)은 0으로 표시
    !>   파일에는 매 반복, 화면에는 print_interval 마다. 정상 문제에서 모두 res_tol 아래면 converged
    subroutine monitor_update(step, t, R, converged)
        integer,  intent(in)  :: step
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: R(:,:,:,:)
        logical,  intent(out) :: converged

        real(dp) :: res(nvar), res_n(nvar)
        integer  :: ic, jc, kc, k
        character(len=*), parameter :: fmt = '(i10,es11.3,2x,5es11.3,3x,5es10.2)'

        res = 0.0_dp
        !$omp parallel do collapse(2) private(ic) reduction(+:res)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    res = res + R(:,ic,jc,kc)**2
                end do
            end do
        end do
        !$omp end parallel do
        res = sqrt(res/real(nx*ny*nz, dp))

        if (step <= 5) res0 = max(res0, res)
        do k = 1,nvar
            if (res0(k) > 1.0e-12_dp*maxval(res0)) then
                res_n(k) = res(k)/res0(k)
            else
                res_n(k) = 0.0_dp
            end if
        end do

        write(ures,'(i8,11es14.6)') step, t, res_n, res

        if (step == 1 .or. mod(step, print_interval) == 0) then
            if (mod(nprint, 20) == 0) then
                print '(23x,a55,3x,a50)', '<------ normalized (L2 / max of first 5 iters) ------>', &
                                          '<----------- L2 norm sqrt(sum R^2/N) ----------->'
                print '(a10,a11,2x,5a11,3x,5a10)', 'Iteration', 'Time', &
                    'cont', 'x-mom', 'y-mom', 'z-mom', 'energy', 'cont', 'x-mom', 'y-mom', 'z-mom', 'energy'
            end if
            print fmt, step, t, res_n, res
            nprint = nprint + 1
        end if

        converged = steady .and. res_tol > 0.0_dp .and. step > 5 .and. maxval(res_n) < res_tol
        if (converged) then
            print fmt, step, t, res_n, res
            print '(a,es9.2)', ' converged: all residuals <', res_tol
        end if
    end subroutine monitor_update

    !> 발산 검사: 내부 셀에 NaN 이나 0 이하 밀도가 있으면 알리고 멈춤
    subroutine check_solution(Q, step, t)
        real(dp), intent(in) :: Q(:,1-ng:,1-ng:,1-ng:)
        integer,  intent(in) :: step
        real(dp), intent(in) :: t

        logical :: bad
        integer :: ic, jc, kc

        bad = .false.
        !$omp parallel do collapse(2) private(ic) reduction(.or.:bad)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    bad = bad .or. any(ieee_is_nan(Q(:,ic,jc,kc))) .or. Q(1,ic,jc,kc) <= 0.0_dp
                end do
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
