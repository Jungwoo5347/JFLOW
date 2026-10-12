module solver
    use euler
    use grid
    use case_setup, only: cfl                              ! CFL 수 (입력 파일 &control 에서)

    implicit none

    real(dp), allocatable :: fluxi(:,:,:), fluxj(:,:,:)    ! 면 플럭스 (작업용, residual 전용)

    contains

    subroutine rusanov(QL,QR,fnx,fny,F)

        real(dp), intent(in)  :: QL(nvar), QR(nvar)   
        real(dp), intent(in)  :: fnx, fny      
        real(dp), intent(out) :: F(nvar)   

        real(dp) :: l_max, WL(nvar), WR(nvar), FL(nvar), FR(nvar), unL, unR, cL, cR

        call c2p(QL,WL)
        call c2p(QR,WR)
        call sos(WL, cL)
        call sos(WR, cR)
        
        unL = WL(2)*fnx+WL(3)*fny
        unR = WR(2)*fnx+WR(3)*fny

        FL(1) = WL(1)*unL
        FL(2) = WL(1)*WL(2)*unL + WL(4)*fnx
        FL(3) = WL(1)*WL(3)*unL + WL(4)*fny
        FL(4) = (QL(4)+WL(4))*unL

        FR(1) = WR(1)*unR
        FR(2) = WR(1)*WR(2)*unR + WR(4)*fnx
        FR(3) = WR(1)*WR(3)*unR + WR(4)*fny
        FR(4) = (QR(4)+WR(4))*unR

        l_max = max(abs(unL)+cL, abs(unR)+cR)

        F = 0.5_dp*(FL+FR) - 0.5_dp* l_max*(QR-QL)
        

    end subroutine rusanov


    subroutine residual(Q,R)
        real(dp), intent(inout) :: Q(:,0:,0:)
        real(dp), intent(out) :: R(nvar,0:nx+1,0:ny+1)

        integer :: fi, fj, ic, jc
        real(dp) :: F(nvar)

        if (.not. allocated(fluxi)) then
            allocate(fluxi(nvar, 0:nx, 1:ny), fluxj(nvar, 1:nx, 0:ny))
        end if

        R = 0.0_dp


        ! OpenMP: 바깥 루프(줄)를 스레드에 나눔. F 는 스레드마다 따로 (private)
        !   fluxi, fluxj, R 은 스레드마다 다른 위치에 쓰므로 공유해도 됨
        !$omp parallel private(fi, fj, ic, jc, F)

        ! i-면 루프: fi = 면 번호 (셀 fi 와 fi+1 사이), jc = 셀 줄
        !$omp do
        do jc = 1,ny
            do fi = 0,nx
                call rusanov(Q(:,fi,jc),Q(:,fi+1,jc),si_x(fi,jc), si_y(fi,jc), F)
                fluxi(:,fi,jc) = F*si_len(fi,jc)
            end do
        end do
        !$omp end do nowait

        ! j-면 루프: fj = 면 번호 (셀 fj 와 fj+1 사이), ic = 셀 칸
        !$omp do
        do fj = 0,ny
            do ic = 1,nx
                call rusanov(Q(:,ic,fj), Q(:,ic,fj+1), sj_x(ic,fj), sj_y(ic,fj), F)
                fluxj(:,ic,fj) = F*sj_len(ic,fj)
            end do
        end do
        !$omp end do
        ! (여기서 모든 스레드가 기다림: R 은 fluxi, fluxj 가 전부 채워진 뒤 계산)

        !$omp do
        do jc = 1, ny
            do ic = 1, nx
                R(:,ic,jc) = -( fluxi(:,ic,jc) - fluxi(:,ic-1,jc)&
                            + fluxj(:,ic,jc) - fluxj(:,ic,jc-1) ) / area(ic,jc)
            end do
        end do
        !$omp end do

        !$omp end parallel


    end subroutine residual


    !> CFL 조건으로 dt 계산
    !>   셀마다 dt = cfl * A / (lam_i + lam_j)
    !>   lam_i = |u·Si| + c|Si|,  Si = 셀 양쪽 i-면 벡터(법선 × 길이)의 평균  (j도 동일)
    !>   직교 격자에서는 dt = cfl / ((|u|+c)/dx + (|v|+c)/dy) 와 같음
    subroutine compute_dt(Q, dt)
        real(dp), intent(in)  :: Q(:,0:,0:)
        real(dp), intent(out) :: dt

        integer  :: ic, jc
        real(dp) :: W(nvar), c, six, siy, sjx, sjy, lam_i, lam_j

        dt = huge(1.0_dp)

        ! OpenMP: 스레드마다 자기 최솟값을 구한 뒤 마지막에 합침 (reduction)
        !$omp parallel do private(ic, W, c, six, siy, sjx, sjy, lam_i, lam_j) reduction(min:dt)
        do jc = 1,ny
            do ic = 1,nx
                call c2p(Q(:,ic,jc), W)
                call sos(W, c)

                six = 0.5_dp*(si_x(ic-1,jc)*si_len(ic-1,jc) + si_x(ic,jc)*si_len(ic,jc))
                siy = 0.5_dp*(si_y(ic-1,jc)*si_len(ic-1,jc) + si_y(ic,jc)*si_len(ic,jc))
                sjx = 0.5_dp*(sj_x(ic,jc-1)*sj_len(ic,jc-1) + sj_x(ic,jc)*sj_len(ic,jc))
                sjy = 0.5_dp*(sj_y(ic,jc-1)*sj_len(ic,jc-1) + sj_y(ic,jc)*sj_len(ic,jc))

                lam_i = abs(W(2)*six + W(3)*siy) + c*sqrt(six**2 + siy**2)
                lam_j = abs(W(2)*sjx + W(3)*sjy) + c*sqrt(sjx**2 + sjy**2)

                dt = min(dt, cfl*area(ic,jc)/(lam_i + lam_j))
            end do
        end do
        !$omp end parallel do

    end subroutine compute_dt

end module solver