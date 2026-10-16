module solver
    use euler
    use grid
    use case_setup, only: cfl,sigma,kappa,b, flux_id, FLUX_RUSANOV, FLUX_ROE

    implicit none

    real(dp), allocatable :: fluxi(:,:,:,:), fluxj(:,:,:,:), fluxk(:,:,:,:)    ! 면 플럭스 (작업용, residual 전용)

    contains

    !> 면 법선 n 방향 물리 플럭스 F(Q).n
    pure subroutine flux_n(Qc, W, n, F)
        real(dp), intent(in)  :: Qc(nvar), W(nvar), n(3)
        real(dp), intent(out) :: F(nvar)
        real(dp) :: un

        un = W(2)*n(1) + W(3)*n(2) + W(4)*n(3)
        F(1)   = W(1)*un
        F(2:4) = W(1)*W(2:4)*un + W(5)*n
        F(5)   = (Qc(5) + W(5))*un
    end subroutine flux_n

    subroutine rusanov(QL,QR,n,F)

        real(dp), intent(in)  :: QL(nvar), QR(nvar)
        real(dp), intent(in)  :: n(3)
        real(dp), intent(out) :: F(nvar)

        real(dp) :: l_max, WL(nvar), WR(nvar), FL(nvar), FR(nvar), cL, cR

        call c2p(QL,WL)
        call c2p(QR,WR)
        call sos(WL, cL)
        call sos(WR, cR)

        call flux_n(QL, WL, n, FL)
        call flux_n(QR, WR, n, FR)

        l_max = max(abs(dot_product(WL(2:4),n))+cL, abs(dot_product(WR(2:4),n))+cR)

        F = 0.5_dp*(FL+FR) - 0.5_dp* l_max*(QR-QL)

    end subroutine rusanov

    !> 법선 n 에 수직인 단위 벡터 두 개 (t1, t2, n 이 직교 기저)
    !>   n 과 가장 덜 나란한 좌표축과 외적해서 만듦 (2D 의 접선 (-ny, nx) 에 해당)
    pure subroutine tangents(n, t1, t2)
        real(dp), intent(in)  :: n(3)
        real(dp), intent(out) :: t1(3), t2(3)

        if (abs(n(1)) <= abs(n(2)) .and. abs(n(1)) <= abs(n(3))) then
            t1 = [0.0_dp, n(3), -n(2)]            ! n x ex
        else if (abs(n(2)) <= abs(n(3))) then
            t1 = [-n(3), 0.0_dp, n(1)]            ! n x ey
        else
            t1 = [n(2), -n(1), 0.0_dp]            ! n x ez
        end if
        t1 = t1/norm2(t1)
        t2 = cross(n, t1)
    end subroutine tangents

    !> MUSCL 재구성 (특성 변수), 인자 순서: 가운데 Q (면 왼쪽), 그 왼쪽 Qm1, 면 오른쪽 Q1, 그 오른쪽 Q2
    subroutine MUSCL(Q,Qm1,Q1,Q2,n,QL,QR)
        real(dp), intent(in) :: Q(:),Qm1(:),Q1(:),Q2(:),n(3)
        real(dp), intent(out) :: QL(:), QR(:)

        real(dp) :: W(nvar), Wm1(nvar), W1(nvar), W2(nvar),WL(nvar),WR(nvar)
        real(dp) :: deltaWmh(nvar),deltaWph(nvar),deltaWpph(nvar)
        real(dp) :: am(nvar), a0(nvar), ap(nvar), sL(nvar), sR(nvar), dL(nvar), dR(nvar)
        real(dp) :: rhof, cf, c0, c1, t1(3), t2(3)

        call c2p(Q,W)
        call c2p(Qm1,Wm1)
        call c2p(Q1,W1)
        call c2p(Q2,W2)

        !delta 계산
        deltaWmh  = W  - Wm1
        deltaWph  = W1 - W
        deltaWpph = W2 - W1

        call sos(W,c0)
        call sos(W1,c1)
        rhof = 0.5_dp*(W(1)+W1(1))
        cf   = 0.5_dp*(c0+c1)
        call tangents(n, t1, t2)

        call to_char(deltaWmh, rhof,cf,n,t1,t2,am)
        call to_char(deltaWph, rhof,cf,n,t1,t2,a0)
        call to_char(deltaWpph,rhof,cf,n,t1,t2,ap)

        sL = (1-kappa)*minmod(am,b*a0)+(1+kappa)*minmod(a0,b*am)
        sR = (1-kappa)*minmod(ap,b*a0)+(1+kappa)*minmod(a0,b*ap)

        call from_char(sL,rhof,cf,n,t1,t2,dL)
        call from_char(sR,rhof,cf,n,t1,t2,dR)

        WL = W +dL*(sigma/4.0_dp)
        WR = W1-dR*(sigma/4.0_dp)

        if (WL(1) <= 0.0_dp .or. WL(5) <= 0.0_dp .or. WR(1) <= 0.0_dp .or. WR(5) <= 0.0_dp) then
            WL = W
            WR = W1
        end if

        call p2c(WL,QL)
        call p2c(WR,QR)

    end subroutine MUSCL

    !> 원시 변수 차이 dW -> 특성 변수 a  (파동: un-c, 엔트로피, 전단 t1, 전단 t2, un+c)
    subroutine to_char(dW,rho,c,n,t1,t2,a)
        real(dp), intent(in)  :: dW(nvar), rho, c, n(3), t1(3), t2(3)
        real(dp), intent(out) :: a(nvar)
        real(dp) :: dun

        dun = dot_product(dW(2:4), n)

        a(1) = (dW(5) - rho*c*dun)/(2.0_dp*c**2)
        a(2) =  dW(1) - dW(5)/c**2
        a(3) =  rho*dot_product(dW(2:4), t1)
        a(4) =  rho*dot_product(dW(2:4), t2)
        a(5) = (dW(5) + rho*c*dun)/(2.0_dp*c**2)

    end subroutine to_char

    subroutine from_char(a,rho,c,n,t1,t2,dW)
        real(dp), intent(in)  :: a(nvar), rho, c, n(3), t1(3), t2(3)
        real(dp), intent(out) :: dW(nvar)
        real(dp) :: dun

        dun = c/rho*(a(5) - a(1))

        dW(1)   = a(1) + a(2) + a(5)
        dW(2:4) = dun*n + (a(3)*t1 + a(4)*t2)/rho
        dW(5)   = c**2*(a(1) + a(5))

    end subroutine from_char

    !> Roe 플럭스 (3D)
    !>   전단파 두 개는 접선 기저 없이 한 번에:  rho~ (du - n dun) 로 씀 (2D 의 a3 r3 와 같은 것)
    subroutine Roe(QL,QR,n,F)
        real(dp), intent(in) :: QL(:), QR(:), n(3)
        real(dp), intent(out) :: F(:)

        real(dp) :: WL(nvar), WR(nvar), FL(nvar), FR(nvar)
        real(dp) :: rho_tild, vel_tild(3), H_tild, c_tild, q2_tild, Rw, HL, HR, unL, unR, un_tild
        real(dp) :: dp_, drho, du(3), dun, dut(3)
        real(dp) :: a1, a2, a5, l1, l2, l5, r1(nvar), r2(nvar), r5(nvar), rs(nvar)

        call c2p(QL,WL)
        call c2p(QR,WR)

        HL = (QL(5) + WL(5))/WL(1)
        HR = (QR(5) + WR(5))/WR(1)

        ! Roe 평균
        Rw = sqrt(WR(1)/WL(1))
        rho_tild = sqrt(WL(1)*WR(1))
        vel_tild = (WL(2:4)+Rw*WR(2:4))/(1+Rw)
        H_tild   = (HL + Rw*HR)/(1+Rw)
        q2_tild  = dot_product(vel_tild, vel_tild)
        c_tild   = sqrt((gam-1)*(H_tild-0.5_dp*q2_tild))

        !법선 속도
        unL = dot_product(WL(2:4), n)
        unR = dot_product(WR(2:4), n)
        un_tild = dot_product(vel_tild, n)

        drho = WR(1) - WL(1)
        dp_  = WR(5) - WL(5)
        du   = WR(2:4) - WL(2:4)
        dun  = unR - unL
        dut  = du - dun*n

        a1 = (dp_ - rho_tild*c_tild*dun)/(2*c_tild**2)
        a2 = drho - dp_/c_tild**2
        a5 = (dp_ + rho_tild*c_tild*dun)/(2*c_tild**2)

        l1 = un_tild-c_tild
        l2 = un_tild
        l5 = un_tild+c_tild

        r1(1) = 1; r1(2:4) = vel_tild - c_tild*n; r1(5) = H_tild - c_tild*un_tild
        r2(1) = 1; r2(2:4) = vel_tild;            r2(5) = 0.5_dp*q2_tild
        r5(1) = 1; r5(2:4) = vel_tild + c_tild*n; r5(5) = H_tild + c_tild*un_tild
        rs(1) = 0; rs(2:4) = rho_tild*dut;        rs(5) = rho_tild*dot_product(vel_tild, dut)

        call flux_n(QL, WL, n, FL)
        call flux_n(QR, WR, n, FR)

        F = 0.5_dp*(FL+FR) - 0.5_dp*(abs(l1)*a1*r1 + abs(l2)*(a2*r2 + rs) + abs(l5)*a5*r5)

    end subroutine Roe

    !> 리만 솔버 고르기 (flow.in 의 flux = rusanov | roe)
    subroutine riemann(QL, QR, n, F)
        real(dp), intent(in)  :: QL(nvar), QR(nvar), n(3)
        real(dp), intent(out) :: F(nvar)
        select case (flux_id)
        case (FLUX_ROE)
            call Roe(QL, QR, n, F)
        case default
            call rusanov(QL, QR, n, F)
        end select
    end subroutine riemann

    !> minmod 리미터: 부호가 같으면 절댓값이 작은 쪽, 다르면 (극값) 0
    elemental function minmod(a,b) result(m)
        real(dp),intent(in) :: a,b
        real(dp) :: m

        if (a*b <= 0.0_dp) then
            m = 0.0_dp
        else
            m = sign(min(abs(a),abs(b)), a)
        end if
    end function  minmod

    subroutine residual(Q,R)
        real(dp), intent(inout) :: Q(:,1-ng:,1-ng:,1-ng:)
        real(dp), intent(out) :: R(nvar,nx,ny,nz)

        integer :: fi, fj, fk, ic, jc, kc
        real(dp) :: F(nvar), QL(nvar), QR(nvar)

        if (.not. allocated(fluxi)) then
            allocate(fluxi(nvar,0:nx,ny,nz), fluxj(nvar,nx,0:ny,nz), fluxk(nvar,nx,ny,0:nz))
        end if

        ! OpenMP: 바깥 두 루프를 묶어 스레드에 나눔 (collapse(2): 한 방향 셀이 적어도 고르게 나뉨)
        !   F, QL, QR 은 스레드마다 따로 (private). fluxi/j/k, R 은 스레드마다 다른 위치에 씀
        !$omp parallel private(fi, fj, fk, ic, jc, kc, F, QL, QR)

        ! i-면: fi = 면 번호 (셀 fi 와 fi+1 사이)
        !   MUSCL 인자 순서: 가운데(fi), 왼쪽(fi-1), 오른쪽(fi+1), 그 오른쪽(fi+2)
        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do fi = 0,nx
                    call MUSCL(Q(:,fi,jc,kc), Q(:,fi-1,jc,kc), Q(:,fi+1,jc,kc), Q(:,fi+2,jc,kc), si(:,fi,jc,kc), QL, QR)
                    call riemann(QL, QR, si(:,fi,jc,kc), F)
                    fluxi(:,fi,jc,kc) = F*si_a(fi,jc,kc)
                end do
            end do
        end do
        !$omp end do nowait

        ! j-면
        !$omp do collapse(2)
        do kc = 1,nz
            do fj = 0,ny
                do ic = 1,nx
                    call MUSCL(Q(:,ic,fj,kc), Q(:,ic,fj-1,kc), Q(:,ic,fj+1,kc), Q(:,ic,fj+2,kc), sj(:,ic,fj,kc), QL, QR)
                    call riemann(QL, QR, sj(:,ic,fj,kc), F)
                    fluxj(:,ic,fj,kc) = F*sj_a(ic,fj,kc)
                end do
            end do
        end do
        !$omp end do nowait

        ! k-면
        !$omp do collapse(2)
        do fk = 0,nz
            do jc = 1,ny
                do ic = 1,nx
                    call MUSCL(Q(:,ic,jc,fk), Q(:,ic,jc,fk-1), Q(:,ic,jc,fk+1), Q(:,ic,jc,fk+2), sk(:,ic,jc,fk), QL, QR)
                    call riemann(QL, QR, sk(:,ic,jc,fk), F)
                    fluxk(:,ic,jc,fk) = F*sk_a(ic,jc,fk)
                end do
            end do
        end do
        !$omp end do
        ! (여기서 모든 스레드가 기다림: R 은 플럭스가 전부 채워진 뒤 계산)

        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    R(:,ic,jc,kc) = -( fluxi(:,ic,jc,kc) - fluxi(:,ic-1,jc,kc) &
                                     + fluxj(:,ic,jc,kc) - fluxj(:,ic,jc-1,kc) &
                                     + fluxk(:,ic,jc,kc) - fluxk(:,ic,jc,kc-1) ) / vol(ic,jc,kc)
                end do
            end do
        end do
        !$omp end do

        !$omp end parallel

    end subroutine residual

    !> CFL 조건으로 dt 계산
    !>   셀마다 dt = cfl * V / (lam_i + lam_j + lam_k)
    !>   lam_i = |u.Si| + c|Si|,  Si = 셀 양쪽 i-면 벡터(법선 x 넓이)의 평균  (j, k 도 같음)
    !>   dtl(ic,jc,kc) : 셀마다의 dt (정상 문제의 국소 시간 간격),  dt : 그 최솟값 (비정상 문제의 전체 dt)
    subroutine compute_dt(Q, dtl, dt)
        real(dp), intent(in)  :: Q(:,1-ng:,1-ng:,1-ng:)
        real(dp), intent(out) :: dtl(:,:,:)
        real(dp), intent(out) :: dt

        integer  :: ic, jc, kc
        real(dp) :: W(nvar), c, s_i(3), s_j(3), s_k(3), lam

        dt = huge(1.0_dp)

        !$omp parallel do collapse(2) private(ic, W, c, s_i, s_j, s_k, lam) reduction(min:dt)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    call c2p(Q(:,ic,jc,kc), W)
                    call sos(W, c)

                    s_i = 0.5_dp*(si(:,ic-1,jc,kc)*si_a(ic-1,jc,kc) + si(:,ic,jc,kc)*si_a(ic,jc,kc))
                    s_j = 0.5_dp*(sj(:,ic,jc-1,kc)*sj_a(ic,jc-1,kc) + sj(:,ic,jc,kc)*sj_a(ic,jc,kc))
                    s_k = 0.5_dp*(sk(:,ic,jc,kc-1)*sk_a(ic,jc,kc-1) + sk(:,ic,jc,kc)*sk_a(ic,jc,kc))

                    lam = abs(dot_product(W(2:4), s_i)) + c*norm2(s_i) &
                        + abs(dot_product(W(2:4), s_j)) + c*norm2(s_j) &
                        + abs(dot_product(W(2:4), s_k)) + c*norm2(s_k)

                    dtl(ic,jc,kc) = cfl*vol(ic,jc,kc)/lam
                    dt = min(dt, dtl(ic,jc,kc))
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine compute_dt

end module solver
