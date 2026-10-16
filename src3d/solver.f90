module solver
    use euler
    use grid
    use case_setup, only: cfl,sigma,kappa,b, flux_id, FLUX_RUSANOV, FLUX_ROE, entropy_fix, entropy_fix_all

    implicit none

    ! 면 플럭스 작업 버퍼 (residual 전용): i, j, k 방향을 차례로 같은 버퍼에 담음 (모든 블록이 같이 씀)
    real(dp), allocatable, target :: fbuf(:)

    contains

    !> 면 법선 n 방향 물리 플럭스 F(W).n   (W = rho, u, v, w, p)
    pure subroutine flux_n(W, n, F)
        real(dp), intent(in)  :: W(nvar), n(3)
        real(dp), intent(out) :: F(nvar)
        real(dp) :: un

        un = W(2)*n(1) + W(3)*n(2) + W(4)*n(3)
        F(1)   = W(1)*un
        F(2:4) = W(1)*W(2:4)*un + W(5)*n
        F(5)   = (gam/(gam - 1.0_dp)*W(5) + 0.5_dp*W(1)*(W(2)**2 + W(3)**2 + W(4)**2))*un    ! (E + p) un
    end subroutine flux_n

    subroutine rusanov(WL,WR,n,F)

        real(dp), intent(in)  :: WL(nvar), WR(nvar)
        real(dp), intent(in)  :: n(3)
        real(dp), intent(out) :: F(nvar)

        real(dp) :: l_max, QL(nvar), QR(nvar), FL(nvar), FR(nvar), cL, cR

        call sos(WL, cL)
        call sos(WR, cR)
        call p2c(WL, QL)
        call p2c(WR, QR)

        call flux_n(WL, n, FL)
        call flux_n(WR, n, FR)

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

    !> MUSCL 재구성 (특성 변수), 원시 변수로 받고 원시 변수로 돌려줌
    !>   인자 순서: 가운데 W (면 왼쪽), 그 왼쪽 Wm1, 면 오른쪽 W1, 그 오른쪽 W2  ->  면 양쪽 값 WL, WR
    subroutine MUSCL(W,Wm1,W1,W2,n,WL,WR)
        real(dp), intent(in)  :: W(:),Wm1(:),W1(:),W2(:),n(3)
        real(dp), intent(out) :: WL(:), WR(:)

        real(dp) :: deltaWmh(nvar),deltaWph(nvar),deltaWpph(nvar)
        real(dp) :: am(nvar), a0(nvar), ap(nvar), sL(nvar), sR(nvar), dL(nvar), dR(nvar)
        real(dp) :: rhof, cf, c0, c1, t1(3), t2(3)

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

    !> Roe 플럭스 (3D), 면 양쪽 원시 변수 WL, WR
    !>   전단파 두 개는 접선 기저 없이 한 번에:  rho~ (du - n dun) 로 씀 (2D 의 a3 r3 와 같은 것)
    subroutine Roe(WL,WR,n,F)
        real(dp), intent(in) :: WL(:), WR(:), n(3)
        real(dp), intent(out) :: F(:)

        real(dp) :: FL(nvar), FR(nvar)
        real(dp) :: rho_tild, vel_tild(3), H_tild, c_tild, q2_tild, Rw, HL, HR, unL, unR, un_tild
        real(dp) :: dp_, drho, du(3), dun, dut(3)
        real(dp) :: a1, a2, a5, l1, l2, l5, r1(nvar), r2(nvar), r5(nvar), rs(nvar)
        real(dp) :: al1, al2, al5, delta

        ! 전 엔탈피 H = (E + p)/rho = g/(g-1) p/rho + |u|^2/2
        HL = gam/(gam - 1.0_dp)*WL(5)/WL(1) + 0.5_dp*(WL(2)**2 + WL(3)**2 + WL(4)**2)
        HR = gam/(gam - 1.0_dp)*WR(5)/WR(1) + 0.5_dp*(WR(2)**2 + WR(3)**2 + WR(4)**2)

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

        ! 고유값 절댓값 (소산 크기). 엔트로피 보정 (Harten): |lam| < delta 이면 (lam^2 + delta^2)/(2 delta) 로
        !   0 근처에서 소산이 사라지지 않게 함. 음속점 (un -+ c = 0) 의 팽창 충격파 방지, all 이면 정체점 (un = 0) 도
        al1 = abs(l1);  al2 = abs(l2);  al5 = abs(l5)
        if (entropy_fix > 0.0_dp) then
            delta = entropy_fix*c_tild
            if (al1 < delta) al1 = (l1**2 + delta**2)/(2.0_dp*delta)
            if (al5 < delta) al5 = (l5**2 + delta**2)/(2.0_dp*delta)
            if (entropy_fix_all .and. al2 < delta) al2 = (l2**2 + delta**2)/(2.0_dp*delta)
        end if

        r1(1) = 1; r1(2:4) = vel_tild - c_tild*n; r1(5) = H_tild - c_tild*un_tild
        r2(1) = 1; r2(2:4) = vel_tild;            r2(5) = 0.5_dp*q2_tild
        r5(1) = 1; r5(2:4) = vel_tild + c_tild*n; r5(5) = H_tild + c_tild*un_tild
        rs(1) = 0; rs(2:4) = rho_tild*dut;        rs(5) = rho_tild*dot_product(vel_tild, dut)

        call flux_n(WL, n, FL)
        call flux_n(WR, n, FR)

        F = 0.5_dp*(FL+FR) - 0.5_dp*(al1*a1*r1 + al2*(a2*r2 + rs) + al5*a5*r5)

    end subroutine Roe

    !> 리만 솔버 고르기 (flow.in 의 flux = rusanov | roe), 면 양쪽 원시 변수
    subroutine riemann(WL, WR, n, F)
        real(dp), intent(in)  :: WL(nvar), WR(nvar), n(3)
        real(dp), intent(out) :: F(nvar)
        select case (flux_id)
        case (FLUX_ROE)
            call Roe(WL, WR, n, F)
        case default
            call rusanov(WL, WR, n, F)
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

    !> 블록 bk 의 잔차 R = -(면 플럭스 합)/V  (고스트 셀은 bc() 로 미리 채워 둠)
    !>   1) 원시 변수 W 를 셀마다 한 번만 계산 (내부 + 면 고스트. 모서리 고스트는 안 씀)
    !>   2) i, j, k 방향마다 면 플럭스를 같은 버퍼에 계산 -> R 에 더함 (버퍼 하나라 메모리 1/3)
    subroutine residual(bk, R)
        integer,  intent(in)  :: bk
        real(dp), intent(out) :: R(:,:,:,:)

        integer :: fi, fj, fk, ic, jc, kc, nf, bb
        real(dp) :: F(nvar), WL(nvar), WR(nvar)
        real(dp), pointer, contiguous :: fl(:,:,:,:), Q(:,:,:,:), W(:,:,:,:)

        call use_block(bk)
        Q => blk(bk)%Q
        if (.not. allocated(blk(bk)%W)) allocate(blk(bk)%W, mold=blk(bk)%Q)
        if (.not. allocated(fbuf)) then
            nf = 0
            do bb = 1,nblk
                associate (m => blk(bb))
                    nf = max(nf, (m%nx+1)*m%ny*m%nz, m%nx*(m%ny+1)*m%nz, m%nx*m%ny*(m%nz+1))
                end associate
            end do
            allocate(fbuf(nvar*nf))
        end if

        W => blk(bk)%W

        ! OpenMP: 바깥 두 루프를 묶어 스레드에 나눔 (collapse(2): 한 방향 셀이 적어도 고르게 나뉨)
        !   F, WL, WR 은 스레드마다 따로 (private). 버퍼, R 은 스레드마다 다른 위치에 씀
        !$omp parallel private(fi, fj, fk, ic, jc, kc, F, WL, WR)

        ! ---- 원시 변수: 내부 셀, 그리고 i / j / k 방향 고스트 (각 방향 면 루프가 쓰는 것만) ----
        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1-ng,nx+ng
                    call c2p(Q(:,ic,jc,kc), W(:,ic,jc,kc))
                end do
            end do
        end do
        !$omp end do nowait
        !$omp do collapse(2)
        do kc = 1,nz
            do ic = 1,nx
                do jc = 1-ng,0
                    call c2p(Q(:,ic,jc,kc), W(:,ic,jc,kc))
                    call c2p(Q(:,ic,ny+ng+jc,kc), W(:,ic,ny+ng+jc,kc))
                end do
            end do
        end do
        !$omp end do nowait
        !$omp do collapse(2)
        do jc = 1,ny
            do ic = 1,nx
                do kc = 1-ng,0
                    call c2p(Q(:,ic,jc,kc), W(:,ic,jc,kc))
                    call c2p(Q(:,ic,jc,nz+ng+kc), W(:,ic,jc,nz+ng+kc))
                end do
            end do
        end do
        !$omp end do

        ! ---- i-면: fi = 면 번호 (셀 fi 와 fi+1 사이) ----
        !   MUSCL 인자 순서: 가운데(fi), 왼쪽(fi-1), 오른쪽(fi+1), 그 오른쪽(fi+2)
        !$omp single
        fl(1:nvar,0:nx,1:ny,1:nz) => fbuf(1:nvar*(nx+1)*ny*nz)
        !$omp end single
        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do fi = 0,nx
                    call MUSCL(W(:,fi,jc,kc), W(:,fi-1,jc,kc), W(:,fi+1,jc,kc), W(:,fi+2,jc,kc), si(:,fi,jc,kc), WL, WR)
                    call riemann(WL, WR, si(:,fi,jc,kc), F)
                    fl(:,fi,jc,kc) = F*si_a(fi,jc,kc)
                end do
            end do
        end do
        !$omp end do
        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    R(:,ic,jc,kc) = fl(:,ic,jc,kc) - fl(:,ic-1,jc,kc)
                end do
            end do
        end do
        !$omp end do

        ! ---- j-면 ----
        !$omp single
        fl(1:nvar,1:nx,0:ny,1:nz) => fbuf(1:nvar*nx*(ny+1)*nz)
        !$omp end single
        !$omp do collapse(2)
        do kc = 1,nz
            do fj = 0,ny
                do ic = 1,nx
                    call MUSCL(W(:,ic,fj,kc), W(:,ic,fj-1,kc), W(:,ic,fj+1,kc), W(:,ic,fj+2,kc), sj(:,ic,fj,kc), WL, WR)
                    call riemann(WL, WR, sj(:,ic,fj,kc), F)
                    fl(:,ic,fj,kc) = F*sj_a(ic,fj,kc)
                end do
            end do
        end do
        !$omp end do
        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    R(:,ic,jc,kc) = R(:,ic,jc,kc) + (fl(:,ic,jc,kc) - fl(:,ic,jc-1,kc))
                end do
            end do
        end do
        !$omp end do

        ! ---- k-면, 그리고 R = -(합)/V ----
        !$omp single
        fl(1:nvar,1:nx,1:ny,0:nz) => fbuf(1:nvar*nx*ny*(nz+1))
        !$omp end single
        !$omp do collapse(2)
        do fk = 0,nz
            do jc = 1,ny
                do ic = 1,nx
                    call MUSCL(W(:,ic,jc,fk), W(:,ic,jc,fk-1), W(:,ic,jc,fk+1), W(:,ic,jc,fk+2), sk(:,ic,jc,fk), WL, WR)
                    call riemann(WL, WR, sk(:,ic,jc,fk), F)
                    fl(:,ic,jc,fk) = F*sk_a(ic,jc,fk)
                end do
            end do
        end do
        !$omp end do
        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    R(:,ic,jc,kc) = -(R(:,ic,jc,kc) + (fl(:,ic,jc,kc) - fl(:,ic,jc,kc-1))) / vol(ic,jc,kc)
                end do
            end do
        end do
        !$omp end do

        !$omp end parallel

    end subroutine residual

    !> 블록 bk 의 CFL 조건 dt
    !>   셀마다 dt = cfl * V / (lam_i + lam_j + lam_k)
    !>   lam_i = |u.Si| + c|Si|,  Si = 셀 양쪽 i-면 벡터(법선 x 넓이)의 평균  (j, k 도 같음)
    !>   blk(bk)%dtl : 셀마다의 dt (정상 문제의 국소 시간 간격),  dt : 그 최솟값
    subroutine compute_dt(bk, dt)
        integer,  intent(in)  :: bk
        real(dp), intent(out) :: dt

        integer  :: ic, jc, kc
        real(dp) :: W(nvar), c, s_i(3), s_j(3), s_k(3), lam
        real(dp), pointer, contiguous :: Q(:,:,:,:), dtl(:,:,:)

        call use_block(bk)
        if (.not. allocated(blk(bk)%dtl)) allocate(blk(bk)%dtl(nx,ny,nz))
        dt = huge(1.0_dp)

        Q => blk(bk)%Q;  dtl => blk(bk)%dtl
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
