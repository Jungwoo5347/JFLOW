module solver
    use euler
    use grid
    use case_setup, only: cfl,sigma,kappa,b, flux_id, FLUX_RUSANOV, FLUX_ROE, &
                          viscous, prandtl, sutherland, S_suth, visc_coef                            ! CFL 수 (입력 파일 &control 에서)

    implicit none

    real(dp), allocatable :: fluxi(:,:,:), fluxj(:,:,:)    ! 면 플럭스 (작업용, residual 전용)
    real(dp), allocatable :: fluxvi(:,:,:), fluxvj(:,:,:)  ! 면 점성 플럭스 F_v.n (단위 길이당, M_inf/Re 포함, viscous 일 때만)

    ! 점성항 작업 배열 (셀, 고스트 첫 층 포함 0:nx+1, 0:ny+1): vis_stress 가 채우고 cal_tau_four 가 읽음
    real(dp), allocatable :: Wc(:,:,:), Tpr(:,:)                 ! 셀 원시변수 (rho, u, v, p), 온도 T = gam p / rho
    real(dp), allocatable :: dudx(:,:), dudy(:,:), dvdx(:,:), dvdy(:,:), dTdx(:,:), dTdy(:,:)   ! 셀 기울기

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

    subroutine MUSCL(Q,Qm1,Q1,Q2,fnx,fny,QL,QR)
        real(dp), intent(in) :: Q(:),Qm1(:),Q1(:),Q2(:),fnx,fny
        real(dp), intent(out) :: QL(:), QR(:)

        real(dp) :: W(nvar), Wm1(nvar), W1(nvar), W2(nvar),WL(nvar),WR(nvar)
        real(dp) :: deltaWmh(nvar),deltaWph(nvar),deltaWpph(nvar)
        real(dp) :: am(nvar), a0(nvar), ap(nvar), sL(nvar), sR(nvar), dL(nvar), dR(nvar)
        real(dp) :: rhof, cf, c0, c1

        call c2p(Q,W)
        call c2p(Qm1,Wm1)
        call c2p(Q1,W1)
        call c2p(Q2,W2)

        !delta 계산
        call delta(W,WM1,deltaWmh)
        call delta(W1,W,deltaWph)
        call delta(W2,W1,deltaWpph)

        call sos(W,c0)
        call sos(W1,c1)
        rhof = 0.5_dp*(W(1)+W1(1))
        cf   = 0.5_dp*(c0+c1)

        call to_char(deltaWmh, rhof,cf,fnx,fny,am)
        call to_char(deltaWph, rhof,cf,fnx,fny,a0)
        call to_char(deltaWpph,rhof,cf,fnx,fny,ap)

        sL = (1-kappa)*minmod(am,b*a0)+(1+kappa)*minmod(a0,b*am)
        sR = (1-kappa)*minmod(ap,b*a0)+(1+kappa)*minmod(a0,b*ap)

        call from_char(sL,rhof,cf,fnx,fny,dL)
        call from_char(sR,rhof,cf,fnx,fny,dR)

        WL = W +dL*(sigma/4.0_dp)
        WR = W1-dR*(sigma/4.0_dp)

        if (WL(1) <= 0.0_dp .or. WL(4) <= 0.0_dp .or. WR(1) <= 0.0_dp .or. WR(4) <= 0.0_dp) then
            WL = W
            WR = W1
        end if

        call p2c(WL,QL)
        call p2c(WR,QR)


    end subroutine MUSCL

    subroutine to_char(dW,rho,c,nx,ny,a)
        real(dp), intent(in)  :: dW(nvar), rho, c, nx, ny
        real(dp), intent(out) :: a(nvar)
        real(dp) :: dun, dut

        dun =  dW(2)*nx + dW(3)*ny
        dut = -dW(2)*ny + dW(3)*nx

        a(1) = (dW(4) - rho*c*dun)/(2.0_dp*c**2)
        a(2) =  dW(1) - dW(4)/c**2
        a(3) =  rho*dut
        a(4) = (dW(4) + rho*c*dun)/(2.0_dp*c**2)

    end subroutine to_char

    subroutine from_char(a,rho,c,nx,ny,dW)
        real(dp), intent(in)  :: a(nvar), rho, c, nx, ny
        real(dp), intent(out) :: dW(nvar)
        real(dp) :: dun, dut

        dun = c/rho*(a(4) - a(1))
        dut = a(3)/rho

        dW(1) = a(1) + a(2) + a(4)
        dW(2) = dun*nx - dut*ny
        dW(3) = dun*ny + dut*nx
        dW(4) = c**2*(a(1) + a(4))

    end subroutine from_char

    subroutine Roe(QL,QR,fnx,fny,F)
        real(dp), intent(in) :: QL(:), QR(:), fnx, fny
        real(dp), intent(out) :: F(:)

        real(dp) :: WL(nvar), WR(nvar), FL(nvar), FR(nvar)
        real(dp) :: rho_tild, u_tild, v_tild, H_tild, c_tild, Rw, HL, HR, unL,unR, utL,utR,un_tild,ut_tild
        real(dp) :: a1, a2, a3, a4, l1, l2, l3, l4, r1(nvar),r2(nvar),r3(nvar),r4(nvar)

        call c2p(QL,WL)
        call c2p(QR,WR)

        HL = (QL(4) + WL(4))/WL(1)
        HR = (QR(4) + WR(4))/WR(1)

        ! rho, u, v, p 순서
        Rw = sqrt(WR(1)/WL(1))
        rho_tild = sqrt(WL(1)*WR(1))
        u_tild = (WL(2)+Rw*WR(2))/(1+Rw)
        v_tild = (WL(3)+Rw*WR(3))/(1+Rw)
        H_tild = (HL + Rw*HR)/(1+Rw)
        c_tild = sqrt((gam-1)*(H_tild-0.5_dp*(u_tild**2+v_tild**2)))

        !법선 / 접선 속도
        unL = WL(2)*fnx + WL(3)*fny        
        utL = -WL(2)*fny + WL(3)*fnx
        unR = WR(2)*fnx + WR(3)*fny
        utR = -WR(2)*fny + WR(3)*fnx

        un_tild = u_tild*fnx+v_tild*fny
        ut_tild = -u_tild*fny+v_tild*fnx

        a1 = ((WR(4)-WL(4))-rho_tild*c_tild*(unR-unL))/(2*c_tild**2)
        a2 = (WR(1)-WL(1))-(WR(4)-WL(4))/c_tild**2
        a3 = rho_tild*(utR-utL)
        a4 = ((WR(4)-WL(4))+rho_tild*c_tild*(unR-unL))/(2*c_tild**2)
        
        l1 = un_tild-c_tild
        l2 = un_tild
        l3 = un_tild
        l4 = un_tild+c_tild

        r1(1) = 1; r1(2) = u_tild-c_tild*fnx; r1(3) = v_tild-c_tild*fny; r1(4) = H_tild-c_tild*un_tild
        r2(1) = 1; r2(2) = u_tild; r2(3) = v_tild; r2(4) = 0.5_dp*(u_tild**2+v_tild**2)
        r3(1) = 0; r3(2) = -fny; r3(3) = fnx; r3(4) = ut_tild
        r4(1) = 1; r4(2) = u_tild+c_tild*fnx; r4(3) = v_tild+c_tild*fny; r4(4) = H_tild+c_tild*un_tild

        FL(1) = WL(1)*unL
        FL(2) = WL(1)*WL(2)*unL + WL(4)*fnx
        FL(3) = WL(1)*WL(3)*unL + WL(4)*fny
        FL(4) = (QL(4)+WL(4))*unL

        FR(1) = WR(1)*unR
        FR(2) = WR(1)*WR(2)*unR + WR(4)*fnx
        FR(3) = WR(1)*WR(3)*unR + WR(4)*fny
        FR(4) = (QR(4)+WR(4))*unR

        F = 0.5_dp*(FL+FR) - 0.5_dp*(abs(l1)*a1*r1 + abs(l2)*a2*r2 + abs(l3)*a3*r3 + abs(l4)*a4*r4)


    end subroutine Roe

    !> 면 하나의 점성 플럭스 F_v.n (단위 길이당): 셀 L = (iL,jL) 과 R = (iR,jR) 사이 면
    !>   tau, q 는 cal_tau_four 결과, (fnx, fny) 는 면 단위 법선 (L -> R)
    !>   F_v.n = (M_inf/Re) * ( 0,  txx nx + txy ny,  txy nx + tyy ny,  u (txx nx + txy ny) + v (txy nx + tyy ny) - q.n )
    subroutine vis_flux(iL,jL, iR,jR, txxf,txyf,tyyf,qxf,qyf, fnx,fny, Fv)
        integer,  intent(in)  :: iL, jL, iR, jR
        real(dp), intent(in)  :: txxf, txyf, tyyf, qxf, qyf, fnx, fny
        real(dp), intent(out) :: Fv(nvar)

        real(dp) :: u, v

        ! 면 속도 = 양옆 셀 평균 (점착 벽이면 고스트가 -u 라서 0)
        u = 0.5_dp*(Wc(2,iL,jL) + Wc(2,iR,jR))
        v = 0.5_dp*(Wc(3,iL,jL) + Wc(3,iR,jR))

        Fv(1) = 0.0_dp
        Fv(2) = txxf*fnx + txyf*fny
        Fv(3) = txyf*fnx + tyyf*fny
        Fv(4) = u*Fv(2) + v*Fv(3) - (qxf*fnx + qyf*fny)
        Fv    = visc_coef*Fv

    end subroutine vis_flux


    !> 모든 면의 점성 플럭스 F_v.n (단위 길이당)
    !>   셀 원시변수, 온도 -> 셀 기울기 (Green-Gauss) -> 면마다 tau, q (cal_tau_four) -> F_v.n (vis_flux)
    !>   Q 의 고스트 첫 층은 bc() 로 채워져 있어야 함
    subroutine vis_stress(Q, fvi, fvj)
        real(dp), intent(in)  :: Q(nvar,1-ng:nx+ng,1-ng:ny+ng)
        real(dp), intent(out) :: fvi(nvar,0:nx,1:ny)       ! i-면
        real(dp), intent(out) :: fvj(nvar,1:nx,0:ny)       ! j-면

        real(dp) :: txxf, txyf, tyyf, qxf, qyf
        integer :: i,j

        if (.not. allocated(Wc)) then
            allocate(Wc(nvar,0:nx+1,0:ny+1), Tpr(0:nx+1,0:ny+1))
            allocate(dudx(0:nx+1,0:ny+1), dudy(0:nx+1,0:ny+1), dvdx(0:nx+1,0:ny+1), dvdy(0:nx+1,0:ny+1))
            allocate(dTdx(0:nx+1,0:ny+1), dTdy(0:nx+1,0:ny+1))
            Wc = 0.0_dp;  Tpr = 0.0_dp            ! 모서리 고스트 (안 씀) 에 쓰레기값이 남지 않게
        end if

        ! 셀 원시변수, 온도 (고스트 첫 층 포함, 모서리 제외: 모서리 Q 는 0 이라 c2p 에서 0 으로 나눔)
        do j = 0,ny+1
            do i = 0,nx+1
                if ((i == 0 .or. i == nx+1) .and. (j == 0 .or. j == ny+1)) cycle
                call c2p(Q(:,i,j), Wc(:,i,j))
                Tpr(i,j) = gam*Wc(4,i,j)/Wc(1,i,j)
            end do
        end do

        ! 셀 기울기 (Green-Gauss)
        call grad(Wc(2,:,:), dudx, dudy)
        call grad(Wc(3,:,:), dvdx, dvdy)
        call grad(Tpr,       dTdx, dTdy)

        ! i-면 (i,j): 셀 (i,j) | (i+1,j)
        !$omp parallel do private(i, txxf, txyf, tyyf, qxf, qyf)
        do j = 1,ny
            do i = 0,nx
                call cal_tau_four(i,j, i+1,j, txxf, txyf, tyyf, qxf, qyf)
                call vis_flux(i,j, i+1,j, txxf, txyf, tyyf, qxf, qyf, si_x(i,j), si_y(i,j), fvi(:,i,j))
            end do
        end do
        !$omp end parallel do

        ! j-면 (i,j): 셀 (i,j) | (i,j+1)
        !$omp parallel do private(i, txxf, txyf, tyyf, qxf, qyf)
        do j = 0,ny
            do i = 1,nx
                call cal_tau_four(i,j, i,j+1, txxf, txyf, tyyf, qxf, qyf)
                call vis_flux(i,j, i,j+1, txxf, txyf, tyyf, qxf, qyf, sj_x(i,j), sj_y(i,j), fvj(:,i,j))
            end do
        end do
        !$omp end parallel do

    end subroutine vis_stress

    !> 무차원 점성계수 mu* = mu / mu_inf (자유류 T = 1 에서 1): Sutherland, 또는 상수 1 (viscosity_law)
    elemental real(dp) function mu_star(T)
        real(dp), intent(in) :: T
        if (sutherland) then
            mu_star = T**1.5_dp * (1.0_dp + S_suth)/(T + S_suth)
        else
            mu_star = 1.0_dp
        end if
    end function mu_star

    !> 면 하나의 tau, q: 셀 L = (iL,jL) 과 R = (iR,jR) 사이 면 (i-면, j-면 공통)
    !>   i-면 (i,j): L = (i,j), R = (i+1,j)  /  j-면 (i,j): L = (i,j), R = (i,j+1)
    !>   셀 값, 셀 기울기는 vis_stress 가 채운 모듈 배열 (Wc, Tpr, dudx ...), 셀 중심은 grid 의 xcg, ycg
    subroutine cal_tau_four(iL,jL, iR,jR, txxf,txyf,tyyf,qxf,qyf)
        integer,  intent(in)  :: iL, jL, iR, jR
        real(dp), intent(out) :: txxf,txyf,tyyf,qxf,qyf

        real(dp) :: TL, TR, Tf, muf
        real(dp) :: dist, tx, ty, corr
        real(dp) :: dudx_f, dudy_f, dvdx_f, dvdy_f, dTdx_f, dTdy_f    ! 면에서의 기울기 (평균 + 보정)

        TL = Tpr(iL,jL)
        TR = Tpr(iR,jR)

        ! 면 점성계수 mu* (자유류에서 1): 면 온도로 Sutherland, 또는 상수
        Tf = 0.5_dp*(TL + TR)
        muf = mu_star(Tf)

        ! 두 셀 중심을 잇는 단위 벡터 (tx, ty) 와 거리
        dist = sqrt((xcg(iR,jR) - xcg(iL,jL))**2 + (ycg(iR,jR) - ycg(iL,jL))**2)
        tx = (xcg(iR,jR) - xcg(iL,jL))/dist
        ty = (ycg(iR,jR) - ycg(iL,jL))/dist

        ! 면 기울기 = 양옆 평균 + 보정 (t 방향 성분만 두 셀 값의 직접 차분으로 바꿈)
        dudx_f = 0.5_dp*(dudx(iL,jL) + dudx(iR,jR))
        dudy_f = 0.5_dp*(dudy(iL,jL) + dudy(iR,jR))
        corr   = (Wc(2,iR,jR) - Wc(2,iL,jL))/dist - (dudx_f*tx + dudy_f*ty)
        dudx_f = dudx_f + corr*tx
        dudy_f = dudy_f + corr*ty

        dvdx_f = 0.5_dp*(dvdx(iL,jL) + dvdx(iR,jR))
        dvdy_f = 0.5_dp*(dvdy(iL,jL) + dvdy(iR,jR))
        corr   = (Wc(3,iR,jR) - Wc(3,iL,jL))/dist - (dvdx_f*tx + dvdy_f*ty)
        dvdx_f = dvdx_f + corr*tx
        dvdy_f = dvdy_f + corr*ty

        dTdx_f = 0.5_dp*(dTdx(iL,jL) + dTdx(iR,jR))
        dTdy_f = 0.5_dp*(dTdy(iL,jL) + dTdy(iR,jR))
        corr   = (TR - TL)/dist - (dTdx_f*tx + dTdy_f*ty)
        dTdx_f = dTdx_f + corr*tx
        dTdy_f = dTdy_f + corr*ty

        ! 응력 (Stokes 가정), 열유속 q = -mu*/((gam-1) Pr) grad T
        !   M_inf/Re (visc_coef) 는 여기서 곱하지 않음: 점성 플럭스를 만들 때 한 번만 곱할 것
        txxf = muf*(4.0_dp/3.0_dp*dudx_f - 2.0_dp/3.0_dp*dvdy_f)
        tyyf = muf*(4.0_dp/3.0_dp*dvdy_f - 2.0_dp/3.0_dp*dudx_f)
        txyf = muf*(dudy_f + dvdx_f)
        qxf  = -muf/((gam - 1.0_dp)*prandtl)*dTdx_f
        qyf  = -muf/((gam - 1.0_dp)*prandtl)*dTdy_f

    end subroutine cal_tau_four

    subroutine grad(phi,phix,phiy)
        real(dp), intent(in) :: phi(0:,0:)              ! (0:nx+1, 0:ny+1), 띄엄띄엄한 단면 (Wc(2,:,:)) 도 복사 없이 받음
        real(dp), intent(out) :: phix(0:nx+1,0:ny+1), phiy(0:nx+1,0:ny+1)
        
        real(dp) :: pE, pW, pS, pN
        integer :: i,j


        !$omp parallel do private(i, pE, pW, pN, pS)
        do j = 1,ny
            do i = 1,nx
                pE = 0.5_dp*(phi(i,j) + phi(i+1,j))
                pW = 0.5_dp*(phi(i-1,j) + phi(i,j))
                pN = 0.5_dp*(phi(i,j) + phi(i,j+1))
                pS = 0.5_dp*(phi(i,j-1) + phi(i,j))

                phix(i,j) = ( pE*si_x(i,j)*si_len(i,j) - pW*si_x(i-1,j)*si_len(i-1,j) &
                            + pN*sj_x(i,j)*sj_len(i,j) - pS*sj_x(i,j-1)*sj_len(i,j-1) ) / area(i,j)
                phiy(i,j) = ( pE*si_y(i,j)*si_len(i,j) - pW*si_y(i-1,j)*si_len(i-1,j) &
                            + pN*sj_y(i,j)*sj_len(i,j) - pS*sj_y(i,j-1)*sj_len(i,j-1) ) / area(i,j)
            end do
        end do
        !$omp end parallel do


        phix(0,:) = 0.0_dp;  phix(nx+1,:) = 0.0_dp;  phix(:,0) = 0.0_dp;  phix(:,ny+1) = 0.0_dp
        phiy(0,:) = 0.0_dp;  phiy(nx+1,:) = 0.0_dp;  phiy(:,0) = 0.0_dp;  phiy(:,ny+1) = 0.0_dp
        phix(0,   1:ny) = phix(1, 1:ny);  phix(nx+1,1:ny) = phix(nx,1:ny)
        phiy(0,   1:ny) = phiy(1, 1:ny);  phiy(nx+1,1:ny) = phiy(nx,1:ny)
        phix(1:nx,0   ) = phix(1:nx,1 );  phix(1:nx,ny+1) = phix(1:nx,ny)
        phiy(1:nx,0   ) = phiy(1:nx,1 );  phiy(1:nx,ny+1) = phiy(1:nx,ny)

    end subroutine grad


    !> 리만 솔버 고르기 (flow.in 의 flux = rusanov | roe)
    !>   새 리만 솔버를 만들면 case_setup 에 FLUX_ 번호와 이름을 추가하고 여기에 case 한 줄 추가
    subroutine riemann(QL, QR, fnx, fny, F)
        real(dp), intent(in)  :: QL(nvar), QR(nvar), fnx, fny
        real(dp), intent(out) :: F(nvar)
        select case (flux_id)
        case (FLUX_ROE)
            call Roe(QL, QR, fnx, fny, F)
        case default
            call rusanov(QL, QR, fnx, fny, F)
        end select
    end subroutine riemann

    !> minmod 리미터: 부호가 같으면 절댓값이 작은 쪽, 다르면 (극값) 0
    elemental function minmod(a,b) result(m)
        real(dp),intent(in) :: a,b
        real(dp) :: m

        if (a*b <= 0.0_dp) then
            m = 0.0_dp                          
        else
            m = sign(min(abs(a),abs(b)), a)     !
        end if
    end function  minmod

    subroutine delta(uL,uLL,du)
        real(dp),intent(in) :: uL(:),uLL(:)
        real(dp), intent(out) :: du(:)

        du = uL-uLL

    end subroutine

    subroutine residual(Q,R)
        real(dp), intent(inout) :: Q(:,1-ng:,1-ng:)
        real(dp), intent(out) :: R(nvar,0:nx+1,0:ny+1)

        integer :: fi, fj, ic, jc
        real(dp) :: F(nvar), QL(nvar), QR(nvar)

        if (.not. allocated(fluxi)) then
            allocate(fluxi(nvar, 0:nx, 1:ny), fluxj(nvar, 1:nx, 0:ny))
            if (viscous) allocate(fluxvi(nvar, 0:nx, 1:ny), fluxvj(nvar, 1:nx, 0:ny))
        end if

        R = 0.0_dp

        ! 점성 플럭스 (모든 면): 면 플럭스 = (비점성 - 점성) * 면 길이
        if (viscous) call vis_stress(Q, fluxvi, fluxvj)

        !$omp parallel private(fi, fj, ic, jc, F, QL, QR)


        !$omp do
        do jc = 1,ny
            do fi = 0,nx
                call MUSCL(Q(:,fi,jc), Q(:,fi-1,jc), Q(:,fi+1,jc), Q(:,fi+2,jc), si_x(fi,jc), si_y(fi,jc), QL, QR)
                call riemann(QL, QR, si_x(fi,jc), si_y(fi,jc), F)
                if (viscous) F = F - fluxvi(:,fi,jc)
                fluxi(:,fi,jc) = F*si_len(fi,jc)
            end do
        end do
        !$omp end do nowait

        !$omp do
        do fj = 0,ny
            do ic = 1,nx
                call MUSCL(Q(:,ic,fj), Q(:,ic,fj-1), Q(:,ic,fj+1), Q(:,ic,fj+2), sj_x(ic,fj), sj_y(ic,fj), QL, QR)
                call riemann(QL, QR, sj_x(ic,fj), sj_y(ic,fj), F)
                if (viscous) F = F - fluxvj(:,ic,fj)
                fluxj(:,ic,fj) = F*sj_len(ic,fj)
            end do
        end do
        !$omp end do

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
    !>   dtl(ic,jc) : 셀마다의 dt (정상 문제의 국소 시간 간격),  dt : 그 최솟값 (비정상 문제의 전체 dt)
    subroutine compute_dt(Q, dtl, dt)
        real(dp), intent(in)  :: Q(:,1-ng:,1-ng:)
        real(dp), intent(out) :: dtl(:,:)
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

                dtl(ic,jc) = cfl*area(ic,jc)/(lam_i + lam_j)
                dt = min(dt, dtl(ic,jc))
            end do
        end do
        !$omp end parallel do

    end subroutine compute_dt

end module solver