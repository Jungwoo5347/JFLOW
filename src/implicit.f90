!> LU-SGS 암시적 반복 (Yoon & Jameson), 2D  (flow.in 의 time_scheme = implicit, implicit_solver = lusgs)
!>
!>   풀 식 (셀 I, 면적 A):  A dQ/dtau = A R*(Q),   R* = R (정상)  또는  R - (a0 Q + a1 Q^n + a2 Q^(n-1))/dt (dual time)
!>   역방향 오일러 + 면 플럭스 선형화 (스펙트럼 반경 분리):
!>     F_f(Q+dQ).S - F_f(Q).S  ~  1/2 (dF_I + dF_J).S + 1/2 r_f (dQ_I - dQ_J)
!>     r_f = w (|u.n| + c) |S|  +  2 r_v,   r_v = (M/Re) mu*/rho max(4/3, gam/Pr) |S| / d   (점성일 때, d = 두 셀 중심 거리)
!>   셀 I 에 대해 (닫힌 셀이라 sum_f dF_I.S = 0):
!>     D_I dQ_I + sum_J 1/2 (dF_J.S_out - r_f dQ_J) = A R*_I,     D_I = A/dtau_I + A diag_add + 1/2 sum_f r_f   (스칼라)
!>     dF_J = F(Q_J + dQ_J) - F(Q_J)   (비점성 플럭스 차이 그대로, 야코비안 행렬 없음)
!>     diag_add = a0/dt (dual time 의 시간 미분 항을 암시적으로), 정상이면 0
!>   L (i-1, j-1 이웃) / U (i+1, j+1 이웃) 로 나눠 두 번 스윕:
!>     앞으로 : dQ*_I = [ A R*_I - sum_{J in L} 1/2 (dF*_J.S - r_f dQ*_J) ] / D_I
!>     뒤로   : dQ_I  = dQ*_I - [ sum_{J in U} 1/2 (dF_J.S - r_f dQ_J) ] / D_I
!>   우변 R* 은 explicit 과 같은 2차 잔차 -> 수렴한 해는 같음. 좌변은 1차 근사라 수렴 속도만 바꿈
!>
!>   병렬: 격자를 tile x tile 셀 타일로 나누고 타일마다 독립적으로 앞/뒤 스윕
!>     (타일 경계 건너편의 dQ 는 0, 격자 경계와 같은 처리 -> 스레드끼리 기다릴 일이 없음)
!>     타일 크기가 고정이라 스레드 수와 상관없이 결과가 같음
!>   고스트 셀의 dQ 는 0 (경계조건은 매 반복 bc() 로 명시적으로)
!>
!>   다른 암시적 방법 (예: DADI) 을 넣을 때: 같은 인자 (Q, Rs, dtl, diag_add) 의 서브루틴을 만들고
!>   main 의 implicit_step 에서 implicit_solver 에 따라 고르면 됨
module implicit
    use euler
    use grid
    use solver, only: mu_star
    use case_setup, only: lusgs_omega, viscous, prandtl, visc_coef
    implicit none
    private
    public :: lusgs

    integer, parameter :: tile = 32        ! 타일 한 변의 셀 수

    integer :: ntile = 0
    integer,  allocatable :: lo(:,:), hi(:,:)              ! 타일 t 의 셀 범위 (i, j)
    real(dp), allocatable :: Wp(:,:,:)                     ! 원시 변수 (고스트 첫 층 포함)
    real(dp), allocatable :: ri(:,:), rj(:,:), dg(:,:)     ! 면마다 r_f (i-면 (0:nx,ny), j-면 (nx,0:ny)), 셀마다 D_I
    real(dp), allocatable :: dQ(:,:,:)

contains

    !> 작업 배열과 타일 목록 (처음 한 번)
    subroutine setup_lu()
        integer :: nti, ntj, ti, tj, t
        if (ntile > 0) return
        allocate(Wp(nvar,0:nx+1,0:ny+1), ri(0:nx,ny), rj(nx,0:ny), dg(nx,ny), dQ(nvar,nx,ny))
        Wp = 0.0_dp
        nti = (nx + tile - 1)/tile
        ntj = (ny + tile - 1)/tile
        ntile = nti*ntj
        allocate(lo(2,ntile), hi(2,ntile))
        t = 0
        do tj = 1,ntj
            do ti = 1,nti
                t = t + 1
                lo(:,t) = [(ti-1)*tile + 1, (tj-1)*tile + 1]
                hi(:,t) = [min(ti*tile, nx), min(tj*tile, ny)]
            end do
        end do
    end subroutine setup_lu

    !> LU-SGS 한 반복: dQ 를 구해 Q (내부 셀) 에 더함
    !>   Q        : 보존 변수 (고스트 첫 층은 bc() 로 채워져 있어야 함)
    !>   Rs       : 우변 dQ/dtau (정상: residual 의 R,  dual time: R*)
    !>   dtl      : 셀마다 가상 시간 간격 dtau
    !>   diag_add : 대각항에 더할 값 / 면적 (dual time: a0/dt,  정상: 0)
    subroutine lusgs(Q, Rs, dtl, diag_add)
        real(dp), intent(inout) :: Q(:,1-ng:,1-ng:)
        real(dp), intent(in)    :: Rs(:,0:,0:)
        real(dp), intent(in)    :: dtl(:,:)
        real(dp), intent(in)    :: diag_add

        integer  :: t, ic, jc, fi, fj
        real(dp) :: rhs(nvar)

        call setup_lu()

        !$omp parallel private(t, ic, jc, fi, fj, rhs)

        ! ---- 원시 변수 (고스트 첫 층 포함, 모서리 제외) ----
        !$omp do
        do jc = 0,ny+1
            do ic = 0,nx+1
                if ((ic == 0 .or. ic == nx+1) .and. (jc == 0 .or. jc == ny+1)) cycle
                call c2p(Q(:,ic,jc), Wp(:,ic,jc))
            end do
        end do
        !$omp end do

        ! ---- 면마다 r_f ----
        !$omp do
        do jc = 1,ny
            do fi = 0,nx
                ri(fi,jc) = face_rad(fi,jc, fi+1,jc, si_x(fi,jc), si_y(fi,jc), si_len(fi,jc))
            end do
        end do
        !$omp end do nowait
        !$omp do
        do fj = 0,ny
            do ic = 1,nx
                rj(ic,fj) = face_rad(ic,fj, ic,fj+1, sj_x(ic,fj), sj_y(ic,fj), sj_len(ic,fj))
            end do
        end do
        !$omp end do

        ! ---- 대각항 D_I = A/dtau + A diag_add + 1/2 sum_4면 r_f ----
        !$omp do
        do jc = 1,ny
            do ic = 1,nx
                dg(ic,jc) = area(ic,jc)*(1.0_dp/dtl(ic,jc) + diag_add) &
                          + 0.5_dp*(ri(ic-1,jc) + ri(ic,jc) + rj(ic,jc-1) + rj(ic,jc))
            end do
        end do
        !$omp end do

        ! ---- 타일마다 앞/뒤 스윕 (타일끼리 독립) ----
        !$omp do schedule(dynamic,1)
        do t = 1,ntile
            ! 앞으로 (i, j 오름차순): 아래 이웃 (타일 안) 의 dQ* 사용
            do jc = lo(2,t),hi(2,t)
                do ic = lo(1,t),hi(1,t)
                    rhs = area(ic,jc)*Rs(:,ic,jc)
                    if (ic > lo(1,t)) call nb(ic-1, jc, -si_x(ic-1,jc), -si_y(ic-1,jc), si_len(ic-1,jc), ri(ic-1,jc), rhs)
                    if (jc > lo(2,t)) call nb(ic, jc-1, -sj_x(ic,jc-1), -sj_y(ic,jc-1), sj_len(ic,jc-1), rj(ic,jc-1), rhs)
                    dQ(:,ic,jc) = rhs/dg(ic,jc)
                end do
            end do
            ! 뒤로 (내림차순): 위 이웃 (타일 안) 의 최종 dQ 사용
            do jc = hi(2,t),lo(2,t),-1
                do ic = hi(1,t),lo(1,t),-1
                    rhs = 0.0_dp
                    if (ic < hi(1,t)) call nb(ic+1, jc, si_x(ic,jc), si_y(ic,jc), si_len(ic,jc), ri(ic,jc), rhs)
                    if (jc < hi(2,t)) call nb(ic, jc+1, sj_x(ic,jc), sj_y(ic,jc), sj_len(ic,jc), rj(ic,jc), rhs)
                    dQ(:,ic,jc) = dQ(:,ic,jc) + rhs/dg(ic,jc)
                end do
            end do
        end do
        !$omp end do

        ! ---- 갱신 (모든 타일의 dQ 가 끝난 뒤: nb 가 Q^n 을 읽으므로) ----
        !$omp do
        do jc = 1,ny
            do ic = 1,nx
                Q(:,ic,jc) = Q(:,ic,jc) + dQ(:,ic,jc)
            end do
        end do
        !$omp end do

        !$omp end parallel

    contains

        !> 이웃 셀 J = (i,j) 의 몫:  rhs -= 1/2 (dF_J.S_out - r_f dQ_J)
        !>   (fnx, fny) : I -> J 방향 단위 법선,  len : 면 길이,  rf : 그 면의 r_f
        subroutine nb(i, j, fnx, fny, len, rf, rhs)
            integer,  intent(in)    :: i, j
            real(dp), intent(in)    :: fnx, fny, len, rf
            real(dp), intent(inout) :: rhs(nvar)
            real(dp) :: Wn(nvar), F0(nvar), F1(nvar)
            call c2p(Q(:,i,j) + dQ(:,i,j), Wn)
            call euler_flux(Wp(:,i,j), fnx, fny, F0)
            call euler_flux(Wn, fnx, fny, F1)
            rhs = rhs - 0.5_dp*((F1 - F0)*len - rf*dQ(:,i,j))
        end subroutine nb

    end subroutine lusgs

    !> 면 r_f = w (|u.n| + c) |S|  (+ 2 r_v 점성):  셀 L, R 의 원시 변수 평균으로
    real(dp) function face_rad(iL, jL, iR, jR, fnx, fny, len) result(rf)
        integer,  intent(in) :: iL, jL, iR, jR
        real(dp), intent(in) :: fnx, fny, len
        real(dp) :: Wf(nvar), Tf, dist

        Wf = 0.5_dp*(Wp(:,iL,jL) + Wp(:,iR,jR))
        rf = lusgs_omega*(abs(Wf(2)*fnx + Wf(3)*fny) + sqrt(gam*Wf(4)/Wf(1)))*len
        if (viscous) then
            Tf   = 0.5_dp*(gam*Wp(4,iL,jL)/Wp(1,iL,jL) + gam*Wp(4,iR,jR)/Wp(1,iR,jR))
            dist = sqrt((xcg(iR,jR) - xcg(iL,jL))**2 + (ycg(iR,jR) - ycg(iL,jL))**2)
            rf   = rf + 2.0_dp*visc_coef*mu_star(Tf)/Wf(1)*max(4.0_dp/3.0_dp, gam/prandtl)*len/dist
        end if
    end function face_rad

    !> 비점성 플럭스 F.n (단위 법선, 넓이 1 당):  W = (rho, u, v, p)
    pure subroutine euler_flux(W, fnx, fny, F)
        real(dp), intent(in)  :: W(nvar), fnx, fny
        real(dp), intent(out) :: F(nvar)
        real(dp) :: un, E
        un = W(2)*fnx + W(3)*fny
        E  = W(4)/(gam - 1.0_dp) + 0.5_dp*W(1)*(W(2)**2 + W(3)**2)
        F(1) = W(1)*un
        F(2) = W(1)*W(2)*un + W(4)*fnx
        F(3) = W(1)*W(3)*un + W(4)*fny
        F(4) = (E + W(4))*un
    end subroutine euler_flux

end module implicit
