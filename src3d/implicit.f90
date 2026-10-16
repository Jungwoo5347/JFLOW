!> LU-SGS 암시적 시간 적분 (Yoon & Jameson), 정상 문제 전용  (flow.in 의 time_scheme = implicit)
!>
!>   역방향 오일러 + 면 플럭스 선형화 (Jameson 의 스펙트럼 반경 분리):
!>     F_f(Q+dQ).S - F_f(Q).S  ~  1/2 (dF_I + dF_J).S + 1/2 w r_f (dQ_I - dQ_J),   r_f = (|u.n| + c) |S|  (면 평균)
!>   셀 I 에 대해 (닫힌 셀이라 sum_f dF_I.S = 0):
!>     D_I dQ_I + sum_J 1/2 (dF_J.S_out - w r_f dQ_J) = V R_I,      D_I = V/dt_I + w/2 sum_f r_f   (스칼라)
!>     dF_J = F(Q_J + dQ_J) - F(Q_J)  (플럭스 차이 그대로, 야코비안 행렬 없음)
!>   L (i-1, j-1, k-1 이웃) / U (i+1, j+1, k+1 이웃) 로 나눠 두 번 스윕:
!>     앞으로 : dQ*_I = [ V R_I - sum_{J in L} 1/2 (dF*_J.S - w r_f dQ*_J) ] / D_I
!>     뒤로   : dQ_I  = dQ*_I - [ sum_{J in U} 1/2 (dF_J.S - w r_f dQ_J) ] / D_I
!>   우변 R 은 explicit 과 같은 2차 잔차 -> 수렴한 해는 같음. 좌변은 1차 근사라 수렴 속도만 바꿈
!>
!>   병렬: 블록의 큰 두 방향을 tile 셀 크기의 타일로 나누고, 타일마다 독립적으로 앞/뒤 스윕
!>     (타일 경계 건너편의 dQ 는 0, 블록 경계와 같은 처리 -> 스레드끼리 기다릴 일이 없음)
!>     타일 크기가 고정이라 스레드 수와 상관없이 결과가 같음
!>   블록 경계, 고스트 셀의 dQ 는 0 (경계조건은 매 반복 bc() 로 명시적으로)
!>   w r_f (면마다) 와 D_I 는 반복마다 한 번 미리 계산해 둠
module implicit
    use euler
    use grid
    use solver, only: flux_n
    use case_setup, only: lusgs_omega
    implicit none
    private
    public :: lusgs

    integer, parameter :: tile = 16        ! 타일 한 변의 셀 수 (큰 두 방향)

    ! 블록마다 작업 배열
    !   lo(:,t), hi(:,t) : 타일 t 의 셀 범위 (i, j, k)
    !   ri, rj, rk : 면마다 w r_f (i-면 (0:nx,ny,nz), j-면, k-면),   dg : 셀마다 D_I
    type lu_t
        integer :: ntile = 0
        integer,  allocatable :: lo(:,:), hi(:,:)
        real(dp), allocatable :: ri(:,:,:), rj(:,:,:), rk(:,:,:), dg(:,:,:)
    end type
    type(lu_t), allocatable, target :: lu(:)

contains

    !> 블록 b 의 작업 배열과 타일 목록 (처음 한 번)
    subroutine setup_lu(b)
        integer, intent(in) :: b
        integer :: n(3), d1, d2, d3, nt(3), t1, t2, t
        integer :: lo(3), hi(3)

        if (.not. allocated(lu)) allocate(lu(nblk))
        if (lu(b)%ntile > 0) return
        n = [blk(b)%nx, blk(b)%ny, blk(b)%nz]
        associate (L => lu(b))
            allocate(L%ri(0:n(1),n(2),n(3)), L%rj(n(1),0:n(2),n(3)), L%rk(n(1),n(2),0:n(3)), L%dg(n(1),n(2),n(3)))

            ! 셀이 가장 많은 두 방향 d1, d2 를 타일로 나눔 (d3 방향은 통째로)
            d3 = minloc(n, dim=1)
            d1 = merge(2, 1, d3 == 1);  d2 = 6 - d1 - d3
            nt = 1
            nt(d1) = (n(d1) + tile - 1)/tile
            nt(d2) = (n(d2) + tile - 1)/tile
            L%ntile = nt(d1)*nt(d2)
            allocate(L%lo(3,L%ntile), L%hi(3,L%ntile))
            t = 0
            do t2 = 1,nt(d2)
                do t1 = 1,nt(d1)
                    t = t + 1
                    lo = 1;  hi = n
                    lo(d1) = (t1-1)*tile + 1;  hi(d1) = min(t1*tile, n(d1))
                    lo(d2) = (t2-1)*tile + 1;  hi(d2) = min(t2*tile, n(d2))
                    L%lo(:,t) = lo;  L%hi(:,t) = hi
                end do
            end do
        end associate
    end subroutine setup_lu

    !> 블록 b 한 반복: dQ 를 구해 Q 에 더함
    !>   입력: blk(b)%R (Q^n 의 잔차 dQ/dt), blk(b)%W (Q^n 원시 변수, 고스트 포함), blk(b)%dtl (셀마다 dt)
    !>   blk(b)%Rt 를 dQ 저장에 씀
    subroutine lusgs(b)
        integer, intent(in) :: b
        integer  :: t, ic, jc, kc, fi, fj, fk, l(3), h(3)
        real(dp) :: rhs(nvar)
        real(dp), pointer, contiguous :: Q(:,:,:,:), W(:,:,:,:), R(:,:,:,:), dQ(:,:,:,:), dtl(:,:,:)
        real(dp), pointer, contiguous :: ri(:,:,:), rj(:,:,:), rk(:,:,:), dg(:,:,:)

        call setup_lu(b)
        call use_block(b)
        Q => blk(b)%Q;  W => blk(b)%W;  R => blk(b)%R;  dQ => blk(b)%Rt;  dtl => blk(b)%dtl
        ri => lu(b)%ri;  rj => lu(b)%rj;  rk => lu(b)%rk;  dg => lu(b)%dg

        !$omp parallel private(t, ic, jc, kc, fi, fj, fk, l, h, rhs)

        ! ---- 면마다 w r_f (양쪽 셀 원시 변수 평균, 경계 면은 고스트와 평균) ----
        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do fi = 0,nx
                    ri(fi,jc,kc) = lusgs_omega*rad(W(:,fi,jc,kc), W(:,fi+1,jc,kc), si(:,fi,jc,kc))*si_a(fi,jc,kc)
                end do
            end do
        end do
        !$omp end do nowait
        !$omp do collapse(2)
        do kc = 1,nz
            do fj = 0,ny
                do ic = 1,nx
                    rj(ic,fj,kc) = lusgs_omega*rad(W(:,ic,fj,kc), W(:,ic,fj+1,kc), sj(:,ic,fj,kc))*sj_a(ic,fj,kc)
                end do
            end do
        end do
        !$omp end do nowait
        !$omp do collapse(2)
        do fk = 0,nz
            do jc = 1,ny
                do ic = 1,nx
                    rk(ic,jc,fk) = lusgs_omega*rad(W(:,ic,jc,fk), W(:,ic,jc,fk+1), sk(:,ic,jc,fk))*sk_a(ic,jc,fk)
                end do
            end do
        end do
        !$omp end do

        ! ---- 대각항 D_I = V/dt + 1/2 sum_6면 w r_f ----
        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    dg(ic,jc,kc) = vol(ic,jc,kc)/dtl(ic,jc,kc) + 0.5_dp*( ri(ic-1,jc,kc) + ri(ic,jc,kc) &
                                 + rj(ic,jc-1,kc) + rj(ic,jc,kc) + rk(ic,jc,kc-1) + rk(ic,jc,kc) )
                end do
            end do
        end do
        !$omp end do

        ! ---- 타일마다 앞/뒤 스윕 + 갱신 (타일끼리 독립) ----
        !$omp do schedule(dynamic,1)
        do t = 1,lu(b)%ntile
            l = lu(b)%lo(:,t);  h = lu(b)%hi(:,t)

            ! 앞으로 (i, j, k 오름차순): 아래 이웃 (타일 안) 의 dQ* 사용
            do kc = l(3),h(3)
                do jc = l(2),h(2)
                    do ic = l(1),h(1)
                        rhs = vol(ic,jc,kc)*R(:,ic,jc,kc)
                        if (ic > l(1)) call nb(ic-1, jc, kc, -si(:,ic-1,jc,kc), si_a(ic-1,jc,kc), ri(ic-1,jc,kc), rhs)
                        if (jc > l(2)) call nb(ic, jc-1, kc, -sj(:,ic,jc-1,kc), sj_a(ic,jc-1,kc), rj(ic,jc-1,kc), rhs)
                        if (kc > l(3)) call nb(ic, jc, kc-1, -sk(:,ic,jc,kc-1), sk_a(ic,jc,kc-1), rk(ic,jc,kc-1), rhs)
                        dQ(:,ic,jc,kc) = rhs/dg(ic,jc,kc)
                    end do
                end do
            end do

            ! 뒤로 (내림차순): 위 이웃 (타일 안) 의 최종 dQ 사용
            do kc = h(3),l(3),-1
                do jc = h(2),l(2),-1
                    do ic = h(1),l(1),-1
                        rhs = 0.0_dp
                        if (ic < h(1)) call nb(ic+1, jc, kc, si(:,ic,jc,kc), si_a(ic,jc,kc), ri(ic,jc,kc), rhs)
                        if (jc < h(2)) call nb(ic, jc+1, kc, sj(:,ic,jc,kc), sj_a(ic,jc,kc), rj(ic,jc,kc), rhs)
                        if (kc < h(3)) call nb(ic, jc, kc+1, sk(:,ic,jc,kc), sk_a(ic,jc,kc), rk(ic,jc,kc), rhs)
                        dQ(:,ic,jc,kc) = dQ(:,ic,jc,kc) + rhs/dg(ic,jc,kc)
                    end do
                end do
            end do
        end do
        !$omp end do

        ! ---- 갱신 (모든 타일의 dQ 가 끝난 뒤: nb 가 Q^n 을 읽으므로) ----
        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    Q(:,ic,jc,kc) = Q(:,ic,jc,kc) + dQ(:,ic,jc,kc)
                end do
            end do
        end do
        !$omp end do

        !$omp end parallel

    contains

        !> 이웃 셀 J 의 몫:  rhs -= 1/2 (dF_J.S_out - w r_f dQ_J)
        !>   n : I -> J 방향 단위 법선,  a : 면 넓이,  wr : 그 면의 w r_f
        subroutine nb(i, j, k, n, a, wr, rhs)
            integer,  intent(in)    :: i, j, k
            real(dp), intent(in)    :: n(3), a, wr
            real(dp), intent(inout) :: rhs(nvar)
            real(dp) :: Wn(nvar), F0(nvar), F1(nvar), dqj(nvar)
            dqj = dQ(:,i,j,k)
            call c2p(Q(:,i,j,k) + dqj, Wn)
            call flux_n(W(:,i,j,k), n, F0)
            call flux_n(Wn, n, F1)
            rhs = rhs - 0.5_dp*((F1 - F0)*a - wr*dqj)
        end subroutine nb

    end subroutine lusgs

    !> 면 스펙트럼 반경 (넓이 1 당):  |u.n| + c,  양쪽 원시 변수 평균으로
    pure real(dp) function rad(WI, WJ, n)
        real(dp), intent(in) :: WI(nvar), WJ(nvar), n(3)
        real(dp) :: Wf(nvar)
        Wf  = 0.5_dp*(WI + WJ)
        rad = abs(Wf(2)*n(1) + Wf(3)*n(2) + Wf(4)*n(3)) + sqrt(gam*Wf(5)/Wf(1))
    end function rad

end module implicit
