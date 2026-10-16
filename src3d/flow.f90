module flow
    use euler
    use grid
    use case_setup
    implicit none

    real(dp), allocatable :: Q(:,:,:,:)

    ! 경계조건은 경계 면 하나하나마다 따로 저장 (한 변을 여러 구간으로 나눌 수 있게)
    !   side : 1 imin, 2 imax, 3 jmin, 4 jmax, 5 kmin, 6 kmax
    !   (fa, fb) : 그 변 위의 면 번호 두 개 (1 부터).  변마다 (fa, fb) 가 가리키는 방향:
    !            imin/imax : (jc, kc),   jmin/jmax : (ic, kc),   kmin/kmax : (ic, jc)
    !   bc_face(fa, fb, side)  : 경계조건 코드 (case_setup 의 BC_...)
    !   bc_state(fa, fb, side) : 입구/원방일 때 고스트에 넣을 상태 번호 (case_setup 의 states)
    !   cut_side, cut_p, cut_q : BC_cut 면이 맞닿은 짝 면 (wake cut, O-grid 이음매)
    ! 값은 grid.in 의 bc-1-f 줄로 setup_bc 에서 채움
    integer, allocatable :: bc_face(:,:,:), bc_state(:,:,:)
    integer, allocatable :: cut_side(:,:,:), cut_p(:,:,:), cut_q(:,:,:)


    contains

    !> 초기조건: 전체를 initial 상태로 채운 뒤, region 상자 안의 셀은 그 상태로 (뒤에 쓴 region 우선)
    subroutine init_flow()

        integer :: ic, jc, kc, r
        real(dp) :: W(nvar)

        allocate(Q(nvar,1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
        Q = 0.0_dp                       ! 고스트(모서리 포함)에 쓰레기값이 남지 않게

        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    W = states(init_state)%W
                    do r = 1,nregions
                        associate (g => regions(r))
                            if (xc(ic,jc,kc) >= g%x0 .and. xc(ic,jc,kc) < g%x1 .and. &
                                yc(ic,jc,kc) >= g%y0 .and. yc(ic,jc,kc) < g%y1 .and. &
                                zc(ic,jc,kc) >= g%z0 .and. zc(ic,jc,kc) < g%z1) then
                                W = states(g%istate)%W
                            end if
                        end associate
                    end do
                    call p2c(W,Q(:,ic,jc,kc))
                end do
            end do
        end do

    end subroutine init_flow

    !> 경계조건 저장 배열 할당 (처음 한 번). 처음엔 전부 '지정 안 됨'
    subroutine alloc_bc()
        integer :: na, nb
        if (allocated(bc_face)) return
        na = max(nx, ny);  nb = max(ny, nz)
        allocate(bc_face(na,nb,nside), bc_state(na,nb,nside))
        allocate(cut_side(na,nb,nside), cut_p(na,nb,nside), cut_q(na,nb,nside))
        bc_face  = BC_unset
        bc_state = 0
        cut_side = 0
        cut_p    = 0
        cut_q    = 0
    end subroutine alloc_bc

    !> 변 side 의 면 개수 두 방향 (fa 방향 na, fb 방향 nb)
    subroutine nface(side, na, nb)
        integer, intent(in)  :: side
        integer, intent(out) :: na, nb
        select case (side)
        case (1,2);    na = ny;  nb = nz
        case (3,4);    na = nx;  nb = nz
        case default;  na = nx;  nb = ny
        end select
    end subroutine nface

    !> 경계 면 (side, fa, fb) 에서 안쪽으로 m 번째 내부 셀 (m = 1 이 경계에 붙은 셀)
    !>   격자가 m 칸보다 얇으면 가장 안쪽 셀 (예: Sod 처럼 한 방향 셀이 1 개)
    subroutine inner_cell(side, fa, fb, m, ic, jc, kc)
        integer, intent(in)  :: side, fa, fb, m
        integer, intent(out) :: ic, jc, kc
        select case (side)
        case (1);  ic = min(m, nx);       jc = fa;  kc = fb
        case (2);  ic = max(nx+1-m, 1);   jc = fa;  kc = fb
        case (3);  ic = fa;  jc = min(m, ny);       kc = fb
        case (4);  ic = fa;  jc = max(ny+1-m, 1);   kc = fb
        case (5);  ic = fa;  jc = fb;  kc = min(m, nz)
        case default;  ic = fa;  jc = fb;  kc = max(nz+1-m, 1)
        end select
    end subroutine inner_cell

    !> 경계 면 (side, fa, fb) 에서 바깥으로 m 번째 고스트 셀 (m = 1 이 경계에 붙은 고스트)
    subroutine ghost_cell(side, fa, fb, m, ic, jc, kc)
        integer, intent(in)  :: side, fa, fb, m
        integer, intent(out) :: ic, jc, kc
        select case (side)
        case (1);  ic = 1-m;    jc = fa;  kc = fb
        case (2);  ic = nx+m;   jc = fa;  kc = fb
        case (3);  ic = fa;  jc = 1-m;    kc = fb
        case (4);  ic = fa;  jc = ny+m;   kc = fb
        case (5);  ic = fa;  jc = fb;  kc = 1-m
        case default;  ic = fa;  jc = fb;  kc = nz+m
        end select
    end subroutine ghost_cell

    !> 경계 면 (side, fa, fb) 의 네 꼭짓점 좌표 pts(:, 1:4)
    subroutine face_corners(side, fa, fb, pts)
        integer,  intent(in)  :: side, fa, fb
        real(dp), intent(out) :: pts(3,4)
        integer :: f
        select case (side)
        case (1,2)
            f = merge(0, nx, side == 1)
            pts(:,1) = node(f,fa-1,fb-1);  pts(:,2) = node(f,fa,fb-1);  pts(:,3) = node(f,fa,fb);  pts(:,4) = node(f,fa-1,fb)
        case (3,4)
            f = merge(0, ny, side == 3)
            pts(:,1) = node(fa-1,f,fb-1);  pts(:,2) = node(fa,f,fb-1);  pts(:,3) = node(fa,f,fb);  pts(:,4) = node(fa-1,f,fb)
        case default
            f = merge(0, nz, side == 5)
            pts(:,1) = node(fa-1,fb-1,f);  pts(:,2) = node(fa,fb-1,f);  pts(:,3) = node(fa,fb,f);  pts(:,4) = node(fa-1,fb,f)
        end select
    end subroutine face_corners

    !> 경계 면 (side, fa, fb) 의 단위 법선, 영역 바깥 방향
    !>   (격자 법선은 +i, +j, +k 방향이라 imin, jmin, kmin 에서는 부호를 뒤집음)
    function out_normal(side, fa, fb) result(n)
        integer, intent(in) :: side, fa, fb
        real(dp) :: n(3)
        select case (side)
        case (1);  n = -si(:,0,fa,fb)
        case (2);  n =  si(:,nx,fa,fb)
        case (3);  n = -sj(:,fa,0,fb)
        case (4);  n =  sj(:,fa,ny,fb)
        case (5);  n = -sk(:,fa,fb,0)
        case default;  n = sk(:,fa,fb,nz)
        end select
    end function out_normal

    !> grid.in 의 bc 줄들을 경계 면에 적용 (뒤에 쓴 줄이 우선)
    !>   노드 a1 ~ a2 사이의 면 = 면 번호 a1 ~ a2-1 (b 도 같음)
    !>   지정 안 된 면은 기하학적으로 맞닿은 다른 지정 안 된 면을 찾아 자동 연결 (BC_cut)
    !>   그래도 남는 면이 있으면 어디인지 알려주고 멈춤
    subroutine setup_bc()
        integer :: k, side, na, nb, a1, a2, b1, b2
        logical :: any_cut

        call alloc_bc()

        do k = 1,npatches
            side = patches(k)%side
            call nface(side, na, nb)
            a1 = patches(k)%a1;  a2 = patches(k)%a2
            b1 = patches(k)%b1;  b2 = patches(k)%b2
            if (a2 == 0) a2 = na + 1
            if (b2 == 0) b2 = nb + 1
            if (a1 < 1 .or. a2 > na+1 .or. a1 >= a2 .or. b1 < 1 .or. b2 > nb+1 .or. b1 >= b2) then
                write(*,'(a,i0,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') ' Error: grid.in bc line #', k, &
                    ' (face ', side_name(side), ') nodes ', a1, '-', a2, ' ; ', b1, '-', b2, &
                    ' out of range 1-', na+1, ' ; 1-', nb+1, ' (need a1 < a2, b1 < b2)'
                error stop
            end if
            bc_face (a1:a2-1, b1:b2-1, side) = patches(k)%btype
            bc_state(a1:a2-1, b1:b2-1, side) = patches(k)%istate
        end do

        call connect_cuts()
        any_cut = report_faces(BC_cut, 'auto-connected (cut)')
        if (report_faces(BC_unset, 'Error: no boundary condition and no matching face on')) error stop
    end subroutine setup_bc

    !> 지정 안 된 경계 면끼리 네 꼭짓점이 겹치면 (순서 무관) 서로 연결
    !>   C-grid wake cut 은 jmin 의 앞쪽과 뒤쪽이, O-grid 이음매는 imin 과 imax 가 겹침
    !>   면이 많아서 (3D) 모든 짝을 비교하지 않고, 면 중심 x 로 정렬한 뒤 x 가 가까운 면끼리만 비교
    subroutine connect_cuts()
        integer  :: side, na, nb, fa, fb, n, a, c, i1, i2
        integer,  allocatable :: fs(:), fp(:), fq(:), ord(:)
        real(dp), allocatable :: pts(:,:,:), cen(:,:), tol(:)
        real(dp) :: tmax

        ! 지정 안 된 면 목록
        n = 0
        do side = 1,nside
            call nface(side, na, nb)
            n = n + count(bc_face(1:na,1:nb,side) == BC_unset)
        end do
        allocate(fs(n), fp(n), fq(n), pts(3,4,n), cen(3,n), tol(n))
        n = 0
        do side = 1,nside
            call nface(side, na, nb)
            do fb = 1,nb
                do fa = 1,na
                    if (bc_face(fa,fb,side) /= BC_unset) cycle
                    n = n + 1
                    fs(n) = side;  fp(n) = fa;  fq(n) = fb
                    call face_corners(side, fa, fb, pts(:,:,n))
                    cen(:,n) = 0.25_dp*sum(pts(:,:,n), dim=2)
                    tol(n)   = 1.0e-6_dp*max(norm2(pts(:,3,n) - pts(:,1,n)), norm2(pts(:,4,n) - pts(:,2,n)))
                end do
            end do
        end do
        if (n == 0) return
        tmax = maxval(tol(1:n))

        allocate(ord(n))
        ord = [(a, a = 1,n)]
        call sort_index(cen(1,1:n), ord)

        do a = 1,n
            i1 = ord(a)
            if (bc_face(fp(i1),fq(i1),fs(i1)) /= BC_unset) cycle
            do c = a+1,n
                i2 = ord(c)
                if (cen(1,i2) - cen(1,i1) > tmax) exit
                if (bc_face(fp(i2),fq(i2),fs(i2)) /= BC_unset) cycle
                if (norm2(cen(:,i2) - cen(:,i1)) > tol(i1)) cycle
                if (.not. same_corners(pts(:,:,i1), pts(:,:,i2), tol(i1))) cycle
                bc_face (fp(i1),fq(i1),fs(i1)) = BC_cut
                cut_side(fp(i1),fq(i1),fs(i1)) = fs(i2)
                cut_p   (fp(i1),fq(i1),fs(i1)) = fp(i2)
                cut_q   (fp(i1),fq(i1),fs(i1)) = fq(i2)
                bc_face (fp(i2),fq(i2),fs(i2)) = BC_cut
                cut_side(fp(i2),fq(i2),fs(i2)) = fs(i1)
                cut_p   (fp(i2),fq(i2),fs(i2)) = fp(i1)
                cut_q   (fp(i2),fq(i2),fs(i2)) = fq(i1)
                exit
            end do
        end do

    contains

        !> A 의 네 꼭짓점이 모두 B 의 어떤 꼭짓점과 tol 안에서 겹치면 .true.
        logical function same_corners(pa, pb, tol)
            real(dp), intent(in) :: pa(3,4), pb(3,4), tol
            integer :: i, j
            logical :: hit
            same_corners = .false.
            do i = 1,4
                hit = .false.
                do j = 1,4
                    if (norm2(pa(:,i) - pb(:,j)) < tol) hit = .true.
                end do
                if (.not. hit) return
            end do
            same_corners = .true.
        end function same_corners

    end subroutine connect_cuts

    !> key 값이 작은 순서로 idx 를 정렬 (병합 정렬, 같은 값은 원래 순서 유지)
    subroutine sort_index(key, idx)
        real(dp), intent(in)    :: key(:)
        integer,  intent(inout) :: idx(:)
        integer, allocatable :: tmp(:)
        integer :: w, lo, mid, hi, i, j, k, n
        n = size(idx)
        allocate(tmp(n))
        w = 1
        do while (w < n)
            lo = 1
            do while (lo <= n - w)
                mid = lo + w - 1
                hi  = min(lo + 2*w - 1, n)
                i = lo;  j = mid + 1;  k = lo
                do while (i <= mid .and. j <= hi)
                    if (key(idx(j)) < key(idx(i))) then
                        tmp(k) = idx(j);  j = j + 1
                    else
                        tmp(k) = idx(i);  i = i + 1
                    end if
                    k = k + 1
                end do
                tmp(k:k+mid-i) = idx(i:mid);   k = k + mid - i + 1
                tmp(k:k+hi-j)  = idx(j:hi)
                idx(lo:hi) = tmp(lo:hi)
                lo = lo + 2*w
            end do
            w = 2*w
        end do
    end subroutine sort_index

    !> 경계조건이 code 인 면을 변마다 개수와 노드 범위 (감싸는 상자) 로 출력. 하나라도 있으면 .true.
    !>   cut 은 짝 변별로 나눠서 출력
    logical function report_faces(code, label) result(found)
        integer,          intent(in) :: code
        character(len=*), intent(in) :: label
        integer :: side, s2, na, nb
        logical, allocatable :: mask(:,:)
        found = .false.
        do side = 1,nside
            call nface(side, na, nb)
            if (code == BC_cut) then
                do s2 = 1,nside
                    mask = bc_face(1:na,1:nb,side) == BC_cut .and. cut_side(1:na,1:nb,side) == s2
                    if (.not. any(mask)) cycle
                    call show(side, mask, '  <->  '//side_name(s2))
                    found = .true.
                end do
            else
                mask = bc_face(1:na,1:nb,side) == code
                if (.not. any(mask)) cycle
                call show(side, mask, '')
                found = .true.
            end if
        end do

    contains

        subroutine show(side, mask, tail)
            integer,          intent(in) :: side
            logical,          intent(in) :: mask(:,:)
            character(len=*), intent(in) :: tail
            integer :: fa, fb, a0, a1, b0, b1
            a0 = huge(1);  a1 = 0;  b0 = huge(1);  b1 = 0
            do fb = 1,size(mask,2)
                do fa = 1,size(mask,1)
                    if (.not. mask(fa,fb)) cycle
                    a0 = min(a0, fa);  a1 = max(a1, fa);  b0 = min(b0, fb);  b1 = max(b1, fb)
                end do
            end do
            write(*,'(1x,a,1x,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,a)') label, side_name(side), ', nodes ', a0, '-', a1+1, &
                ' ; ', b0, '-', b1+1, '  (', count(mask), ' faces)', tail
        end subroutine show

    end function report_faces

    !> 고스트 셀 채우기 (매 스텝, 잔차 계산 전에 호출)
    !>   경계 면마다 고스트 ng 층을 채움. m 번째 고스트는 경계를 사이에 두고 m 번째 내부 셀과 마주봄
    !>     벽/대칭 : m 번째 내부 셀을 거울 반사 (경계 면 법선 기준)
    !>     cut     : 짝 면 건너편의 m 번째 내부 셀
    !>     출구    : 모든 층에 경계에 붙은 내부 셀 (기울기 0 외삽)
    !>     원방    : 모든 층에 경계에 붙은 내부 셀로 계산한 같은 경계 상태
    !>     입구    : 모든 층에 지정 상태
    subroutine bc()
        integer :: side, na, nb, fa, fb, m

        call alloc_bc()

        do side = 1,nside
            call nface(side, na, nb)
            !$omp parallel do private(fa, m)
            do fb = 1,nb
                do fa = 1,na
                    do m = 1,ng
                        call fill(side, fa, fb, m)
                    end do
                end do
            end do
            !$omp end parallel do
        end do

    contains

        !> 경계 면 (side, fa, fb) 의 m 번째 고스트 셀 채우기
        subroutine fill(side, fa, fb, m)
            integer, intent(in) :: side, fa, fb, m
            integer  :: ig, jg, kg, ii, ji, ki, k
            real(dp) :: Wb(nvar)
            call ghost_cell(side, fa, fb, m, ig, jg, kg)
            k = bc_state(fa,fb,side)
            if (k > 0) then
                Wb = states(k)%W
            else
                Wb = 0.0_dp
            end if
            select case (bc_face(fa,fb,side))
            case (BC_cut)
                ! 맞닿은 짝 면 건너편의 m 번째 내부 셀 값을 그대로
                call inner_cell(cut_side(fa,fb,side), cut_p(fa,fb,side), cut_q(fa,fb,side), m, ii, ji, ki)
                Q(:,ig,jg,kg) = Q(:,ii,ji,ki)
            case (BC_outlet, BC_farfield)
                call inner_cell(side, fa, fb, 1, ii, ji, ki)
                call ghost(bc_face(fa,fb,side), Wb, Q(:,ii,ji,ki), out_normal(side,fa,fb), Q(:,ig,jg,kg))
            case default
                call inner_cell(side, fa, fb, m, ii, ji, ki)
                call ghost(bc_face(fa,fb,side), Wb, Q(:,ii,ji,ki), out_normal(side,fa,fb), Q(:,ig,jg,kg))
            end select
        end subroutine fill

    end subroutine bc

    !> 경계 면 하나의 고스트 셀 값
    !>   Wb  : 입구/원방 상태 (rho, u, v, w, fa)
    !>   qin : 경계에 붙은 내부 셀,  n : 경계 면 단위 법선 (영역 바깥 방향),  qg : 고스트 셀
    subroutine ghost(btype, Wb, qin, n, qg)
        integer,  intent(in)  :: btype
        real(dp), intent(in)  :: Wb(nvar), qin(nvar), n(3)
        real(dp), intent(out) :: qg(nvar)

        real(dp) :: mn
        real(dp) :: Wi(nvar), Wg(nvar)
        real(dp) :: uni, ci, si_, unf, cf, sf, uti(3), utf(3), utb(3)
        real(dp) :: Rp, Rm, unb, cb, sb, rhob

        select case (btype)

        case (BC_inlet)
            ! 지정한 상태 고정 (초음속 입구)
            call p2c(Wb, qg)

        case (BC_farfield)
            ! 리만 불변량 원방 조건 (n 은 바깥 방향, un > 0 이 유출)
            !   R+ = un + 2c/(g-1) 은 속도 un + c,  R- = un - 2c/(g-1) 은 속도 un - c 로 이동
            !   나가는 파동은 안쪽 셀(i)에서, 들어오는 파동은 자유류(f)에서 가져옴
            !   엔트로피 s = fa/rho^g 와 접선 속도 (벡터) 는 유입이면 자유류, 유출이면 안쪽 셀
            call c2p(qin, Wi)
            uni = dot_product(Wi(2:4), n)
            uti = Wi(2:4) - uni*n
            ci  = sqrt(gam*Wi(5)/Wi(1))
            si_ = Wi(5)/Wi(1)**gam
            unf = dot_product(Wb(2:4), n)
            utf = Wb(2:4) - unf*n
            cf  = sqrt(gam*Wb(5)/Wb(1))
            sf  = Wb(5)/Wb(1)**gam

            if (unf <= -cf) then
                Wg = Wb                                    ! 초음속 유입: 모두 자유류
            else if (uni >= ci) then
                Wg = Wi                                    ! 초음속 유출: 모두 안쪽
            else
                Rp  = uni + 2.0_dp*ci/(gam - 1.0_dp)
                Rm  = unf - 2.0_dp*cf/(gam - 1.0_dp)
                unb = 0.5_dp*(Rp + Rm)
                cb  = 0.25_dp*(gam - 1.0_dp)*(Rp - Rm)
                if (unb < 0.0_dp) then
                    sb = sf;   utb = utf
                else
                    sb = si_;  utb = uti
                end if
                rhob    = (cb**2/(gam*sb))**(1.0_dp/(gam - 1.0_dp))
                Wg(1)   = rhob
                Wg(2:4) = unb*n + utb
                Wg(5)   = rhob*cb**2/gam
            end if
            call p2c(Wg, qg)

        case (BC_outlet)
            ! 내부 셀 그대로 복사
            qg = qin

        case (BC_wall, BC_symmetry)
            ! 밀도, 에너지는 그대로, 운동량의 법선 성분만 반전 (비점성이라 대칭면도 같음)
            mn = dot_product(qin(2:4), n)
            qg(1)   = qin(1)
            qg(2:4) = qin(2:4) - 2.0_dp*mn*n
            qg(5)   = qin(5)

        case default
            write(*,*) 'Unknown BC type:', btype
            error stop

        end select

    end subroutine ghost

end module flow
