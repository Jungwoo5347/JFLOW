module flow
    use euler
    use grid
    use case_setup
    implicit none

    real(dp), allocatable :: Q(:,:,:)

    ! 경계조건은 경계 면 하나하나마다 따로 저장 (한 변을 여러 구간으로 나눌 수 있게)
    !   side : 1 = imin (고스트 ic=0,-1..), 2 = imax (ic=nx+1,nx+2..), 3 = jmin (jc=0,-1..), 4 = jmax (jc=ny+1,ny+2..)
    !   k    : 그 변 위의 면 번호.  imin/imax 는 k = jc (1:ny),  jmin/jmax 는 k = ic (1:nx)
    !          면 k 는 (1부터 센) 노드 k 와 k+1 사이
    !   bc_face(k, side) : 경계조건 코드 (case_setup 의 BC_...)
    !   W_bc(:, k, side) : 입구/원방일 때 고스트에 넣을 상태 (rho, u, v, p)
    !   cut_side, cut_k  : BC_cut 면이 맞닿은 짝 면 (wake cut, O-grid 이음매)
    ! 값은 grid.in 의 bc-1-f 줄로 setup_bc 에서 채움
    integer,  allocatable :: bc_face(:,:), cut_side(:,:), cut_k(:,:)
    real(dp), allocatable :: W_bc(:,:,:)


    contains

    !> 초기조건: 전체를 &initial 상태로 채운 뒤, &region 상자 안의 셀은 그 상태로 (뒤에 쓴 region 우선)
    subroutine init_flow()

        integer :: ic, jc, r
        real(dp) :: W(nvar)

        allocate(Q(nvar,1-ng:nx+ng,1-ng:ny+ng))
        Q = 0.0_dp                       ! 고스트(모서리 포함)에 쓰레기값이 남지 않게

        do jc = 1,ny
            do ic = 1,nx
                W = states(init_state)%W
                do r = 1,nregions
                    if (xc(ic,jc) >= regions(r)%x0 .and. xc(ic,jc) < regions(r)%x1 .and. &
                        yc(ic,jc) >= regions(r)%y0 .and. yc(ic,jc) < regions(r)%y1) then
                        W = states(regions(r)%istate)%W
                    end if
                end do
                call p2c(W,Q(:,ic,jc))
            end do
        end do

    end subroutine init_flow

    !> 경계조건 저장 배열 할당 (처음 한 번). 처음엔 전부 '지정 안 됨'
    subroutine alloc_bc()
        if (allocated(bc_face)) return
        allocate(bc_face(max(nx,ny), 4), W_bc(nvar, max(nx,ny), 4))
        allocate(cut_side(max(nx,ny), 4), cut_k(max(nx,ny), 4))
        bc_face  = BC_unset
        W_bc     = 0.0_dp
        cut_side = 0
        cut_k    = 0
    end subroutine alloc_bc

    !> 변 side 의 면 개수 (imin/imax: ny, jmin/jmax: nx)
    integer function nface(side)
        integer, intent(in) :: side
        if (side <= 2) then
            nface = ny
        else
            nface = nx
        end if
    end function nface

    !> 경계 면 (side, k) 에서 안쪽으로 m 번째 내부 셀 (m = 1 이 경계에 붙은 셀)
    !>   격자가 m 칸보다 얇으면 가장 안쪽 셀 (예: Sod 처럼 한 방향 셀이 1 개)
    subroutine inner_cell(side, k, m, ic, jc)
        integer, intent(in)  :: side, k, m
        integer, intent(out) :: ic, jc
        select case (side)
        case (1);  ic = min(m, nx);       jc = k
        case (2);  ic = max(nx+1-m, 1);   jc = k
        case (3);  ic = k;   jc = min(m, ny)
        case default;  ic = k;   jc = max(ny+1-m, 1)
        end select
    end subroutine inner_cell

    !> 경계 면 (side, k) 에서 바깥으로 m 번째 고스트 셀 (m = 1 이 경계에 붙은 고스트)
    subroutine ghost_cell(side, k, m, ic, jc)
        integer, intent(in)  :: side, k, m
        integer, intent(out) :: ic, jc
        select case (side)
        case (1);  ic = 1-m;    jc = k
        case (2);  ic = nx+m;   jc = k
        case (3);  ic = k;   jc = 1-m
        case default;  ic = k;   jc = ny+m
        end select
    end subroutine ghost_cell

    !> 경계 면 (side, k) 의 양 끝 노드 좌표
    subroutine face_ends(side, k, xa, ya, xb, yb)
        integer,  intent(in)  :: side, k
        real(dp), intent(out) :: xa, ya, xb, yb
        integer :: ia, ja, ib, jb
        select case (side)
        case (1);  ia = 0;   ja = k-1;  ib = 0;   jb = k
        case (2);  ia = nx;  ja = k-1;  ib = nx;  jb = k
        case (3);  ia = k-1; ja = 0;    ib = k;   jb = 0
        case default;  ia = k-1; ja = ny;   ib = k;   jb = ny
        end select
        xa = x(ia,ja);  ya = y(ia,ja);  xb = x(ib,jb);  yb = y(ib,jb)
    end subroutine face_ends

    !> grid.in 의 bc 줄들을 경계 면에 적용 (뒤에 쓴 줄이 우선)
    !>   노드 n1 ~ n2 사이의 면 = 면 번호 n1 ~ n2-1
    !>   지정 안 된 면은 기하학적으로 맞닿은 다른 지정 안 된 면을 찾아 자동 연결 (BC_cut)
    !>   그래도 남는 면이 있으면 어디인지 알려주고 멈춤
    subroutine setup_bc()
        integer :: p, side, nnode, n1, n2, k
        logical :: any_cut

        call alloc_bc()

        do p = 1,npatches
            side  = patches(p)%side
            nnode = nface(side) + 1
            n1 = patches(p)%n1
            n2 = patches(p)%n2
            if (n2 == 0) n2 = nnode
            if (n1 < 1 .or. n2 > nnode .or. n1 >= n2) then
                write(*,'(a,i0,a,a,a,i0,a,i0,a,i0,a)') ' Error: grid.in bc line #', p, ' (face ', side_name(side), &
                    ') nodes ', n1, '-', n2, ' out of range 1-', nnode, ' (need n1 < n2)'
                error stop
            end if
            do k = n1, n2-1
                bc_face(k,side) = patches(p)%btype
                if (patches(p)%istate > 0) W_bc(:,k,side) = states(patches(p)%istate)%W
            end do
        end do

        call connect_cuts()
        any_cut = report_faces(BC_cut, 'auto-connected (cut)')
        if (report_faces(BC_unset, 'Error: no boundary condition and no matching face on')) error stop
    end subroutine setup_bc

    !> 지정 안 된 경계 면끼리 양 끝 노드가 겹치면 (순서 무관) 서로 연결
    !>   C-grid wake cut 은 jmin 의 앞쪽과 뒤쪽이, O-grid 이음매는 imin 과 imax 가 겹침
    subroutine connect_cuts()
        integer  :: s1, k1, s2, k2
        real(dp) :: ax, ay, bx, by, cx, cy, dx, dy, tol

        do s1 = 1,4
            do k1 = 1,nface(s1)
                if (bc_face(k1,s1) /= BC_unset) cycle
                call face_ends(s1, k1, ax, ay, bx, by)
                tol = 1.0e-6_dp*hypot(bx - ax, by - ay)
                search: do s2 = 1,4
                    do k2 = 1,nface(s2)
                        if (s2 == s1 .and. k2 == k1) cycle
                        if (bc_face(k2,s2) /= BC_unset) cycle
                        call face_ends(s2, k2, cx, cy, dx, dy)
                        if ((hypot(ax-cx, ay-cy) < tol .and. hypot(bx-dx, by-dy) < tol) .or. &
                            (hypot(ax-dx, ay-dy) < tol .and. hypot(bx-cx, by-cy) < tol)) then
                            bc_face(k1,s1) = BC_cut;  cut_side(k1,s1) = s2;  cut_k(k1,s1) = k2
                            bc_face(k2,s2) = BC_cut;  cut_side(k2,s2) = s1;  cut_k(k2,s2) = k1
                            exit search
                        end if
                    end do
                end do search
            end do
        end do
    end subroutine connect_cuts

    !> 경계조건이 code 인 면을 연속 구간마다 노드 번호로 출력. 하나라도 있으면 .true.
    logical function report_faces(code, label) result(found)
        integer,          intent(in) :: code
        character(len=*), intent(in) :: label
        integer :: side, k, k0
        found = .false.
        do side = 1,4
            k = 1
            do while (k <= nface(side))
                if (bc_face(k,side) == code) then
                    k0 = k
                    do while (k < nface(side))
                        if (bc_face(k+1,side) /= code) exit
                        k = k + 1
                    end do
                    if (code == BC_cut) then
                        write(*,'(1x,a,1x,a,a,i0,a,i0,a,a,a,i0,a,i0)') label, side_name(side), ' nodes ', k0, '-', k+1, &
                            '  <->  ', side_name(cut_side(k0,side)), ' nodes ', &
                            min(cut_k(k0,side), cut_k(k,side)), '-', max(cut_k(k0,side), cut_k(k,side)) + 1
                    else
                        write(*,'(1x,a,1x,a,a,i0,a,i0)') label, side_name(side), ', nodes ', k0, '-', k+1
                    end if
                    found = .true.
                end if
                k = k + 1
            end do
        end do
    end function report_faces

    !> 고스트 셀 채우기 (매 스텝, 잔차 계산 전에 호출)
    !>   경계 면마다 고스트 ng 층을 채움. m 번째 고스트는 경계를 사이에 두고 m 번째 내부 셀과 마주봄
    !>     벽/대칭 : m 번째 내부 셀을 거울 반사 (경계 면 법선 기준)
    !>     cut     : 짝 면 건너편의 m 번째 내부 셀
    !>     출구    : 모든 층에 경계에 붙은 내부 셀 (기울기 0 외삽)
    !>     입구    : 모든 층에 지정 상태
    subroutine bc()
        integer :: ic, jc, m

        call alloc_bc()

        do m = 1,ng
            ! imin, imax 고스트 : 경계 i-면 fi=0, fi=nx 의 법선 사용
            do jc = 1,ny
                call fill(1, jc, m, si_x(0,jc),  si_y(0,jc))
                call fill(2, jc, m, si_x(nx,jc), si_y(nx,jc))
            end do

            ! jmin, jmax 고스트 : 경계 j-면 fj=0, fj=ny 의 법선 사용
            do ic = 1,nx
                call fill(3, ic, m, sj_x(ic,0),  sj_y(ic,0))
                call fill(4, ic, m, sj_x(ic,ny), sj_y(ic,ny))
            end do
        end do

    contains

        !> 경계 면 (side, k) 의 m 번째 고스트 셀 채우기
        subroutine fill(side, k, m, fnx, fny)
            integer,  intent(in) :: side, k, m
            real(dp), intent(in) :: fnx, fny
            integer :: ig, jg, ii, ji, ip, jp
            call ghost_cell(side, k, m, ig, jg)
            select case (bc_face(k,side))
            case (BC_cut)
                ! 맞닿은 짝 면 건너편의 m 번째 내부 셀 값을 그대로
                call inner_cell(cut_side(k,side), cut_k(k,side), m, ip, jp)
                Q(:,ig,jg) = Q(:,ip,jp)
            case (BC_outlet)
                call inner_cell(side, k, 1, ii, ji)
                call ghost(bc_face(k,side), W_bc(:,k,side), Q(:,ii,ji), fnx, fny, Q(:,ig,jg))
            case default
                call inner_cell(side, k, m, ii, ji)
                call ghost(bc_face(k,side), W_bc(:,k,side), Q(:,ii,ji), fnx, fny, Q(:,ig,jg))
            end select
        end subroutine fill

    end subroutine bc

    !> 경계 면 하나의 고스트 셀 값
    !>   Wb  : 입구/원방 상태 (rho, u, v, p)
    !>   qin : 경계에 붙은 내부 셀,  (fnx, fny) : 경계 면 단위 법선,  qg : 고스트 셀
    subroutine ghost(btype, Wb, qin, fnx, fny, qg)
        integer,  intent(in)  :: btype
        real(dp), intent(in)  :: Wb(nvar), qin(nvar), fnx, fny
        real(dp), intent(out) :: qg(nvar)

        real(dp) :: mn

        select case (btype)

        case (BC_inlet, BC_farfield)
            ! 지정한 상태 고정 (초음속 입구 기준, 원방은 나중에 특성 조건으로 교체)
            call p2c(Wb, qg)

        case (BC_outlet)
            ! 내부 셀 그대로 복사
            qg = qin

        case (BC_wall, BC_symmetry)
            ! 밀도, 에너지는 그대로, 운동량의 법선 성분만 반전 (비점성이라 대칭면도 같음)
            mn = qin(2)*fnx + qin(3)*fny
            qg(1) = qin(1)
            qg(2) = qin(2) - 2.0_dp*mn*fnx
            qg(3) = qin(3) - 2.0_dp*mn*fny
            qg(4) = qin(4)

        case default
            write(*,*) 'Unknown BC type:', btype
            error stop

        end select

    end subroutine ghost

end module flow