module grid
    use euler
    implicit none


    ! 인덱스 이름 규칙
    !   셀   : ic, jc, kc   (내부 1:nx, 1:ny, 1:nz / 고스트 ng 층)
    !   노드 : in, jn, kn   (0:nx, 0:ny, 0:nz)
    !   면   : fi  (i-면 번호 0:nx, 셀 fi 와 fi+1 사이),  fj, fk 도 같음

    ! 고스트 셀 층 수: MUSCL 이 면 하나에 양쪽 2 셀씩 쓰므로 2
    !   Q 는 Q(nvar, 1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng),  모서리 고스트 (두 방향 이상 바깥) 는 안 씀
    integer, parameter :: ng = 2

    integer  :: nx, ny, nz                              ! 셀 개수
    real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(dp), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)   ! 노드 좌표 x(in,jn,kn) (0:nx, 0:ny, 0:nz)

    ! 셀 (1:nx, 1:ny, 1:nz) : 셀 (ic,jc,kc) 는 노드 (ic-1:ic, jc-1:jc, kc-1:kc) 8 개로 둘러싸인 육면체
    real(dp), allocatable :: xc(:,:,:), yc(:,:,:), zc(:,:,:)  ! 셀 중심 (8 꼭짓점 평균)
    real(dp), allocatable :: vol(:,:,:)                       ! 셀 부피

    ! 면 법선은 단위 벡터 (3, ...) 와 넓이를 따로 저장
    ! i-면 (fi,jc,kc) (0:nx, 1:ny, 1:nz) : 셀 fi 와 fi+1 사이, 법선은 +i 방향
    real(dp), allocatable :: si(:,:,:,:), si_a(:,:,:)
    ! j-면 (ic,fj,kc) (1:nx, 0:ny, 1:nz) : 법선은 +j 방향
    real(dp), allocatable :: sj(:,:,:,:), sj_a(:,:,:)
    ! k-면 (ic,jc,fk) (1:nx, 1:ny, 0:nz) : 법선은 +k 방향
    real(dp), allocatable :: sk(:,:,:,:), sk_a(:,:,:)

contains

    !> 균일 직교 격자 생성
    !>   nx_in x ny_in x nz_in 셀, 영역 [x0, x1] x [y0, y1] x [z0, z1]
    subroutine init_grid(nx_in, ny_in, nz_in, x0, x1, y0, y1, z0, z1)
        integer,  intent(in) :: nx_in, ny_in, nz_in
        real(dp), intent(in) :: x0, x1, y0, y1, z0, z1
        integer :: in, jn, kn

        nx = nx_in;  ny = ny_in;  nz = nz_in
        allocate(x(0:nx,0:ny,0:nz), y(0:nx,0:ny,0:nz), z(0:nx,0:ny,0:nz))
        do kn = 0,nz
            do jn = 0,ny
                do in = 0,nx
                    x(in,jn,kn) = x0 + (x1 - x0)*in/nx
                    y(in,jn,kn) = y0 + (y1 - y0)*jn/ny
                    z(in,jn,kn) = z0 + (z1 - z0)*kn/nz
                end do
            end do
        end do

        call compute_metrics()
    end subroutine init_grid

    !> PLOT3D 격자 읽기 (단일 블록)
    !>   form = 1 : ASCII,  form = 0 : Fortran unformatted (바이너리, little endian)
    !>   파일 모양은 자동으로 알아냄 (블록 수 줄 유무, 2D/3D, 정밀도, iblank)
    !>   3D 파일 : 그대로.   2D 파일 : ext_nk >= 2 이면 z = 0 ~ ext_span 으로 밀어냄 (아니면 에러)
    subroutine read_plot3d(fname, form, ext_nk, ext_span)
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: form, ext_nk
        real(dp),         intent(in) :: ext_span
        real(dp), allocatable :: x2(:,:), y2(:,:)
        integer :: kn

        if (form == 1) then
            call read_plot3d_ascii(fname, x2, y2)
        else
            call read_plot3d_binary(fname, x2, y2)
        end if

        if (allocated(x2)) then
            ! 2D 파일
            if (ext_nk < 2) call grid_fail(fname, '2D PLOT3D grid: add  extrude-1 = nk span  to grid.in to make it 3D')
            nx = size(x2,1) - 1;  ny = size(x2,2) - 1;  nz = ext_nk - 1
            allocate(x(0:nx,0:ny,0:nz), y(0:nx,0:ny,0:nz), z(0:nx,0:ny,0:nz))
            do kn = 0,nz
                x(:,:,kn) = x2
                y(:,:,kn) = y2
                z(:,:,kn) = ext_span*kn/nz
            end do
            write(*,'(a,i0,a)') ' note: 2D PLOT3D grid extruded in z to nk = ', ext_nk, ' nodes'
        else if (ext_nk > 0) then
            call grid_fail(fname, 'extrude-1 is only for 2D grid files, this one is 3D')
        end if

        call compute_metrics()
    end subroutine read_plot3d

    !> ASCII PLOT3D.  3D 면 모듈 배열 x, y, z 를 채우고, 2D 면 x2, y2 에 담아 돌려줌
    subroutine read_plot3d_ascii(fname, x2, y2)
        character(len=*), intent(in) :: fname
        real(dp), allocatable, intent(out) :: x2(:,:), y2(:,:)
        integer :: u, ios, nt, nblk, in, jn, kn, dims(3), ni, nj, nk
        character(len=1024) :: line

        open(newunit=u, file=fname, status='old', action='read', iostat=ios)
        if (ios /= 0) call grid_fail(fname, 'cannot open grid file')

        ! 첫 줄 숫자 개수: 1 = 블록 수 줄, 2 = ni nj (2D), 3 = ni nj nk (3D)
        read(u,'(a)',iostat=ios) line
        if (ios /= 0) call grid_fail(fname, 'empty file')
        nt = count_tokens(line)
        if (nt == 1) then
            read(line,*) nblk
            if (nblk /= 1) call grid_fail(fname, 'multi-block grid (nblock > 1) is not supported yet')
            read(u,'(a)',iostat=ios) line
            if (ios /= 0) call grid_fail(fname, 'cannot read block size line')
            nt = count_tokens(line)
        end if
        if (nt /= 2 .and. nt /= 3) call grid_fail(fname, 'cannot read grid size (expected ni nj or ni nj nk)')
        dims = 1
        read(line,*,iostat=ios) dims(1:nt)
        if (ios /= 0) call grid_fail(fname, 'cannot read grid size (expected ni nj or ni nj nk)')
        ni = dims(1);  nj = dims(2);  nk = dims(3)
        if (ni < 2 .or. nj < 2 .or. nk < 1) call grid_fail(fname, 'bad grid size')
        if (nt == 3 .and. nk < 2) call grid_fail(fname, '3D grid with nk = 1: write it as a 2D file and use extrude-1')

        ! x, y (, z) 를 한 read 문으로 (x 가 줄 중간에서 끝나고 y 가 같은 줄에서 시작하는 파일이 많음)
        if (nt == 2) then
            allocate(x2(0:ni-1,0:nj-1), y2(0:ni-1,0:nj-1))
            read(u,*,iostat=ios) ((x2(in,jn), in=0,ni-1), jn=0,nj-1), ((y2(in,jn), in=0,ni-1), jn=0,nj-1)
        else
            nx = ni - 1;  ny = nj - 1;  nz = nk - 1
            allocate(x(0:nx,0:ny,0:nz), y(0:nx,0:ny,0:nz), z(0:nx,0:ny,0:nz))
            read(u,*,iostat=ios) (((x(in,jn,kn), in=0,nx), jn=0,ny), kn=0,nz), &
                                 (((y(in,jn,kn), in=0,nx), jn=0,ny), kn=0,nz), &
                                 (((z(in,jn,kn), in=0,nx), jn=0,ny), kn=0,nz)
        end if
        if (ios /= 0) call grid_fail(fname, 'file ended early or has bad data (is form-1 = 1 right?)')
        close(u)
    end subroutine read_plot3d_ascii

    !> Fortran unformatted 순차 파일: 각 레코드 앞뒤에 4바이트 길이 표시가 붙음
    !>   레코드 길이로 블록 수 줄 / 2D, 3D / 정밀도 / iblank 를 알아냄
    subroutine read_plot3d_binary(fname, x2, y2)
        use, intrinsic :: iso_fortran_env, only: int32, real32, real64
        character(len=*), intent(in) :: fname
        real(dp), allocatable, intent(out) :: x2(:,:), y2(:,:)

        integer :: u, ios, ndim, npts, ni, nj, nk, c
        integer(int32) :: m, mend, nblk, d(3)
        integer(int32), allocatable :: iblank(:)
        real(real64),   allocatable :: b8(:)
        real(real32),   allocatable :: b4(:)
        real(dp),       allocatable :: crd(:,:)          ! crd(:, c) : c = 1 x, 2 y, 3 z 전체
        integer :: prec
        logical :: has_iblank

        open(newunit=u, file=fname, access='stream', form='unformatted', status='old', action='read', iostat=ios)
        if (ios /= 0) call grid_fail(fname, 'cannot open grid file')

        ! 첫 레코드: 4바이트 = 블록 수,  8바이트 = ni nj,  12바이트 = ni nj nk
        read(u,iostat=ios) m
        if (ios /= 0) call grid_fail(fname, 'empty file')
        if (m == 4) then
            read(u,iostat=ios) nblk, mend
            if (ios /= 0 .or. mend /= 4) call grid_fail(fname, 'broken record (block count)')
            if (nblk /= 1) call grid_fail(fname, 'multi-block grid (nblock > 1) is not supported yet')
            read(u,iostat=ios) m
            if (ios /= 0) call grid_fail(fname, 'cannot read block size record')
        end if
        if (m /= 8 .and. m /= 12) then
            if (any(byteswap(m) == [4, 8, 12])) then
                call grid_fail(fname, 'looks like a big-endian file; please export little-endian')
            end if
            call grid_fail(fname, 'not a Fortran unformatted PLOT3D file (is form-1 = 0 right? ascii is 1)')
        end if
        ndim = m/4
        d = 1
        read(u,iostat=ios) d(1:ndim), mend
        if (ios /= 0 .or. mend /= m) call grid_fail(fname, 'broken record (grid size)')
        ni = d(1);  nj = d(2);  nk = d(3)
        if (ni < 2 .or. nj < 2 .or. nk < 1) call grid_fail(fname, 'bad grid size')
        if (ndim == 3 .and. nk < 2) call grid_fail(fname, '3D grid with nk = 1: write it as a 2D file and use extrude-1')

        ! 좌표 레코드 길이로 정밀도와 iblank 판단
        npts = ni*nj*nk
        read(u,iostat=ios) m
        if (ios /= 0) call grid_fail(fname, 'cannot read coordinate record')
        if (m == ndim*npts*8) then
            prec = 8;  has_iblank = .false.
        else if (m == ndim*npts*4) then
            prec = 4;  has_iblank = .false.
        else if (m == ndim*npts*8 + npts*4) then
            prec = 8;  has_iblank = .true.
        else if (m == ndim*npts*4 + npts*4) then
            prec = 4;  has_iblank = .true.
        else
            call grid_fail(fname, 'coordinate record size does not match ni*nj(*nk)')
        end if

        ! x 전체, y 전체 (, z 전체)
        allocate(crd(npts, ndim))
        if (prec == 8) then
            allocate(b8(npts))
            do c = 1,ndim
                if (ios == 0) read(u,iostat=ios) b8
                crd(:,c) = b8
            end do
        else
            allocate(b4(npts))
            do c = 1,ndim
                if (ios == 0) read(u,iostat=ios) b4
                crd(:,c) = real(b4, dp)
            end do
        end if
        if (ios == 0 .and. has_iblank) then
            allocate(iblank(npts))
            read(u,iostat=ios) iblank
        end if
        if (ios == 0) read(u,iostat=ios) mend
        if (ios /= 0 .or. mend /= m) call grid_fail(fname, 'file ended early or broken coordinate record')
        close(u)

        if (ndim == 2) then
            allocate(x2(0:ni-1,0:nj-1), y2(0:ni-1,0:nj-1))
            x2 = reshape(crd(:,1), [ni, nj])
            y2 = reshape(crd(:,2), [ni, nj])
        else
            nx = ni - 1;  ny = nj - 1;  nz = nk - 1
            allocate(x(0:nx,0:ny,0:nz), y(0:nx,0:ny,0:nz), z(0:nx,0:ny,0:nz))
            x = reshape(crd(:,1), [ni, nj, nk])
            y = reshape(crd(:,2), [ni, nj, nk])
            z = reshape(crd(:,3), [ni, nj, nk])
        end if

    contains

        elemental integer(int32) function byteswap(v)
            integer(int32), intent(in) :: v
            integer(int32) :: k
            byteswap = 0
            do k = 0,3
                call mvbits(v, 8*k, 8, byteswap, 8*(3-k))
            end do
        end function byteswap

    end subroutine read_plot3d_binary

    integer function count_tokens(s)
        character(len=*), intent(in) :: s
        integer :: k
        logical :: in_tok
        count_tokens = 0;  in_tok = .false.
        do k = 1,len_trim(s)
            if (s(k:k) == ' ' .or. s(k:k) == ',' .or. s(k:k) == achar(9)) then
                in_tok = .false.
            else if (.not. in_tok) then
                in_tok = .true.
                count_tokens = count_tokens + 1
            end if
        end do
    end function count_tokens

    subroutine grid_fail(fname, text)
        character(len=*), intent(in) :: fname, text
        write(*,'(a,a,a,a)') ' Error: ', trim(fname), ': ', text
        error stop
    end subroutine grid_fail

    !> 노드 (in,jn,kn) 좌표
    pure function node(in, jn, kn) result(p)
        integer, intent(in) :: in, jn, kn
        real(dp) :: p(3)
        p = [x(in,jn,kn), y(in,jn,kn), z(in,jn,kn)]
    end function node

    !> 사각형 면 (꼭짓점 a -> b -> c -> d 를 한 바퀴 도는 순서) 의 면 벡터 (법선 x 넓이)
    !>   = 두 대각선 외적의 절반. 네 점이 한 평면에 없어도 (bilinear 면) 정확하고,
    !>   닫힌 셀의 여섯 면 벡터 합이 정확히 0 이 됨 (오른손 법칙: a->b->c->d 를 감는 방향)
    pure function face_vec(a, b, c, d) result(s)
        real(dp), intent(in) :: a(3), b(3), c(3), d(3)
        real(dp) :: s(3)
        s = 0.5_dp*cross(c - a, d - b)
    end function face_vec

    !> i-면 fi, j-면 fj, k-면 fk 의 면 벡터 (+i, +j, +k 방향, 오른손 격자일 때)
    pure function svec_i(fi, jc, kc) result(s)
        integer, intent(in) :: fi, jc, kc
        real(dp) :: s(3)
        s = face_vec(node(fi,jc-1,kc-1), node(fi,jc,kc-1), node(fi,jc,kc), node(fi,jc-1,kc))
    end function svec_i

    pure function svec_j(ic, fj, kc) result(s)
        integer, intent(in) :: ic, fj, kc
        real(dp) :: s(3)
        s = face_vec(node(ic-1,fj,kc-1), node(ic-1,fj,kc), node(ic,fj,kc), node(ic,fj,kc-1))
    end function svec_j

    pure function svec_k(ic, jc, fk) result(s)
        integer, intent(in) :: ic, jc, fk
        real(dp) :: s(3)
        s = face_vec(node(ic-1,jc-1,fk), node(ic,jc-1,fk), node(ic,jc,fk), node(ic-1,jc,fk))
    end function svec_k

    !> 면 중심 (네 꼭짓점 평균)
    pure function fcen_i(fi, jc, kc) result(p)
        integer, intent(in) :: fi, jc, kc
        real(dp) :: p(3)
        p = 0.25_dp*(node(fi,jc-1,kc-1) + node(fi,jc,kc-1) + node(fi,jc,kc) + node(fi,jc-1,kc))
    end function fcen_i

    pure function fcen_j(ic, fj, kc) result(p)
        integer, intent(in) :: ic, fj, kc
        real(dp) :: p(3)
        p = 0.25_dp*(node(ic-1,fj,kc-1) + node(ic-1,fj,kc) + node(ic,fj,kc) + node(ic,fj,kc-1))
    end function fcen_j

    pure function fcen_k(ic, jc, fk) result(p)
        integer, intent(in) :: ic, jc, fk
        real(dp) :: p(3)
        p = 0.25_dp*(node(ic-1,jc-1,fk) + node(ic,jc-1,fk) + node(ic,jc,fk) + node(ic-1,jc,fk))
    end function fcen_k

    !> 노드 좌표만으로 셀 부피, 셀 중심, 면 법선/넓이 계산
    !>   부피 : 발산 정리  V = 1/3 * sum_면 (면 중심 - 셀 중심) . (바깥 면 벡터)
    subroutine compute_metrics()
        integer  :: ic, jc, kc, fi, fj, fk
        real(dp) :: s(3), c(3)

        xmin = minval(x);  xmax = maxval(x)
        ymin = minval(y);  ymax = maxval(y)
        zmin = minval(z);  zmax = maxval(z)

        allocate(xc(nx,ny,nz), yc(nx,ny,nz), zc(nx,ny,nz), vol(nx,ny,nz))
        allocate(si(3,0:nx,ny,nz), si_a(0:nx,ny,nz))
        allocate(sj(3,nx,0:ny,nz), sj_a(nx,0:ny,nz))
        allocate(sk(3,nx,ny,0:nz), sk_a(nx,ny,0:nz))

        !$omp parallel private(ic, jc, kc, fi, fj, fk, s, c)
        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do fi = 0,nx
                    s = svec_i(fi, jc, kc)
                    si_a(fi,jc,kc)   = norm2(s)
                    si(:,fi,jc,kc)   = s/max(si_a(fi,jc,kc), tiny(1.0_dp))    ! 넓이 0 (특이선) 면은 법선 0
                end do
            end do
        end do
        !$omp end do nowait
        !$omp do collapse(2)
        do kc = 1,nz
            do fj = 0,ny
                do ic = 1,nx
                    s = svec_j(ic, fj, kc)
                    sj_a(ic,fj,kc)   = norm2(s)
                    sj(:,ic,fj,kc)   = s/max(sj_a(ic,fj,kc), tiny(1.0_dp))    ! 넓이 0 (특이선) 면은 법선 0
                end do
            end do
        end do
        !$omp end do nowait
        !$omp do collapse(2)
        do fk = 0,nz
            do jc = 1,ny
                do ic = 1,nx
                    s = svec_k(ic, jc, fk)
                    sk_a(ic,jc,fk)   = norm2(s)
                    sk(:,ic,jc,fk)   = s/max(sk_a(ic,jc,fk), tiny(1.0_dp))    ! 넓이 0 (특이선) 면은 법선 0
                end do
            end do
        end do
        !$omp end do nowait

        !$omp do collapse(2)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    c = 0.125_dp*( node(ic-1,jc-1,kc-1) + node(ic,jc-1,kc-1) + node(ic,jc,kc-1) + node(ic-1,jc,kc-1) &
                                 + node(ic-1,jc-1,kc  ) + node(ic,jc-1,kc  ) + node(ic,jc,kc  ) + node(ic-1,jc,kc  ) )
                    xc(ic,jc,kc) = c(1);  yc(ic,jc,kc) = c(2);  zc(ic,jc,kc) = c(3)
                    vol(ic,jc,kc) = ( dot_product(fcen_i(ic,  jc,kc) - c,  svec_i(ic,  jc,kc)) &
                                    - dot_product(fcen_i(ic-1,jc,kc) - c,  svec_i(ic-1,jc,kc)) &
                                    + dot_product(fcen_j(ic,jc,  kc) - c,  svec_j(ic,jc,  kc)) &
                                    - dot_product(fcen_j(ic,jc-1,kc) - c,  svec_j(ic,jc-1,kc)) &
                                    + dot_product(fcen_k(ic,jc,kc  ) - c,  svec_k(ic,jc,kc  )) &
                                    - dot_product(fcen_k(ic,jc,kc-1) - c,  svec_k(ic,jc,kc-1)) ) / 3.0_dp
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine compute_metrics

    !> 격자 점검: 부피 부호, 전체 부피, 셀마다 면 벡터 합(= 0이어야 함)
    subroutine check_grid()
        integer  :: ic, jc, kc
        real(dp) :: cs(3), max_closure

        if (any(vol <= 0.0_dp)) then
            write(*,'(a,es12.4,a,i0,a,i0)') ' Error: non-positive cell volume, min =', minval(vol), &
                ',  cells with V <= 0: ', count(vol <= 0.0_dp), ' of ', size(vol)
            if (all(vol < 0.0_dp)) write(*,'(a)') '   the grid is left-handed: reverse one index direction (i, j or k)'
            error stop
        end if

        ! 셀의 여섯 면: +i 면 ic, -i 면 ic-1, ...
        max_closure = 0.0_dp
        !$omp parallel do collapse(2) private(ic, cs) reduction(max:max_closure)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    cs = si(:,ic,jc,kc)*si_a(ic,jc,kc) - si(:,ic-1,jc,kc)*si_a(ic-1,jc,kc) &
                       + sj(:,ic,jc,kc)*sj_a(ic,jc,kc) - sj(:,ic,jc-1,kc)*sj_a(ic,jc-1,kc) &
                       + sk(:,ic,jc,kc)*sk_a(ic,jc,kc) - sk(:,ic,jc,kc-1)*sk_a(ic,jc,kc-1)
                    max_closure = max(max_closure, maxval(abs(cs)))
                end do
            end do
        end do
        !$omp end parallel do

        write(*,'(a,i0,a,i0,a,i0,a,i0,a)') ' grid       : ', nx, ' x ', ny, ' x ', nz, '  (', nx*ny*nz, ' cells)'
        write(*,'(a,6es11.3)')   ' bounds     : ', xmin, xmax, ymin, ymax, zmin, zmax
        write(*,'(a,es12.4)')    ' total vol  : ', sum(vol)
        write(*,'(a,2es12.4)')   ' min/max V  : ', minval(vol), maxval(vol)
        write(*,'(a,es12.4)')    ' closure    : ', max_closure
    end subroutine check_grid

end module grid
