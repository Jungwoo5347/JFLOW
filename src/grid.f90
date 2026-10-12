module grid
    use euler
    implicit none
    

    ! 인덱스 이름 규칙
    !   셀   : ic, jc   (내부 1:nx, 1:ny / 고스트 ng 층: 1-ng:0, nx+1:nx+ng, 1-ng:0, ny+1:ny+ng)
    !   노드 : in, jn   (0:nx, 0:ny)
    !   면   : fi       (i-면 번호 0:nx, 셀 fi 와 fi+1 사이)
    !          fj       (j-면 번호 0:ny, 셀 fj 와 fj+1 사이)

    ! 고스트 셀 층 수: 면 하나에 양쪽 2 셀씩 (fi-1, fi | fi+1, fi+2) 쓰는 MUSCL 까지 되도록 2
    !   Q 는 Q(nvar, 1-ng:nx+ng, 1-ng:ny+ng),  모서리 고스트 (i, j 둘 다 바깥) 는 안 씀
    integer, parameter :: ng = 2

    integer  :: nx, ny                          ! 셀 개수
    real(dp) :: xmin, xmax, ymin, ymax
    real(dp), allocatable :: x(:,:), y(:,:)     ! 노드 좌표 x(in,jn) (0:nx, 0:ny)

    ! 셀 (1:nx, 1:ny) : 셀 (ic,jc)는 노드 (ic-1,jc-1), (ic,jc-1), (ic,jc), (ic-1,jc)로 둘러싸임
    real(dp), allocatable :: xc(:,:), yc(:,:)   ! 셀 중심
    real(dp), allocatable :: area(:,:)          ! 셀 면적

    ! i-면 (fi, jc) (0:nx, 1:ny) : 셀 fi 와 fi+1 사이, 법선은 +i 방향
    real(dp), allocatable :: si_x(:,:), si_y(:,:), si_len(:,:)
    ! j-면 (ic, fj) (1:nx, 0:ny) : 셀 fj 와 fj+1 사이, 법선은 +j 방향
    real(dp), allocatable :: sj_x(:,:), sj_y(:,:), sj_len(:,:)

contains

    !> 균일 직교 격자 생성. dx, dy는 여기서만 쓰고 버림
    !>   nx_in x ny_in 셀, 영역 [x0, x1] x [y0, y1]
    subroutine init_grid(nx_in, ny_in, x0, x1, y0, y1)
        integer,  intent(in) :: nx_in, ny_in
        real(dp), intent(in) :: x0, x1, y0, y1
        integer  :: in, jn
        real(dp) :: dx, dy

        nx = nx_in;  ny = ny_in
        xmin = x0;  xmax = x1
        ymin = y0;  ymax = y1

        allocate(x(0:nx, 0:ny), y(0:nx, 0:ny))

        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny
        do jn = 0, ny
            do in = 0, nx
                x(in,jn) = xmin + in*dx
                y(in,jn) = ymin + jn*dy
            end do
        end do

        call compute_metrics()
    end subroutine init_grid

    !> PLOT3D 격자 읽기 (단일 블록)
    !>   form = 1 : ASCII,  form = 0 : Fortran unformatted (바이너리, little endian)
    !>   파일 모양은 자동으로 알아냄
    !>     - 맨 앞 블록 수 줄 (nblock) 있음 / 없음
    !>     - 2D (ni nj, x y) / 3D (ni nj nk, x y z)  : 3D 면 k=1 면만 사용
    !>     - (바이너리) 단정밀도 / 배정밀도, iblank 있음 / 없음
    !>   노드 ni x nj  ->  셀 nx = ni-1, ny = nj-1   (i가 가장 빠르게 변함)
    subroutine read_plot3d(fname, form)
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: form
        integer :: ni, nj, nk

        if (form == 1) then
            call read_plot3d_ascii(fname, ni, nj, nk)
        else
            call read_plot3d_binary(fname, ni, nj, nk)
        end if

        if (nk > 1) write(*,'(a,i0,a)') ' note: 3D PLOT3D grid with nk = ', nk, ', using the k = 1 plane'

        xmin = minval(x);  xmax = maxval(x)
        ymin = minval(y);  ymax = maxval(y)

        call compute_metrics()
    end subroutine read_plot3d

    subroutine read_plot3d_ascii(fname, ni, nj, nk)
        character(len=*), intent(in)  :: fname
        integer,          intent(out) :: ni, nj, nk
        integer :: u, ios, nt, nblk, in, jn, kn, dims(3)
        character(len=1024) :: line
        real(dp), allocatable :: x3(:,:,:), y3(:,:,:), z3(:,:,:)

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

        nx = ni - 1
        ny = nj - 1
        allocate(x(0:nx, 0:ny), y(0:nx, 0:ny))

        ! x와 y(와 z)를 한 read 문으로 읽어야 함: read 문은 매번 새 줄에서 시작하는데,
        ! x가 줄 중간에서 끝나고 y가 같은 줄에서 시작하는 파일이 많음
        if (nt == 2) then
            read(u,*,iostat=ios) ((x(in,jn), in=0,nx), jn=0,ny), ((y(in,jn), in=0,nx), jn=0,ny)
        else
            allocate(x3(0:nx,0:ny,nk), y3(0:nx,0:ny,nk), z3(0:nx,0:ny,nk))
            read(u,*,iostat=ios) (((x3(in,jn,kn), in=0,nx), jn=0,ny), kn=1,nk), &
                                 (((y3(in,jn,kn), in=0,nx), jn=0,ny), kn=1,nk), &
                                 (((z3(in,jn,kn), in=0,nx), jn=0,ny), kn=1,nk)
            if (ios == 0) then
                x = x3(:,:,1);  y = y3(:,:,1)
            end if
        end if
        if (ios /= 0) call grid_fail(fname, 'file ended early or has bad data (is form-1 = 1 right?)')
        close(u)
    end subroutine read_plot3d_ascii

    !> Fortran unformatted 순차 파일: 각 레코드 앞뒤에 4바이트 길이 표시가 붙음
    !>   레코드 길이로 블록 수 줄 / 2D, 3D / 정밀도 / iblank 를 알아냄
    subroutine read_plot3d_binary(fname, ni, nj, nk)
        use, intrinsic :: iso_fortran_env, only: int32, real32, real64
        character(len=*), intent(in)  :: fname
        integer,          intent(out) :: ni, nj, nk

        integer :: u, ios, ndim, npts, nplane
        integer(int32) :: m, mend, nblk, d(3)
        integer(int32), allocatable :: iblank(:)
        real(real64),   allocatable :: b8(:)
        real(real32),   allocatable :: b4(:)
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

        ! 좌표 레코드 길이로 정밀도와 iblank 판단
        npts   = ni*nj*nk
        nplane = ni*nj
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

        nx = ni - 1
        ny = nj - 1
        allocate(x(0:nx, 0:ny), y(0:nx, 0:ny))

        ! x 전체, y 전체 (z 전체), 각각 k=1 면이 맨 앞 ni*nj 개
        if (prec == 8) then
            allocate(b8(npts))
            read(u,iostat=ios) b8;  x = reshape(b8(1:nplane), [ni, nj])
            if (ios == 0) read(u,iostat=ios) b8;  y = reshape(b8(1:nplane), [ni, nj])
            if (ios == 0 .and. ndim == 3) read(u,iostat=ios) b8
        else
            allocate(b4(npts))
            read(u,iostat=ios) b4;  x = reshape(real(b4(1:nplane), dp), [ni, nj])
            if (ios == 0) read(u,iostat=ios) b4;  y = reshape(real(b4(1:nplane), dp), [ni, nj])
            if (ios == 0 .and. ndim == 3) read(u,iostat=ios) b4
        end if
        if (ios == 0 .and. has_iblank) then
            allocate(iblank(npts))
            read(u,iostat=ios) iblank
        end if
        if (ios == 0) read(u,iostat=ios) mend
        if (ios /= 0 .or. mend /= m) call grid_fail(fname, 'file ended early or broken coordinate record')
        close(u)

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

    !> 노드 좌표 x, y만으로 셀 면적, 셀 중심, 면 법선/길이 계산
    !> (격자를 어디서 가져왔든 x, y만 채운 뒤 이걸 부르면 됨)
    subroutine compute_metrics()
        integer  :: ic, jc, fi, fj
        real(dp) :: ex, ey, d1x, d1y, d2x, d2y

        allocate(xc(nx,ny), yc(nx,ny), area(nx,ny))
        allocate(si_x(0:nx,ny), si_y(0:nx,ny), si_len(0:nx,ny))
        allocate(sj_x(nx,0:ny), sj_y(nx,0:ny), sj_len(nx,0:ny))

        ! 셀 면적 (두 대각선 외적의 절반), 셀 중심 (네 꼭짓점 평균)
        !   셀 (ic,jc)의 꼭짓점 노드: (ic-1,jc-1), (ic,jc-1), (ic,jc), (ic-1,jc)
        do jc = 1, ny
            do ic = 1, nx
                d1x = x(ic,jc)   - x(ic-1,jc-1);  d1y = y(ic,jc)   - y(ic-1,jc-1)
                d2x = x(ic-1,jc) - x(ic,jc-1);    d2y = y(ic-1,jc) - y(ic,jc-1)
                area(ic,jc) = 0.5_dp * (d1x*d2y - d1y*d2x)

                xc(ic,jc) = 0.25_dp * (x(ic-1,jc-1) + x(ic,jc-1) + x(ic,jc) + x(ic-1,jc))
                yc(ic,jc) = 0.25_dp * (y(ic-1,jc-1) + y(ic,jc-1) + y(ic,jc) + y(ic-1,jc))
            end do
        end do

        ! i-면 (fi,jc): 노드 (fi,jc-1) -> (fi,jc) 변 벡터 (ex,ey), 법선 = (ey, -ex)
        do jc = 1, ny
            do fi = 0, nx
                ex = x(fi,jc) - x(fi,jc-1)
                ey = y(fi,jc) - y(fi,jc-1)
                si_len(fi,jc) = sqrt(ex**2 + ey**2)
                si_x(fi,jc)   =  ey / si_len(fi,jc)
                si_y(fi,jc)   = -ex / si_len(fi,jc)
            end do
        end do

        ! j-면 (ic,fj): 노드 (ic-1,fj) -> (ic,fj) 변 벡터 (ex,ey), 법선 = (-ey, ex)
        do fj = 0, ny
            do ic = 1, nx
                ex = x(ic,fj) - x(ic-1,fj)
                ey = y(ic,fj) - y(ic-1,fj)
                sj_len(ic,fj) = sqrt(ex**2 + ey**2)
                sj_x(ic,fj)   = -ey / sj_len(ic,fj)
                sj_y(ic,fj)   =  ex / sj_len(ic,fj)
            end do
        end do
    end subroutine compute_metrics

    !> 격자 점검: 면적 부호, 전체 면적, 셀마다 면 벡터 합(= 0이어야 함)
    subroutine check_grid()
        integer  :: ic, jc
        real(dp) :: cx, cy, max_closure

        if (any(area <= 0.0_dp)) then
            write(*,*) 'Error: non-positive cell area, min =', minval(area)
            error stop
        end if

        ! 셀 (ic,jc)의 네 면: 오른쪽 i-면 ic, 왼쪽 i-면 ic-1, 위 j-면 jc, 아래 j-면 jc-1
        max_closure = 0.0_dp
        do jc = 1, ny
            do ic = 1, nx
                cx = si_x(ic,jc)*si_len(ic,jc) - si_x(ic-1,jc)*si_len(ic-1,jc) &
                   + sj_x(ic,jc)*sj_len(ic,jc) - sj_x(ic,jc-1)*sj_len(ic,jc-1)
                cy = si_y(ic,jc)*si_len(ic,jc) - si_y(ic-1,jc)*si_len(ic-1,jc) &
                   + sj_y(ic,jc)*sj_len(ic,jc) - sj_y(ic,jc-1)*sj_len(ic,jc-1)
                max_closure = max(max_closure, abs(cx), abs(cy))
            end do
        end do

        write(*,'(a,i0,a,i0)')   ' grid       : ', nx, ' x ', ny
        write(*,'(a,es12.4)')    ' total area : ', sum(area)
        write(*,'(a,2es12.4)')   ' min/max A  : ', minval(area), maxval(area)
        write(*,'(a,es12.4)')    ' closure    : ', max_closure
    end subroutine check_grid

end module grid
