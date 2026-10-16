module grid
    use euler
    implicit none


    ! 인덱스 이름 규칙 (블록마다)
    !   셀   : ic, jc, kc   (내부 1:nx, 1:ny, 1:nz / 고스트 ng 층)
    !   노드 : in, jn, kn   (0:nx, 0:ny, 0:nz)
    !   면   : fi  (i-면 번호 0:nx, 셀 fi 와 fi+1 사이),  fj, fk 도 같음

    ! 고스트 셀 층 수: MUSCL 이 면 하나에 양쪽 2 셀씩 쓰므로 2
    !   Q 는 Q(nvar, 1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng),  모서리 고스트 (두 방향 이상 바깥) 는 안 씀
    integer, parameter :: ng = 2

    ! ---- 블록 하나 (격자 + 유동 배열) ----
    !   셀 (ic,jc,kc) 는 노드 (ic-1:ic, jc-1:jc, kc-1:kc) 8 개로 둘러싸인 육면체
    !   면 법선은 단위 벡터 (3, ...) 와 넓이를 따로 저장
    !     i-면 (fi,jc,kc) (0:nx, 1:ny, 1:nz) : 셀 fi 와 fi+1 사이, 법선은 +i 방향  (j, k 도 같음)
    type block_t
        integer :: nx = 0, ny = 0, nz = 0                       ! 셀 개수
        real(dp), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)    ! 노드 좌표 (0:nx, 0:ny, 0:nz)
        real(dp), allocatable :: xc(:,:,:), yc(:,:,:), zc(:,:,:) ! 셀 중심 (8 꼭짓점 평균)
        real(dp), allocatable :: vol(:,:,:)                      ! 셀 부피
        real(dp), allocatable :: si(:,:,:,:), si_a(:,:,:)
        real(dp), allocatable :: sj(:,:,:,:), sj_a(:,:,:)
        real(dp), allocatable :: sk(:,:,:,:), sk_a(:,:,:)
        ! 유동 (flow, solver, main 에서 채움)
        real(dp), allocatable :: Q(:,:,:,:)                      ! 보존 변수 (고스트 포함)
        real(dp), allocatable :: Q0(:,:,:,:)                     ! RK 시작 값 (내부)
        real(dp), allocatable :: W(:,:,:,:)                      ! 원시 변수 (고스트 포함, residual 에서 계산)
        real(dp), allocatable :: R(:,:,:,:), Rt(:,:,:,:)         ! 잔차 (RK 첫 단계 / 나머지 단계)
        real(dp), allocatable :: dtl(:,:,:)                      ! 셀마다의 dt
    end type

    integer :: nblk = 0
    type(block_t), allocatable, target :: blk(:)

    ! ---- 현재 블록 (use_block 으로 고름): 셀/면 루프는 이 이름들을 그대로 씀 ----
    integer :: nx, ny, nz
    real(dp), pointer, contiguous :: x(:,:,:) => null(), y(:,:,:) => null(), z(:,:,:) => null()
    real(dp), pointer, contiguous :: xc(:,:,:) => null(), yc(:,:,:) => null(), zc(:,:,:) => null()
    real(dp), pointer, contiguous :: vol(:,:,:) => null()
    real(dp), pointer, contiguous :: si(:,:,:,:) => null(), si_a(:,:,:) => null()
    real(dp), pointer, contiguous :: sj(:,:,:,:) => null(), sj_a(:,:,:) => null()
    real(dp), pointer, contiguous :: sk(:,:,:,:) => null(), sk_a(:,:,:) => null()

    ! PLOT3D 파일에서 읽은 블록 하나 (2D 면 nk = 1, z 없음)
    type rawblk_t
        integer :: ni, nj, nk, ndim
        real(dp), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
    end type

contains

    !> 블록 b 를 현재 블록으로 (nx, ny, nz 와 격자 포인터)
    subroutine use_block(b)
        integer, intent(in) :: b
        nx = blk(b)%nx;  ny = blk(b)%ny;  nz = blk(b)%nz
        x  => blk(b)%x;   y  => blk(b)%y;   z  => blk(b)%z
        if (allocated(blk(b)%vol)) then
            xc => blk(b)%xc;  yc => blk(b)%yc;  zc => blk(b)%zc;  vol => blk(b)%vol
            si => blk(b)%si;  si_a => blk(b)%si_a
            sj => blk(b)%sj;  sj_a => blk(b)%sj_a
            sk => blk(b)%sk;  sk_a => blk(b)%sk_a
        end if
    end subroutine use_block

    !> grid.in 의 block-n 줄들로 블록 만들기 (격자 읽기 + 기하 계산)
    !>   block 줄이 하나뿐이고 파일에 블록이 여러 개면 전부 사용, 아니면 줄 하나 = 블록 하나
    subroutine build_blocks(specs, nspecs, nblock_in)
        use case_setup, only: gridspec_t
        type(gridspec_t), intent(in) :: specs(:)
        integer,          intent(in) :: nspecs, nblock_in
        type(rawblk_t), allocatable :: rb(:)
        integer :: s, b, m, nfile

        ! 블록 수 세기
        nblk = 0
        do s = 1,nspecs
            if (specs(s)%kind == 'cartesian') then
                nblk = nblk + 1
            else
                nfile = plot3d_nblocks(trim(specs(s)%file), specs(s)%form)
                if (specs(s)%index > nfile) then
                    write(*,'(a,a,a,i0,a,i0)') ' Error: ', trim(specs(s)%file), ' has ', nfile, &
                        ' block(s), block-line asks for block ', specs(s)%index
                    error stop
                end if
                if (specs(s)%index > 0 .or. nfile == 1) then
                    nblk = nblk + 1
                else if (nspecs == 1) then
                    nblk = nblk + nfile
                else
                    write(*,'(a,a,a,i0,a)') ' Error: ', trim(specs(s)%file), ' has ', nfile, &
                        ' blocks: write  block-n = file m  for each block (or use a single block-1 line)'
                    error stop
                end if
            end if
        end do
        if (nblock_in > 0 .and. nblock_in /= nblk) then
            write(*,'(a,i0,a,i0,a)') ' Error: grid.in nblock = ', nblock_in, ' but the grid has ', nblk, ' block(s)'
            error stop
        end if
        allocate(blk(nblk))

        ! 읽기
        b = 0
        do s = 1,nspecs
            associate (sp => specs(s))
                if (sp%kind == 'cartesian') then
                    b = b + 1
                    call cartesian_block(b, sp%ni-1, sp%nj-1, sp%nk-1, sp%x0, sp%x1, sp%y0, sp%y1, sp%z0, sp%z1)
                else
                    call read_plot3d(trim(sp%file), sp%form, rb)
                    do m = 1,size(rb)
                        if (sp%index > 0 .and. m /= sp%index) cycle
                        b = b + 1
                        call take_block(b, rb(m), sp%extrude_nk, sp%extrude_span, trim(sp%file))
                    end do
                    deallocate(rb)
                end if
            end associate
        end do

        do b = 1,nblk
            call compute_metrics(b)
        end do
    end subroutine build_blocks

    !> 균일 직교 격자 블록:  nx_in x ny_in x nz_in 셀, 영역 [x0, x1] x [y0, y1] x [z0, z1]
    subroutine cartesian_block(b, nx_in, ny_in, nz_in, x0, x1, y0, y1, z0, z1)
        integer,  intent(in) :: b, nx_in, ny_in, nz_in
        real(dp), intent(in) :: x0, x1, y0, y1, z0, z1
        integer :: in, jn, kn

        associate (bk => blk(b))
            bk%nx = nx_in;  bk%ny = ny_in;  bk%nz = nz_in
            allocate(bk%x(0:nx_in,0:ny_in,0:nz_in), bk%y(0:nx_in,0:ny_in,0:nz_in), bk%z(0:nx_in,0:ny_in,0:nz_in))
            do kn = 0,nz_in
                do jn = 0,ny_in
                    do in = 0,nx_in
                        bk%x(in,jn,kn) = x0 + (x1 - x0)*in/nx_in
                        bk%y(in,jn,kn) = y0 + (y1 - y0)*jn/ny_in
                        bk%z(in,jn,kn) = z0 + (z1 - z0)*kn/nz_in
                    end do
                end do
            end do
        end associate
    end subroutine cartesian_block

    !> 파일에서 읽은 블록을 블록 b 로: 3D 는 그대로, 2D 는 ext_nk >= 2 이면 z = 0 ~ ext_span 으로 밀어냄
    subroutine take_block(b, r, ext_nk, ext_span, fname)
        integer,          intent(in)    :: b, ext_nk
        type(rawblk_t),   intent(inout) :: r
        real(dp),         intent(in)    :: ext_span
        character(len=*), intent(in)    :: fname
        integer :: kn

        associate (bk => blk(b))
            if (r%ndim == 2) then
                if (ext_nk < 2) call grid_fail(fname, '2D PLOT3D grid: add  extrude-n = nk span  to grid.in to make it 3D')
                bk%nx = r%ni - 1;  bk%ny = r%nj - 1;  bk%nz = ext_nk - 1
                allocate(bk%x(0:bk%nx,0:bk%ny,0:bk%nz), bk%y(0:bk%nx,0:bk%ny,0:bk%nz), bk%z(0:bk%nx,0:bk%ny,0:bk%nz))
                do kn = 0,bk%nz
                    bk%x(:,:,kn) = r%x(:,:,0)
                    bk%y(:,:,kn) = r%y(:,:,0)
                    bk%z(:,:,kn) = ext_span*kn/bk%nz
                end do
                write(*,'(a,i0,a,i0,a)') ' note: block ', b, ': 2D PLOT3D grid extruded in z to nk = ', ext_nk, ' nodes'
            else
                if (ext_nk > 0) call grid_fail(fname, 'extrude-n is only for 2D grid files, this one is 3D')
                bk%nx = r%ni - 1;  bk%ny = r%nj - 1;  bk%nz = r%nk - 1
                call move_alloc(r%x, bk%x);  call move_alloc(r%y, bk%y);  call move_alloc(r%z, bk%z)
            end if
        end associate
    end subroutine take_block

    !> PLOT3D 파일의 블록 수 (머리만 읽음)
    integer function plot3d_nblocks(fname, form) result(nb)
        use, intrinsic :: iso_fortran_env, only: int32
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: form
        integer :: u, ios
        integer(int32) :: m, nb32
        character(len=1024) :: line

        if (form == 1) then
            open(newunit=u, file=fname, status='old', action='read', iostat=ios)
            if (ios /= 0) call grid_fail(fname, 'cannot open grid file')
            read(u,'(a)',iostat=ios) line
            if (ios /= 0) call grid_fail(fname, 'empty file')
            if (count_tokens(line) == 1) then
                read(line,*) nb
            else
                nb = 1
            end if
        else
            open(newunit=u, file=fname, access='stream', form='unformatted', status='old', action='read', iostat=ios)
            if (ios /= 0) call grid_fail(fname, 'cannot open grid file')
            read(u,iostat=ios) m
            if (ios /= 0) call grid_fail(fname, 'empty file')
            if (m == 4) then
                read(u,iostat=ios) nb32
                nb = nb32
            else
                nb = 1
            end if
        end if
        close(u)
        if (nb < 1) call grid_fail(fname, 'bad block count')
    end function plot3d_nblocks

    !> PLOT3D 격자 읽기 (단일/다중 블록)
    !>   form = 1 : ASCII,  form = 0 : Fortran unformatted (바이너리, little endian)
    !>   파일 모양은 자동으로 알아냄 (블록 수 줄 유무, 2D/3D, 정밀도, iblank)
    subroutine read_plot3d(fname, form, rb)
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: form
        type(rawblk_t), allocatable, intent(out) :: rb(:)
        if (form == 1) then
            call read_plot3d_ascii(fname, rb)
        else
            call read_plot3d_binary(fname, rb)
        end if
    end subroutine read_plot3d

    subroutine read_plot3d_ascii(fname, rb)
        character(len=*), intent(in) :: fname
        type(rawblk_t), allocatable, intent(out) :: rb(:)
        integer :: u, ios, nt, nb, ndim, m, in, jn, kn
        integer, allocatable :: dims(:,:)
        logical :: has_nb
        character(len=1024) :: line

        open(newunit=u, file=fname, status='old', action='read', iostat=ios)
        if (ios /= 0) call grid_fail(fname, 'cannot open grid file')

        ! 첫 줄 숫자 개수: 1 = 블록 수 줄, 2 = ni nj (2D), 3 = ni nj nk (3D)
        read(u,'(a)',iostat=ios) line
        if (ios /= 0) call grid_fail(fname, 'empty file')
        nt = count_tokens(line)
        has_nb = (nt == 1)
        nb = 1
        if (has_nb) then
            read(line,*) nb
            read(u,'(a)',iostat=ios) line
            if (ios /= 0) call grid_fail(fname, 'cannot read block size line')
            nt = count_tokens(line)
        end if
        ! 크기 줄: 블록마다 한 줄 (2 또는 3 개) 이거나, 한 줄에 전부 (2*nb 또는 3*nb 개)
        if (nt == 2 .or. nt == 3) then
            ndim = nt
        else if (nb > 1 .and. nt == 3*nb) then
            ndim = 3
        else if (nb > 1 .and. nt == 2*nb) then
            ndim = 2
        else
            call grid_fail(fname, 'cannot read grid size (expected ni nj or ni nj nk per block)')
        end if

        ! 처음부터 다시: 크기를 목록 읽기로 (줄바꿈 위치 상관없음)
        rewind(u)
        if (has_nb) read(u,*)
        allocate(dims(3,nb));  dims = 1
        read(u,*,iostat=ios) (dims(1:ndim,m), m=1,nb)
        if (ios /= 0) call grid_fail(fname, 'cannot read grid size')
        if (any(dims(1:2,:) < 2) .or. any(dims(3,:) < 1)) call grid_fail(fname, 'bad grid size')
        if (ndim == 3 .and. any(dims(3,:) < 2)) &
            call grid_fail(fname, '3D grid with nk = 1: write it as a 2D file and use extrude-n')

        ! 블록마다 x, y (, z) 를 한 read 문으로 (x 가 줄 중간에서 끝나고 y 가 같은 줄에서 시작하는 파일이 많음)
        allocate(rb(nb))
        do m = 1,nb
            associate (r => rb(m))
                r%ni = dims(1,m);  r%nj = dims(2,m);  r%nk = dims(3,m);  r%ndim = ndim
                allocate(r%x(0:r%ni-1,0:r%nj-1,0:r%nk-1), r%y(0:r%ni-1,0:r%nj-1,0:r%nk-1))
                if (ndim == 2) then
                    read(u,*,iostat=ios) ((r%x(in,jn,0), in=0,r%ni-1), jn=0,r%nj-1), &
                                         ((r%y(in,jn,0), in=0,r%ni-1), jn=0,r%nj-1)
                else
                    allocate(r%z(0:r%ni-1,0:r%nj-1,0:r%nk-1))
                    read(u,*,iostat=ios) (((r%x(in,jn,kn), in=0,r%ni-1), jn=0,r%nj-1), kn=0,r%nk-1), &
                                         (((r%y(in,jn,kn), in=0,r%ni-1), jn=0,r%nj-1), kn=0,r%nk-1), &
                                         (((r%z(in,jn,kn), in=0,r%ni-1), jn=0,r%nj-1), kn=0,r%nk-1)
                end if
            end associate
            if (ios /= 0) call grid_fail(fname, 'file ended early or has bad data (is form = 1 right?)')
        end do
        close(u)
    end subroutine read_plot3d_ascii

    !> Fortran unformatted 순차 파일: 각 레코드 앞뒤에 4바이트 길이 표시가 붙음
    !>   레코드 길이로 블록 수 줄 / 2D, 3D / 정밀도 / iblank 를 알아냄 (블록마다 좌표 레코드 하나)
    subroutine read_plot3d_binary(fname, rb)
        use, intrinsic :: iso_fortran_env, only: int32, real32, real64
        character(len=*), intent(in) :: fname
        type(rawblk_t), allocatable, intent(out) :: rb(:)

        integer :: u, ios, ndim, npts, nb, mb, c, prec
        integer(int32) :: m, mend, nb32
        integer(int32), allocatable :: d(:), iblank(:)
        real(real64),   allocatable :: b8(:)
        real(real32),   allocatable :: b4(:)
        real(dp),       allocatable :: crd(:,:)          ! crd(:, c) : c = 1 x, 2 y, 3 z 전체
        logical :: has_iblank

        open(newunit=u, file=fname, access='stream', form='unformatted', status='old', action='read', iostat=ios)
        if (ios /= 0) call grid_fail(fname, 'cannot open grid file')

        ! 첫 레코드: 4바이트 = 블록 수,  아니면 바로 크기 (8 = ni nj, 12 = ni nj nk)
        read(u,iostat=ios) m
        if (ios /= 0) call grid_fail(fname, 'empty file')
        nb = 1
        if (m == 4) then
            read(u,iostat=ios) nb32, mend
            if (ios /= 0 .or. mend /= 4) call grid_fail(fname, 'broken record (block count)')
            nb = nb32
            if (nb < 1) call grid_fail(fname, 'bad block count')
            read(u,iostat=ios) m
            if (ios /= 0) call grid_fail(fname, 'cannot read block size record')
        end if
        if (m == 8*nb) then
            ndim = 2
        else if (m == 12*nb) then
            ndim = 3
        else
            if (any(byteswap(m) == [4, 8, 12])) then
                call grid_fail(fname, 'looks like a big-endian file; please export little-endian')
            end if
            call grid_fail(fname, 'not a Fortran unformatted PLOT3D file (is form = 0 right? ascii is 1)')
        end if
        allocate(d(ndim*nb))
        read(u,iostat=ios) d, mend
        if (ios /= 0 .or. mend /= m) call grid_fail(fname, 'broken record (grid size)')

        allocate(rb(nb))
        do mb = 1,nb
            associate (r => rb(mb))
                r%ndim = ndim
                r%ni = d(ndim*(mb-1)+1);  r%nj = d(ndim*(mb-1)+2);  r%nk = 1
                if (ndim == 3) r%nk = d(ndim*(mb-1)+3)
                if (r%ni < 2 .or. r%nj < 2 .or. r%nk < 1) call grid_fail(fname, 'bad grid size')
                if (ndim == 3 .and. r%nk < 2) &
                    call grid_fail(fname, '3D grid with nk = 1: write it as a 2D file and use extrude-n')

                ! 좌표 레코드 길이로 정밀도와 iblank 판단
                npts = r%ni*r%nj*r%nk
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
                    deallocate(b8)
                else
                    allocate(b4(npts))
                    do c = 1,ndim
                        if (ios == 0) read(u,iostat=ios) b4
                        crd(:,c) = real(b4, dp)
                    end do
                    deallocate(b4)
                end if
                if (ios == 0 .and. has_iblank) then
                    allocate(iblank(npts))
                    read(u,iostat=ios) iblank
                    deallocate(iblank)
                end if
                if (ios == 0) read(u,iostat=ios) mend
                if (ios /= 0 .or. mend /= m) call grid_fail(fname, 'file ended early or broken coordinate record')

                allocate(r%x(0:r%ni-1,0:r%nj-1,0:r%nk-1), r%y(0:r%ni-1,0:r%nj-1,0:r%nk-1))
                r%x = reshape(crd(:,1), [r%ni, r%nj, r%nk])
                r%y = reshape(crd(:,2), [r%ni, r%nj, r%nk])
                if (ndim == 3) then
                    allocate(r%z(0:r%ni-1,0:r%nj-1,0:r%nk-1))
                    r%z = reshape(crd(:,3), [r%ni, r%nj, r%nk])
                end if
                deallocate(crd)
            end associate
        end do
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

    !> 현재 블록의 노드 (in,jn,kn) 좌표
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

    !> 블록 b 의 노드 좌표만으로 셀 부피, 셀 중심, 면 법선/넓이 계산
    !>   부피 : 발산 정리  V = 1/3 * sum_면 (면 중심 - 셀 중심) . (바깥 면 벡터)
    subroutine compute_metrics(b)
        integer, intent(in) :: b
        integer  :: ic, jc, kc, fi, fj, fk, mx, my, mz
        real(dp) :: s(3), c(3)

        mx = blk(b)%nx;  my = blk(b)%ny;  mz = blk(b)%nz
        associate (bk => blk(b))
            allocate(bk%xc(mx,my,mz), bk%yc(mx,my,mz), bk%zc(mx,my,mz), bk%vol(mx,my,mz))
            allocate(bk%si(3,0:mx,my,mz), bk%si_a(0:mx,my,mz))
            allocate(bk%sj(3,mx,0:my,mz), bk%sj_a(mx,0:my,mz))
            allocate(bk%sk(3,mx,my,0:mz), bk%sk_a(mx,my,0:mz))
        end associate
        call use_block(b)

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

    !> 격자 점검 (블록마다): 부피 부호, 셀마다 면 벡터 합(= 0이어야 함). 마지막에 전체 합계
    subroutine check_grid()
        integer  :: ic, jc, kc, b
        integer(8) :: ncell
        real(dp) :: cs(3), max_closure, vtot, vmin, vmax, bmin(3), bmax(3)
        logical  :: bad

        bad = .false.;  ncell = 0;  vtot = 0;  vmin = huge(1.0_dp);  vmax = 0
        bmin = huge(1.0_dp);  bmax = -huge(1.0_dp)
        do b = 1,nblk
            call use_block(b)
            if (any(vol <= 0.0_dp)) then
                write(*,'(a,i0,a,es12.4,a,i0,a,i0)') ' Error: block ', b, ': non-positive cell volume, min =', minval(vol), &
                    ',  cells with V <= 0: ', count(vol <= 0.0_dp), ' of ', size(vol)
                if (all(vol < 0.0_dp)) write(*,'(a)') '   the block is left-handed: reverse one index direction (i, j or k)'
                bad = .true.
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

            write(*,'(a,i5,a,i0,a,i0,a,i0,a,i0,a,es10.3,a,es10.3)') ' block', b, ' : ', nx, ' x ', ny, ' x ', nz, &
                '  (', nx*ny*nz, ' cells)   min V ', minval(vol), '   closure ', max_closure
            ncell = ncell + int(nx,8)*ny*nz
            vtot  = vtot + sum(vol);  vmin = min(vmin, minval(vol));  vmax = max(vmax, maxval(vol))
            bmin  = min(bmin, [minval(x), minval(y), minval(z)])
            bmax  = max(bmax, [maxval(x), maxval(y), maxval(z)])
        end do
        if (bad) error stop

        write(*,'(a,i0,a,i0,a)')  ' grid       : ', nblk, ' block(s), ', ncell, ' cells'
        write(*,'(a,6es11.3)')   ' bounds     : ', bmin(1), bmax(1), bmin(2), bmax(2), bmin(3), bmax(3)
        write(*,'(a,es12.4)')    ' total vol  : ', vtot
        write(*,'(a,2es12.4)')   ' min/max V  : ', vmin, vmax
    end subroutine check_grid

end module grid
