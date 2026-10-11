module grid
    use euler
    implicit none
    

    ! 인덱스 이름 규칙
    !   셀   : ic, jc   (내부 1:nx, 1:ny / 고스트 0, nx+1, 0, ny+1)
    !   노드 : in, jn   (0:nx, 0:ny)
    !   면   : fi       (i-면 번호 0:nx, 셀 fi 와 fi+1 사이)
    !          fj       (j-면 번호 0:ny, 셀 fj 와 fj+1 사이)

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

    !> 직교 격자 생성 (테스트용). dx, dy는 여기서만 쓰고 버림
    subroutine init_grid()
        integer  :: in, jn
        real(dp) :: dx, dy

        nx = 100;  ny = 20
        xmin = 0.0_dp;  xmax = 1.0_dp
        ymin = 0.0_dp;  ymax = 0.2_dp

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

    !> PLOT3D 격자 읽기 (2D, 단일 블록, 포맷(ASCII), 배정밀도)
    !>   파일 구조:  ni nj  /  x 전부  /  y 전부   (i가 가장 빠르게 변함)
    !>   노드 ni x nj  ->  셀 nx = ni-1, ny = nj-1
    subroutine read_plot3d(fname)
        character(len=*), intent(in) :: fname
        integer :: u, ios, ni, nj, in, jn

        open(newunit=u, file=fname, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: cannot open grid file ', trim(fname)
            error stop
        end if

        read(u,*,iostat=ios) ni, nj
        if (ios /= 0) then
            write(*,*) 'Error: cannot read grid size from ', trim(fname)
            error stop
        end if
        nx = ni - 1
        ny = nj - 1

        ! x와 y를 한 read 문으로 읽어야 함: read 문은 매번 새 줄에서 시작하는데,
        ! x가 줄 중간에서 끝나고 y가 같은 줄에서 시작하는 파일이 많음
        allocate(x(0:nx, 0:ny), y(0:nx, 0:ny))
        read(u,*,iostat=ios) ((x(in,jn), in=0,nx), jn=0,ny), ((y(in,jn), in=0,nx), jn=0,ny)
        if (ios /= 0) then
            write(*,*) 'Error: grid file ended early or has bad data: ', trim(fname)
            error stop
        end if
        close(u)

        xmin = minval(x);  xmax = maxval(x)
        ymin = minval(y);  ymax = maxval(y)

        call compute_metrics()
    end subroutine read_plot3d

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
