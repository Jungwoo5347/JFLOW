module grid
    use euler
    implicit none
    

    integer  :: nx, ny                          ! 셀 개수
    real(dp) :: xmin, xmax, ymin, ymax
    real(dp), allocatable :: x(:,:), y(:,:)     ! 노드 좌표 (0:nx, 0:ny)

    ! 셀 (1:nx, 1:ny) : 셀 (i,j)는 노드 (i-1,j-1), (i,j-1), (i,j), (i-1,j)로 둘러싸임
    real(dp), allocatable :: xc(:,:), yc(:,:)   ! 셀 중심
    real(dp), allocatable :: area(:,:)          ! 셀 면적

    ! i-면 (0:nx, 1:ny) : 면 i는 셀 i와 i+1 사이, 법선은 +i 방향
    real(dp), allocatable :: si_x(:,:), si_y(:,:), si_len(:,:)
    ! j-면 (1:nx, 0:ny) : 면 j는 셀 j와 j+1 사이, 법선은 +j 방향
    real(dp), allocatable :: sj_x(:,:), sj_y(:,:), sj_len(:,:)

contains

    !> 직교 격자 생성 (테스트용). dx, dy는 여기서만 쓰고 버림
    subroutine init_grid()
        integer  :: i, j
        real(dp) :: dx, dy

        nx = 100;  ny = 20
        xmin = 0.0_dp;  xmax = 1.0_dp
        ymin = 0.0_dp;  ymax = 0.2_dp

        allocate(x(0:nx, 0:ny), y(0:nx, 0:ny))

        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny
        do j = 0, ny
            do i = 0, nx
                x(i,j) = xmin + i*dx
                y(i,j) = ymin + j*dy
            end do
        end do

        call compute_metrics()
    end subroutine init_grid

    !> 노드 좌표 x, y만으로 셀 면적, 셀 중심, 면 법선/길이 계산
    !> (격자를 어디서 가져왔든 x, y만 채운 뒤 이걸 부르면 됨)
    subroutine compute_metrics()
        integer  :: i, j
        real(dp) :: ex, ey, d1x, d1y, d2x, d2y

        allocate(xc(nx,ny), yc(nx,ny), area(nx,ny))
        allocate(si_x(0:nx,ny), si_y(0:nx,ny), si_len(0:nx,ny))
        allocate(sj_x(nx,0:ny), sj_y(nx,0:ny), sj_len(nx,0:ny))

        ! 셀 면적 (두 대각선 외적의 절반), 셀 중심 (네 꼭짓점 평균)
        do j = 1, ny
            do i = 1, nx
                d1x = x(i,j)   - x(i-1,j-1);  d1y = y(i,j)   - y(i-1,j-1)
                d2x = x(i-1,j) - x(i,j-1);    d2y = y(i-1,j) - y(i,j-1)
                area(i,j) = 0.5_dp * (d1x*d2y - d1y*d2x)

                xc(i,j) = 0.25_dp * (x(i-1,j-1) + x(i,j-1) + x(i,j) + x(i-1,j))
                yc(i,j) = 0.25_dp * (y(i-1,j-1) + y(i,j-1) + y(i,j) + y(i-1,j))
            end do
        end do

        ! i-면: 노드 (i,j-1) -> (i,j) 변 벡터 (ex,ey), 법선 = (ey, -ex)
        do j = 1, ny
            do i = 0, nx
                ex = x(i,j) - x(i,j-1)
                ey = y(i,j) - y(i,j-1)
                si_len(i,j) = sqrt(ex**2 + ey**2)
                si_x(i,j)   =  ey / si_len(i,j)
                si_y(i,j)   = -ex / si_len(i,j)
            end do
        end do

        ! j-면: 노드 (i-1,j) -> (i,j) 변 벡터 (ex,ey), 법선 = (-ey, ex)
        do j = 0, ny
            do i = 1, nx
                ex = x(i,j) - x(i-1,j)
                ey = y(i,j) - y(i-1,j)
                sj_len(i,j) = sqrt(ex**2 + ey**2)
                sj_x(i,j)   = -ey / sj_len(i,j)
                sj_y(i,j)   =  ex / sj_len(i,j)
            end do
        end do
    end subroutine compute_metrics

    !> 격자 점검: 면적 부호, 전체 면적, 셀마다 면 벡터 합(= 0이어야 함)
    subroutine check_grid()
        integer  :: i, j
        real(dp) :: cx, cy, max_closure

        if (any(area <= 0.0_dp)) then
            write(*,*) 'Error: non-positive cell area, min =', minval(area)
            error stop
        end if

        max_closure = 0.0_dp
        do j = 1, ny
            do i = 1, nx
                cx = si_x(i,j)*si_len(i,j) - si_x(i-1,j)*si_len(i-1,j) &
                   + sj_x(i,j)*sj_len(i,j) - sj_x(i,j-1)*sj_len(i,j-1)
                cy = si_y(i,j)*si_len(i,j) - si_y(i-1,j)*si_len(i-1,j) &
                   + sj_y(i,j)*sj_len(i,j) - sj_y(i,j-1)*sj_len(i,j-1)
                max_closure = max(max_closure, abs(cx), abs(cy))
            end do
        end do

        write(*,'(a,i0,a,i0)')   ' grid       : ', nx, ' x ', ny
        write(*,'(a,es12.4)')    ' total area : ', sum(area)
        write(*,'(a,2es12.4)')   ' min/max A  : ', minval(area), maxval(area)
        write(*,'(a,es12.4)')    ' closure    : ', max_closure
    end subroutine check_grid

end module grid
