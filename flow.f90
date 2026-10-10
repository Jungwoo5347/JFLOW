module flow
    use euler
    use grid
    use init_cond
    implicit none

    real(dp), allocatable :: Q(:,:,:)
    integer, parameter :: BC_inlet = 1, BC_outlet = 2, BC_wall = 3, BC_farfield = 4

    ! 경계 번호: 1 = 왼쪽(i=1), 2 = 오른쪽(i=nx), 3 = 아래(j=1), 4 = 위(j=ny)
    integer :: bc_type(4) = [BC_outlet, BC_outlet, BC_outlet, BC_outlet]


    contains

    subroutine init_flow()

        integer :: i,j
        real(dp) :: W(nvar)

        allocate(Q(nvar,0:nx+1,0:ny+1))

        W(1) = rho_init
        W(2) = u_init
        W(3) = v_init
        W(4) = p_init

        do j = 1,ny
            do i = 1,nx
                call p2c(W,Q(:,i,j))
            end do
        end do

    end subroutine init_flow

    !> 고스트 셀 채우기 (매 스텝, 잔차 계산 전에 호출)
    subroutine bc()
        integer :: i,j

        ! 왼쪽 (i=0), 오른쪽 (i=nx+1) : i-면 법선 사용
        do j = 1,ny
            call ghost(bc_type(1), Q(:,1,j),  si_x(0,j),  si_y(0,j),  Q(:,0,j))
            call ghost(bc_type(2), Q(:,nx,j), si_x(nx,j), si_y(nx,j), Q(:,nx+1,j))
        end do

        ! 아래 (j=0), 위 (j=ny+1) : j-면 법선 사용
        do i = 1,nx
            call ghost(bc_type(3), Q(:,i,1),  sj_x(i,0),  sj_y(i,0),  Q(:,i,0))
            call ghost(bc_type(4), Q(:,i,ny), sj_x(i,ny), sj_y(i,ny), Q(:,i,ny+1))
        end do

    end subroutine bc

    !> 경계 면 하나의 고스트 셀 값
    !>   qin : 경계에 붙은 내부 셀,  (fnx, fny) : 경계 면 단위 법선,  qg : 고스트 셀
    subroutine ghost(btype, qin, fnx, fny, qg)
        integer,  intent(in)  :: btype
        real(dp), intent(in)  :: qin(nvar), fnx, fny
        real(dp), intent(out) :: qg(nvar)

        real(dp) :: mn, W_inf(nvar)

        select case (btype)

        case (BC_inlet, BC_farfield)
            ! 자유류 값 고정 (초음속 입구 기준, 원방은 나중에 특성 조건으로 교체)
            W_inf = [rho_init, u_init, v_init, p_init]
            call p2c(W_inf, qg)

        case (BC_outlet)
            ! 내부 셀 그대로 복사
            qg = qin

        case (BC_wall)
            ! 밀도, 에너지는 그대로, 운동량의 법선 성분만 반전
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