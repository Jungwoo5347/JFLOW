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

        integer :: ic,jc
        real(dp) :: W(nvar)

        allocate(Q(nvar,0:nx+1,0:ny+1))
        Q = 0.0_dp                       ! 고스트(모서리 포함)에 쓰레기값이 남지 않게

        do jc = 1,ny
            do ic = 1,nx
                select case (case_name)
                case ('sod')
                    if (xc(ic,jc) < x_mid) then
                        W = W_left
                    else
                        W = W_right
                    end if
                case default
                    ! 전체를 자유류로 시작
                    W = [rho_init, u_init, v_init, p_init]
                end select
                call p2c(W,Q(:,ic,jc))
            end do
        end do

    end subroutine init_flow

    !> 고스트 셀 채우기 (매 스텝, 잔차 계산 전에 호출)
    subroutine bc()
        integer :: ic,jc

        ! 왼쪽 고스트 (ic=0), 오른쪽 고스트 (ic=nx+1) : 경계 i-면 fi=0, fi=nx 의 법선 사용
        do jc = 1,ny
            call ghost(bc_type(1), Q(:,1,jc),  si_x(0,jc),  si_y(0,jc),  Q(:,0,jc))
            call ghost(bc_type(2), Q(:,nx,jc), si_x(nx,jc), si_y(nx,jc), Q(:,nx+1,jc))
        end do

        ! 아래 고스트 (jc=0), 위 고스트 (jc=ny+1) : 경계 j-면 fj=0, fj=ny 의 법선 사용
        do ic = 1,nx
            call ghost(bc_type(3), Q(:,ic,1),  sj_x(ic,0),  sj_y(ic,0),  Q(:,ic,0))
            call ghost(bc_type(4), Q(:,ic,ny), sj_x(ic,ny), sj_y(ic,ny), Q(:,ic,ny+1))
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