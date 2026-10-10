module euler
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: gam = 1.4_dp
    integer, parameter :: nvar = 4

    contains

    ! W = (rho, u, v, p),  Q = (rho, rho*u, rho*v, E)

    subroutine p2c(W, Q)
        real(dp), intent(in)  :: W(:)
        real(dp), intent(out) :: Q(:)

        Q(1) = W(1)
        Q(2) = W(1)*W(2)
        Q(3) = W(1)*W(3)
        Q(4) = W(4)/(gam - 1.0_dp) + 0.5_dp*W(1)*(W(2)**2 + W(3)**2)

    end subroutine p2c


    subroutine c2p(Q,W)
        real(dp), intent(in) :: Q(:)
        real(dp), intent(inout) :: W(:)

        W(1) = Q(1)
        W(2) = Q(2)/Q(1)
        W(3) = Q(3)/Q(1)
        W(4) = (gam - 1.0_dp)*(Q(4) - 0.5_dp*W(1)*(W(2)**2 + W(3)**2))

    end subroutine c2p

    subroutine sos(W,c)
        real(dp), intent(in) :: W(:)
        real(dp), intent(out) :: c

        c = sqrt(gam*W(4)/W(1))

    end subroutine sos





end module euler