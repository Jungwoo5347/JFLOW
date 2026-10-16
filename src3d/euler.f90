module euler
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp) :: gam = 1.4_dp                    ! 비열비 (flow.in 의 gamma 로 바뀜)
    integer, parameter :: nvar = 5

    contains

    ! W = (rho, u, v, w, p),  Q = (rho, rho*u, rho*v, rho*w, E)

    subroutine p2c(W, Q)
        real(dp), intent(in)  :: W(:)
        real(dp), intent(out) :: Q(:)

        Q(1) = W(1)
        Q(2) = W(1)*W(2)
        Q(3) = W(1)*W(3)
        Q(4) = W(1)*W(4)
        Q(5) = W(5)/(gam - 1.0_dp) + 0.5_dp*W(1)*(W(2)**2 + W(3)**2 + W(4)**2)

    end subroutine p2c


    subroutine c2p(Q,W)
        real(dp), intent(in) :: Q(:)
        real(dp), intent(inout) :: W(:)

        W(1) = Q(1)
        W(2) = Q(2)/Q(1)
        W(3) = Q(3)/Q(1)
        W(4) = Q(4)/Q(1)
        W(5) = (gam - 1.0_dp)*(Q(5) - 0.5_dp*W(1)*(W(2)**2 + W(3)**2 + W(4)**2))

    end subroutine c2p

    subroutine sos(W,c)
        real(dp), intent(in) :: W(:)
        real(dp), intent(out) :: c

        c = sqrt(gam*W(5)/W(1))

    end subroutine sos

    !> 외적 a x b
    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross

end module euler
