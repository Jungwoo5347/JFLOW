program main
    use grid
    use euler

    implicit none

    real(dp) :: W(nvar), Q(nvar), W2(nvar), c

    call init_grid()
    call check_grid()

    W = [1.0_dp, 0.5_dp, 0.2_dp, 2.0_dp]
    call p2c(W, Q)
    call c2p(Q, W2)
    call sos(W, c)
    print *, 'W  =', W
    print *, 'W2 =', W2      ! W와 같아야 함
    print *, 'c  =', c       ! 1.6733 (= sqrt(1.4*2/1))





end program main