program main
    use grid
    use euler
    use flow

    implicit none

    call init_grid()
    call init_flow()
    call bc()

    print*, Q(:,1,1)

  





end program main