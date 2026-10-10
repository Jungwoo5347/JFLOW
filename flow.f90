module flow
    use euler
    use grid
    implicit none

    real(dp), allocatable :: Q(:,:,:)
    

    contains

    subroutine init_flow()

        integer :: i,j
        real(dp) :: W(nvar)

        allocate(Q(nvar,0:nx+1,0:ny+1))

        call c2p(Q(:,i,j),W)




    end subroutine init_flow

end module flow