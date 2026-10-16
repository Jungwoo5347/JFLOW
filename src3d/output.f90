module output
    use euler
    use grid
    implicit none

    contains

    !> Tecplot ASCII 파일 (ParaView 에서도 열림)
    !>   노드 좌표 x, y, z + 셀 중심 값 rho, u, v, w, p, Mach
    !>   3D 는 파일이 커서 한 줄에 6 개씩 씀 (BLOCK 형식은 줄바꿈 위치 상관없음)
    subroutine write_tecplot(fname, Q, t)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: Q(:,1-ng:,1-ng:,1-ng:)
        real(dp), intent(in) :: t

        integer  :: u, ic, jc, kc, k
        real(dp) :: W(nvar), c
        real(dp), allocatable :: val(:,:,:,:)

        ! 셀 중심 값: k = 1 rho, 2 u, 3 v, 4 w, 5 p, 6 Mach
        allocate(val(nx,ny,nz,6))
        !$omp parallel do collapse(2) private(ic, W, c)
        do kc = 1,nz
            do jc = 1,ny
                do ic = 1,nx
                    call c2p(Q(:,ic,jc,kc), W)
                    call sos(W, c)
                    val(ic,jc,kc,1:5) = W
                    val(ic,jc,kc,6)   = norm2(W(2:4))/c
                end do
            end do
        end do
        !$omp end parallel do

        open(newunit=u, file=fname, status='replace', action='write')
        write(u,'(a)') 'TITLE = "JFLOW3D"'
        write(u,'(a)') 'VARIABLES = "x" "y" "z" "rho" "u" "v" "w" "p" "Mach"'
        write(u,'(a,es14.6,a,i0,a,i0,a,i0,a)') 'ZONE T="t=', t, '", I=', nx+1, ', J=', ny+1, ', K=', nz+1, &
            ', DATAPACKING=BLOCK, VARLOCATION=([4-9]=CELLCENTERED)'

        ! 노드 좌표 (i 가 가장 빠르게)
        write(u,'(6es16.8)') x
        write(u,'(6es16.8)') y
        write(u,'(6es16.8)') z
        do k = 1,6
            write(u,'(6es16.8)') val(:,:,:,k)
        end do

        close(u)
    end subroutine write_tecplot

    !> i 방향 셀 줄 (jc, kc) 하나를 따라 xc, yc, zc, rho, u, v, w, p 출력 (1D 그래프용)
    subroutine write_line(fname, Q, jc, kc)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: Q(:,1-ng:,1-ng:,1-ng:)
        integer,  intent(in) :: jc, kc

        integer  :: u, ic
        real(dp) :: W(nvar)

        open(newunit=u, file=fname, status='replace', action='write')
        write(u,'(a)') '#  xc  yc  zc  rho  u  v  w  p'
        do ic = 1,nx
            call c2p(Q(:,ic,jc,kc), W)
            write(u,'(8es16.8)') xc(ic,jc,kc), yc(ic,jc,kc), zc(ic,jc,kc), W
        end do
        close(u)
    end subroutine write_line

end module output
