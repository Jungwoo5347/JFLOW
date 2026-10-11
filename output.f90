module output
    use euler
    use grid
    implicit none

    contains

    !> Tecplot ASCII 파일 (ParaView에서도 열림)
    !>   노드 좌표 x, y + 셀 중심 값 rho, u, v, p, Mach
    !>   곡선 격자도 그대로 그려짐
    subroutine write_tecplot(fname, Q, t)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: Q(:,0:,0:)
        real(dp), intent(in) :: t

        integer  :: u, in, jn, ic, jc, k
        real(dp) :: W(nvar), c, val

        open(newunit=u, file=fname, status='replace', action='write')
        write(u,'(a)') 'TITLE = "JFLOW"'
        write(u,'(a)') 'VARIABLES = "x" "y" "rho" "u" "v" "p" "Mach"'
        write(u,'(a,es14.6,a,i0,a,i0,a)') 'ZONE T="t=', t, '", I=', nx+1, ', J=', ny+1, &
            ', DATAPACKING=BLOCK, VARLOCATION=([3-7]=CELLCENTERED)'

        ! 노드 좌표
        do jn = 0,ny
            do in = 0,nx
                write(u,'(es16.8)') x(in,jn)
            end do
        end do
        do jn = 0,ny
            do in = 0,nx
                write(u,'(es16.8)') y(in,jn)
            end do
        end do

        ! 셀 중심 값: k = 1 rho, 2 u, 3 v, 4 p, 5 Mach
        do k = 1,5
            do jc = 1,ny
                do ic = 1,nx
                    call c2p(Q(:,ic,jc), W)
                    if (k <= 4) then
                        val = W(k)
                    else
                        call sos(W, c)
                        val = sqrt(W(2)**2 + W(3)**2)/c
                    end if
                    write(u,'(es16.8)') val
                end do
            end do
        end do

        close(u)
    end subroutine write_tecplot

    !> 셀 줄 jc 하나를 따라 xc, rho, u, v, p 출력 (1D 그래프용)
    subroutine write_line(fname, Q, jc)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: Q(:,0:,0:)
        integer,  intent(in) :: jc

        integer  :: u, ic
        real(dp) :: W(nvar)

        open(newunit=u, file=fname, status='replace', action='write')
        write(u,'(a)') '#  xc  rho  u  v  p'
        do ic = 1,nx
            call c2p(Q(:,ic,jc), W)
            write(u,'(5es16.8)') xc(ic,jc), W
        end do
        close(u)
    end subroutine write_line

end module output
