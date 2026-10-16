module output
    use, intrinsic :: iso_fortran_env, only: int32, real32, real64
    use euler
    use grid
    implicit none

    ! 결과 변수: 노드 좌표 x, y, z + 셀 중심 값 rho, u, v, w, p, Mach
    integer, parameter :: nout = 9, ncellvar = 6
    character(len=4), parameter :: var_names(nout) = ['x   ', 'y   ', 'z   ', 'rho ', 'u   ', 'v   ', 'w   ', 'p   ', 'Mach']

    character(len=16), parameter :: restart_magic = 'JFLOW3D-RESTART1'

    contains

    !> 블록 b 의 셀 중심 결과 값 val(nx,ny,nz,6) : rho, u, v, w, p, Mach
    subroutine cell_values(b, val)
        integer,  intent(in) :: b
        real(dp), allocatable, intent(out) :: val(:,:,:,:)
        integer  :: ic, jc, kc
        real(dp) :: W(nvar), c
        real(dp), pointer, contiguous :: Q(:,:,:,:)

        Q => blk(b)%Q
        allocate(val(blk(b)%nx, blk(b)%ny, blk(b)%nz, ncellvar))
        !$omp parallel do collapse(2) private(ic, W, c)
        do kc = 1,blk(b)%nz
            do jc = 1,blk(b)%ny
                do ic = 1,blk(b)%nx
                    call c2p(Q(:,ic,jc,kc), W)
                    call sos(W, c)
                    val(ic,jc,kc,1:5) = W
                    val(ic,jc,kc,6)   = norm2(W(2:4))/c
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine cell_values

    !> Tecplot ASCII 파일 (ParaView 에서도 열림), 블록 하나 = zone 하나
    !>   3D 는 파일이 커서 한 줄에 6 개씩 씀 (BLOCK 형식은 줄바꿈 위치 상관없음)
    !>   지수는 항상 E 와 세 자리 (es17.8e3): 1e-100 보다 작은 값이 E 없이 "1.0-109" 로 찍히는 것 방지
    subroutine write_tecplot(fname, t)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: t

        integer  :: u, k, b
        real(dp), allocatable :: val(:,:,:,:)

        open(newunit=u, file=fname, status='replace', action='write')
        write(u,'(a)') 'TITLE = "JFLOW3D"'
        write(u,'(a)') 'VARIABLES = "x" "y" "z" "rho" "u" "v" "w" "p" "Mach"'
        do b = 1,nblk
            call cell_values(b, val)
            write(u,'(a,i0,a,i0,a,i0,a,i0,a,es14.6,a)') 'ZONE T="block ', b, '", I=', blk(b)%nx+1, ', J=', blk(b)%ny+1, &
                ', K=', blk(b)%nz+1, ', SOLUTIONTIME=', t, ', DATAPACKING=BLOCK, VARLOCATION=([4-9]=CELLCENTERED)'
            ! 노드 좌표 (i 가 가장 빠르게)
            write(u,'(6es17.8e3)') blk(b)%x
            write(u,'(6es17.8e3)') blk(b)%y
            write(u,'(6es17.8e3)') blk(b)%z
            do k = 1,ncellvar
                write(u,'(6es17.8e3)') val(:,:,:,k)
            end do
        end do
        close(u)
    end subroutine write_tecplot

    !> Tecplot 바이너리 .plt (형식 #!TDV112, 단정밀도), 블록 하나 = zone 하나
    !>   preplot 이 만드는 파일과 같은 구조. ASCII 의 약 1/4 크기, 쓰고 읽기가 훨씬 빠름
    !>   셀 중심 변수는 ordered zone 규칙대로 IMax x JMax x (KMax-1) 로 채워 씀 (i, j 끝 칸은 0)
    subroutine write_plt(fname, t)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: t

        integer :: u, b, k
        real(dp), allocatable :: val(:,:,:,:)
        real(real32), allocatable :: pad(:,:,:)

        open(newunit=u, file=fname, access='stream', form='unformatted', status='replace', action='write')

        ! ---- 머리 ----
        write(u) '#!TDV112'
        write(u) 1_int32, 0_int32                       ! 바이트 순서 확인용 1, 파일 종류 0 (격자 + 해)
        call put_str(u, 'JFLOW3D')
        write(u) int(nout, int32)
        do k = 1,nout
            call put_str(u, trim(var_names(k)))
        end do
        do b = 1,nblk
            write(u) 299.0_real32                         ! zone 표시
            call put_str(u, 'block '//itoa(b))
            write(u) -1_int32, -2_int32                   ! parent zone, strand (-2: 자동, preplot 과 같음)
            write(u) real(t, real64)                      ! solution time
            write(u) -1_int32, 0_int32                    ! (사용 안 함), zone type 0 = ordered
            write(u) 1_int32                              ! 변수 위치 지정함
            write(u) (0_int32, k=1,3), (1_int32, k=1,ncellvar)   ! 0 = 노드, 1 = 셀 중심
            write(u) 0_int32, 0_int32                     ! face neighbor 없음
            write(u) int(blk(b)%nx+1, int32), int(blk(b)%ny+1, int32), int(blk(b)%nz+1, int32)
            write(u) 0_int32                              ! aux data 없음
        end do
        write(u) 357.0_real32                             ! 머리 끝

        ! ---- 데이터 (zone 마다) ----
        do b = 1,nblk
            associate (bk => blk(b))
                call cell_values(b, val)
                write(u) 299.0_real32
                write(u) (1_int32, k=1,nout)              ! 모든 변수 단정밀도
                write(u) 0_int32, 0_int32, -1_int32       ! passive 없음, 변수 공유 없음, 연결 공유 없음
                write(u) real(minval(bk%x),real64), real(maxval(bk%x),real64), &
                         real(minval(bk%y),real64), real(maxval(bk%y),real64), &
                         real(minval(bk%z),real64), real(maxval(bk%z),real64)
                do k = 1,ncellvar
                    write(u) real(minval(val(:,:,:,k)),real64), real(maxval(val(:,:,:,k)),real64)
                end do
                write(u) real(bk%x, real32)
                write(u) real(bk%y, real32)
                write(u) real(bk%z, real32)
                allocate(pad(bk%nx+1, bk%ny+1, bk%nz))
                pad = 0.0_real32
                do k = 1,ncellvar
                    pad(1:bk%nx, 1:bk%ny, :) = real(val(:,:,:,k), real32)
                    write(u) pad
                end do
                deallocate(pad)
            end associate
        end do
        close(u)

    contains

        !> 문자열: 글자마다 int32, 끝에 0
        subroutine put_str(u, s)
            integer,          intent(in) :: u
            character(len=*), intent(in) :: s
            integer :: c
            write(u) (int(iachar(s(c:c)), int32), c=1,len(s)), 0_int32
        end subroutine put_str

        function itoa(i) result(s)
            integer, intent(in) :: i
            character(len=:), allocatable :: s
            character(len=16) :: buf
            write(buf,'(i0)') i
            s = trim(buf)
        end function itoa

    end subroutine write_plt

    !> 블록 b 의 i 방향 셀 줄 (jc, kc) 하나를 따라 xc, yc, zc, rho, u, v, w, p 출력 (1D 그래프용)
    subroutine write_line(fname, b, jc, kc)
        character(len=*), intent(in) :: fname
        integer,  intent(in) :: b, jc, kc

        integer  :: u, ic
        real(dp) :: W(nvar)

        open(newunit=u, file=fname, status='replace', action='write')
        write(u,'(a)') '#  xc  yc  zc  rho  u  v  w  p'
        do ic = 1,blk(b)%nx
            call c2p(blk(b)%Q(:,ic,jc,kc), W)
            write(u,'(8es17.8e3)') blk(b)%xc(ic,jc,kc), blk(b)%yc(ic,jc,kc), blk(b)%zc(ic,jc,kc), W
        end do
        close(u)
    end subroutine write_line

    !> 재시작 파일 쓰기 (배정밀도, 내부 셀 Q 전부 + 반복 수, 시간, 잔차 정규화 기준)
    !>   먼저 .tmp 로 쓰고 이름을 바꿔서, 쓰는 도중에 멈춰도 이전 재시작 파일이 남게 함
    subroutine write_restart(fname, step, t, res0)
        character(len=*), intent(in) :: fname
        integer,  intent(in) :: step
        real(dp), intent(in) :: t, res0(nvar)
        integer :: u, b, ios

        open(newunit=u, file=fname//'.tmp', access='stream', form='unformatted', status='replace', action='write')
        write(u) restart_magic, int(nblk, int32)
        write(u) (int(blk(b)%nx,int32), int(blk(b)%ny,int32), int(blk(b)%nz,int32), b=1,nblk)
        write(u) int(step, int32), real(t, real64), real(res0, real64)
        do b = 1,nblk
            associate (bk => blk(b))
                write(u) real(bk%Q(:,1:bk%nx,1:bk%ny,1:bk%nz), real64)
            end associate
        end do
        close(u)

        ! 이전 파일을 지우고 .tmp 를 제 이름으로
        open(newunit=u, file=fname, status='old', iostat=ios)
        if (ios == 0) close(u, status='delete')
        call rename(fname//'.tmp', fname)
    end subroutine write_restart

    !> 재시작 파일 읽기: 블록 수와 크기가 지금 격자와 같아야 함
    subroutine read_restart(fname, step, t, res0)
        character(len=*), intent(in)  :: fname
        integer,  intent(out) :: step
        real(dp), intent(out) :: t, res0(nvar)
        integer :: u, ios, b
        integer(int32) :: nb, d(3), st
        character(len=16) :: magic

        open(newunit=u, file=fname, access='stream', form='unformatted', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(a,a,a)') ' Error: restart = 1 but cannot open ', trim(fname), ' (run once without restart first)'
            error stop
        end if
        read(u, iostat=ios) magic, nb
        if (ios /= 0 .or. magic /= restart_magic) then
            write(*,'(a,a,a)') ' Error: ', trim(fname), ' is not a JFLOW3D restart file'
            error stop
        end if
        if (nb /= nblk) then
            write(*,'(a,i0,a,i0,a)') ' Error: restart file has ', nb, ' block(s), the grid has ', nblk
            error stop
        end if
        do b = 1,nblk
            read(u) d
            if (any(d /= [blk(b)%nx, blk(b)%ny, blk(b)%nz])) then
                write(*,'(a,i0,a,3(1x,i0),a,3(1x,i0))') ' Error: restart file block ', b, ' cells', d, &
                    ' but the grid has', blk(b)%nx, blk(b)%ny, blk(b)%nz
                error stop
            end if
        end do
        read(u) st, t, res0
        step = st
        do b = 1,nblk
            associate (bk => blk(b))
                read(u, iostat=ios) bk%Q(:,1:bk%nx,1:bk%ny,1:bk%nz)
            end associate
            if (ios /= 0) then
                write(*,'(a,a,a)') ' Error: ', trim(fname), ' ended early'
                error stop
            end if
        end do
        close(u)
    end subroutine read_restart

end module output
