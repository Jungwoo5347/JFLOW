module output
    use, intrinsic :: iso_fortran_env, only: int32, real32, real64
    use euler
    use grid
    use case_setup, only: dim_out, ref_length, rho_ref, c_ref
    implicit none

    character(len=16), parameter :: restart_magic = 'JFLOW2D-RESTART1'

    contains

    !> 결과 파일의 단위 환산 (ref_length 가 있으면 SI, 없으면 1 = 무차원 그대로)
    !>   sL 길이 [m], sR 밀도 [kg/m3], sU 속도 [m/s], sP 압력 [Pa], sT 시간 [s]
    !>   names : 변수 이름 (x y rho u v p Mach, SI 면 단위 붙임)
    subroutine out_scales(sL, sR, sU, sP, sT, names)
        real(dp), intent(out) :: sL, sR, sU, sP, sT
        character(len=12), intent(out) :: names(7)
        if (dim_out) then
            sL = ref_length;  sR = rho_ref;  sU = c_ref;  sP = rho_ref*c_ref**2;  sT = ref_length/c_ref
            names = [character(len=12) :: 'x [m]', 'y [m]', 'rho [kg/m3]', 'u [m/s]', 'v [m/s]', 'p [Pa]', 'Mach']
        else
            sL = 1.0_dp;  sR = 1.0_dp;  sU = 1.0_dp;  sP = 1.0_dp;  sT = 1.0_dp
            names = [character(len=12) :: 'x', 'y', 'rho', 'u', 'v', 'p', 'Mach']
        end if
    end subroutine out_scales

    !> Tecplot ASCII 파일 (ParaView에서도 열림)
    !>   노드 좌표 x, y + 셀 중심 값 rho, u, v, p, Mach
    !>   곡선 격자도 그대로 그려짐.  ref_length 가 있으면 SI 단위 (out_scales)
    subroutine write_tecplot(fname, Q, t_nd)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: Q(:,1-ng:,1-ng:)
        real(dp), intent(in) :: t_nd                 ! 무차원 시간

        integer  :: u, in, jn, ic, jc, k
        real(dp) :: W(nvar), c, val, t, sL, sR, sU, sP, sT, sc(5)
        character(len=12) :: names(7)

        call out_scales(sL, sR, sU, sP, sT, names)
        sc = [sR, sU, sU, sP, 1.0_dp]
        t  = t_nd*sT
        open(newunit=u, file=fname, status='replace', action='write')
        write(u,'(a)') 'TITLE = "JFLOW"'
        write(u,'(a,7(a,a,a))') 'VARIABLES =', (' "', trim(names(k)), '"', k=1,7)
        ! STRANDID, SOLUTIONTIME: 여러 시간의 파일을 함께 열면 Tecplot 이 시간 이력 (애니메이션) 으로 묶음
        write(u,'(a,es14.6,a,i0,a,i0,a,es22.14)') 'ZONE T="t=', t, '", I=', nx+1, ', J=', ny+1, &
            ', DATAPACKING=BLOCK, VARLOCATION=([3-7]=CELLCENTERED), STRANDID=1, SOLUTIONTIME=', t

        ! 노드 좌표
        do jn = 0,ny
            do in = 0,nx
                write(u,'(es16.8)') x(in,jn)*sL
            end do
        end do
        do jn = 0,ny
            do in = 0,nx
                write(u,'(es16.8)') y(in,jn)*sL
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
                    write(u,'(es16.8)') val*sc(k)
                end do
            end do
        end do

        close(u)
    end subroutine write_tecplot

    !> Tecplot 바이너리 (.plt, #!TDV112): write_tecplot 과 같은 내용 (x, y + 셀 중심 rho, u, v, p, Mach), 단정밀도
    !>   ASCII 보다 작고 Tecplot 에서 빨리 열림. tools/read_result.m 으로 MATLAB 에서 읽음
    !>   셀 중심 값은 I x (J-1) 개로 저장 (i 방향 끝 칸은 채움 값 0, 마지막 방향 j 는 채우지 않음: preplot 과 같음)
    subroutine write_plt(fname, Q, t_nd)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: Q(:,1-ng:,1-ng:)
        real(dp), intent(in) :: t_nd                 ! 무차원 시간

        integer, parameter :: nout = 7, ncellvar = 5
        character(len=12) :: var_names(nout)
        integer  :: u, ic, jc, k
        real(dp) :: W(nvar), c, t, sL, sR, sU, sP, sT
        real(dp), allocatable :: val(:,:,:)
        real(real32), allocatable :: pad(:,:)

        call out_scales(sL, sR, sU, sP, sT, var_names)
        t = t_nd*sT

        ! 셀 중심 값: k = 1 rho, 2 u, 3 v, 4 p, 5 Mach  (SI 면 환산)
        allocate(val(nx, ny, ncellvar))
        do jc = 1,ny
            do ic = 1,nx
                call c2p(Q(:,ic,jc), W)
                call sos(W, c)
                val(ic,jc,1:4) = W*[sR, sU, sU, sP]
                val(ic,jc,5)   = sqrt(W(2)**2 + W(3)**2)/c
            end do
        end do

        open(newunit=u, file=fname, access='stream', form='unformatted', status='replace', action='write')

        ! ---- 머리 ----
        write(u) '#!TDV112'
        write(u) 1_int32, 0_int32                       ! 바이트 순서 확인용 1, 파일 종류 0 (격자 + 해)
        call put_str(u, 'JFLOW')
        write(u) int(nout, int32)
        do k = 1,nout
            call put_str(u, trim(var_names(k)))
        end do
        write(u) 299.0_real32                           ! zone 표시
        call put_str(u, 'JFLOW')
        write(u) -1_int32, 1_int32                      ! parent zone, strand 1 (여러 시간의 파일을 함께 열면 시간 이력으로 묶임)
        write(u) real(t, real64)                        ! solution time
        write(u) -1_int32, 0_int32                      ! (사용 안 함), zone type 0 = ordered
        write(u) 1_int32                                ! 변수 위치 지정함
        write(u) (0_int32, k=1,2), (1_int32, k=1,ncellvar)   ! 0 = 노드, 1 = 셀 중심
        write(u) 0_int32, 0_int32                       ! face neighbor 없음
        write(u) int(nx+1, int32), int(ny+1, int32), 1_int32
        write(u) 0_int32                                ! aux data 없음
        write(u) 357.0_real32                           ! 머리 끝

        ! ---- 데이터 ----
        write(u) 299.0_real32
        write(u) (1_int32, k=1,nout)                    ! 모든 변수 단정밀도
        write(u) 0_int32, 0_int32, -1_int32             ! passive 없음, 변수 공유 없음, 연결 공유 없음
        write(u) real(minval(x)*sL,real64), real(maxval(x)*sL,real64), real(minval(y)*sL,real64), real(maxval(y)*sL,real64)
        do k = 1,ncellvar
            write(u) real(minval(val(:,:,k)),real64), real(maxval(val(:,:,k)),real64)
        end do
        write(u) real(x*sL, real32)
        write(u) real(y*sL, real32)
        allocate(pad(nx+1, ny))
        pad = 0.0_real32
        do k = 1,ncellvar
            pad(1:nx, :) = real(val(:,:,k), real32)
            write(u) pad
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

    end subroutine write_plt

    !> 셀 줄 jc 하나를 따라 xc, rho, u, v, p 출력 (1D 그래프용)
    subroutine write_line(fname, Q, jc)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: Q(:,1-ng:,1-ng:)
        integer,  intent(in) :: jc

        integer  :: u, ic
        real(dp) :: W(nvar), sL, sR, sU, sP, sT
        character(len=12) :: names(7)

        call out_scales(sL, sR, sU, sP, sT, names)
        open(newunit=u, file=fname, status='replace', action='write')
        if (dim_out) then
            write(u,'(a)') '#  xc [m]  rho [kg/m3]  u [m/s]  v [m/s]  p [Pa]'
        else
            write(u,'(a)') '#  xc  rho  u  v  p'
        end if
        do ic = 1,nx
            call c2p(Q(:,ic,jc), W)
            write(u,'(5es16.8)') xc(ic,jc)*sL, W*[sR, sU, sU, sP]
        end do
        close(u)
    end subroutine write_line


    !> 재시작 파일 쓰기 (배정밀도, 내부 셀 Q 전부 + 반복 수, 시간, 잔차 정규화 기준)
    !>   먼저 .tmp 로 쓰고 이름을 바꿔서, 쓰는 도중에 멈춰도 이전 재시작 파일이 남게 함
    !>   dual time (BDF2) 이면 Qold (Q^(n-1), 내부 셀), dt_old (지난 스텝 dt), iter (안쪽 반복 누적 수) 도 뒤에 붙여 씀
    !>   -> 이어서 계산해도 같은 결과
    subroutine write_restart(fname, Q, step, t, res0, Qold, dt_old, iter)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: Q(:,1-ng:,1-ng:)
        integer,  intent(in) :: step
        real(dp), intent(in) :: t, res0(nvar)
        real(dp), intent(in), optional :: Qold(:,:,:), dt_old
        integer,  intent(in), optional :: iter
        integer :: u, ios

        open(newunit=u, file=fname//'.tmp', access='stream', form='unformatted', status='replace', action='write')
        write(u) restart_magic, int(nx, int32), int(ny, int32)
        write(u) int(step, int32), real(t, real64), real(res0, real64)
        write(u) real(Q(:,1:nx,1:ny), real64)
        if (present(Qold)) write(u) 2_int32, int(iter, int32), real(dt_old, real64), real(Qold, real64)
        close(u)

        ! 이전 파일을 지우고 .tmp 를 제 이름으로
        open(newunit=u, file=fname, status='old', iostat=ios)
        if (ios == 0) close(u, status='delete')
        call rename(fname//'.tmp', fname)
    end subroutine write_restart

    !> 재시작 파일 읽기: 셀 수가 지금 격자와 같아야 함
    !>   Qold, dt_old, have_old, iter (선택): dual time 의 Q^(n-1), 지난 dt, 안쪽 반복 누적 수. 파일에 없으면 have_old = .false.
    subroutine read_restart(fname, Q, step, t, res0, Qold, dt_old, have_old, iter)
        character(len=*), intent(in)    :: fname
        real(dp),         intent(inout) :: Q(:,1-ng:,1-ng:)
        integer,          intent(out)   :: step
        real(dp),         intent(out)   :: t, res0(nvar)
        real(dp), intent(inout), optional :: Qold(:,:,:), dt_old
        logical,  intent(out),   optional :: have_old
        integer,  intent(inout), optional :: iter
        integer :: u, ios
        integer(int32) :: d(2), st, flag, it
        character(len=16) :: magic

        open(newunit=u, file=fname, access='stream', form='unformatted', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(a,a,a)') ' Error: restart = 1 but cannot open ', trim(fname), ' (run once without restart first)'
            error stop
        end if
        read(u, iostat=ios) magic, d
        if (ios /= 0 .or. magic /= restart_magic) then
            write(*,'(a,a,a)') ' Error: ', trim(fname), ' is not a JFLOW 2D restart file'
            error stop
        end if
        if (any(d /= [nx, ny])) then
            write(*,'(a,2(1x,i0),a,2(1x,i0))') ' Error: restart file cells', d, ' but the grid has', nx, ny
            error stop
        end if
        read(u) st, t, res0
        step = st
        read(u, iostat=ios) Q(:,1:nx,1:ny)
        if (ios /= 0) then
            write(*,'(a,a,a)') ' Error: ', trim(fname), ' ended early'
            error stop
        end if
        if (present(have_old)) then
            have_old = .false.
            read(u, iostat=ios) flag
            if (ios == 0 .and. flag == 1) then             ! (iter 없는 형식)
                read(u, iostat=ios) dt_old, Qold
                have_old = (ios == 0)
            else if (ios == 0 .and. flag == 2) then
                read(u, iostat=ios) it, dt_old, Qold
                have_old = (ios == 0)
                if (have_old .and. present(iter)) iter = it
            end if
        end if
        close(u)
    end subroutine read_restart

end module output
