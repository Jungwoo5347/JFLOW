!> 케이스 입력 읽기 (KFLOW 스타일) - 3D
!>
!> 실행:  .\main3d.exe cases3d\sod        (인자 없으면 현재 폴더)
!>   케이스 폴더 안의 두 파일을 읽음
!>     grid.in : 격자 파일과 경계조건 (KFLOW grid.in 형식)
!>     flow.in : 계산 제어, 기체, 유동 상태, 초기조건
!>   결과도 케이스 폴더에 저장:  <output>_result.plt (또는 .dat), _line.dat, _residual.dat, _restart.bin
!>
!> 공통 형식 (KFLOW):  $key = value   ( # 뒤는 주석, 빈 줄 무시, 대소문자 무시, $ 는 생략 가능 )
!>
!> ---- grid.in ----------------------------------------------------------------
!>   dim      = 3                             3D 솔버는 3 만 받음 (2D 는 main.exe)
!>   nblock   = 2                             블록 수 (생략하면 block-n 줄과 격자 파일에서 셈)
!>   block-n  = mesh.x                        블록 n 의 PLOT3D 격자 파일 (케이스 폴더 기준 경로)
!>   block-n  = mesh.x m                      다중 블록 파일의 m 번째 블록을 블록 n 으로
!>   block-1  = mesh.x                        block 줄이 이것 하나뿐이고 파일에 블록이 여러 개면 전부 순서대로 사용
!>   block-n  = cartesian ni nj nk xmin xmax ymin ymax zmin zmax   균일 직교 격자 생성 (노드 수)
!>   form-n   = 0                             0: unformatted (바이너리), 1: ascii   (없으면 form-1, 그것도 없으면 0)
!>   extrude-n = nk span                      2D PLOT3D 블록을 z 방향 (0 ~ span) 으로 nk 노드만큼 밀어내 3D 로
!>                                            (없으면 extrude-1. 2D 결과와 비교할 때: 면 5, 6 을 대칭면 300 으로)
!>   bc-n-f   = code                          블록 n, 면 f 전체
!>   bc-n-f   = code a1 a2 b1 b2              면 f 의 노드 구간 (1부터 센 노드 번호, 면 위 두 방향)
!>   bc-n-f   = code [a1 a2 b1 b2] state      입구/원방이면 flow.in 의 상태 이름 (생략 시 freestream)
!>     면 번호 f : 1 = imin, 2 = imax, 3 = jmin, 4 = jmax, 5 = kmin, 6 = kmax
!>     범위의 두 방향 (a, b) : 면 1,2 는 (j, k),  면 3,4 는 (i, k),  면 5,6 은 (i, j)
!>     같은 면을 여러 줄로 나눠 쓸 수 있음, 뒤에 쓴 줄이 우선
!>     경계조건이 없는 면은 네 꼭짓점이 겹치는 다른 경계 면 (같은 블록이든 다른 블록이든) 을 자동으로 찾아 연결
!>     (블록 사이 면, C-grid wake cut, O-grid 이음매). 짝이 없으면 에러
!>   output   = name                          결과 파일 이름 머리
!>   decompose-n, ncpu-n, diswall             KFLOW 병렬/난류용 항목: 읽기만 하고 사용 안 함
!>
!>   경계조건 코드
!>     110 : 벽 (비점성 미끄럼 벽)
!>     200 : 원방 (리만 불변량, 기본 freestream)
!>     210 : 초음속 입구 (상태 고정)
!>     220 : 초음속 출구 (내부 값 복사)
!>     300 : 대칭면
!>
!> ---- flow.in ----------------------------------------------------------------
!>   2D (src/case_setup.f90) 와 같은 항목에, 다른 점만:
!>   state-name     = rho u v w p              속도 성분 3 개 (숫자 5 개)
!>   state-name     = mach M alpha a           또는 마하수와 받음각 [deg] 으로 (키워드, 순서 상관없음)
!>                    [beta b] [rho r] [p p]     옆미끄럼각 [deg] (기본 0), 밀도, 압력 (기본 rho = 1, p = 1/gamma: 음속 1)
!>                    speed V                    마하수 대신 속도 크기
!>                    속도: u = V cos a cos b,  (lift_axis 성분) = V sin a cos b,  (나머지 축 성분) = V sin b
!>   lift_axis      = y                        받음각 a 가 기울이는 양력 방향 축: y (기본) | z
!>   region-k       = name x0 x1 y0 y1 z0 z1   상자 안의 셀만 그 상태로 (* 는 무한대, k 순서대로 = 뒤가 우선)
!>   line_block     = 1                        <output>_line.dat 으로 뽑을 i 방향 셀 줄의 블록 번호
!>   line_row       = 1                        그 줄의 j 셀 번호
!>   line_k         = 1                        그 줄의 k 셀 번호
!>   entropy_fix    = 0                        Roe 엔트로피 보정 (Harten) 계수: delta = entropy_fix * c~ 에 대해
!>                                             |lam| < delta 이면 |lam| -> (lam^2 + delta^2)/(2 delta).  0 = 끔, 보통 0.05 ~ 0.2
!>   entropy_fix_waves = acoustic              보정할 고유값: acoustic (un -+ c, 음속점 팽창 충격파 방지) | all (un 도: 정체점)
!>   time_scheme    = explicit                 시간 적분: explicit (SSP-RK3) | implicit (LU-SGS, 정상 문제만; cfl 을 크게, 예 5 ~ 50)
!>   lusgs_omega    = 1.0                      LU-SGS 대각 계수 w (1 ~ 2, 클수록 안정, 수렴은 느려짐)
!>   output_format  = plt                      결과 파일: plt (Tecplot 바이너리) | ascii (Tecplot .dat) | both
!>   save_interval  = 0                        N 회마다 결과와 재시작 파일 저장 (0: 끝날 때만)
!>   restart        = 0                        1: <output>_restart.bin 에서 이어서 계산
!>                                             (반복 수와 시간은 이어서 셈: max_iter, t_end 는 전체 기준)
!>   계산 중 케이스 폴더에 STOP 이라는 파일을 만들면, 그 반복을 마치고 결과를 저장한 뒤 끝냄
module case_setup
    use euler
    implicit none
    private

    ! ---- 경계조건 코드 (KFLOW 와 맞출 때는 여기 숫자만 바꾸면 됨) ----
    integer, parameter, public :: BC_unset    = 0
    integer, parameter, public :: BC_wall     = 110
    integer, parameter, public :: BC_farfield = 200
    integer, parameter, public :: BC_inlet    = 210
    integer, parameter, public :: BC_outlet   = 220
    integer, parameter, public :: BC_symmetry = 300
    integer, parameter, public :: BC_cut      = -1      ! 내부용: 자동으로 찾은 맞닿은 면 (블록 사이, wake cut 등)

    integer, parameter, public :: nside = 6             ! 경계 변 수: imin, imax, jmin, jmax, kmin, kmax
    integer, parameter, public :: max_states = 50, max_regions = 50, max_patches = 2000, max_specs = 500

    ! ---- 케이스 ----
    character(len=256), public :: case_dir = '.'
    character(len=512), public :: output_name = 'jflow'

    ! ---- 계산 제어 (flow.in) ----
    logical,  public :: steady = .true.
    logical,  public :: local_dt = .true.            ! 정상 문제에서 국소 시간 간격 사용
    real(dp), public :: t_end = 1.0_dp
    integer,  public :: max_steps = 100000
    real(dp), public :: res_tol = 1.0e-6_dp
    real(dp), public :: cfl = 0.4_dp
    integer,  public :: print_interval = 10
    integer,  public :: line_block = 1, line_row = 1, line_k = 1
    real(dp), public :: sigma = 1.0_dp, kappa = 1.0_dp/3.0_dp, b = 4.0_dp
    real(dp), public :: entropy_fix = 0.0_dp         ! Roe 엔트로피 보정 계수 (0: 끔)
    integer,  public :: lift_axis = 2                 ! 받음각의 양력 방향 축 (2: y, 3: z)
    logical,  public :: entropy_fix_all = .false.    ! .true.: 모든 고유값, .false.: 음향파만
    logical,  public :: implicit_ts = .false.        ! time_scheme = implicit (LU-SGS)
    real(dp), public :: lusgs_omega = 1.0_dp         ! LU-SGS 대각 계수
    logical,  public :: restart = .false.            ! 재시작 파일에서 이어서 계산
    integer,  public :: save_interval = 0            ! 중간 저장 간격 (0: 끝날 때만)
    logical,  public :: out_plt = .true., out_ascii = .false.   ! 결과 파일 형식 (output_format)

    ! 리만 솔버 선택 (flow.in 의 flux)
    integer, parameter, public :: FLUX_RUSANOV = 1, FLUX_ROE = 2
    integer,           public :: flux_id   = FLUX_RUSANOV
    character(len=16), public :: flux_name = 'rusanov'

    ! ---- 격자 (grid.in): block-n 줄 하나 = spec 하나 ----
    type, public :: gridspec_t
        character(len=16)  :: kind = 'plot3d'       ! 'plot3d' | 'cartesian'
        character(len=512) :: file = ''
        integer            :: form = 0              ! 0: unformatted, 1: ascii
        integer            :: index = 0             ! 파일 안 블록 번호 (0: 지정 안 함)
        integer            :: ni = 0, nj = 0, nk = 0
        real(dp)           :: x0 = 0, x1 = 1, y0 = 0, y1 = 1, z0 = 0, z1 = 1
        integer            :: extrude_nk = 0        ! 2D 격자 밀어내기 노드 수 (0: 안 함)
        real(dp)           :: extrude_span = 1.0_dp
    end type
    integer,          public :: nspecs = 0
    integer,          public :: nblock_in = 0       ! grid.in 의 nblock (0: 지정 안 함)
    type(gridspec_t), public :: specs(max_specs)

    ! ---- 유동 상태, 초기조건 ----
    type, public :: state_t
        character(len=64) :: name
        real(dp)          :: W(nvar)            ! rho, u, v, w, p
    end type
    integer,       public :: nstates = 0
    type(state_t), public :: states(max_states)

    type, public :: region_t
        integer  :: istate
        real(dp) :: x0, x1, y0, y1, z0, z1
    end type
    integer,        public :: init_state = 0
    integer,        public :: nregions = 0
    type(region_t), public :: regions(max_regions)

    ! ---- 경계조건 (grid.in 의 bc-n-f 한 줄 = patch 하나) ----
    type, public :: patch_t
        integer :: blk                          ! 블록 번호
        integer :: side                         ! 1 imin, 2 imax, 3 jmin, 4 jmax, 5 kmin, 6 kmax
        integer :: a1, a2, b1, b2               ! 노드 번호 구간 (면 위 두 방향, a2/b2 = 0 이면 면 끝까지)
        integer :: btype
        integer :: istate                       ! 입구/원방 상태 (0 이면 없음)
    end type
    integer,       public :: npatches = 0
    type(patch_t), public :: patches(max_patches)

    public :: read_case, print_case, side_name, bc_name, case_path

    ! ---- key = value 저장 ----
    integer, parameter :: max_kv = 5000
    type :: kv_t
        character(len=64)  :: key
        character(len=512) :: val
        integer            :: line
        logical            :: used
    end type
    integer :: nkv
    type(kv_t) :: kv(max_kv)
    character(len=512) :: cur_file

contains

    !> 케이스 폴더의 grid.in, flow.in 읽기
    subroutine read_case(dir)
        character(len=*), intent(in) :: dir

        case_dir = dir
        if (len_trim(case_dir) == 0) case_dir = '.'
        call read_flow_in(case_path('flow.in'))
        call read_grid_in(case_path('grid.in'))
    end subroutine read_case

    !> 케이스 폴더 기준 경로 (절대 경로는 그대로)
    function case_path(fname) result(path)
        character(len=*), intent(in) :: fname
        character(len=512) :: path
        if (is_absolute(fname)) then
            path = fname
        else
            path = trim(case_dir)//'/'//trim(fname)
        end if
    end function case_path

    ! ======================================================================
    !  flow.in
    ! ======================================================================
    subroutine read_flow_in(fname)
        character(len=*), intent(in) :: fname
        integer :: k, m, nt
        character(len=64) :: tok(16)
        character(len=16) :: fmt
        real(dp) :: lim(6)
        integer, allocatable :: order(:)
        integer :: nreg_keys

        call load_kv(fname)

        steady         = get_int('steady', 1) /= 0
        local_dt       = get_int('local_dt', 1) /= 0
        max_steps      = get_int('max_iter', max_steps)
        res_tol        = get_real('res_tol', res_tol)
        t_end          = get_real('t_end', t_end)
        cfl            = get_real('cfl', cfl)
        gam            = get_real('gamma', gam)
        sigma          = get_real('sigma', sigma)
        kappa          = get_real('kappa', kappa)
        b              = get_real('b', b)
        flux_name      = lower(adjustl(get_str('flux', flux_name)))
        select case (flux_name)
        case ('rusanov');  flux_id = FLUX_RUSANOV
        case ('roe');      flux_id = FLUX_ROE
        case default;      call fail_kv('flux = '//trim(flux_name)//' : unknown (use rusanov or roe)')
        end select
        if (b < 1.0_dp) call fail_kv('b must be >= 1')
        if (kappa < 1.0_dp) then
            if (b > (3.0_dp - kappa)/(1.0_dp - kappa) + 1.0e-12_dp) call fail_kv('b must be <= (3-kappa)/(1-kappa)')
        end if
        entropy_fix    = get_real('entropy_fix', entropy_fix)
        if (entropy_fix < 0.0_dp) call fail_kv('entropy_fix must be >= 0')
        fmt            = lower(adjustl(get_str('entropy_fix_waves', 'acoustic')))
        select case (fmt)
        case ('acoustic');  entropy_fix_all = .false.
        case ('all');       entropy_fix_all = .true.
        case default;       call fail_kv('entropy_fix_waves = '//trim(fmt)//' : unknown (use acoustic or all)')
        end select
        fmt            = lower(adjustl(get_str('time_scheme', 'explicit')))
        select case (fmt)
        case ('explicit', 'rk3');  implicit_ts = .false.
        case ('implicit', 'lusgs', 'lu-sgs');  implicit_ts = .true.
        case default;  call fail_kv('time_scheme = '//trim(fmt)//' : unknown (use explicit or implicit)')
        end select
        if (implicit_ts .and. .not. steady) call fail_kv('time_scheme = implicit (LU-SGS) is for steady problems only (steady = 1)')
        lusgs_omega    = get_real('lusgs_omega', lusgs_omega)
        if (lusgs_omega < 1.0_dp) call fail_kv('lusgs_omega must be >= 1')
        print_interval = get_int('print_interval', print_interval)
        line_block     = get_int('line_block', line_block)
        line_row       = get_int('line_row', line_row)
        line_k         = get_int('line_k', line_k)
        restart        = get_int('restart', 0) /= 0
        save_interval  = get_int('save_interval', save_interval)
        fmt            = lower(adjustl(get_str('output_format', 'plt')))
        select case (fmt)
        case ('plt');    out_plt = .true.;   out_ascii = .false.
        case ('ascii');  out_plt = .false.;  out_ascii = .true.
        case ('both');   out_plt = .true.;   out_ascii = .true.
        case default;    call fail_kv('output_format = '//trim(fmt)//' : unknown (use plt, ascii or both)')
        end select
        if (cfl <= 0) call fail_kv('cfl must be > 0')
        if (print_interval < 1) call fail_kv('print_interval must be >= 1')
        if (save_interval < 0) call fail_kv('save_interval must be >= 0')

        ! 받음각 기준 축: 양력 방향 (y 또는 z), 나머지 축이 옆 방향
        fmt = lower(adjustl(get_str('lift_axis', 'y')))
        select case (fmt)
        case ('y');  lift_axis = 2
        case ('z');  lift_axis = 3
        case default;  call fail_kv('lift_axis = '//trim(fmt)//' : use y or z')
        end select

        ! state-<name> = rho u v w p                          (숫자 5 개)
        !   또는         mach M | speed V,  alpha a [deg], beta b [deg], rho r, p p   (키워드, 순서 상관없음)
        do k = 1,nkv
            if (index(kv(k)%key, 'state-') /= 1) cycle
            kv(k)%used = .true.
            call split_tokens(kv(k)%val, tok, nt)
            if (nstates == max_states) call fail_line(k, 'too many states')
            nstates = nstates + 1
            states(nstates)%name = kv(k)%key(7:)
            if (is_number(tok(1))) then
                if (nt /= 5) call fail_line(k, 'state needs 5 numbers: rho u v w p  (or: mach M alpha a ...)')
                do m = 1,5
                    states(nstates)%W(m) = to_real(tok(m), k)
                end do
            else
                states(nstates)%W = state_from_angles(k, tok, nt)
            end if
            if (states(nstates)%W(1) <= 0 .or. states(nstates)%W(5) <= 0) &
                call fail_line(k, 'rho and p must be > 0')
        end do
        if (nstates == 0) call fail_kv('no state-<name> lines (rho u v w p, or mach M alpha a)')

        init_state = need_state(get_str('initial', ''), 'initial')

        ! region-<k> = name x0 x1 y0 y1 z0 z1   (k 순서대로)
        nreg_keys = 0
        allocate(order(nkv))
        do k = 1,nkv
            if (index(kv(k)%key, 'region-') == 1) then
                nreg_keys = nreg_keys + 1
                order(nreg_keys) = k
            end if
        end do
        call sort_by_suffix(order(1:nreg_keys), 8)
        do k = 1,nreg_keys
            associate (e => kv(order(k)))
                e%used = .true.
                call split_tokens(e%val, tok, nt)
                if (nt /= 7) call fail_line(order(k), 'region needs: state x0 x1 y0 y1 z0 z1  (* = unbounded)')
                do m = 1,6
                    lim(m) = to_bound(tok(m+1), merge(-1, 1, mod(m,2) == 1), order(k))
                end do
                if (nregions == max_regions) call fail_line(order(k), 'too many regions')
                nregions = nregions + 1
                regions(nregions) = region_t(need_state(tok(1), trim(e%key)), lim(1), lim(2), lim(3), lim(4), lim(5), lim(6))
            end associate
        end do

        call warn_unused([character(len=16) :: ''])
    end subroutine read_flow_in

    !> state 줄의 키워드 형식 -> W = (rho, u, v, w, p)
    !>   mach M | speed V (둘 중 하나),  alpha a, beta b [deg] (기본 0),  rho r, p p (기본 rho = 1, p = 1/gamma)
    function state_from_angles(k, tok, nt) result(W)
        integer,          intent(in) :: k, nt
        character(len=*), intent(in) :: tok(:)
        real(dp) :: W(nvar)
        real(dp) :: mach, speed, alpha, beta, rho, p, vmag, deg
        integer  :: m, side
        mach = -1;  speed = -1;  alpha = 0;  beta = 0;  rho = 1.0_dp;  p = 1.0_dp/gam
        if (mod(nt, 2) /= 0) call fail_line(k, 'expected keyword value pairs: mach M alpha a [beta b] [rho r] [p p]')
        do m = 1, nt, 2
            select case (lower(tok(m)))
            case ('mach', 'm');               mach  = to_real(tok(m+1), k)
            case ('speed', 'v');              speed = to_real(tok(m+1), k)
            case ('alpha', 'aoa');            alpha = to_real(tok(m+1), k)
            case ('beta', 'sideslip');        beta  = to_real(tok(m+1), k)
            case ('rho');                      rho   = to_real(tok(m+1), k)
            case ('p');                        p     = to_real(tok(m+1), k)
            case default
                call fail_line(k, 'unknown keyword '''//trim(tok(m))//''' (use mach, speed, alpha, beta, rho, p)')
            end select
        end do
        if ((mach >= 0) .eqv. (speed >= 0)) call fail_line(k, 'give exactly one of mach or speed')
        if (rho <= 0 .or. p <= 0) call fail_line(k, 'rho and p must be > 0')
        if (mach >= 0) then
            vmag = mach*sqrt(gam*p/rho)
        else
            vmag = speed
        end if
        deg  = acos(-1.0_dp)/180.0_dp
        side = 5 - lift_axis                                 ! y 가 양력이면 옆은 z, z 면 y
        W = 0.0_dp
        W(1) = rho;  W(5) = p
        W(2)           = vmag*cos(alpha*deg)*cos(beta*deg)
        W(1+lift_axis) = vmag*sin(alpha*deg)*cos(beta*deg)
        W(1+side)      = vmag*sin(beta*deg)
    end function state_from_angles

    ! ======================================================================
    !  grid.in
    ! ======================================================================
    subroutine read_grid_in(fname)
        character(len=*), intent(in) :: fname
        integer :: k, n, nt, face, code, istate, nnum, r(4), nmax, p
        character(len=512) :: fval
        character(len=64) :: tok(16)
        character(len=16) :: key
        integer, allocatable :: kblk(:)

        call load_kv(fname)

        if (get_int('dim', 3) /= 3) call fail_kv('main3d.exe needs dim = 3 (use main.exe for 2D cases)')
        nblock_in   = get_int('nblock', 0)
        output_name = get_str('output', 'jflow')

        ! block-<n> 줄들: n 은 1, 2, ... 빠짐없이
        nmax = 0
        do k = 1,nkv
            if (index(kv(k)%key, 'block-') /= 1) cycle
            n = to_int(kv(k)%key(7:), k)
            if (n < 1 .or. n > max_specs) call fail_line(k, 'block number out of range')
            nmax = max(nmax, n)
        end do
        if (nmax == 0) call fail_kv('block-1 = <grid file> is missing')
        allocate(kblk(nmax));  kblk = 0
        do k = 1,nkv
            if (index(kv(k)%key, 'block-') /= 1) cycle
            kv(k)%used = .true.
            kblk(to_int(kv(k)%key(7:), k)) = k
        end do
        if (any(kblk == 0)) call fail_kv('block-n lines must be numbered 1, 2, ... without gaps')
        nspecs = nmax

        do n = 1,nspecs
            k = kblk(n)
            call split_tokens(kv(k)%val, tok, nt)
            associate (s => specs(n))
                if (lower(tok(1)) == 'cartesian') then
                    if (nt /= 10) call fail_line(k, 'cartesian needs: ni nj nk xmin xmax ymin ymax zmin zmax')
                    s%kind = 'cartesian'
                    s%ni = to_int(tok(2), k);   s%nj = to_int(tok(3), k);   s%nk = to_int(tok(4), k)
                    s%x0 = to_real(tok(5), k);  s%x1 = to_real(tok(6), k)
                    s%y0 = to_real(tok(7), k);  s%y1 = to_real(tok(8), k)
                    s%z0 = to_real(tok(9), k);  s%z1 = to_real(tok(10), k)
                    if (s%ni < 2 .or. s%nj < 2 .or. s%nk < 2) call fail_line(k, 'ni, nj, nk (node counts) must be >= 2')
                else
                    ! 값 전체가 파일 경로 (공백 포함 가능), 끝에 숫자가 하나 더 있으면 파일 안 블록 번호
                    s%kind = 'plot3d'
                    fval = kv(k)%val
                    if (nt >= 2 .and. verify(trim(tok(nt)), '0123456789') == 0) then
                        s%index = to_int(tok(nt), k)
                        if (s%index < 1) call fail_line(k, 'block index in file must be >= 1')
                        fval = fval(:index(trim(fval), trim(tok(nt)), back=.true.) - 1)
                    end if
                    s%file = case_path(trim(fval))
                    write(key,'(a,i0)') 'form-', n
                    s%form = get_int(trim(key), get_int('form-1', 0))
                    if (s%form /= 0 .and. s%form /= 1) call fail_kv(trim(key)//' must be 0 (unformatted) or 1 (ascii)')
                end if

                ! extrude-n = nk span   (2D PLOT3D 블록을 3D 로, 없으면 extrude-1)
                write(key,'(a,i0)') 'extrude-', n
                p = find_key(trim(key))
                if (p == 0) p = find_key('extrude-1')
                if (p > 0 .and. s%kind == 'plot3d') then
                    kv(p)%used = .true.
                    call split_tokens(kv(p)%val, tok, nt)
                    if (nt /= 2) call fail_line(p, 'extrude needs: nk span')
                    s%extrude_nk   = to_int(tok(1), p)
                    s%extrude_span = to_real(tok(2), p)
                    if (s%extrude_nk < 2 .or. s%extrude_span <= 0.0_dp) call fail_line(p, 'need nk >= 2 and span > 0')
                end if
            end associate
        end do
        do k = 1,nkv                                       ! form-n, extrude-n 은 위에서 씀
            if (index(kv(k)%key, 'form-') == 1 .or. index(kv(k)%key, 'extrude-') == 1) kv(k)%used = .true.
        end do

        ! bc-<n>-<face> = code [a1 a2 b1 b2] [state]   (파일에 쓴 순서대로 적용 = 뒤가 우선)
        !   블록 번호가 실제 블록 수 안에 있는지는 격자를 읽은 뒤 (flow 의 setup_bc) 에서 검사
        do k = 1,nkv
            if (index(kv(k)%key, 'bc-') /= 1) cycle
            kv(k)%used = .true.
            p = index(kv(k)%key(4:), '-') + 3                ! 블록 번호 뒤의 '-' 위치
            if (p <= 4 .or. len_trim(kv(k)%key) /= p + 1) call fail_line(k, 'expected bc-<block>-<face>, face = 1..6')
            n    = to_int(kv(k)%key(4:p-1), k)
            face = to_int(kv(k)%key(p+1:p+1), k)
            if (n < 1) call fail_line(k, 'block number must be >= 1')
            if (face < 1 .or. face > 6) &
                call fail_line(k, 'face must be 1 (imin), 2 (imax), 3 (jmin), 4 (jmax), 5 (kmin) or 6 (kmax)')

            call split_tokens(kv(k)%val, tok, nt)
            if (nt < 1) call fail_line(k, 'boundary condition code is missing')
            code = to_int(tok(1), k)
            if (bc_name(code) == '?') call fail_line(k, 'unknown boundary condition code '//trim(tok(1)))

            ! 코드 뒤: 숫자 0개 또는 4개 (노드 범위 두 방향), 그다음 선택적으로 상태 이름
            nnum = 0
            do while (1 + nnum + 1 <= nt)
                if (.not. is_number(tok(2+nnum))) exit
                nnum = nnum + 1
            end do
            r = [1, 0, 1, 0]
            if (nnum == 4) then
                r = [to_int(tok(2), k), to_int(tok(3), k), to_int(tok(4), k), to_int(tok(5), k)]
            else if (nnum /= 0) then
                call fail_line(k, 'range must be four node numbers: a1 a2 b1 b2')
            end if
            if (nt > 1 + nnum + 1) call fail_line(k, 'too many items')

            istate = 0
            if (code == BC_inlet .or. code == BC_farfield) then
                if (nt == 1 + nnum + 1) then
                    istate = need_state(tok(nt), trim(kv(k)%key))
                else
                    istate = need_state('freestream', trim(kv(k)%key)//' (default state)')
                end if
            else if (nt == 1 + nnum + 1) then
                call fail_line(k, 'a state name is only used with codes 200 and 210')
            end if

            if (npatches == max_patches) call fail_line(k, 'too many bc lines')
            npatches = npatches + 1
            patches(npatches) = patch_t(n, face, r(1), r(2), r(3), r(4), code, istate)
        end do

        ! KFLOW 에서 쓰는 항목: 읽기만 함
        call warn_unused([character(len=16) :: 'decompose-', 'ncpu-', 'diswall'])
    end subroutine read_grid_in

    !> 읽은 설정 요약
    subroutine print_case()
        integer :: k
        character(len=24) :: range

        write(*,'(a)')   ' ================================================================'
        write(*,'(a,a,a)') '  case     : ', trim(case_dir), '   (3D)'
        do k = 1,nspecs
            associate (s => specs(k))
                if (s%kind == 'cartesian') then
                    write(*,'(a,i3,a,i0,a,i0,a,i0,a)') '  block-', k, ': cartesian ', s%ni, ' x ', s%nj, ' x ', s%nk, ' nodes'
                else
                    write(*,'(a,i3,a,a,a)', advance='no') '  block-', k, ': ', trim(s%file), &
                        merge('  (unformatted)', '  (ascii)      ', s%form == 0)
                    if (s%index > 0) write(*,'(a,i0)', advance='no') '  block ', s%index
                    if (s%extrude_nk > 0) write(*,'(a,i0,a,es9.2)', advance='no') '  extrude nk = ', s%extrude_nk, &
                        ', span = ', s%extrude_span
                    write(*,*)
                end if
            end associate
        end do
        if (steady) then
            if (res_tol > 0.0_dp) then
                write(*,'(a,es9.2,a,i0)') '  steady   : res_tol = ', res_tol, ',  max_iter = ', max_steps
            else
                write(*,'(a,i0)') '  steady   : res_tol = off (run to max_iter),  max_iter = ', max_steps
            end if
            write(*,'(a,a)') '  time step: ', merge('local (each cell its own dt)', 'global (min dt)             ', local_dt)
        else
            write(*,'(a,es12.4)')     '  unsteady : t_end = ', t_end
        end if
        write(*,'(a,f6.3,a,f6.3)')    '  CFL      : ', cfl, '    gamma = ', gam
        if (implicit_ts) then
            write(*,'(a,f5.2,a)')     '  time int : implicit LU-SGS  (omega = ', lusgs_omega, ')'
        else
            write(*,'(a)')            '  time int : explicit SSP-RK3'
        end if
        write(*,'(a,f6.3,a,f7.4,a,f6.3)') '  MUSCL    : sigma = ', sigma, '    kappa = ', kappa, '    b = ', b
        if (flux_id == FLUX_ROE .and. entropy_fix > 0.0_dp) then
            write(*,'(a,a,a,f6.3,a,a,a)') '  flux     : ', trim(flux_name), '   entropy fix delta = ', entropy_fix, &
                ' c  (', merge('all waves     ', 'acoustic waves', entropy_fix_all), ')'
        else
            write(*,'(a,a)')          '  flux     : ', trim(flux_name)
        end if
        write(*,'(a)') '  states   :         rho           u           v           w           p'
        do k = 1,nstates
            write(*,'(4x,a12,5es12.4)') states(k)%name, states(k)%W
        end do
        write(*,'(a,a,a,i0,a)') '  initial  : ', trim(states(init_state)%name), '  (+ ', nregions, ' region(s))'
        if (restart) write(*,'(a,a)') '  restart  : ', trim(case_path(trim(output_name)//'_restart.bin'))
        write(*,'(a)') '  bc       : block face        nodes (a ; b)           code  type       state'
        do k = 1,npatches
            associate (p => patches(k))
                if (p%a2 == 0 .and. p%b2 == 0) then
                    range = 'all'
                else
                    write(range,'(i0,a,a,a,i0,a,a)') p%a1, '-', trim(endnum(p%a2)), ' ; ', p%b1, '-', trim(endnum(p%b2))
                end if
                write(*,'(13x,i5,2x,i1,1x,a4,2x,a22,1x,i4,2x,a9)', advance='no') p%blk, p%side, side_name(p%side), &
                    range, p%btype, bc_name(p%btype)
                if (p%istate > 0) then
                    write(*,'(2x,a)') trim(states(p%istate)%name)
                else
                    write(*,*)
                end if
            end associate
        end do
        write(*,'(a,a,a,a,a)') '  output   : ', trim(case_path(trim(output_name)//'_*')), &
            merge('  (result: plt', '  (result:    ', out_plt), merge(' ascii)', ')      ', out_ascii)
        if (save_interval > 0) write(*,'(a,i0,a)') '  save     : every ', save_interval, ' iterations (result + restart)'
        write(*,'(a)')   ' ================================================================'

    contains

        function endnum(n) result(s)
            integer, intent(in) :: n
            character(len=8) :: s
            if (n == 0) then
                s = 'end'
            else
                write(s,'(i0)') n
            end if
        end function endnum

    end subroutine print_case

    function side_name(id) result(s)
        integer, intent(in) :: id
        character(len=4) :: s
        character(len=4), parameter :: names(6) = ['imin', 'imax', 'jmin', 'jmax', 'kmin', 'kmax']
        s = names(id)
    end function side_name

    !> 경계조건 코드 -> 이름  (모르는 코드는 '?')
    function bc_name(code) result(s)
        integer, intent(in) :: code
        character(len=9) :: s
        select case (code)
        case (BC_wall);      s = 'wall'
        case (BC_farfield);  s = 'farfield'
        case (BC_inlet);     s = 'inlet'
        case (BC_outlet);    s = 'outlet'
        case (BC_symmetry);  s = 'symmetry'
        case (BC_cut);       s = 'cut(auto)'
        case (BC_unset);     s = 'unset'
        case default;        s = '?'
        end select
    end function bc_name

    ! ======================================================================
    !  key = value 파일 읽기와 도우미
    ! ======================================================================
    subroutine load_kv(fname)
        character(len=*), intent(in) :: fname
        integer :: u, ios, p, c, ln
        character(len=1024) :: line

        cur_file = fname
        nkv = 0
        open(newunit=u, file=fname, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: cannot open ', trim(fname)
            error stop
        end if
        ln = 0
        do
            read(u,'(a)',iostat=ios) line
            if (ios /= 0) exit
            ln = ln + 1
            if (ln == 1 .and. line(1:3) == char(239)//char(187)//char(191)) line = line(4:)   ! UTF-8 BOM (메모장 저장)
            line = replace_tabs(line)
            c = scan(line, '#!')                      ! 주석 자르기
            if (c > 0) line = line(:c-1)
            line = adjustl(line)
            if (line(1:1) == '$') line = line(2:)     ! KFLOW 형식: 입력 줄은 $ 로 시작
            if (len_trim(line) == 0) cycle
            p = index(line, '=')
            if (p == 0) then
                write(*,'(a,a,a,i0,a)') ' Error: ', trim(fname), ' line ', ln, ': expected  key = value'
                write(*,'(a,a)') '   ', trim(line)
                error stop
            end if
            if (nkv == max_kv) then
                write(*,*) 'Error: too many lines in ', trim(fname)
                error stop
            end if
            nkv = nkv + 1
            kv(nkv)%key  = lower(adjustl(line(:p-1)))
            kv(nkv)%val  = adjustl(line(p+1:))
            kv(nkv)%line = ln
            kv(nkv)%used = .false.
        end do
        close(u)
    end subroutine load_kv

    !> 같은 key 가 여러 번 있으면 마지막 것 (bc-, state-, region- 은 따로 처리)
    integer function find_key(key)
        character(len=*), intent(in) :: key
        integer :: k
        find_key = 0
        do k = 1,nkv
            if (kv(k)%key == key) find_key = k
        end do
    end function find_key

    function get_str(key, default) result(s)
        character(len=*), intent(in) :: key, default
        character(len=512) :: s
        integer :: k
        k = find_key(key)
        if (k == 0) then
            s = default
        else
            kv(k)%used = .true.
            s = kv(k)%val
        end if
    end function get_str

    integer function get_int(key, default)
        character(len=*), intent(in) :: key
        integer,          intent(in) :: default
        integer :: k
        k = find_key(key)
        get_int = default
        if (k > 0) then
            kv(k)%used = .true.
            get_int = to_int(kv(k)%val, k)
        end if
    end function get_int

    real(dp) function get_real(key, default)
        character(len=*), intent(in) :: key
        real(dp),         intent(in) :: default
        integer :: k
        k = find_key(key)
        get_real = default
        if (k > 0) then
            kv(k)%used = .true.
            get_real = to_real(kv(k)%val, k)
        end if
    end function get_real

    integer function to_int(s, k)
        character(len=*), intent(in) :: s
        integer,          intent(in) :: k
        integer :: ios
        read(s,*,iostat=ios) to_int
        if (ios /= 0) call fail_line(k, 'not an integer: '//trim(s))
    end function to_int

    real(dp) function to_real(s, k)
        character(len=*), intent(in) :: s
        integer,          intent(in) :: k
        integer :: ios
        read(s,*,iostat=ios) to_real
        if (ios /= 0) call fail_line(k, 'not a number: '//trim(s))
    end function to_real

    !> region 경계: * 는 무한대 (sgn = -1 이면 -무한대)
    real(dp) function to_bound(s, sgn, k)
        character(len=*), intent(in) :: s
        integer,          intent(in) :: sgn, k
        if (trim(s) == '*') then
            to_bound = sgn*huge(1.0_dp)
        else
            to_bound = to_real(s, k)
        end if
    end function to_bound

    logical function is_number(s)
        character(len=*), intent(in) :: s
        real(dp) :: r
        integer :: ios
        read(s,*,iostat=ios) r
        is_number = (ios == 0) .and. scan(trim(s), '0123456789') > 0
    end function is_number

    integer function need_state(sname, where)
        character(len=*), intent(in) :: sname, where
        integer :: k
        need_state = 0
        do k = 1,nstates
            if (states(k)%name == lower(sname)) need_state = k
        end do
        if (need_state == 0) then
            write(*,'(a,a,a,a,a,a)') ' Error: ', trim(where), ': unknown state ''', trim(sname), &
                ''' (define it in flow.in as state-', trim(lower(sname))//' = rho u v w p)'
            error stop
        end if
    end function need_state

    !> 공백, 쉼표로 나누기
    subroutine split_tokens(s, tok, nt)
        character(len=*), intent(in)  :: s
        character(len=*), intent(out) :: tok(:)
        integer,          intent(out) :: nt
        integer :: i, i0, n
        logical :: in_tok
        character(len=1) :: ch
        nt = 0;  in_tok = .false.;  i0 = 1
        n = len_trim(s)
        do i = 1,n+1
            if (i <= n) then
                ch = s(i:i)
            else
                ch = ' '
            end if
            if (ch == ' ' .or. ch == ',') then
                if (in_tok) then
                    if (nt < size(tok)) then
                        nt = nt + 1
                        tok(nt) = s(i0:i-1)
                    end if
                    in_tok = .false.
                end if
            else if (.not. in_tok) then
                in_tok = .true.;  i0 = i
            end if
        end do
    end subroutine split_tokens

    !> region-<k> 를 k 숫자 순서로 정렬 (prefix 길이 = off-1)
    subroutine sort_by_suffix(idx, off)
        integer, intent(inout) :: idx(:)
        integer, intent(in)    :: off
        integer :: a, b, t
        do a = 2,size(idx)
            do b = a, 2, -1
                if (to_int(kv(idx(b))%key(off:), idx(b)) < to_int(kv(idx(b-1))%key(off:), idx(b-1))) then
                    t = idx(b);  idx(b) = idx(b-1);  idx(b-1) = t
                else
                    exit
                end if
            end do
        end do
    end subroutine sort_by_suffix

    !> 쓰지 않은 key 경고 (오타 잡기). ignore 로 시작하는 key 는 조용히 넘어감
    subroutine warn_unused(ignore)
        character(len=*), intent(in) :: ignore(:)
        integer :: k, m
        logical :: skip
        do k = 1,nkv
            if (kv(k)%used) cycle
            skip = .false.
            do m = 1,size(ignore)
                if (len_trim(ignore(m)) > 0) then
                    if (index(kv(k)%key, trim(ignore(m))) == 1) skip = .true.
                end if
            end do
            if (.not. skip) write(*,'(a,a,a,i0,a,a,a)') ' Warning: ', trim(cur_file), ' line ', kv(k)%line, &
                ': unknown key ''', trim(kv(k)%key), ''' ignored'
        end do
    end subroutine warn_unused

    subroutine fail_line(k, text)
        integer,          intent(in) :: k
        character(len=*), intent(in) :: text
        write(*,'(a,a,a,i0,a,a)') ' Error: ', trim(cur_file), ' line ', kv(k)%line, ': ', text
        write(*,'(a,a,a,a)') '   ', trim(kv(k)%key), ' = ', trim(kv(k)%val)
        error stop
    end subroutine fail_line

    subroutine fail_kv(text)
        character(len=*), intent(in) :: text
        write(*,'(a,a,a,a)') ' Error: ', trim(cur_file), ': ', text
        error stop
    end subroutine fail_kv

    logical function is_absolute(p)
        character(len=*), intent(in) :: p
        is_absolute = .false.
        if (len_trim(p) == 0) return
        if (p(1:1) == '/' .or. p(1:1) == '\') is_absolute = .true.
        if (len_trim(p) >= 2) then
            if (p(2:2) == ':') is_absolute = .true.
        end if
    end function is_absolute

    pure function replace_tabs(s) result(t)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: t
        integer :: k
        t = s
        do k = 1,len(s)
            if (t(k:k) == achar(9)) t(k:k) = ' '
        end do
    end function replace_tabs

    pure function lower(s) result(t)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: t
        integer :: k, c
        t = s
        do k = 1,len(s)
            c = iachar(s(k:k))
            if (c >= iachar('A') .and. c <= iachar('Z')) t(k:k) = achar(c + 32)
        end do
    end function lower

end module case_setup
