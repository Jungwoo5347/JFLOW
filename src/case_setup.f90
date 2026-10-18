!> 케이스 입력 읽기 (KFLOW 스타일)
!>
!> 실행:  .\main.exe cases\wedge        (인자 없으면 현재 폴더)
!>   케이스 폴더 안의 두 파일을 읽음
!>     grid.in : 격자 파일과 경계조건 (KFLOW grid.in 형식)
!>     flow.in : 계산 제어, 기체, 유동 상태, 초기조건
!>   결과도 케이스 폴더에 저장:  <output>_result.dat, <output>_line.dat, <output>_residual.dat
!>
!> 공통 형식 (KFLOW):  $key = value   ( # 뒤는 주석, 빈 줄 무시, 대소문자 무시, $ 는 생략 가능 )
!>
!> ---- grid.in ----------------------------------------------------------------
!>   dim      = 2
!>   nblock   = 1
!>   block-1  = mesh.x                        PLOT3D 격자 파일 (케이스 폴더 기준 경로)
!>   block-1  = cartesian ni nj xmin xmax ymin ymax     또는 균일 직교 격자 생성 (ni, nj: 노드 수)
!>   form-1   = 0                             0: unformatted (바이너리), 1: ascii
!>   bc-1-f   = code                          블록 1, 면 f 전체
!>   bc-1-f   = code n1 n2                    면 f 의 노드 n1 ~ n2 구간 (1부터 센 노드 번호)
!>   bc-1-f   = code [n1 n2] state            입구/원방이면 flow.in 의 상태 이름 (생략 시 freestream)
!>     면 번호 f : 1 = imin, 2 = imax, 3 = jmin, 4 = jmax
!>     범위는 그 면을 따라가는 방향의 노드 번호 (면 1,2 는 j, 면 3,4 는 i)
!>     같은 면을 여러 줄로 나눠 쓸 수 있음 (bc-1-3 = ... 를 여러 번), 뒤에 쓴 줄이 우선
!>     경계조건이 없는 면은 기하학적으로 맞닿은 다른 경계 면을 자동으로 찾아 연결
!>     (C-grid wake cut, O-grid 이음매). 짝이 없으면 에러
!>   output   = name                          결과 파일 이름 머리
!>   decompose-1, ncpu-1, diswall             KFLOW 병렬/난류용 항목: 읽기만 하고 사용 안 함
!>
!>   경계조건 코드 (bc_code_table 한 곳에서 정의)
!>     110 : 벽 (비점성 미끄럼 벽)
!>     200 : 원방 (상태 고정, 기본 freestream)  * 특성 조건은 나중에
!>     210 : 초음속 입구 (상태 고정)
!>     220 : 초음속 출구 (내부 값 복사)
!>     300 : 대칭면
!>
!> ---- flow.in ----------------------------------------------------------------
!>   steady         = 1          1: 정상 (수렴할 때까지), 0: 비정상 (t_end 까지)
!>   max_iter       = 50000      정상: 최대 반복
!>   res_tol        = 1.0e-6     정상: 정규화 잔차 수렴 기준 (0 이하: 판정 안 함, max_iter 까지)
!>   local_dt       = 1          정상: 1 이면 셀마다 자기 dt (국소 시간 간격, 수렴 빠름), 0 이면 전체 최솟값
!>   t_end          = 0.2        비정상: 끝 시간 (무차원. ref_length 가 있으면 초 [s])
!>   cfl            = 0.4
!>   gamma          = 1.4
!>   sigma          = 1          MUSCL 재구성 on/off (0: 1차, 1: 2차 이상)
!>   kappa          = 0.3333     MUSCL kappa (-1 풍상 2차, 0 Fromm, 1/3 3차, 1 중심)
!>   b              = 1          minmod 압축 계수, 1 <= b <= (3-kappa)/(1-kappa)  (1: 가장 안정, 클수록 날카로움)
!>   flux           = rusanov    리만 솔버: rusanov | roe
!>   time_scheme    = explicit   시간 적분: explicit (SSP-RK3) | implicit
!>                                 정상 (steady = 1): LU-SGS 반복 (cfl 을 크게, 예 5 ~ 50)
!>                                 비정상 (steady = 0): dual time stepping, 실제 시간은 BDF2 (dt 고정),
!>                                   스텝마다 가상 시간 LU-SGS 반복 (국소 dtau, cfl 은 가상 시간용)
!>   implicit_solver = lusgs     암시적 반복 방법 (지금은 lusgs 만)
!>   lusgs_omega    = 1.0        LU-SGS 대각 계수 w (1 ~ 2, 클수록 안정, 수렴은 느려짐)
!>   dt             = 0.1        dual time: 실제 시간 간격 (무차원. ref_length 가 있으면 초 [s]. 마지막 스텝만 t_end 에 맞춰 짧아짐)
!>   inner_iter     = 30         dual time: 스텝마다 안쪽 최대 반복 수
!>   inner_tol      = 1.0e-3     dual time: 안쪽 잔차 (L2) 가 첫 반복의 이 비율 아래로 내려가면 다음 스텝
!>   print_interval = 10         화면 출력 간격 (잔차 파일에는 매 반복). 출력할 때 반복당 시간, 남은 시간도 표시
!>   print_inner    = 1          dual time: print_interval 스텝 안에서 몇 안쪽 반복마다 화면에 찍을지
!>                               (1: 전부, 10: 10, 20, ... 번째, 0: 마지막만. 마지막 안쪽 반복은 항상 찍음)
!>   line_row       = 1          <output>_line.dat 으로 뽑을 셀 줄 (j 셀 번호)
!>   output_format  = plt        결과 파일: plt (Tecplot 바이너리) | ascii (Tecplot .dat) | both
!>   save_interval  = 0          N 회마다 결과와 재시작 파일 저장 (0: 끝날 때만)
!>   save_series    = 0          1: 중간 저장 결과를 res/<output>_result_00012000.plt 처럼 케이스 폴더의 res/ 에 반복 수를 붙여 따로 저장
!>                               (비정상 문제의 시간 이력. 재시작 파일은 하나만 덮어씀. 파일마다 solution time 이 들어가서
!>                                Tecplot 에서 여러 파일을 함께 열면 애니메이션이 됨)
!>   restart        = 0          1: <output>_restart.bin 에서 이어서 계산
!>                               (반복 수와 시간은 이어서 셈: max_iter, t_end 는 전체 기준)
!>   viscous        = 0          1: 층류 Navier-Stokes (점성항)
!>   reynolds       = 1.0e6      Re = rho U L / mu  (U = state-freestream 의 속도, L = 격자 단위 1)
!>   prandtl        = 0.72
!>   viscosity_law  = sutherland 점성계수: sutherland | constant (mu* = 1)
!>   t_inf          = 288.15     자유류 온도 [K] (Sutherland, 차원 출력용: 무차원 T = gamma p / rho = 1 인 온도)
!>   ref_length     = 0.01       [m] 격자 길이 1 의 실제 길이. 있으면 결과 파일 (result, line) 을 SI 단위로 씀
!>                               (x, y [m], rho [kg/m3], u, v [m/s], p [Pa], solution time [s], Mach 그대로)
!>                               화면, 잔차 파일의 시간도 초. 입력 t_end, dt 도 초로 받음
!>                               없으면 모두 무차원. 계산 자체와 재시작 파일은 항상 무차원
!>                               기준 밀도: viscous = 1 이면 reynolds 로 정해짐, 아니면 p_inf 로
!>   p_inf          = 101325     [Pa] 비점성 + ref_length 일 때 p* = 1/gamma 인 압력 (점성이면 쓰지 않음)
!>   state-name     = rho u v p  이름 붙인 유동 상태 (여러 개)
!>   initial        = name       초기조건 기본 상태
!>   region-k       = name x0 x1 y0 y1    상자 x0 <= xc < x1, y0 <= yc < y1 안의 셀만 그 상태로
!>                                         (* 는 무한대, 번호 k 순서대로 적용 = 뒤가 우선)
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
    integer, parameter, public :: BC_wall_ad  = 120
    integer, parameter, public :: BC_wall_iso = 130
    integer, parameter, public :: BC_cut      = -1      ! 내부용: 자동으로 찾은 맞닿은 면 (wake cut 등)

    integer, parameter, public :: max_states = 50, max_regions = 50, max_patches = 200

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
    integer,  public :: print_inner = 1              ! dual time: 화면에 찍을 안쪽 반복 간격 (0: 마지막만)
    integer,  public :: line_row = 1
    logical,  public :: restart = .false.            ! 재시작 파일에서 이어서 계산
    logical,  public :: out_plt = .true., out_ascii = .false.   ! 결과 파일 형식 (output_format)
    integer,  public :: save_interval = 0            ! 중간 저장 간격 (0: 끝날 때만)
    logical,  public :: save_series = .false.        ! 중간 저장 결과를 반복 수 붙인 새 파일로 (비정상 시간 이력)
    real(dp), public :: sigma = 1.0_dp, kappa = 1.0_dp/3.0_dp, b = 4.0_dp

    ! 리만 솔버 선택 (flow.in 의 flux)
    integer, parameter, public :: FLUX_RUSANOV = 1, FLUX_ROE = 2
    integer,           public :: flux_id   = FLUX_RUSANOV

    ! 시간 적분 (flow.in 의 time_scheme)
    logical,           public :: implicit_ts = .false.        ! time_scheme = implicit
    character(len=16), public :: implicit_solver = 'lusgs'    ! 암시적 반복 방법 (나중에 dadi 등)
    real(dp),          public :: lusgs_omega = 1.0_dp         ! LU-SGS 대각 계수
    real(dp),          public :: dt_phys    = 0.0_dp          ! dual time: 실제 시간 간격
    integer,           public :: inner_iter = 30              ! dual time: 안쪽 (가상 시간) 최대 반복
    real(dp),          public :: inner_tol  = 1.0e-3_dp       ! dual time: 안쪽 잔차가 첫 반복의 이 비율 아래면 다음 스텝
    character(len=16), public :: flux_name = 'rusanov'

    ! ---- 점성 (flow.in) ----
    logical,  public :: viscous    = .false.
    real(dp), public :: reynolds   = 1.0e6_dp      ! rho U L / mu, U = freestream 속도, L = 격자 단위
    real(dp), public :: prandtl    = 0.72_dp
    real(dp), public :: t_inf      = 288.15_dp     ! [K], Sutherland 용
    logical,  public :: sutherland = .true.        ! .false. 이면 mu* = 1 (상수)
    real(dp), public :: S_suth     = 0.0_dp        ! 110.4 / t_inf  (읽은 뒤 계산)
    real(dp), public :: visc_coef  = 0.0_dp        ! M_inf / Re     (읽은 뒤 계산)

    ! ---- 차원 출력 (flow.in 의 ref_length 가 있을 때): 결과 파일을 SI 단위로 ----
    real(dp), parameter, public :: R_gas = 287.05_dp     ! 공기 기체 상수 [J/(kg K)]
    logical,  public :: dim_out    = .false.
    real(dp), public :: ref_length = 0.0_dp        ! 격자 길이 1 의 실제 길이 [m]
    real(dp), public :: p_inf      = 101325.0_dp   ! 비점성일 때 p* = 1/gam 인 압력 [Pa]
    real(dp), public :: c_ref = 1.0_dp, rho_ref = 1.0_dp, mu_ref = 0.0_dp   ! 음속 [m/s], 밀도 [kg/m3], 점성계수 [Pa s] (T* = 1)
    real(dp),         public :: t_unit  = 1.0_dp       ! 표시용 시간 환산 (SI: ref_length / c_ref [s],  무차원: 1)
    character(len=8), public :: t_label = 'Time'       ! 표시용 시간 이름 (SI: 'Time [s]')

    ! ---- 격자 (grid.in) ----
    character(len=16),  public :: grid_kind = 'plot3d'      ! 'plot3d' | 'cartesian'
    character(len=512), public :: grid_file = ''
    integer,            public :: grid_form = 0             ! 0: unformatted, 1: ascii
    integer,            public :: grid_ni = 0, grid_nj = 0
    real(dp),           public :: grid_x0 = 0, grid_x1 = 1, grid_y0 = 0, grid_y1 = 1

    ! ---- 유동 상태, 초기조건 ----
    type, public :: state_t
        character(len=64) :: name
        real(dp)          :: W(nvar)            ! rho, u, v, p
    end type
    integer,       public :: nstates = 0
    type(state_t), public :: states(max_states)

    type, public :: region_t
        integer  :: istate
        real(dp) :: x0, x1, y0, y1
    end type
    integer,        public :: init_state = 0
    integer,        public :: nregions = 0
    type(region_t), public :: regions(max_regions)

    ! ---- 경계조건 (grid.in 의 bc-1-f 한 줄 = patch 하나) ----
    type, public :: patch_t
        integer :: side                         ! 1 imin, 2 imax, 3 jmin, 4 jmax
        integer :: n1, n2                       ! 노드 번호 구간 (n2 = 0 이면 면 끝까지)
        integer :: btype
        integer :: istate                       ! 입구/원방 상태 (0 이면 없음)
    end type
    integer,       public :: npatches = 0
    type(patch_t), public :: patches(max_patches)

    public :: read_case, print_case, side_name, bc_name, case_path

    ! ---- key = value 저장 ----
    integer, parameter :: max_kv = 500
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
        integer :: k, nt
        character(len=64) :: tok(16)
        real(dp) :: lim(4), minf
        integer, allocatable :: order(:)
        integer :: nreg_keys
        character(len=16) :: fmt

        call load_kv(fname)

        steady        = get_int('steady', 1) /= 0
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

        ! 시간 적분: explicit (SSP-RK3) | implicit (정상: LU-SGS 반복,  비정상: dual time = BDF2 + 안쪽 LU-SGS 반복)
        fmt            = lower(adjustl(get_str('time_scheme', 'explicit')))
        select case (fmt)
        case ('explicit', 'rk3');  implicit_ts = .false.
        case ('implicit');         implicit_ts = .true.
        case default;  call fail_kv('time_scheme = '//trim(fmt)//' : unknown (use explicit or implicit)')
        end select
        implicit_solver = lower(adjustl(get_str('implicit_solver', implicit_solver)))
        if (implicit_solver /= 'lusgs') call fail_kv('implicit_solver = '//trim(implicit_solver)//' : unknown (use lusgs)')
        lusgs_omega    = get_real('lusgs_omega', lusgs_omega)
        if (lusgs_omega < 1.0_dp) call fail_kv('lusgs_omega must be >= 1')
        dt_phys        = get_real('dt', dt_phys)
        inner_iter     = get_int('inner_iter', inner_iter)
        inner_tol      = get_real('inner_tol', inner_tol)
        if (implicit_ts .and. .not. steady) then
            if (dt_phys <= 0.0_dp) call fail_kv('time_scheme = implicit with steady = 0 needs dt > 0 (physical time step)')
            if (inner_iter < 1) call fail_kv('inner_iter must be >= 1')
        end if

        if (b < 1.0_dp) call fail_kv('b must be >= 1')
        if (kappa < 1.0_dp) then
            if (b > (3.0_dp - kappa)/(1.0_dp - kappa) + 1.0e-12_dp) call fail_kv('b must be <= (3-kappa)/(1-kappa)')
        end if
        print_interval = get_int('print_interval', print_interval)
        print_inner    = get_int('print_inner', print_inner)
        line_row       = get_int('line_row', line_row)
        if (cfl <= 0) call fail_kv('cfl must be > 0')
        restart        = get_int('restart', 0) /= 0
        save_interval  = get_int('save_interval', save_interval)
        if (print_interval < 1) call fail_kv('print_interval must be >= 1')
        if (print_inner < 0) call fail_kv('print_inner must be >= 0 (0: last inner iteration only)')
        if (save_interval < 0) call fail_kv('save_interval must be >= 0')
        save_series    = get_int('save_series', 0) /= 0
        if (save_series .and. save_interval == 0) call fail_kv('save_series = 1 needs save_interval > 0')
        fmt            = lower(adjustl(get_str('output_format', 'plt')))
        select case (fmt)
        case ('plt');    out_plt = .true.;   out_ascii = .false.
        case ('ascii');  out_plt = .false.;  out_ascii = .true.
        case ('both');   out_plt = .true.;   out_ascii = .true.
        case default;    call fail_kv('output_format = '//trim(fmt)//' : unknown (use plt, ascii or both)')
        end select

        ! state-<name> = rho u v p
        do k = 1,nkv
            if (index(kv(k)%key, 'state-') /= 1) cycle
            kv(k)%used = .true.
            call split_tokens(kv(k)%val, tok, nt)
            if (nt /= 4) call fail_line(k, 'state needs 4 numbers: rho u v p')
            if (nstates == max_states) call fail_line(k, 'too many states')
            nstates = nstates + 1
            states(nstates)%name = kv(k)%key(7:)
            states(nstates)%W = [to_real(tok(1), k), to_real(tok(2), k), to_real(tok(3), k), to_real(tok(4), k)]
            if (states(nstates)%W(1) <= 0 .or. states(nstates)%W(4) <= 0) &
                call fail_line(k, 'rho and p must be > 0')
        end do
        if (nstates == 0) call fail_kv('no state-<name> = rho u v p lines')

        init_state = need_state(get_str('initial', ''), 'initial')

        ! 점성: state-freestream 을 읽은 뒤라야 M_inf 를 구할 수 있음
        viscous  = get_int('viscous', 0) /= 0
        reynolds = get_real('reynolds', reynolds)
        prandtl  = get_real('prandtl', prandtl)
        t_inf    = get_real('t_inf', t_inf)
        select case (lower(adjustl(get_str('viscosity_law', 'sutherland'))))
        case ('sutherland');  sutherland = .true.
        case ('constant');    sutherland = .false.
        case default;         call fail_kv('viscosity_law must be sutherland or constant')
        end select
        if (viscous) then
            k = need_state('freestream', 'viscous (reference state for Re)')
            minf = sqrt(states(k)%W(2)**2 + states(k)%W(3)**2) / sqrt(gam*states(k)%W(4)/states(k)%W(1))
            if (minf <= 0.0_dp) call fail_kv('viscous: state-freestream must have a nonzero speed (defines Re)')
            if (reynolds <= 0.0_dp) call fail_kv('reynolds must be > 0')
            if (prandtl <= 0.0_dp) call fail_kv('prandtl must be > 0')
            if (t_inf <= 0.0_dp) call fail_kv('t_inf must be > 0 (kelvin)')
            visc_coef = minf/reynolds
            S_suth    = 110.4_dp/t_inf
        end if

        ! 차원 출력 (ref_length 가 있을 때만): 무차원 기준 -> SI 기준값
        !   T* = 1 <-> t_inf,  c_ref = sqrt(gam R t_inf),  길이 1 <-> ref_length
        !   rho_ref: 점성이면 Re 로 (visc_coef = mu_ref / (rho_ref c_ref L)),  비점성이면 p_inf (= rho_ref R t_inf, p* = 1/gam 인 압력)
        ref_length = get_real('ref_length', 0.0_dp)
        p_inf      = get_real('p_inf', p_inf)
        dim_out    = ref_length > 0.0_dp
        if (find_key('ref_length') > 0 .and. .not. dim_out) call fail_kv('ref_length must be > 0 [m]')
        if (dim_out) then
            if (t_inf <= 0.0_dp) call fail_kv('t_inf must be > 0 (kelvin)')
            c_ref  = sqrt(gam*R_gas*t_inf)
            t_unit = ref_length/c_ref
            t_label = 'Time [s]'
            ! 시간 입력 (t_end, dt) 은 초 [s] 로 받음 -> 내부 무차원으로
            t_end   = t_end/t_unit
            dt_phys = dt_phys/t_unit
            mu_ref = 1.716e-5_dp*(t_inf/273.15_dp)**1.5_dp*(273.15_dp + 110.4_dp)/(t_inf + 110.4_dp)   ! Sutherland [Pa s]
            if (viscous) then
                rho_ref = mu_ref/(visc_coef*c_ref*ref_length)
                if (find_key('p_inf') > 0) write(*,'(a)') &
                    ' Warning: flow.in: p_inf is not used with viscous = 1 (density follows from reynolds and ref_length)'
            else
                if (p_inf <= 0.0_dp) call fail_kv('p_inf must be > 0 [Pa]')
                rho_ref = p_inf/(R_gas*t_inf)
            end if
        end if

        ! region-<k> = name x0 x1 y0 y1   (k 순서대로)
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
                if (nt /= 5) call fail_line(order(k), 'region needs: state x0 x1 y0 y1  (* = unbounded)')
                lim(1) = to_bound(tok(2), -1, order(k));  lim(2) = to_bound(tok(3), 1, order(k))
                lim(3) = to_bound(tok(4), -1, order(k));  lim(4) = to_bound(tok(5), 1, order(k))
                if (nregions == max_regions) call fail_line(order(k), 'too many regions')
                nregions = nregions + 1
                regions(nregions) = region_t(need_state(tok(1), trim(e%key)), lim(1), lim(2), lim(3), lim(4))
            end associate
        end do

        call warn_unused([character(len=16) :: ''])
    end subroutine read_flow_in

    ! ======================================================================
    !  grid.in
    ! ======================================================================
    subroutine read_grid_in(fname)
        character(len=*), intent(in) :: fname
        integer :: k, nt, face, n1, n2, code, istate, nnum
        character(len=64) :: tok(16)
        character(len=512) :: blk

        call load_kv(fname)

        if (get_int('dim', 2) /= 2) call fail_kv('only dim = 2 is supported')
        if (get_int('nblock', 1) /= 1) call fail_kv('only nblock = 1 is supported for now')
        output_name = get_str('output', 'jflow')

        ! block-1 = file   또는   block-1 = cartesian ni nj xmin xmax ymin ymax
        blk = get_str('block-1', '')
        if (len_trim(blk) == 0) call fail_kv('block-1 = <grid file> is missing')
        call split_tokens(blk, tok, nt)
        if (lower(tok(1)) == 'cartesian') then
            k = find_key('block-1')
            if (nt /= 7) call fail_line(k, 'cartesian needs: ni nj xmin xmax ymin ymax')
            grid_kind = 'cartesian'
            grid_ni = to_int(tok(2), k);   grid_nj = to_int(tok(3), k)
            grid_x0 = to_real(tok(4), k);  grid_x1 = to_real(tok(5), k)
            grid_y0 = to_real(tok(6), k);  grid_y1 = to_real(tok(7), k)
            if (grid_ni < 2 .or. grid_nj < 2) call fail_line(k, 'ni, nj (node counts) must be >= 2')
        else
            grid_kind = 'plot3d'
            grid_file = case_path(trim(blk))
            grid_form = get_int('form-1', 0)
            if (grid_form /= 0 .and. grid_form /= 1) call fail_kv('form-1 must be 0 (unformatted) or 1 (ascii)')
        end if

        ! bc-1-<face> = code [n1 n2] [state]   (파일에 쓴 순서대로 적용 = 뒤가 우선)
        do k = 1,nkv
            if (index(kv(k)%key, 'bc-') /= 1) cycle
            kv(k)%used = .true.
            if (index(kv(k)%key, 'bc-1-') /= 1 .or. len_trim(kv(k)%key) /= 6) &
                call fail_line(k, 'expected bc-1-<face>, face = 1..4')
            face = to_int(kv(k)%key(6:6), k)
            if (face < 1 .or. face > 4) call fail_line(k, 'face must be 1 (imin), 2 (imax), 3 (jmin) or 4 (jmax)')

            call split_tokens(kv(k)%val, tok, nt)
            if (nt < 1) call fail_line(k, 'boundary condition code is missing')
            code = to_int(tok(1), k)
            if (bc_name(code) == '?') call fail_line(k, 'unknown boundary condition code '//trim(tok(1)))

            ! 코드 뒤: 숫자 0개 또는 2개 (노드 범위), 그다음 선택적으로 상태 이름
            nnum = 0
            do while (1 + nnum + 1 <= nt)
                if (.not. is_number(tok(2+nnum))) exit
                nnum = nnum + 1
            end do
            n1 = 1;  n2 = 0
            if (nnum == 2) then
                n1 = to_int(tok(2), k);  n2 = to_int(tok(3), k)
            else if (nnum /= 0) then
                call fail_line(k, 'range must be two node numbers: n1 n2')
            end if
            if (nt > 1 + nnum + 1) call fail_line(k, 'too many items')

            istate = 0
            if (code == BC_inlet .or. code == BC_farfield .or. code == BC_wall_iso) then
                if (nt == 1 + nnum + 1) then
                    istate = need_state(tok(nt), trim(kv(k)%key))
                else if (code == BC_wall_iso) then
                    call fail_line(k, 'code 130 needs a wall state name (wall velocity u v, T = gamma p / rho)')
                else
                    istate = need_state('freestream', trim(kv(k)%key)//' (default state)')
                end if
            else if (nt == 1 + nnum + 1) then
                call fail_line(k, 'a state name is only used with codes 130, 200 and 210')
            end if

            if (npatches == max_patches) call fail_line(k, 'too many bc lines')
            npatches = npatches + 1
            patches(npatches) = patch_t(face, n1, n2, code, istate)
        end do

        ! KFLOW 에서 쓰는 항목: 읽기만 함
        call warn_unused([character(len=16) :: 'decompose-', 'ncpu-', 'diswall'])
    end subroutine read_grid_in

    !> 읽은 설정 요약
    subroutine print_case()
        integer :: k
        character(len=16) :: range

        write(*,'(a)')   ' ================================================================'
        write(*,'(a,a)') '  case     : ', trim(case_dir)
        if (grid_kind == 'cartesian') then
            write(*,'(a,i0,a,i0,a)') '  grid     : cartesian ', grid_ni, ' x ', grid_nj, ' nodes'
        else
            write(*,'(a,a,a,a)') '  grid     : ', trim(grid_file), merge('  (unformatted)', '  (ascii)      ', grid_form == 0)
        end if
        if (steady) then
            if (res_tol > 0.0_dp) then
                write(*,'(a,es9.2,a,i0)') '  steady   : res_tol = ', res_tol, ',  max_iter = ', max_steps
            else
                write(*,'(a,i0)') '  steady   : res_tol = off (run to max_iter),  max_iter = ', max_steps
            end if
            write(*,'(a,a)') '  time step: ', merge('local (each cell its own dt)', 'global (min dt)             ', local_dt)
        else
            write(*,'(a,es12.4,a)')   '  unsteady : t_end = ', t_end*t_unit, trim(merge(' s', '  ', dim_out))
        end if
        write(*,'(a,f6.3,a,f6.3)')    '  CFL      : ', cfl, '    gamma = ', gam
        write(*,'(a,f6.3,a,f7.4,a,f6.3)') '  MUSCL    : sigma = ', sigma, '    kappa = ', kappa, '    b = ', b
        write(*,'(a,a)')              '  flux     : ', trim(flux_name)
        if (implicit_ts .and. steady) then
            write(*,'(a,f5.2,a)') '  time int : implicit LU-SGS  (omega = ', lusgs_omega, ')'
        else if (implicit_ts) then
            write(*,'(a,es10.3,2a,i0,a,es9.2,a,f5.2,a)') '  time int : dual time BDF2, dt = ', dt_phys*t_unit, &
                trim(merge(' s', '  ', dim_out)), &
                ',  inner LU-SGS <= ', inner_iter, ' its to ', inner_tol, '  (omega = ', lusgs_omega, ')'
        else
            write(*,'(a)')        '  time int : explicit SSP-RK3'
        end if
        write(*,'(a)') '  states   :         rho           u           v           p'
        do k = 1,nstates
            write(*,'(4x,a12,4es12.4)') states(k)%name, states(k)%W
        end do
        write(*,'(a,a,a,i0,a)') '  initial  : ', trim(states(init_state)%name), '  (+ ', nregions, ' region(s))'
        if (restart) write(*,'(a,a)') '  restart  : ', trim(case_path(trim(output_name)//'_restart.bin'))
        if (viscous) write(*,'(a,es10.3,a,f5.2,a,es10.3,a,a)') '  viscous  : Re = ', reynolds, &
            '  Pr = ', prandtl, '  M/Re = ', visc_coef, '  mu: ', merge('sutherland', 'constant  ', sutherland)
        if (dim_out) then
            write(*,'(a,es10.3,a,f7.2,a,es10.3,a,es10.3,a)') '  SI out   : L = ', ref_length, ' m,  T = ', t_inf, &
                ' K,  rho = ', rho_ref, ' kg/m3,  p = ', rho_ref*R_gas*t_inf, ' Pa  (at rho* = 1, p* = 1/gamma)'
            write(*,'(a,f8.2,a,es10.3,a,es11.4,a)') '             c = ', c_ref, ' m/s,  mu = ', mu_ref, &
                ' Pa s,  time 1 = ', ref_length/c_ref, ' s'
            if (.not. steady) write(*,'(a,es11.4,a,es11.4,a)') '             t_end = ', t_end*t_unit, ' s  (', t_end, ' nondim)'
            if (.not. steady .and. implicit_ts) write(*,'(a,es11.4,a,es11.4,a)') '             dt    = ', dt_phys*t_unit, &
                ' s  (', dt_phys, ' nondim)'
        end if
        write(*,'(a)') '  bc       :  face        nodes      code  type       state'
        do k = 1,npatches
            if (patches(k)%n2 == 0) then
                write(range,'(i0,a)') patches(k)%n1, '-end'
            else
                write(range,'(i0,a,i0)') patches(k)%n1, '-', patches(k)%n2
            end if
            write(*,'(14x,i1,1x,a4,2x,a12,1x,i4,2x,a9)', advance='no') patches(k)%side, side_name(patches(k)%side), &
                range, patches(k)%btype, bc_name(patches(k)%btype)
            if (patches(k)%istate > 0) then
                write(*,'(2x,a)') trim(states(patches(k)%istate)%name)
            else
                write(*,*)
            end if
        end do
        write(*,'(a,a,a,a,a)') '  output   : ', trim(case_path(trim(output_name)//'_*')), &
            merge('  (result: plt', '  (result:    ', out_plt), merge(' ascii)', ')      ', out_ascii)
        if (save_interval > 0) write(*,'(a,i0,a,a)') '  save     : every ', save_interval, ' iterations (result + restart)', &
            merge(', result numbered by iteration (series)', '                                       ', save_series)
        write(*,'(a)')   ' ================================================================'
    end subroutine print_case

    function side_name(id) result(s)
        integer, intent(in) :: id
        character(len=4) :: s
        character(len=4), parameter :: names(4) = ['imin', 'imax', 'jmin', 'jmax']
        s = names(id)
    end function side_name

    !> 경계조건 코드 -> 이름  (모르는 코드는 '?')
    function bc_name(code) result(s)
        integer, intent(in) :: code
        character(len=9) :: s
        select case (code)
        case (BC_wall);      s = 'wall'
        case (BC_wall_ad);   s = 'wall_ad'
        case (BC_wall_iso);  s = 'wall_iso'
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
                ''' (define it in flow.in as state-', trim(lower(sname))//' = rho u v p)'
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
