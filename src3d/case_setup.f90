!> 케이스 입력 읽기 (KFLOW 스타일) - 3D
!>
!> 실행:  .\main3d.exe cases3d\sod        (인자 없으면 현재 폴더)
!>   케이스 폴더 안의 두 파일을 읽음
!>     grid.in : 격자 파일과 경계조건 (KFLOW grid.in 형식)
!>     flow.in : 계산 제어, 기체, 유동 상태, 초기조건
!>   결과도 케이스 폴더에 저장:  <output>_result.dat, <output>_line.dat, <output>_residual.dat
!>
!> 공통 형식 (KFLOW):  $key = value   ( # 뒤는 주석, 빈 줄 무시, 대소문자 무시, $ 는 생략 가능 )
!>
!> ---- grid.in ----------------------------------------------------------------
!>   dim      = 3                             3D 솔버는 3 만 받음 (2D 는 main.exe)
!>   nblock   = 1
!>   block-1  = mesh.x                        PLOT3D 격자 파일 (케이스 폴더 기준 경로)
!>   block-1  = cartesian ni nj nk xmin xmax ymin ymax zmin zmax   또는 균일 직교 격자 생성 (노드 수)
!>   form-1   = 0                             0: unformatted (바이너리), 1: ascii
!>   extrude-1 = nk span                      2D PLOT3D 파일을 z 방향 (0 ~ span) 으로 nk 노드만큼 밀어내 3D 로
!>                                            (2D 결과와 비교할 때: 면 5, 6 을 대칭면 300 으로)
!>   bc-1-f   = code                          블록 1, 면 f 전체
!>   bc-1-f   = code a1 a2 b1 b2              면 f 의 노드 구간 (1부터 센 노드 번호, 면 위 두 방향)
!>   bc-1-f   = code [a1 a2 b1 b2] state      입구/원방이면 flow.in 의 상태 이름 (생략 시 freestream)
!>     면 번호 f : 1 = imin, 2 = imax, 3 = jmin, 4 = jmax, 5 = kmin, 6 = kmax
!>     범위의 두 방향 (a, b) : 면 1,2 는 (j, k),  면 3,4 는 (i, k),  면 5,6 은 (i, j)
!>     같은 면을 여러 줄로 나눠 쓸 수 있음, 뒤에 쓴 줄이 우선
!>     경계조건이 없는 면은 네 꼭짓점이 겹치는 다른 경계 면을 자동으로 찾아 연결
!>     (C-grid wake cut, O-grid 이음매). 짝이 없으면 에러
!>   output   = name                          결과 파일 이름 머리
!>   decompose-1, ncpu-1, diswall             KFLOW 병렬/난류용 항목: 읽기만 하고 사용 안 함
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
!>   region-k       = name x0 x1 y0 y1 z0 z1   상자 안의 셀만 그 상태로 (* 는 무한대, k 순서대로 = 뒤가 우선)
!>   line_row       = 1                        <output>_line.dat 으로 뽑을 i 방향 셀 줄의 j 셀 번호
!>   line_k         = 1                        그 줄의 k 셀 번호
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
    integer, parameter, public :: BC_cut      = -1      ! 내부용: 자동으로 찾은 맞닿은 면 (wake cut 등)

    integer, parameter, public :: nside = 6             ! 경계 변 수: imin, imax, jmin, jmax, kmin, kmax
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
    integer,  public :: line_row = 1, line_k = 1
    real(dp), public :: sigma = 1.0_dp, kappa = 1.0_dp/3.0_dp, b = 4.0_dp

    ! 리만 솔버 선택 (flow.in 의 flux)
    integer, parameter, public :: FLUX_RUSANOV = 1, FLUX_ROE = 2
    integer,           public :: flux_id   = FLUX_RUSANOV
    character(len=16), public :: flux_name = 'rusanov'

    ! ---- 격자 (grid.in) ----
    character(len=16),  public :: grid_kind = 'plot3d'      ! 'plot3d' | 'cartesian'
    character(len=512), public :: grid_file = ''
    integer,            public :: grid_form = 0             ! 0: unformatted, 1: ascii
    integer,            public :: grid_ni = 0, grid_nj = 0, grid_nk = 0
    real(dp),           public :: grid_x0 = 0, grid_x1 = 1, grid_y0 = 0, grid_y1 = 1, grid_z0 = 0, grid_z1 = 1
    integer,            public :: extrude_nk = 0            ! 2D 격자 밀어내기 노드 수 (0: 안 함)
    real(dp),           public :: extrude_span = 1.0_dp

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

    ! ---- 경계조건 (grid.in 의 bc-1-f 한 줄 = patch 하나) ----
    type, public :: patch_t
        integer :: side                         ! 1 imin, 2 imax, 3 jmin, 4 jmax, 5 kmin, 6 kmax
        integer :: a1, a2, b1, b2               ! 노드 번호 구간 (면 위 두 방향, a2/b2 = 0 이면 면 끝까지)
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
        integer :: k, m, nt
        character(len=64) :: tok(16)
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
        print_interval = get_int('print_interval', print_interval)
        line_row       = get_int('line_row', line_row)
        line_k         = get_int('line_k', line_k)
        if (cfl <= 0) call fail_kv('cfl must be > 0')
        if (print_interval < 1) call fail_kv('print_interval must be >= 1')

        ! state-<name> = rho u v w p
        do k = 1,nkv
            if (index(kv(k)%key, 'state-') /= 1) cycle
            kv(k)%used = .true.
            call split_tokens(kv(k)%val, tok, nt)
            if (nt /= 5) call fail_line(k, 'state needs 5 numbers: rho u v w p  (3D)')
            if (nstates == max_states) call fail_line(k, 'too many states')
            nstates = nstates + 1
            states(nstates)%name = kv(k)%key(7:)
            do m = 1,5
                states(nstates)%W(m) = to_real(tok(m), k)
            end do
            if (states(nstates)%W(1) <= 0 .or. states(nstates)%W(5) <= 0) &
                call fail_line(k, 'rho and p must be > 0')
        end do
        if (nstates == 0) call fail_kv('no state-<name> = rho u v w p lines')

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

    ! ======================================================================
    !  grid.in
    ! ======================================================================
    subroutine read_grid_in(fname)
        character(len=*), intent(in) :: fname
        integer :: k, nt, face, code, istate, nnum, r(4)
        character(len=64) :: tok(16)
        character(len=512) :: blk

        call load_kv(fname)

        if (get_int('dim', 3) /= 3) call fail_kv('main3d.exe needs dim = 3 (use main.exe for 2D cases)')
        if (get_int('nblock', 1) /= 1) call fail_kv('only nblock = 1 is supported for now')
        output_name = get_str('output', 'jflow')

        ! block-1 = file   또는   block-1 = cartesian ni nj nk xmin xmax ymin ymax zmin zmax
        blk = get_str('block-1', '')
        if (len_trim(blk) == 0) call fail_kv('block-1 = <grid file> is missing')
        call split_tokens(blk, tok, nt)
        if (lower(tok(1)) == 'cartesian') then
            k = find_key('block-1')
            if (nt /= 10) call fail_line(k, 'cartesian needs: ni nj nk xmin xmax ymin ymax zmin zmax')
            grid_kind = 'cartesian'
            grid_ni = to_int(tok(2), k);   grid_nj = to_int(tok(3), k);   grid_nk = to_int(tok(4), k)
            grid_x0 = to_real(tok(5), k);  grid_x1 = to_real(tok(6), k)
            grid_y0 = to_real(tok(7), k);  grid_y1 = to_real(tok(8), k)
            grid_z0 = to_real(tok(9), k);  grid_z1 = to_real(tok(10), k)
            if (grid_ni < 2 .or. grid_nj < 2 .or. grid_nk < 2) call fail_line(k, 'ni, nj, nk (node counts) must be >= 2')
        else
            grid_kind = 'plot3d'
            grid_file = case_path(trim(blk))
            grid_form = get_int('form-1', 0)
            if (grid_form /= 0 .and. grid_form /= 1) call fail_kv('form-1 must be 0 (unformatted) or 1 (ascii)')
        end if

        ! extrude-1 = nk span   (2D PLOT3D 파일을 3D 로)
        k = find_key('extrude-1')
        if (k > 0) then
            kv(k)%used = .true.
            if (grid_kind /= 'plot3d') call fail_line(k, 'extrude is only for 2D PLOT3D grid files')
            call split_tokens(kv(k)%val, tok, nt)
            if (nt /= 2) call fail_line(k, 'extrude needs: nk span')
            extrude_nk   = to_int(tok(1), k)
            extrude_span = to_real(tok(2), k)
            if (extrude_nk < 2 .or. extrude_span <= 0.0_dp) call fail_line(k, 'need nk >= 2 and span > 0')
        end if

        ! bc-1-<face> = code [a1 a2 b1 b2] [state]   (파일에 쓴 순서대로 적용 = 뒤가 우선)
        do k = 1,nkv
            if (index(kv(k)%key, 'bc-') /= 1) cycle
            kv(k)%used = .true.
            if (index(kv(k)%key, 'bc-1-') /= 1 .or. len_trim(kv(k)%key) /= 6) &
                call fail_line(k, 'expected bc-1-<face>, face = 1..6')
            face = to_int(kv(k)%key(6:6), k)
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
            patches(npatches) = patch_t(face, r(1), r(2), r(3), r(4), code, istate)
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
        if (grid_kind == 'cartesian') then
            write(*,'(a,i0,a,i0,a,i0,a)') '  grid     : cartesian ', grid_ni, ' x ', grid_nj, ' x ', grid_nk, ' nodes'
        else
            write(*,'(a,a,a,a)') '  grid     : ', trim(grid_file), merge('  (unformatted)', '  (ascii)      ', grid_form == 0)
            if (extrude_nk > 0) write(*,'(a,i0,a,es10.3)') '  extrude  : nk = ', extrude_nk, ',  span = ', extrude_span
        end if
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
        write(*,'(a,f6.3,a,f7.4,a,f6.3)') '  MUSCL    : sigma = ', sigma, '    kappa = ', kappa, '    b = ', b
        write(*,'(a,a)')              '  flux     : ', trim(flux_name)
        write(*,'(a)') '  states   :         rho           u           v           w           p'
        do k = 1,nstates
            write(*,'(4x,a12,5es12.4)') states(k)%name, states(k)%W
        end do
        write(*,'(a,a,a,i0,a)') '  initial  : ', trim(states(init_state)%name), '  (+ ', nregions, ' region(s))'
        write(*,'(a)') '  bc       :  face        nodes (a ; b)           code  type       state'
        do k = 1,npatches
            associate (p => patches(k))
                if (p%a2 == 0 .and. p%b2 == 0) then
                    range = 'all'
                else
                    write(range,'(i0,a,a,a,i0,a,a)') p%a1, '-', trim(endnum(p%a2)), ' ; ', p%b1, '-', trim(endnum(p%b2))
                end if
                write(*,'(14x,i1,1x,a4,2x,a22,1x,i4,2x,a9)', advance='no') p%side, side_name(p%side), &
                    range, p%btype, bc_name(p%btype)
                if (p%istate > 0) then
                    write(*,'(2x,a)') trim(states(p%istate)%name)
                else
                    write(*,*)
                end if
            end associate
        end do
        write(*,'(a,a)') '  output   : ', trim(case_path(trim(output_name)//'_*.dat'))
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
    !  key = value 파일 읽기와 도우미  (2D 와 같음)
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
