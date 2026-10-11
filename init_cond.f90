module init_cond
    use euler
    implicit none

    ! 계산할 문제: 'sod'   = Sod 충격관 (직교 격자 생성, t_end = 0.2)
    !              'wedge' = NASA 초음속 쐐기 (grids/wedge.x.fmt, 마하 2.5, 15도, 정상 상태)
    character(len=*), parameter :: case_name = 'wedge'

    ! 자유류 (입구, 원방 경계조건과 wedge 초기조건에서 사용)
    !   무차원: rho = 1, p = 1/gam  ->  음속 c = 1  ->  u 가 곧 마하수
    real(dp), parameter :: M_inf = 2.5_dp
    real(dp), parameter :: rho_init = 1.0_dp, u_init = M_inf, v_init = 0.0_dp, p_init = 1.0_dp/gam

    ! Sod 충격관 초기조건 (무차원, t_end = 0.2)
    !   xc < x_mid : W_left,  나머지 : W_right       W = (rho, u, v, p)
    !   W_left = W_right 로 두면 균일 유동 테스트
    real(dp), parameter :: x_mid = 0.5_dp
    real(dp), parameter :: W_left(nvar)  = [1.0_dp,   0.0_dp, 0.0_dp, 1.0_dp]
    real(dp), parameter :: W_right(nvar) = [0.125_dp, 0.0_dp, 0.0_dp, 0.1_dp]

end module init_cond
