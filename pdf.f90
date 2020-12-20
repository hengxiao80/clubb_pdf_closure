module pdf

  implicit none

  public :: pdf_closure

  private :: &
    svp_ice_flatau, &
    svp_liq_flatau, &
    sat_mixrat_liq, &
    sat_mixrat_ice, &
    w_closure, &
    x_closure, &
    calc_cf_comp

  contains

  elemental function svp_ice_flatau(T, T_freeze_K) result(esat)

    implicit none

    intrinsic :: max
    real(kind=8), intent(in) :: T, T_freeze_K
    real(kind=8), parameter :: min_T_in_C = -90.0
    real(kind=8), parameter, dimension(9) :: a = &
      100.0*(/6.09868993,     0.499320233,     0.184672631e-1, &
              0.402737184e-3, 0.565392987e-5,  0.521693933e-7, &
              0.307839583e-9, 0.105785160e-11, 0.161444444e-14/)
    real(kind=8) :: T_in_C, esat

    T_in_C = T - T_freeze_K
    T_in_C = max(T_in_C, min_T_in_C)
    esat = a(1) + T_in_C*( a(2)+T_in_C*(a(3)+T_in_C*(a(4)+T_in_C &
      *(a(5)+T_in_C*(a(6)+T_in_C*(a(7)+T_in_C*(a(8)+T_in_C*a(9))))))))

    return
  end function svp_ice_flatau

  elemental function svp_liq_flatau(T, T_freeze_K) result(esat)

    implicit none

    intrinsic :: max
    real(kind=8), intent(in) :: T, T_freeze_K
    real(kind=8), parameter :: min_T_in_C = -85.0
    real(kind=8), parameter, dimension(9) :: a = &
      100.0*(/6.11583699,      0.444606896,     0.143177157e-1, &
              0.264224321e-3,  0.299291081e-5,  0.203154182e-7, &
              0.702620698e-10, 0.379534310e-13, -0.321582393e-15/)
    real(kind=8) :: T_in_C, esat

    T_in_C = T - T_freeze_K
    T_in_C = max(T_in_C, min_T_in_C)
    esat = a(1) + T_in_C*( a(2)+T_in_C*(a(3)+T_in_C*(a(4)+T_in_C &
      *(a(5)+T_in_C*(a(6)+T_in_C*(a(7)+T_in_C*(a(8)+T_in_C*a(9))))))))

    return
  end function svp_liq_flatau

  elemental function sat_mixrat_ice(T, p, ep, T_freeze_K) result(mixrat)

    implicit none

    real(kind=8), intent(in) :: T, p, ep, T_freeze_K
    real(kind=8) :: esat, mixrat

    esat = svp_ice_flatau(T, T_freeze_K)
    if (p-esat < 1.0) then
      mixrat = ep
    else
      mixrat = ep*(esat/(p-esat))
    end if

    return
  end function sat_mixrat_ice

  elemental function sat_mixrat_liq(T, p, ep, T_freeze_K) result(mixrat)

    implicit none

    real(kind=8), intent(in) :: T, p, ep, T_freeze_K
    real(kind=8) :: esat, mixrat

    esat = svp_liq_flatau(T, T_freeze_K)
    if (p-esat < 1.0) then
      mixrat = ep
    else
      mixrat = ep*(esat/(p-esat))
    end if

    return
  end function sat_mixrat_liq

  elemental subroutine w_closure &
    ( &
    skw, wp2, sig_w_sqd, uf_max, &                                      ! in
    uf, w_1, w_2, w_1_n, w_2_n, var_w_1, var_w_2 &                      ! out
    )
    implicit none

    intrinsic :: abs, min, max, sqrt
    real(kind=8), intent(in) :: skw, wp2, sig_w_sqd, uf_max
    real(kind=8), intent(out) :: uf, w_1, w_2, w_1_n, w_2_n, var_w_1, var_w_2

    if (abs(skw) <= 1.0e-5) then
      uf = 0.5
    else
      uf = 0.5*( 1.0 - skw/sqrt(4.0*(1.0-sig_w_sqd)**3+skw**2))
      uf = min(max(uf, 1.0-uf_max), uf_max)
    end if
    w_1_n = (((1.0-uf)/uf)*(1.0-sig_w_sqd))**0.5
    w_2_n = -((uf/(1.0-uf))*(1.0-sig_w_sqd))**0.5
    w_1 = 0.0 + (wp2**0.5)*w_1_n ! wm = 0.0
    w_2 = 0.0 + (wp2**0.5)*w_2_n ! wm = 0.0
    var_w_1 = wp2*sig_w_sqd
    var_w_2 = wp2*sig_w_sqd

    return
  end subroutine w_closure

  elemental subroutine x_closure &
    ( &
    xm, xp2, wp2, wpxp, uf, sig_w_sqd, w_1_n, w_2_n, &                 ! in
    width_fac_1, width_fac_2, tol, zero_threshold, &                   ! in
    x_1, x_2, var_x_1, var_x_2, alpha_x &                              ! out
    )
    implicit none

    intrinsic :: min, max, sqrt
    real(kind=8), intent(in) :: &
      xm, xp2, wp2, wpxp, uf, sig_w_sqd, w_1_n, w_2_n, &
      width_fac_1, width_fac_2, tol, zero_threshold
    real(kind=8), intent(out) :: x_1, x_2, var_x_1, var_x_2, alpha_x

    if (xp2 <= tol**2) then
      x_1 = xm
      x_2 = xm
      var_x_1 = 0.0
      var_x_2 = 0.0
      alpha_x = 0.5
    else
      x_1 = xm - (wpxp/sqrt(wp2))/w_2_n
      x_2 = xm - (wpxp/sqrt(wp2))/w_1_n
      alpha_x = 0.5*(1.0-wpxp*wpxp/((1.0-sig_w_sqd)*wp2*xp2))
      alpha_x = max(min(alpha_x, 1.0), zero_threshold)
      var_x_1 = (alpha_x/uf)*width_fac_1*xp2
      var_x_2 = (alpha_x/(1.0-uf))*width_fac_2*xp2
    end if

    return
  end subroutine x_closure

  elemental subroutine calc_cf_comp &
    ( &
    mean_chi, stdev_chi, chi_tol, &
    chi_at_liq_sat, cf, rc &
    )

    implicit none

    intrinsic :: exp, sqrt, erf
    real(kind=8), intent(in) :: mean_chi, stdev_chi, chi_at_liq_sat, chi_tol
    real(kind=8), intent(out) :: cf, rc
    real(kind=8), parameter :: pi = 3.1415926
    real(kind=8) :: zeta

    if (stdev_chi > chi_tol) then
      zeta = (mean_chi-chi_at_liq_sat)/stdev_chi
      cf = 0.5*(1.0+erf(zeta/sqrt(2.0)))
      rc = (mean_chi-chi_at_liq_sat)*cf &
            + stdev_chi*exp(-0.5*(zeta**2))/sqrt(2.0*pi)
    else
      if (mean_chi-chi_at_liq_sat < 0.0) then
        cf = 0.0
        rc = 0.0
      else
        cf = 1.0
        rc = mean_chi - chi_at_liq_sat
      end if
    end if

    return
  end subroutine calc_cf_comp

  elemental subroutine pdf_closure &
    ( &
    wp2, wp3, rtp2, thlp2, rtpthlp, &                                   ! in
    wprtp, wpthlp, thlm, rtm, p, &                                      ! in
    wp4, wp2rtp, wprtp2, wp2thlp, wpthlp2, wprtpthlp, &                 ! out
    wp2rcp, wp2thvp, wprcp, wpthvp, &                                   ! out
    rtprcp, rtpthvp, thlprcp, thlpthvp, rcp2, &                         ! out
    uf, w_1, w_2, var_w_1, var_w_2, &                                   ! out
    w_1_n, w_2_n, sig_w_sqd, &                                          ! out
    rt_1, rt_2, var_rt_1, var_rt_2, &                                   ! out
    thl_1, thl_2, var_thl_1, var_thl_2, rrtthl, &                       ! out
    rcm, rc_1, rc_2, rcmi, rc_1i, rc_2i, &                              ! out
    rsatl_1, rsatl_2, &                                                 ! out
    cf, cf1, cf2, cfi, cf1i, cf2i, &                                    ! out
    alpha_thl, alpha_rt, &                                              ! out
    crt_1, crt_2, cthl_1, cthl_2, &                                     ! out
    chi_1, chi_2, stdev_chi_1, stdev_chi_2, &                           ! out
    stdev_eta_1, stdev_eta_2, &                                         ! out
    covar_chi_eta_1, covar_chi_eta_2, &                                 ! out
    corr_chi_eta_1, corr_chi_eta_2 &                                    ! out
    )

    implicit none

    intrinsic :: max, min, sqrt

    real(kind=8), intent(in) :: wp2, wp3, rtp2, thlp2, rtpthlp, &
                                wprtp, wpthlp,  thlm, rtm, p

    real(kind=8), intent(out) :: &
      wp4, wp2rtp, wprtp2, wp2thlp, wpthlp2, wprtpthlp, &
      wp2rcp, wp2thvp, wprcp, wpthvp, &
      rtprcp, rtpthvp, thlprcp, thlpthvp, rcp2, &
      uf, w_1, w_2, var_w_1, var_w_2, &
      w_1_n, w_2_n, sig_w_sqd, &
      rt_1, rt_2, var_rt_1, var_rt_2, &
      thl_1, thl_2, var_thl_1, var_thl_2, rrtthl, &
      rcm, rc_1, rc_2, rcmi, rc_1i, rc_2i, &
      rsatl_1, rsatl_2, &
      cf, cf1, cf2, cfi, cf1i, cf2i, &
      alpha_thl, alpha_rt, &
      crt_1, crt_2, cthl_1, cthl_2, &
      chi_1, chi_2, stdev_chi_1, stdev_chi_2, &
      stdev_eta_1, stdev_eta_2, &
      covar_chi_eta_1, covar_chi_eta_2, &
      corr_chi_eta_1, corr_chi_eta_2

    logical, parameter :: l_calc_ice_supersat_frac = .True.

    real(kind=8), parameter :: w_tol = 2.0e-2, &
                               thl_tol = 1.0e-2, &
                               rt_tol = 1.0e-4, &
                               chi_tol = 1.0e-8, &
                               zero_threshold = 0.0, &
                               l_v = 2.5e6, &
                               r_d = 287.05, &
                               r_v = 461.5, &
                               c_p = 1004.67, &
                               ep = r_d/r_v, &
                               ep1 = (1.0-ep)/ep, &
                               ep2 = 1.0/ep, &
                               T_freeze_K = 273.15, &
                               chi_at_liq_sat = 0.0, &
                               skw_max = 4.5, &
                              !  beta = 1.0, &
                              !  gamma = 0.25
                               beta = 2.4, &
                               gamma = 0.32
                              !  gamma = 0.64

    real(kind=8) :: exner, thv_ds, skw, width_fac_1, width_fac_2, &
                    tl1, tl2, beta1, beta2, &
                    chi_at_ice_sat1, chi_at_ice_sat2, &
                    rt_at_ice_sat1, rt_at_ice_sat2, &
                    rc_coef, uf_max, &
                    wprxp, wp2rxp, thlprxp, rtprxp

    thv_ds = thlm
    exner = (p*1.0e-5)**(r_d/c_p)
    uf_max = 1.0-0.5*(1.0-(skw_max/sqrt(4.0*(1.0-0.4)**3+skw_max**2)))

    ! ansatz for normalized subplume w variance
    sig_w_sqd = max((wpthlp/(sqrt(wp2*thlp2)+0.01*w_tol*thl_tol))**2, &
                    (wprtp/(sqrt(wp2*rtp2)+0.01*w_tol*rt_tol))**2)
    sig_w_sqd = gamma*(1.0-min(sig_w_sqd, 1.0))

    ! PDF closure for w, theta_l and r_t
    if (wp2 <= w_tol*w_tol) then
      uf = 0.5
      w_1 = 0.0
      w_2 = 0.0
      var_w_1 = 0.0
      var_w_2 = 0.0
      rt_1 = rtm
      rt_2 = rtm
      var_rt_1 = 0.0
      var_rt_2 = 0.0
      alpha_rt = 0.5
      thl_1 = thlm
      thl_2 = thlm
      var_thl_1 = 0.0
      var_thl_2 = 0.0
      alpha_thl = 0.5
      rrtthl = 0.0
    else
      skw = wp3/(wp2**1.5)
      call w_closure(skw, wp2, sig_w_sqd, uf_max, &
                     uf, w_1, w_2, w_1_n, w_2_n, var_w_1, var_w_2)
      width_fac_1 = (2.0/3.0)*beta + 2.0*uf*(1.0-(2.0/3.0)*beta)
      width_fac_2 = 2.0 - width_fac_1
      call x_closure(thlm, thlp2, wp2, wpthlp, uf, sig_w_sqd, w_1_n, w_2_n, &
                     width_fac_1, width_fac_2, thl_tol, zero_threshold, &
                     thl_1, thl_2, var_thl_1, var_thl_2, alpha_thl)
      call x_closure(rtm, rtp2, wp2, wprtp, uf, sig_w_sqd, w_1_n, w_2_n, &
                     width_fac_1, width_fac_2, rt_tol, zero_threshold, &
                     rt_1, rt_2, var_rt_1, var_rt_2, alpha_rt)
      if ((var_rt_1*var_thl_1 > 0.0) .and. (var_rt_2*var_thl_2 > 0.0)) then
        rrtthl = (rtpthlp - uf*(rt_1-rtm)*(thl_1-thlm) &
                  - (1.0-uf)*(rt_2-rtm)*(thl_2-thlm)) &
                 / (uf*sqrt(var_rt_1*var_thl_1) &
                    + (1.0-uf)*sqrt(var_rt_2*var_thl_2))
        rrtthl = min(max(rrtthl, -1.0), 1.0)
      else
        rrtthl = 0.0
      end if
    end if

    ! first batch of parameterized HOMs
    wp2rtp = uf*(w_1**2+var_w_1)*(rt_1-rtm) &
             + (1.0-uf)*(w_2**2+var_w_2)*(rt_2-rtm)
    wp2thlp = uf*(w_1**2+var_w_1)*(thl_1-thlm) &
              + (1.0-uf)*(w_2**2+var_w_2)*(thl_2-thlm)
    wp4 = uf*(3.0*(var_w_1**2)+6.0*(w_1**2)*var_w_1+w_1**4) &
          + (1.0-uf)*(3.0*(var_w_2**2)+6.0*(w_2**2)*var_w_2+w_2**4)
    wprtp2 = uf*w_1*((rt_1-rtm)**2+var_rt_1) &
             + (1.0-uf)*w_2*((rt_2-rtm)**2+var_rt_2)
    wpthlp2 = uf*w_1*((thl_1-thlm)**2+var_thl_1) &
              + (1.0-uf)*w_2*((thl_2-thlm)**2+var_thl_2)
    wprtpthlp = uf*w_1*((rt_1-rtm)*(thl_1-thlm) &
                        + rrtthl*sqrt(var_rt_1*var_thl_1)) &
                + (1.0-uf)*w_2*((rt_2-rtm)*(thl_2-thlm) &
                                + rrtthl*sqrt(var_rt_2*var_thl_2))

    ! PDF closure for potential supersaturation
    tl1 = thl_1*exner
    tl2 = thl_2*exner
    rsatl_1 = sat_mixrat_liq(tl1, p, ep, T_freeze_K)
    rsatl_2 = sat_mixrat_liq(tl2, p, ep, T_freeze_K)
    beta1 = ep*(l_v/(r_d*tl1))*(l_v/(c_p*tl1))
    beta2 = ep*(l_v/(r_d*tl2))*(l_v/(c_p*tl2))
    chi_1 = (rt_1-rsatl_1)/(1.0+beta1*rsatl_1)
    chi_2 = (rt_2-rsatl_2)/(1.0+beta2*rsatl_2)
    crt_1 = 1.0/(1.0+beta1*rsatl_1)
    crt_2 = 1.0/(1.0+beta2*rsatl_2)
    cthl_1 = ((1.0+beta1*rt_1)/((1.0+beta1*rsatl_1)**2)) &
             * (c_p/l_v)*beta1*rsatl_1*exner
    cthl_2 = ((1.0+beta2*rt_2)/((1.0+beta2*rsatl_2)**2)) &
             * (c_p/l_v)*beta2*rsatl_2*exner
    stdev_chi_1 = (crt_1**2)*var_rt_1 &
                  - 2.0*rrtthl*crt_1*cthl_1*sqrt(var_rt_1*var_thl_1) &
                  + (cthl_1**2)*var_thl_1
    stdev_chi_1 = sqrt(max(stdev_chi_1, zero_threshold))
    stdev_chi_2 = (crt_2**2)*var_rt_2 &
                  - 2.0*rrtthl*crt_2*cthl_2*sqrt(var_rt_2*var_thl_2) &
                  + (cthl_2**2)*var_thl_2
    stdev_chi_2 = sqrt(max(stdev_chi_2, zero_threshold))
    stdev_eta_1 = (crt_1**2)*var_rt_1 &
                  + 2.0*rrtthl*crt_1*cthl_1*sqrt(var_rt_1*var_thl_1) &
                  + (cthl_1**2)*var_thl_1
    stdev_eta_1 = sqrt(max(stdev_eta_1, zero_threshold))
    stdev_eta_2 = (crt_2**2)*var_rt_2 &
                  + 2.0*rrtthl*crt_2*cthl_2*sqrt(var_rt_2*var_thl_2) &
                  + (cthl_2**2)*var_thl_2
    stdev_eta_2 = sqrt(max(stdev_eta_2, zero_threshold))
    covar_chi_eta_1 = (crt_1**2)*var_rt_1 - (cthl_1**2)*var_thl_1
    covar_chi_eta_2 = (crt_2**2)*var_rt_2 - (cthl_2**2)*var_thl_2
    if (stdev_chi_1 <= chi_tol) then
      stdev_chi_1 = 0.0
    end if
    if (stdev_chi_2 <= chi_tol) then
      stdev_chi_2 = 0.0
    end if
    if (stdev_chi_1*stdev_eta_1 > 0.0) then
      corr_chi_eta_1 = covar_chi_eta_1/(stdev_chi_1*stdev_eta_1)
    else
      corr_chi_eta_1 = 0.0
    end if
    if (stdev_chi_2*stdev_eta_2 > 0.0) then
      corr_chi_eta_2 = covar_chi_eta_2/(stdev_chi_2*stdev_eta_2)
    else
      corr_chi_eta_2 = 0.0
    end if

    ! cloud water and cloud fraction for each plume
    call calc_cf_comp(chi_1, stdev_chi_1, chi_tol, &
                      chi_at_liq_sat, cf1, rc_1)
    call calc_cf_comp(chi_2, stdev_chi_2, chi_tol, &
                      chi_at_liq_sat, cf2, rc_2)
    cf = min(1.0, max(zero_threshold, uf*cf1+(1.0-uf)*cf2))
    rcm = max(uf*rc_1+(1.0-uf)*rc_2, zero_threshold)

    ! cloud ice and ice supersaturation for each plume
    cf1i = 0.0
    rc_1i = 0.0
    cf2i = 0.0
    rc_2i = 0.0
    cfi = 0.0
    rcmi = 0.0
    if (l_calc_ice_supersat_frac) then
      if (tl1<=T_freeze_K) then
        rt_at_ice_sat1 = sat_mixrat_ice(tl1, p, ep, T_freeze_K)
        chi_at_ice_sat1 = (rt_at_ice_sat1-rsatl_1)/(1.0+beta1*rsatl_1)
      else
        chi_at_ice_sat1 = chi_at_liq_sat
      end if
      if (tl2<=T_freeze_K) then
        rt_at_ice_sat2 = sat_mixrat_ice(tl2, p, ep, T_freeze_K)
        chi_at_ice_sat2 = (rt_at_ice_sat2-rsatl_2)/(1.0+beta2*rsatl_2)
      else
        chi_at_ice_sat2 = chi_at_liq_sat
      end if
      call calc_cf_comp(chi_1, stdev_chi_1, chi_tol, &
                        chi_at_ice_sat1, cf1i, rc_1i)
      call calc_cf_comp(chi_2, stdev_chi_2, chi_tol, &
                        chi_at_ice_sat2, cf2i, rc_2i)
      cfi = min(1.0, max(zero_threshold, uf*cf1i+(1.0-uf)*cf2i))
      rcmi = max(uf*rc_1i+(1.0-uf)*rc_2i, zero_threshold)
    end if

    ! hydrometeor loading for buoyancy term calculation
    wprxp = 0.0
    wp2rxp = 0.0
    thlprxp = 0.0
    rtprxp = 0.0

    ! second batch of HOMs involving cloud water
    ! and virtual potential temperature
    rc_coef = l_v/(exner*c_p) - ep2*thv_ds
    wp2rcp = uf*((w_1**2)+var_w_1)*rc_1 &
             + (1.0-uf)*((w_2**2)+var_w_2)*rc_2 &
              - wp2*(uf*rc_1+(1.0-uf)*rc_2)
    wp2thvp = wp2thlp + ep1*thv_ds*wp2rtp + rc_coef*wp2rcp - thv_ds*wp2rxp
    wprcp = uf*w_1*rc_1 + (1.0-uf)*w_2*rc_2
    wpthvp = wpthlp + ep1*thv_ds*wprtp + rc_coef*wprcp - thv_ds*wprxp
    thlprcp = uf*((thl_1-thlm)*rc_1 &
                  + cf1*(rrtthl*crt_1*sqrt(var_rt_1*var_thl_1) &
                         -cthl_1*var_thl_1)) &
              + (1.0-uf)*((thl_2-thlm)*rc_2 &
                          + cf2*(rrtthl*crt_2*sqrt(var_rt_2*var_thl_2)&
                                 - cthl_2*var_thl_2))
    thlpthvp = thlp2 + ep1*thv_ds*rtpthlp +rc_coef*thlprcp - thv_ds*thlprxp
    rtprcp = uf*((rt_1-rtm)*rc_1 &
                 + cf1*(crt_1*var_rt_1 &
                        - rrtthl*cthl_1*sqrt(var_rt_1*var_thl_1))) &
             + (1.0-uf)*((rt_2-rtm)*rc_2 &
                         + cf2*(crt_2*var_rt_2 &
                                - rrtthl*cthl_2*sqrt(var_rt_2*var_thl_2)))
    rtpthvp = rtpthlp + ep1*thv_ds*rtp2 + rc_coef*rtprcp - thv_ds*rtprxp
    rcp2 = max(uf*(chi_1*rc_1+cf1*(stdev_chi_1**2)) &
               + (1.0-uf)*(chi_2*rc_2+cf2*(stdev_chi_2**2)) - rcm**2, &
               zero_threshold)
    return
  end subroutine pdf_closure
end module pdf
