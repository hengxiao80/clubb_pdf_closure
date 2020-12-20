program main

  use pdf, only : pdf_closure

  implicit none

  ! integer, parameter :: nz = 224, nt = 360 ! LBA
  integer, parameter :: nz = 75, nt = 360 ! BOMEX
  ! integer, parameter :: nz = 128, nt = 1800 ! ARM_97_first_event

  real(kind=8), dimension(nz, nt) :: &
    rtm, thlm, wp2, wp3, wprtp, wpthlp, &
    rtp2, thlp2, rtpthlp, p

  real(kind=8), dimension(nz, nt) :: &
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

  ! open(70, file='../data/input_lba_r1.bin', form='unformatted', &
  !          access='sequential')
  ! open(71, file='../data/output_lba_r1_s1.bin', form='unformatted', &
  !          access='sequential')
  open(70, file='../data/input_bomex_r3.bin', form='unformatted', &
           access='sequential')
  open(71, file='../data/output_bomex_r3.bin', form='unformatted', &
           access='sequential')
  ! open(70, file='../data/input_arm_97_first_event_r5.bin', &
  !          form='unformatted', access='sequential')
  ! open(71, file='../data/output_arm_97_first_event_r5_s1.bin', &
  !          form='unformatted', access='sequential')

  read(70) rtm
  read(70) thlm
  read(70) wp2
  read(70) wp3
  read(70) wprtp
  read(70) wpthlp
  read(70) rtp2
  read(70) thlp2
  read(70) rtpthlp
  read(70) p

  call pdf_closure &
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

  write(71) rtm
  write(71) thlm
  write(71) wp2
  write(71) wp3
  write(71) wprtp
  write(71) wpthlp
  write(71) rtp2
  write(71) thlp2
  write(71) rtpthlp
  write(71) p
  write(71) wp4
  write(71) wp2rtp
  write(71) wprtp2
  write(71) wp2thlp
  write(71) wpthlp2
  write(71) wprtpthlp
  write(71) wp2rcp
  write(71) wp2thvp
  write(71) wprcp
  write(71) wpthvp
  write(71) rtprcp
  write(71) rtpthvp
  write(71) thlprcp
  write(71) thlpthvp
  write(71) rcp2
  write(71) uf
  write(71) w_1
  write(71) w_2
  write(71) var_w_1
  write(71) var_w_2
  write(71) w_1_n
  write(71) w_2_n
  write(71) sig_w_sqd
  write(71) rt_1
  write(71) rt_2
  write(71) var_rt_1
  write(71) var_rt_2
  write(71) thl_1
  write(71) thl_2
  write(71) var_thl_1
  write(71) var_thl_2
  write(71) rrtthl
  write(71) rcm
  write(71) rc_1
  write(71) rc_2
  write(71) rcmi
  write(71) rc_1i
  write(71) rc_2i
  write(71) rsatl_1
  write(71) rsatl_2
  write(71) cf
  write(71) cf1
  write(71) cf2
  write(71) cfi
  write(71) cf1i
  write(71) cf2i
  write(71) alpha_thl
  write(71) alpha_rt
  write(71) crt_1
  write(71) crt_2
  write(71) cthl_1
  write(71) cthl_2
  write(71) chi_1
  write(71) chi_2
  write(71) stdev_chi_1
  write(71) stdev_chi_2
  write(71) stdev_eta_1
  write(71) stdev_eta_2
  write(71) covar_chi_eta_1
  write(71) covar_chi_eta_2
  write(71) corr_chi_eta_1
  write(71) corr_chi_eta_2

  close(71)
  close(70)

end program main
