program flt_haz

    use const_module
    use input_module
    use flt_module
    use GMPE_module

    implicit none

    real(8), allocatable :: Haz(:,:)
    real(8), allocatable :: Haz_Bin (:,:,:,:,:)
    real(8) :: haz_part, p_exceed
    real(8) :: lnSa, Sigma, m_eps
    real(8) :: t0, t1, t2

    call cpu_time(t0)
    call read_input()
    call unit_conversion ()
    call flt_ini()
    call mag_freq_distribution ()

    allocate(Haz(n_inten, n_freq))
    allocate(Haz_Bin(n_eps_bin, n_dist_bin, n_mag_bin, n_inten, n_freq))
    Haz = 0.0
    Haz_bin = 0.0

    call cpu_time(t1)
    write(*,*) new_line('a')//new_line('a')
    write(*,'("Preprocess Time = ",f6.3," seconds.")') t1 - t0
    write(*,*) new_line('a')

    write(fp_rup,*) '      Mw      Rrup    Rjb     Rx'

    do i_mag = 1, n_mag

        Mw = mag_inc(i_mag); rate = rate_inc(i_mag)

        call locate(i_mag_bin, mag_bin,  Mw)

        call rupture_location() ! discrete the fault per rupture length

        p_locS = 1.0d0 / real(n_locS, 8)

        do i_locS = 1, n_locS

            call cal_p_locD_arr()
            call locate_rupture(S1(i_locS),S2(i_locS),rup_coor) !

            do i_locD = 1, n_locD
                p_locD = p_locD_arr(i_locD)

                call cal_coor_d()

                call dist_rup_set(Rrup, Rjb, Rx, coor_D, &
                    rup_top(i_locD), flt_strike_deg, flt_dip_deg, rup_wid)

                call locate(i_dist_bin, dist_bin,  Rrup)

                do i_freq = 1, n_freq
                    Tin = 1.0d0 / frequency(i_freq)

                    call GMPE_interface(m_gmpe_name, Tin, Mw, m_sof, Rrup, Rjb, Rx, &
                        rup_top(i_locD), flt_dip_deg, &
                        Vs30, Z10, gmpe_params, gmpe_opts, &
                        lnSa, Sigma)

                    do i_inten = 1, n_inten

                        m_eps = (log(intensity(i_inten)) - lnSa) / Sigma

                        call locate(i_eps_bin, eps_bin,  m_eps)

                        call prob_exceed (p_exceed, m_eps, m_aleatory_distribution, trunc_level)
                        haz_part = rate * p_locS * p_locD * p_exceed

                        Haz(i_inten, i_freq) = Haz(i_inten, i_freq) + haz_part
                        Haz_bin(i_eps_bin, i_dist_bin, i_mag_bin, i_inten, i_freq) = &
                            Haz_bin(i_eps_bin, i_dist_bin, i_mag_bin, i_inten, i_freq) + haz_part

                    end do ! end n_inten
                end do ! end n_freq

                if (m_wrt_rup .eq. True) write(fp_rup,'(1x,f10.4, 3f8.2, f10.4)') Mw, Rrup, Rjb, Rx, exp(lnSa)
            end do ! end n_locD, loop down dip
        end do ! end n_locS, loop along strike
    end do ! end n_mag, loop for magnitude step

    call cpu_time(t2)
    write(*,'("Hazard Integration Time = ",f8.3," seconds.")') t2 - t1

    call print_haz(haz)
    call print_haz_bin(haz_bin)

    call close_file()

end program flt_haz
