module flt_module

    use const_module
    use input_module
    use utils

    implicit none

    real(8), allocatable :: flt_len_seg(:)
    real(8), allocatable :: flt_az_seg(:)

    real(8), allocatable :: flt_coor(:,:)
    real(8), allocatable :: flt_S_corner(:)
    real(8)              :: site_coor(2)

    real(8), allocatable :: rup_top(:)
    real(8), allocatable :: rup_coor(:,:), coor_D(:,:)
    real(8), allocatable :: S1(:), S2(:)
    real(8), allocatable :: p_locD_arr(:)
    real(8), allocatable :: mag_inc_0(:), rate_inc_0(:) ! lower bound of zero
    real(8), allocatable :: mag_inc(:), rate_inc(:) ! lower bound of Mmin


    real(8) :: flt_area, flt_len, flt_wid

    real(8) :: flt_strike_deg, flt_strike_rad
    real(8) :: step_D, step_S
    real(8) :: step_D_trial, step_S_trial

    real(8) :: rup_len, rup_wid, rup_area
    real(8) :: rup_len_trial, rup_wid_trial, rup_area_trial

    integer :: n_locD, i_locD
    integer :: n_locS, i_locS
    integer :: n_cor, i_seg
    real(8) :: step_D_V, step_D_H
    real(8) :: step_D_Hc, step_D_Hs
    real(8) :: p_locS, p_locD, ftop

    real(8) :: Mw, rate, Rrup, Rjb, Rx

    integer :: i_mag, n_mag, i_mag_bin, i_dist_bin, i_eps_bin
    integer :: i_freq, i_inten
    real(8) :: Tin

contains

    subroutine mag_freq_distribution( )

        real(8) :: beta, Af, b

        real(8) :: c
        real(8) :: c2
        real(8) :: d
        real(8) :: dmag2
        real(8) :: K
        real(8) :: m_local
        real(8) :: M0_Mmax
        real(8) :: Mmin0
        real(8) :: mu
        real(8) :: muAS
        real(8), allocatable :: mag_cum(:), rate_cum (:)
        real(8), allocatable :: rate_inc_0(:), mag_inc_0(:)

        real(8) :: N_Mmin
        ! magnitude frequency distribution for characteristic recurrence relation

        Mmin0 = 0.0 ! Mmin alway set as zero
        b = b_value
        beta = b * log(10.0)
        Af = flt_area * 1.0d10 ! convert km2 to cm**2
        c = 1.5
        d = 16.05
        mu = 3.0d11 ! dyne/cm2

        muAS = mu * Af * slip_rate * 0.1 ! slip rate in mm/yr here

        M0_Mmax = 10**(Mmax * c + d)

        select case ( m_rec_relation) !

            case ( EXPONENTIAL) !
                call mfd_exp()
            case ( CHARACTERISTIC) !
                call mfd_char()
            case ( DELTA) !
                call mfd_delta()

            case default
                stop 'recurrence relation not supported'
        end select

        mag_inc = pack(mag_inc_0, mag_inc_0 > Mmin)
        rate_inc = pack (rate_inc_0, mag_inc_0 > Mmin)
        n_mag = size(mag_inc)

    contains

        subroutine mfd_delta()

            allocate(mag_inc_0(1), rate_inc_0(1))
            mag_inc_0(1) = Mmax
            rate_inc_0(1) = muAS / M0_Mmax

        end subroutine mfd_delta

        !>
        subroutine mfd_exp()

            ! Arguments declarations

            ! Variable declarations

            real(8) :: nu
            real(8) :: tmp
            !

            tmp = exp(-beta * (Mmax - Mmin0))

            nu = muAS * (1.0d0 - tmp) * (c - b) / b / M0_Mmax / tmp

            n_mag = nint(Mmax / mag_step) + 1;

            allocate(mag_cum(n_mag))

            do i_mag = 1,n_mag
                mag_cum (i_mag)= dble(i_mag - 1) * mag_step
            end do


            allocate(rate_cum(n_mag))
            rate_cum = 0.0d0


            do i_mag=1,n_mag

                m_local = mag_cum(i_mag)
                rate_cum(i_mag) = nu * (exp(-beta * (m_local - Mmin0)) - tmp) / (1.0 - tmp)

            end do

            dmag2 = mag_step / 2.0d0

            allocate(mag_inc_0(n_mag - 1))
            do i_mag = 1, (n_mag - 1)
                mag_inc_0(i_mag) = dmag2 + dble(i_mag-1) * mag_step
            end do
            !m2f: mag_inc = mag_inc'

            allocate(rate_inc_0(n_mag - 1))
            rate_inc_0 = rate_cum(1:n_mag-1) - rate_cum(2:n_mag)


        end subroutine mfd_exp

        subroutine mfd_char()
            implicit none

            real(8) :: nu_C
            real(8) :: nu_NC
            real(8) :: tmp_05
            real(8) :: tmp_15
            real(8) :: tmp_c

            tmp_15 = exp(-beta * (Mmax - Mmin0 -1.5))
            tmp_05 = exp(-beta * (Mmax - Mmin0 -0.5))

            c2 = 0.5 * beta * tmp_15 / (1.0 - tmp_05)

            M0_Mmax = 10**(Mmax * c + d)
            tmp_c = 10**(-c/2.0)

            K = b * tmp_c /(c-b) + b * exp(beta) * (1 - tmp_c) / c
            N_Mmin = muAS * (1+c2) * (1-tmp_05) / M0_Mmax / K / tmp_05

            nu_NC = muAS * (1.0 - tmp_05) / M0_Mmax / K / tmp_05
            nu_C = nu_NC * 0.5 * beta * tmp_15 / (1.0 - tmp_05)

            n_mag = nint(Mmax / mag_step) + 1;

            allocate(mag_cum(n_mag))

            do i_mag = 1,n_mag
                mag_cum (i_mag)= dble(i_mag - 1) * mag_step
            end do

            allocate(rate_cum(n_mag))
            rate_cum = 0.0d0

            do i_mag=1,n_mag

                m_local = mag_cum(i_mag)

                if (m_local <= Mmax - 0.5) then
                    rate_cum(i_mag) = nu_C + nu_NC * &
                        ((exp(-beta*(m_local-Mmin0)) - tmp_05)/(1.0 - tmp_05))
                else
                    rate_cum(i_mag) = nu_C * (Mmax - m_local) / 0.5
                end if

            end do


            dmag2 = mag_step / 2.0
            allocate(mag_inc_0(n_mag - 1))
            do i_mag = 1, (n_mag - 1)
                mag_inc_0(i_mag) = dmag2 + dble(i_mag-1) * mag_step
            end do

            allocate(rate_inc_0(n_mag - 1))
            rate_inc_0 = rate_cum(1:n_mag-1) - rate_cum(2:n_mag)

        end subroutine mfd_char

    end subroutine mag_freq_distribution

    subroutine unit_conversion()

        if (m_unit .eq. DEG) then
            call deg2km_model()
        else if(m_unit .eq. KM ) then
            call align_model()
        else
            stop 'wrong unit, check UNIT section in input '
        end if

    end subroutine


    subroutine calDepthProb ()
        implicit none

        real(8) :: focal_depth
        real(8) :: half_wid
        real(8) :: z1
        real(8) :: z2
        real(8) :: z3
        integer(8) :: i_locD


        z1 = Smin
        z2 = Smin + (Smax-Smin) * depth_param
        z3 = Smax


        half_wid = rup_wid * sin(flt_dip_rad) / 2.0

        do i_locD=1,n_locD

            focal_depth = rup_top(i_locD) + half_wid
            if (focal_depth < z1) then
                p_locD_arr(i_locD) = 0.
            else if ( focal_depth < z2 ) then
                p_locD_arr(i_locD) = (focal_depth - z1) / (z2 - z1)
            else if ( focal_depth < z3 ) then
                p_locD_arr(i_locD) = 1. - (focal_depth - z2) / (z3 - z2)
            else if ( focal_depth >= z3 ) then
                p_locD_arr(i_locD) = 0.
            end if

        end do

        p_locD_arr = p_locD_arr / sum(p_locD_arr)

    end subroutine  calDepthProb


    subroutine deg2km_model()

        integer :: i_corner

        allocate(flt_coor(flt_n_corner,2))

        do i_corner=1,flt_n_corner

            call deg2km_simple(flt_coor(i_corner,1), flt_coor(i_corner,2), &
                flt_trace(i_corner,1), flt_trace(i_corner,2), &
                site(1), site(2))
        end do

        site_coor = [ 0.0, 0.0 ]

    end subroutine deg2km_model

    subroutine align_model()
        allocate(flt_coor(flt_n_corner,2))
        flt_coor(:,1) = flt_trace(:,1) - site(1)
        flt_coor(:,2) = flt_trace(:,2) - site(2)
        site_coor = [ 0.0,0.0 ]
    end subroutine  align_model


    subroutine flt_ini()
        real(8) ::  x1, y1, x2, y2, xN, yN
        integer :: i_seg

        flt_n_corner = size(flt_trace) / 2
        flt_n_seg = flt_n_corner - 1

        allocate(flt_len_seg(flt_n_seg))
        allocate(flt_az_seg(flt_n_seg))
        allocate(flt_s_corner(flt_n_corner))

        flt_s_corner(1) = 0.0d0

        do i_seg = 1, flt_n_seg
            y1 = flt_coor(i_seg, 1);
            x1 = flt_coor(i_seg, 2);
            y2 = flt_coor(i_seg + 1, 1);
            x2 = flt_coor(i_seg + 1, 2);

            call delaz2_km (y1,x1,y2,x2, flt_len_seg(i_seg), flt_az_seg(i_seg))
            flt_s_corner(i_seg+1) = flt_s_corner(i_seg) + flt_len_seg(i_seg)
        end do

        flt_len = sum(flt_len_seg)
        flt_wid = (Smax - Smin) / sin(flt_dip_rad)
        flt_area = flt_len * flt_wid

        Y1 = flt_coor(1,1)
        X1 = flt_coor(1,2)
        Yn = flt_coor(flt_n_corner,1)
        Xn = flt_coor(flt_n_corner,2)

        flt_strike_rad = atan2(Xn - X1, Yn - Y1)   ! in radian
        flt_strike_deg = flt_strike_rad * RAD2DEG ! in degree


    end subroutine flt_ini

    subroutine cal_p_locD_arr()
        if (allocated(p_locD_arr)) deallocate (p_locD_arr)

        allocate(p_locD_arr(n_locD))
        p_locD_arr = 0.0d0

        select case(m_DEPTH_distribution)
            case (UNIFORM)
                p_locD_arr = p_locD_arr + 1.0d0 / real(n_locD, 8)
            case (TRIANGULAR)
                call calDepthProb()
            case default
                stop 'unsupported depth distribution'

        end select


    end subroutine

    subroutine cal_coor_d()
        n_cor = size(rup_coor) / 2

        if(allocated(coor_D)) deallocate(coor_D)
        allocate(coor_D(n_cor,2))
        coor_D = 0.0d0

        do i_seg = 1, n_cor

            coor_D(i_seg,1) = rup_coor(i_seg,1) + real((i_locD - 1),8) * step_D_Hc;
            coor_D(i_seg,2) = rup_coor(i_seg,2) + real((i_locD - 1),8) * step_D_Hs;
        end do ! end n_cor


    end subroutine


    subroutine rupture_location ()

        integer :: i_locD, i_locS

        step_S_trial = strike_step
        step_D_trial = dip_step

        call Mw2Arup()

        if (rup_area_trial > flt_area) then ! break the entire fault

            n_locS = 1; n_locD = 1

            step_S = 0.0d0; step_D = 0.0d0

            if (allocated(S1)) deallocate(S1)
            allocate(S1(n_locS))

            if (allocated(S2)) deallocate(S2)
            allocate(S2(n_locS))

            S1(1) = 0.0d0;  S2(1) = flt_len

            if (allocated(rup_top)) deallocate(rup_top)
            allocate(rup_top(n_locD))

            rup_top(1) = Smin
            rup_wid = flt_wid
            rup_len = flt_len
            rup_area = rup_area_trial

        else if (rup_wid_trial > flt_wid) then ! break the whole width

            n_locD = 1
            step_D = 0

            if (allocated(rup_top)) deallocate(rup_top)
            allocate(rup_top(n_locD))

            rup_wid = flt_wid
            rup_top(1) = Smin
            rup_len = rup_area_trial / rup_wid

            if (rup_len > flt_len) then

                rup_len = flt_len
                n_locS = 1
                step_S = 0

                if (allocated(S1)) deallocate(S1)
                allocate(S1(n_locS))

                if (allocated(S2)) deallocate(S2)
                allocate(S2(n_locS))

                S1(1) = 0.0d0;  S2(1) = flt_len

            else

                n_locS = floor((flt_len - rup_len) / step_S_trial) + 1
                step_S = (flt_len - rup_len) / n_locS

                n_locS = n_locS + 1

                !m2f: S1 = zeros(n_locS,n_locD)
                if (allocated(S1)) deallocate(S1)
                allocate(S1(n_locS))
                S1 = 0.0d0

                !m2f: S2 = zeros(n_locS,n_locD)
                if (allocated(S2)) deallocate(S2)
                allocate(S2(n_locS))
                S2 = 0.0d0

                do i_locS=1,n_locS
                    S1(i_locS) = (i_locS - 1)* step_S
                    S2(i_locS) = S1(i_locS) + rup_len
                end do

            end if

        else

            rup_wid = rup_wid_trial
            rup_len = rup_len_trial

            n_locS = floor((flt_len - rup_len) / step_S_trial) + 1
            step_S = (flt_len - rup_len) / n_locS

            n_locD = floor((flt_wid - rup_wid) / step_D_trial) + 1
            step_D = (flt_wid - rup_wid) / n_locD

            step_D_V = step_D * sin(flt_dip_rad)


            n_locS = n_locS + 1
            n_locD = n_locD + 1

            !m2f: S1 = zeros(n_locS,1)
            if (allocated(S1)) deallocate(S1)
            allocate(S1(n_locS))
            S1 = 0.0d0

            !m2f: S2 = zeros(n_locS,1)
            if (allocated(S2)) deallocate(S2)
            allocate(S2(n_locS))
            S2 = 0.0d0


            do i_locS=1,n_locS

                S1(i_locS) = (i_locS - 1)* step_S
                S2(i_locS) = S1(i_locS) + rup_len

            end do

            !m2f: rup_top = zeros(n_locD,1)
            if (allocated(rup_top)) deallocate(rup_top)
            allocate(rup_top(n_locD))

            rup_top = 0.0d0

            do i_locD=1,n_locD
                rup_top(i_locD) = Smin + (i_locD - 1) * step_D_V
            end do

        end if

        if (S2(n_locS) > flt_len) then
            S2(n_locS) = flt_len * 0.99999
        end if

        step_D_V = step_D * sin(flt_dip_rad)
        step_D_H = step_D * cos(flt_dip_rad)

        step_D_Hc = step_D_H * cos(flt_strike_rad + PI / 2.0d0)
        step_D_Hs = step_D_H * sin(flt_strike_rad + PI / 2.0d0)

    end subroutine

    subroutine locate_rupture( S1_local, S2_local , rup_coor)

        ! Arguments declarations
        real(8), allocatable, intent(out) :: rup_coor(:,:)    ! m2f:check dim(Ntot
        real(8), intent(in) :: S1_local, S2_local
        ! Variable declarations
        integer :: N1, N2, Ntot , i
        real(8) :: X1, X2, Y1, Y2          !


        ! S1_local and S2_local are the starting and ending points of the rupture, measured as
        ! distance to the first corner of the fault
        ! rup_coor is the fault surface trace between S1_local and S2_local

        N1 = count(S1_local >= flt_S_corner)
        Y1 = flt_coor(N1,1) + (S1_local - flt_S_corner(N1)) * cos(flt_az_seg(N1))
        X1 = flt_coor(N1,2) + (S1_local - flt_S_corner(N1)) * sin(flt_az_seg(N1))

        if (S2_local > flt_S_corner(flt_n_corner)) then
            stop 'S2_local cannot be greater than fault length'
        end if

        N2 = count(S2_local > flt_S_corner)
        Y2 = flt_coor(N2,1) + (S2_local - flt_S_corner(N2)) * cos(flt_az_seg(N2))
        X2 = flt_coor(N2,2) + (S2_local - flt_S_corner(N2)) * sin(flt_az_seg(N2))

        Ntot = N2 - N1 + 2

        !m2f: rup_coor = zeros(Ntot, 2)
        if (allocated(rup_coor)) deallocate(rup_coor)
        allocate(rup_coor(Ntot, 2))


        if (N1 == N2) then
            rup_coor(1,:) = [Y1, X1]
            rup_coor(2,:) = [Y2, X2]
        else
            rup_coor(1,:) = [ Y1,X1 ]
            rup_coor(Ntot,:) = [ Y2,X2 ]

            do i=2,(Ntot - 1)
                rup_coor(i,1) = flt_coor(N1+i-1,1)
                rup_coor(i,2) = flt_coor(N1+i-1,2)
            end do

        end if

    end subroutine  locate_rupture

    subroutine Mw2Arup()


        ! Variable declarations
        real(8) :: AR2
        real(8) :: A_rup

        ! AR is aspect ratio (length to width ratio)

        select case ( m_scaling) !

            case ( PEER) !
                A_rup = 10**(Mw - 4.0)
            case (CEUS)
                A_rup = 10**(Mw - 4.36)
            case (POINT)
                A_rup = 1.0d-6

            case ( WC94) !

                select case ( m_SOF) !

                    case (SS) !
                        A_rup = 10**(0.9 * Mw - 3.42)
                    case ( RV) !
                        A_rup = 10**(0.98 * Mw - 3.99)
                    case (NM) !
                        A_rup = 10**(0.82 * Mw - 2.87)
                    case (NA) !
                        A_rup = 10**(0.91 * Mw - 3.49)
                    case default
                        stop 'wrong source mechanism '
                end select


            case default
                stop 'wrong scaling model '

        end select

        AR2 = sqrt(Aspect_Ratio)
        rup_len_trial = sqrt(A_rup)*AR2
        rup_wid_trial = sqrt(A_rup)/AR2
        rup_area_trial = A_rup

    end subroutine  Mw2Arup

end module flt_module
