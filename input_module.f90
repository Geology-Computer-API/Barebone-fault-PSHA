module input_module
    use const_module
    implicit none

    integer :: fp_inp
    integer :: fp_log
    integer :: fp_haz
    integer :: fp_dag
    integer :: fp_rup

    integer :: ppos

    logical :: inp_exist

    character(130) :: fnm_inp, arg
    character(130) :: fnm_log
    character(130) :: fnm_haz
    character(130) :: fnm_dag
    character(130) :: fnm_rup

    integer :: eastat, iost
    character(130) :: line, wrt_fmt, str_tmp, gmpe_name
    character(3) :: ext_log = 'log', ext_haz = 'haz', ext_dag = 'dag', ext_rup = 'rup'

    real(8) :: site(2)
    real(8), allocatable :: frequency(:), intensity(:)
    real(8), allocatable :: mag_bin(:), dist_bin(:), eps_bin(:)
    real(8), allocatable :: flt_trace(:,:)
    real(8), allocatable :: gmpe_params(:)
    integer, allocatable :: gmpe_opts(:)

    real(8) :: temp2(500,2), temp1(500)
    integer :: temp_int(500), tmp_int
    real(8) :: tmp1, tmp2
    integer :: numvalues

    real(8) :: slip_rate, b_value, trunc_level
    real(8) :: VS30, Smin, Smax, flt_dip_deg, flt_dip_rad, aspect_ratio
    real(8) :: Z10, Z25
    real(8) :: strike_step, dip_step, depth_param, mag_step
    real(8) :: Mmin, Mmax

    integer :: m_SOF, m_SCALING, m_REC_RELATION, m_UNIT
    integer :: m_SIGMA_TRUNC, m_GMPE_Name, m_DEPTH_distribution
    integer :: m_aleatory_distribution, m_WRT_RUP

    integer :: n_freq, n_inten, n_mag_bin, n_dist_bin, n_eps_bin
    integer :: flt_n_corner, flt_n_seg

contains

    subroutine read_input()

        if (command_argument_count() .ne. 1) then
            stop 'usage: flt_haz inputfile'
        end if

        call get_command_argument(1, arg)
        fnm_inp = trim(adjustl(arg))

        inquire(file = fnm_inp, exist = inp_exist)

        if (inp_exist) then

            ppos = scan(fnm_inp,".", BACK= .true.)

            if (ppos > 0) then
                fnm_log = fnm_inp(1:ppos) // ext_log
                fnm_haz = fnm_inp(1:ppos) // ext_haz
                fnm_dag = fnm_inp(1:ppos) // ext_dag
                fnm_rup = fnm_inp(1:ppos) // ext_rup
            else
                fnm_log = fnm_inp // '.' // ext_log
                fnm_haz = fnm_inp // '.' // ext_haz
                fnm_dag = fnm_inp // '.' // ext_dag
                fnm_rup = fnm_inp // '.' // ext_rup

            end if


            open(newunit = fp_inp,file=fnm_inp,status='old',action='read')
            open(newunit = fp_log,file=fnm_log,status='replace',action='write')
            open(newunit = fp_haz,file=fnm_haz,status='replace',action='write')
            open(newunit = fp_dag,file=fnm_dag,status='replace',action='write')
            open(newunit = fp_rup,file=fnm_rup,status='replace',action='write')

            ReadLoop: do
                read(fp_inp,*,iostat=iost) line !read data line by line

                if (iost < 0) exit !end of file reached

                if (index(line, 'SITE').EQ.1) then
                    call read_site()

                else if (index(line, 'FREQUENCY').EQ.1) then
                    call read_frequency()

                else if (index(line, 'FAULT_TRACE').EQ.1) then
                    call read_fault_trace()

                else if (index(line, 'REC_RELATION').EQ.1) then
                    call read_rec_relation()

                else if (index(line, 'SLIP_RATE').EQ.1) then
                    call read_slip_rate()

                else if (index(line, 'B_VALUE').EQ.1) then
                    call read_b_value()

                else if (index(line, 'SOF').EQ.1) then
                    call read_SOF()

                else if (index(line, 'UNIT').EQ.1) then
                    call read_UNIT()

                else if (index(line, 'ALEATORY_DISTRIBUTION').EQ.1) then
                    call read_aleatory_distribution()

                else if (index(line, 'TRUNC_LEVEL').EQ.1) then
                    call read_trunc_level()

                else if (index(line, 'SCALING_MODEL').EQ.1) then
                    call read_scaling_model()

                else if (index(line, 'DIP').EQ.1) then
                    call read_dip()

                else if (index(line, 'GMPE_NAME').EQ.1) then
                    call read_gmpe_name()

                else if (index(line, 'VS30').EQ.1) then
                    call read_VS30()

                else if (index(line, 'SEISMOGENIC_DEPTH').EQ.1) then
                    call read_seismogenic_depth()

                else if (index(line, 'ASPECT_RATIO').EQ.1) then
                    call read_aspect_ratio()


                else if (index(line, 'INTENSITY').EQ.1) then
                    call read_intensity()

                else if (index(line, 'MAG_BIN').EQ.1) then
                    call read_mag_bin()

                else if (index(line, 'DIST_BIN').EQ.1) then
                    call read_dist_bin()

                else if (index(line, 'EPS_BIN').EQ.1) then
                    call read_eps_bin()

                else if (index(line, 'MAG_STEP') .eq. 1) then
                    call read_mag_step()

                else if (index(line, 'MAG_RANGE') .eq. 1) then
                    call read_mag_range()

                else if (index(line, 'STRIKE_DIP_STEP') .eq. 1) then
                    call read_strike_dip_step()

                else if (index(line, 'DEPTH_DISTRIBUTION') .eq. 1) then
                    call read_depth_distribution()

                else if (index(line, 'DEPTH_PARAM') .eq. 1) then
                    call read_depth_param()

                else if (index(line, 'Z10') .eq. 1) then
                    call read_Z10()

                else if (index(line, 'Z25') .eq. 1) then
                    call read_Z25()

                else if (index(line, 'GMPE_PARAMS') .eq. 1) then
                    call read_gmpe_params()

                else if (index(line, 'GMPE_OPTS') .eq. 1) then
                    call read_gmpe_opts()

                else if (index(line, 'WRT_RUP') .eq. 1) then
                    call read_wrt_rup()

                end if

            end do ReadLoop

        else
            stop 'can not open the inputfile, please check'

        end if

        close(fp_inp)

        n_freq = size(frequency)
        n_inten = size(intensity)
        n_mag_bin = size(mag_bin)
        n_dist_bin = size(dist_bin)
        n_eps_bin = size(eps_bin)
        flt_n_corner = size(flt_trace) / 2
        flt_n_seg = flt_n_corner - 1


    end subroutine

    subroutine read_site()
        read(fp_inp, *) site(1), site(2)
        wrt_fmt = '("site coordinate is ", 2f10.3)'
        write(fp_log,wrt_fmt) site(1), site(2)
        write(fp_log,*) new_line('A')

    end subroutine read_site

    subroutine read_frequency()
        write(fp_log,*) "FREQUENCY found"
        numvalues = 0
        !line = 'null'
        InnerLoop: do
            read(fp_inp,*,iostat=eastat) tmp1
            if (eastat .ne. 0) exit InnerLoop

            numvalues = numvalues+1
            temp1(numvalues) = tmp1

        end do InnerLoop
        write(fp_log,*) "number of frequency points = ", numvalues
        write(fp_log,*) new_line('a')
        allocate(frequency(numvalues))
        frequency (1:numvalues) = temp1(1:numvalues)

    end subroutine read_frequency

    subroutine read_fault_trace()
        write(fp_log,*) "fault_trace found"
        numvalues = 0
        !line = 'null'
        InnerLoop: do
            read(fp_inp,*,iostat=eastat) tmp1, tmp2
            if (eastat .ne. 0) exit InnerLoop

            numvalues = numvalues+1
            temp2(numvalues,1) = tmp1
            temp2(numvalues,2) = tmp2

        end do InnerLoop
        write(fp_log,*) "number of fault trace points = ", numvalues
        write(fp_log,*) new_line('a')
        allocate(flt_trace(numvalues,2))
        flt_trace (1:numvalues,:) = temp2(1:numvalues,:)

    end subroutine read_fault_trace

    subroutine read_rec_relation()
        write(fp_log,*) "recurrence relation found"//new_line('a')
        read(fp_inp, *) str_tmp
        write(fp_log,*) "recurrence relation is " // str_tmp
        if (str_tmp .eq. 'EXP') then
            m_rec_relation = EXPONENTIAL

        else if (str_tmp .eq. 'CHAR') then
            m_rec_relation = CHARACTERISTIC

        else if (str_tmp .eq. 'DELTA') then
            m_REC_RELATION = DELTA

        else
            write(fp_log,*) 'wrong input in rec_relation'
            stop
        end if



    end subroutine read_rec_relation

    subroutine read_slip_rate()
        write(fp_log,*) 'slip rate found'//new_line('a')
        read(fp_inp,*) slip_rate
    end subroutine read_slip_rate

    subroutine read_b_value()
        write(fp_log,*) 'b-value found'//new_line('a')
        read(fp_inp,*) b_value
    end subroutine

    subroutine read_sof ()
        write(fp_log,*) 'style of faulting found in input'//new_line('a')
        read(fp_inp,*) str_tmp
        if(str_tmp .eq. 'SS') then
            m_sof = SS

        else if (str_tmp .eq. 'RV') then
            m_sof = RV

        else if (str_tmp .eq. 'NM') then
            m_sof = NM

        else if (str_tmp .eq. 'NA') then
            m_sof = NA

        else
            write(fp_log,*) 'wrong input in style of faulting'
            stop

        end if
    end subroutine read_sof

    subroutine read_unit()
        write(fp_log,*) 'unit found in input '//new_line('a')
        read(fp_inp,*) str_tmp
        if(str_tmp .eq. 'DEG') then
            m_unit = DEG
        else if (str_tmp .eq. 'KM') then
            m_unit = KM
        else
            write(fp_log,*) 'wrong input in unit'
            stop
        end if
    end subroutine read_unit

    subroutine read_wrt_rup()
        write(fp_log,*) 'wrt_rup found in input '//new_line('a')
        read(fp_inp,*) str_tmp
        if(str_tmp .eq. 'NO') then
            m_wrt_rup = FALSE
        else if (str_tmp .eq. 'YES') then
            m_wrt_rup = TRUE
        else
            write(fp_log,*) 'wrong input in wrt_rup'
            stop
        end if
    end subroutine read_wrt_rup


    subroutine read_aleatory_distribution()
        write(fp_log,*) 'aleatory ditribution found in input'//new_line('a')
        read(fp_inp,*) str_tmp
        if (str_tmp .eq. 'NORMAL') then
            m_aleatory_distribution = NORMAL

        else if (str_tmp .eq. 'TRUNC_NORMAL') then
            m_aleatory_distribution = TRUNC_NORMAL

        else if (str_tmp .eq. 'HEAVISIDE') then
            m_aleatory_distribution = HEAVISIDE
        else
            write(fp_log,*) 'wrong input in aleatory distribution'
            stop 'wrong input in aleatory distribution'

        end if

    end subroutine read_aleatory_distribution

    subroutine read_trunc_level()
        write (fp_log,*) 'trunc_level found in input'//new_line('a')
        read(fp_inp,*) trunc_level

    end subroutine read_trunc_level

    subroutine read_scaling_model
        write(fp_log,*) 'scaling model found in input' // new_line('a')
        read(fp_inp,*) str_tmp
        if (str_tmp .eq. 'WC94') then
            m_scaling = WC94
        else if (str_tmp .eq. 'PEER') then
            m_scaling = PEER
        else if (str_tmp .eq. 'CEUS') then
            m_scaling = CEUS
        else if (str_tmp .eq. 'POINT') then
            m_scaling = POINT
        else
            write(fp_log, *) 'wrong input in scaling model'
            write(*,*) 'wrong input in scaling model'
            stop
        end if

    end subroutine read_scaling_model

    subroutine read_dip()
        write(fp_log,*) 'dip angle found in input' // new_line('a')
        read(fp_inp,*) flt_dip_deg
        flt_dip_rad = flt_dip_deg * DEG2RAD
    end subroutine read_dip

    subroutine read_gmpe_name()
        write(fp_log,*) 'gmpe_name found in input' // new_line('a')
        read(fp_inp,*) str_tmp
        if (str_tmp .eq. 'SADIGH97') then
            m_gmpe_name = SADIGH97
        else if (str_tmp .eq. 'CY14') then
            m_gmpe_name = CY14
        else
            write(fp_log,*) 'wrong gmpe_name in input'
            stop
        end if

    end subroutine read_gmpe_name

    subroutine read_vs30()
        write(fp_log,*) 'VS30 found in input' // new_line('a')
        read(fp_inp,*) Vs30

    end subroutine read_vs30

    subroutine read_z10()
        write(fp_log,*) 'Z10 found in input' // new_line('a')
        read(fp_inp,*) z10

    end subroutine read_z10

    subroutine read_z25()
        write(fp_log,*) 'Z25 found in input' // new_line('a')
        read(fp_inp,*) z25

    end subroutine read_z25

    subroutine read_seismogenic_depth()
        write(fp_log,*) 'seismogenic_depth found in input' // new_line('a')
        read(fp_inp,*) Smin, Smax

    end subroutine read_seismogenic_depth

    subroutine read_depth_distribution()
        write(fp_log,*) 'depth ditribution found in input' // new_line('a')
        read(fp_inp,*) str_tmp

        if (str_tmp .eq. 'UNIFORM') then
            m_DEPTH_distribution = UNIFORM
        else if (str_tmp .eq. 'TRIANGULAR') then
            m_DEPTH_distribution = TRIANGULAR
        else
            write(fp_log,*) 'wrong input in depth distribution'
            stop

        end if

    end subroutine read_depth_distribution

    subroutine read_aspect_ratio()
        write(fp_log,*) 'aspect ratio found in input' // new_line('a')
        read(fp_inp,*) aspect_ratio

    end subroutine read_aspect_ratio

    subroutine read_strike_dip_step ()
        write(fp_log,*) 'strike_dip_step found in input' // new_line('a')
        read(fp_inp,*) strike_step, dip_step
    end subroutine read_strike_dip_step

    subroutine read_mag_range ()
        write(fp_log,*) 'magnitude range found in input' // new_line('a')
        read(fp_inp,*) Mmin, Mmax
    end subroutine read_mag_range

    subroutine read_depth_param()
        write(fp_log,*) 'depth_param found in input' // new_line('a')
        read(fp_inp,*) depth_param
    end subroutine read_depth_param

    subroutine read_mag_step()
        write(fp_log,*) 'mag_step found in input' // new_line('a')
        read(fp_inp,*) mag_step
    end subroutine read_mag_step

    subroutine read_intensity()
        write(fp_log,*) "INTENSITY found"//new_line('A')
        numvalues = 0
        !line = 'null'
        InnerLoop: do
            read(fp_inp,*,iostat=eastat) tmp1
            if (eastat .ne. 0) exit InnerLoop

            numvalues = numvalues+1
            temp1(numvalues) = tmp1

        end do InnerLoop
        write(fp_log,*) "number of intensity points = ", numvalues
        write(fp_log,*) new_line('a')
        allocate(intensity(numvalues))
        intensity (1:numvalues) = temp1(1:numvalues)

    end subroutine read_intensity

    subroutine read_mag_bin()
        write(fp_log,*) "mag_bin found in input"//new_line('A')
        numvalues = 0
        !line = 'null'
        InnerLoop: do
            read(fp_inp,*,iostat=eastat) tmp1
            if (eastat .ne. 0) exit InnerLoop

            numvalues = numvalues+1
            temp1(numvalues) = tmp1

        end do InnerLoop
        write(fp_log,*) "number of mag bins = ", numvalues
        write(fp_log,*) new_line('a')
        allocate(mag_bin(numvalues))
        mag_bin (1:numvalues) = temp1(1:numvalues)

    end subroutine read_mag_bin


    subroutine read_dist_bin()
        write(fp_log,*) "dist_bin found in input"//new_line('A')
        numvalues = 0
        !line = 'null'
        InnerLoop: do
            read(fp_inp,*,iostat=eastat) tmp1
            if (eastat .ne. 0) exit InnerLoop

            numvalues = numvalues+1
            temp1(numvalues) = tmp1

        end do InnerLoop
        write(fp_log,*) "number of dist bins = ", numvalues
        write(fp_log,*) new_line('a')
        allocate(dist_bin(numvalues))
        dist_bin (1:numvalues) = temp1(1:numvalues)

    end subroutine read_dist_bin

    subroutine read_eps_bin()
        write(fp_log,*) "eps_bin found in input"//new_line('A')
        numvalues = 0
        !line = 'null'
        InnerLoop: do
            read(fp_inp,*,iostat=eastat) tmp1
            if (eastat .ne. 0) exit InnerLoop

            numvalues = numvalues+1
            temp1(numvalues) = tmp1

        end do InnerLoop
        write(fp_log,*) "number of eps bins = ", numvalues
        write(fp_log,*) new_line('a')
        allocate(eps_bin(numvalues))
        eps_bin (1:numvalues) = temp1(1:numvalues)

    end subroutine read_eps_bin

    subroutine read_gmpe_params()
        write(fp_log,*) "gmpe_params found in input"//new_line('A')
        numvalues = 0
        !line = 'null'
        InnerLoop: do
            read(fp_inp,*,iostat=eastat) tmp1
            if (eastat .ne. 0) exit InnerLoop

            numvalues = numvalues+1
            temp1(numvalues) = tmp1

        end do InnerLoop
        write(fp_log,*) "number of gmpe parameters = ", numvalues
        write(fp_log,*) new_line('a')
        allocate(gmpe_params(numvalues))
        gmpe_params (1:numvalues) = temp1(1:numvalues)

    end subroutine read_gmpe_params

    subroutine read_gmpe_opts()
        write(fp_log,*) "gmpe_opts found in input"//new_line('A')
        numvalues = 0
        !line = 'null'
        InnerLoop: do
            read(fp_inp,*,iostat=eastat) tmp_int
            if (eastat .ne. 0) exit InnerLoop

            numvalues = numvalues+1
            temp_int(numvalues) = tmp_int

        end do InnerLoop
        write(fp_log,*) "number of gmpe options = ", numvalues
        write(fp_log,*) new_line('a')
        allocate(gmpe_opts(numvalues))
        gmpe_opts(1:numvalues) = temp_int(1:numvalues)

    end subroutine read_gmpe_opts


    subroutine close_file()
        close(fp_haz)
        close(fp_dag)
        close(fp_log)
        close(fp_rup)
    end subroutine

    subroutine print_haz_bin(haz_bin)
        real(8) :: haz_bin(:,:,:,:,:)
        integer :: i_eps, i_dist, i_mag, i_inten, i_freq

        write(fp_dag, *) '    hazard deaggregation '//new_line('a')
        write(fp_dag, *) '  i_eps i_dist  i_mag i_inten  i_freq     haz'

        do i_freq = 1, n_freq
            do i_inten = 1, n_inten
                do i_mag = 1, n_mag_bin
                    do i_dist = 1, n_dist_bin
                        do i_eps = 1, n_eps_bin
                            write(fp_dag, '(5i7,e15.6)') i_eps, i_dist, i_mag, i_inten, i_freq, &
                                Haz_bin(i_eps, i_dist, i_mag, i_inten, i_freq)

                        end do
                    end do
                end do
            end do
        end do
    end subroutine print_haz_bin

    subroutine print_haz(haz)
        real(8) :: haz(:,:)
        integer :: i_freq, i_inten

        write(fp_haz, *) '   Hazard results'

        do i_inten = 1, n_inten
            write(fp_haz,'(1f10.3)', advance = 'NO') intensity(i_inten)
            do i_freq = 1, n_freq
                write(fp_haz,'(1e15.6)', advance = 'NO') Haz(i_inten, i_freq)

            end do
            write(fp_haz,'(A)', advance = 'no') ' '//new_line('a')
        end do

    end subroutine

end module input_module
