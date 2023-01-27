module mod_io
    use mod_types
    use mod_share
    implicit none
    private
    public :: read_input
    ! parameters
    !
    real(dp), parameter :: pi = acos(-1._dp)
    contains
        subroutine read_input()
            ! Read the initial conditions
            integer :: iunit, ierr
            open(newunit = iunit, file = 'input.in', status = 'old', action = 'read', iostat = ierr)
            if (ierr == 0) then
                read(iunit, *) ncell, nvacuum, prog
                read(iunit, *) itmax, tmax
                read(iunit, *) iter0, iter1, iter2, iter3
                read(iunit, *) dti
                read(iunit, *) n0hot
                read(iunit, *) Thmax
                read(iunit, *) lfini, lmax
                read(iunit, *) nb_cons, En_cons
                read(iunit, *) T_MeV, LSS, nLSS
            else
                write(*, *) '*** Error reading the input file *** '
                write(*, ) 'Aborting...'
                stop
            endif
            close(iunit)
            !
            if ((.not.lfini) .and. (nb_cons .or. En_cons)) then
            nb_cons = .false.
            En_cons = .false.
            write(, ) 'lfini=.false. therefore En_cons=.false. and nb_cons=.false.'
            endif
            !
            if ((.not.nb_cons) .and. En_cons) then
            nb_cons = .true.
            write(, ) 'En_cons=.true. therefore nb_cons=.true.'
            endif
            !
            if (nb_cons .and. (n0hot < 1.e-04)) then
            nb_cons = .false.
            write(, ) 'n0hot < 1.e-04 therefore nb_cons=.false.'
            endif
            !
            ntotal = ncell + nvacuum
            if (ntotal > nmax) then
            write(, *) 'ntotal=ncell+nvacuum is bigger than nmax'
            stop
            endif
            !
            nstep = (tmax / dti + 0.5_dp)
            nsort = 1
            if (nstep > 2000) nsort = nstep / 2000
            laststep = .false.
            iline = 1
            if (ncell > 3999) iline = ncell / 2000
            ! Determine the length of the plasma
            if (.not.lfini) then
            time = 0._dp
            dt = dti
            cs2 = (n0hot) / (n0hot / Thmax)
            cs = sqrt(cs2)
            length = 25._dp * cs
            if ((time + dt) > (tmax - 1.d-5)) dt = tmax - time
            cs2old = cs2
            time = time + dt
            cs2 = (n0hot) / (n0hot / Thmax)
            cs = 0.5_dp * (sqrt(cs2old) + sqrt(cs2))
            length = length + cs * dt
            endif
            !
            if (lfini) length = lmax
            ! mesh
            rgauss = 2._dp * length / sqrt(pi)
            ni0 = n0hot
            if (prog == 1._dp) then
                dx0(1) = length / ncell
            else
                dx0(1) = length * (1. - prog) / (1. - prog ** ncell)
            endif
            !
            x0(0) = -length
            x0(1) = x0(0) + dx0(1)
            !
            do i = 0, 1
                niSS(i) = ni0
            end do
            do i = 2, ncell
                if (i == 2) then
                    dx0(2) = (1._dp - prog) * dx0(1)
                else
                    dx0(i) = dx0(i - 1) * prog
                endif
            !
                x0(i) = x0(i - 1) + dx0(i)
            end do
    end subroutine read_input
end module mod_io



