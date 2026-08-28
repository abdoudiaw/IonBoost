module mod_config
    use mod_types, only: dp
    implicit none
    private
    public :: run_config, read_config

    type :: run_config
        ! mesh
        integer  :: ncell        ! ion sheets (sheet index 0:ncell, 0 = rear wall)
        integer  :: nvacuum      ! vacuum-sheath nodes beyond the ion front
        integer  :: ntotal       ! ncell + nvacuum
        real(dp) :: prog         ! mesh progression (<1: cells shrink towards the front)
        logical  :: lfini        ! finite foil (T) or semi-infinite estimate (F)
        real(dp) :: lmax         ! foil thickness when lfini
        real(dp) :: length       ! actual plasma length used to build the mesh
        ! time stepping
        integer  :: itmax
        real(dp) :: tmax, dti
        logical  :: vtt = .true. ! adaptive step: dt = dti / sqrt(ni(0)/ni0)
        ! Poisson-Boltzmann Newton iterations (iter0 for the first 3 steps)
        integer  :: iter0, iter1
        ! plasma
        real(dp) :: n0hot        ! hot-electron reference density
        real(dp) :: Thmax        ! hot-electron temperature
        real(dp) :: n0cold = 0.0_dp  ! cold-electron reference density (0 = single species)
        real(dp) :: Tcmax = 0.0_dp   ! cold-electron temperature
        real(dp) :: ni0          ! unperturbed ion density (= total electron density)
        real(dp) :: nu = 0.0_dp  ! ion viscosity (0 = off)
        ! legacy input fields, read for file-format compatibility but unused:
        ! iter2/iter3 belonged to the energy-conserving Th iteration, T_MeV/LSS/
        ! nLSS to the relativistic/graded-profile variants of the original code
        integer  :: iter2, iter3
        logical  :: nb_cons, En_cons
        real(dp) :: T_MeV, LSS, nLSS
        ! output cadence
        integer  :: nsort = 1    ! history stride
        integer  :: iline = 1    ! profile/spectrum stride
    end type run_config

contains

    subroutine read_config(cfg, filename)
        type(run_config), intent(out) :: cfg
        character(*), intent(in) :: filename
        integer :: iunit, ierr, nstep

        open(newunit=iunit, file=filename, status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
            write(*, *) '*** Error reading the input file ***'
            stop 1
        end if
        read(iunit, *) cfg%ncell, cfg%nvacuum, cfg%prog
        read(iunit, *) cfg%itmax, cfg%tmax
        read(iunit, *) cfg%iter0, cfg%iter1, cfg%iter2, cfg%iter3
        read(iunit, *) cfg%dti
        read(iunit, *) cfg%n0hot
        read(iunit, *) cfg%Thmax
        read(iunit, *) cfg%lfini, cfg%lmax
        read(iunit, *) cfg%nb_cons, cfg%En_cons
        read(iunit, *) cfg%T_MeV, cfg%LSS, cfg%nLSS
        ! optional 10th line: cold electron population "n0cold Tcmax"
        ! (bi-Maxwellian expansion; absent = single species)
        cfg%n0cold = 0.0_dp
        cfg%Tcmax = 0.0_dp
        read(iunit, *, iostat=ierr) cfg%n0cold, cfg%Tcmax
        if (ierr /= 0) then
            cfg%n0cold = 0.0_dp
            cfg%Tcmax = 0.0_dp
        end if
        close(iunit)

        if (cfg%n0cold > 0.0_dp .and. cfg%Tcmax <= 0.0_dp) then
            write(*, *) 'n0cold > 0 requires Tcmax > 0.'
            stop 1
        end if
        cfg%ni0 = cfg%n0hot + cfg%n0cold

        ! plasma-oscillation stability of the leapfrog: omega_p * dt <= 2
        if (sqrt(cfg%ni0) * cfg%dti > 2.0_dp) then
            cfg%dti = 2.0_dp / sqrt(cfg%ni0)
            write(*, *) 'dti reduced for stability: dti =', cfg%dti
        end if

        if (cfg%Thmax <= 0.0_dp) then
            write(*, *) 'Hot temperature must be positive.'
            stop 1
        end if
        if (cfg%nb_cons .or. cfg%En_cons) then
            cfg%nb_cons = .false.
            cfg%En_cons = .false.
            write(*, *) 'isothermal single-species mode: nb_cons = En_cons = .false.'
        end if

        cfg%ntotal = cfg%ncell + cfg%nvacuum

        if (cfg%lfini) then
            cfg%length = cfg%lmax
        else
            ! semi-infinite run: plasma deep enough that the rarefaction never
            ! reaches the rear boundary before tmax
            cfg%length = (25.0_dp + cfg%tmax) * sqrt(cfg%Thmax)
        end if

        nstep = int(cfg%tmax / cfg%dti + 0.5_dp)
        cfg%nsort = 1
        if (nstep > 2000) cfg%nsort = nstep / 2000
        cfg%iline = 1
        if (cfg%ncell > 3999) cfg%iline = cfg%ncell / 2000
    end subroutine read_config

end module mod_config
