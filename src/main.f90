program IonBoost
    ! 1D Lagrangian solver for collisionless plasma expansion into a vacuum
    ! with Boltzmann electrons, after P. Mora, PRL 90, 185002 (2003).
    !
    ! The driver owns the time loop only; the physics lives in the modules:
    !   mod_config        input parameters
    !   mod_grid          Lagrangian node grid (1D seam for a future multi-D grid)
    !   mod_ions          ion species (instantiate more for multi-species runs)
    !   mod_electrons     Boltzmann electron populations (append for two-temperature)
    !   mod_fields        Poisson-Boltzmann Newton solve, sheath, field, energies
    !   mod_linear_solver linear-solve seam (swap Thomas for e.g. an AmgX backend)
    !   mod_pusher        leapfrog sheet dynamics and sheet ordering
    !   mod_diagnostics   output files
    use mod_types, only: dp
    use mod_config, only: run_config, read_config
    use mod_grid, only: grid_t, build_grid, update_ion_spacing
    use mod_ions, only: ion_species_t, init_ions, rebuild_density, ion_kinetic_energy
    use mod_electrons, only: electron_species_t
    use mod_fields, only: fields_t, init_fields, save_previous, solve_potential, &
                          extend_sheath, compute_field, electron_work, &
                          electrostatic_energy, total_electrons
    use mod_linear_solver, only: thomas_solver_t
    use mod_pusher, only: first_half_kick, advance_velocity, sync_full_velocity, &
                          push_positions, sort_sheets
    use mod_diagnostics, only: diagnostics_t, open_histories, write_history, &
                               close_histories, write_profiles, write_spectra
    implicit none

    type(run_config) :: cfg
    type(grid_t) :: grid
    type(ion_species_t) :: ions
    type(electron_species_t), allocatable :: els(:)
    type(fields_t) :: flds
    type(thomas_solver_t) :: solver
    type(diagnostics_t) :: diag

    integer :: i, itime, niter, nsort, nstep
    real(dp) :: time, dt, dtold
    real(dp) :: nti, nthot, En_ion, En_elec
    logical :: laststep

    call read_config(cfg, 'input.in')
    call build_grid(grid, cfg)
    call init_ions(ions, grid, ni0=cfg%n0hot, Zion=1.0_dp, prog=cfg%prog)
    allocate(els(1))
    els(1) = electron_species_t(n0=cfg%n0hot, T=cfg%Thmax)
    call init_fields(flds, grid, els)
    call open_histories(diag)

    time = 0.0_dp
    dt = cfg%dti
    dtold = dt
    laststep = .false.
    nsort = cfg%nsort

    do itime = 1, cfg%itmax
        if (itime <= 3) then
            niter = cfg%iter0
        else
            niter = cfg%iter1
        end if

        if (mod(itime, max(nsort, 1)) == 0) print *, time

        if (itime > 1) call save_previous(flds, grid)

        call update_ion_spacing(grid)
        call rebuild_density(ions, grid)
        call solve_potential(flds, grid, ions, els, solver, niter)
        call extend_sheath(flds, grid, els)

        nti = 0.0_dp
        do i = 1, ions%n
            nti = nti + ions%dens_half(i) * grid%dxt(i)
        end do
        nthot = total_electrons(flds, grid)

        ! node velocities for the electron-work advection term: ion sheets carry
        ! their half-step velocity, sheath nodes move with the front
        flds%unode(0:ions%n) = ions%vint(0:ions%n)
        if (itime > 1) then
            do i = ions%n + 1, grid%ntotal
                flds%unode(i) = (grid%xt(i) - grid%xtold(i)) / max(dt, 1.0e-12_dp)
            end do
            call electron_work(flds, grid, dt)
        else
            call electron_work(flds, grid, 0.0_dp)   ! gradients only, no work yet
        end if

        call compute_field(flds, grid)
        if (itime > 1) call sync_full_velocity(ions, flds, dt)

        En_ion = ion_kinetic_energy(ions)
        En_elec = electrostatic_energy(flds, grid, els(1)%T)

        if (cfg%itmax == 1) exit

        nstep = int(cfg%tmax / max(dt, 1.0e-12_dp))
        nsort = 1
        if (nstep > 2000) nsort = max(nstep / 2000, 1)
        if (mod(itime - 1, nsort) == 0 .or. laststep) then
            call write_history(diag, grid, ions, flds, time, cfg%n0hot, els(1)%T, &
                               nti, nthot, En_ion, En_elec, first_step=(itime == 1))
        end if

        if (time >= cfg%tmax - 1.0e-5_dp) exit
        if (itime == cfg%itmax) exit

        dtold = dt
        dt = cfg%dti
        if (cfg%vtt) dt = cfg%dti / sqrt(max(ions%dens(0) / cfg%n0hot, 1.0e-12_dp))
        if ((time + dt) > (cfg%tmax - 1.0e-5_dp)) then
            dt = cfg%tmax - time
            laststep = .true.
        end if
        time = time + dt

        if (itime == 1) then
            call first_half_kick(ions, flds, dt)
        else
            call advance_velocity(ions, grid, flds, solver, cfg%nu, dtold, dt)
        end if
        call push_positions(ions, grid, dt)
        call sort_sheets(ions, grid, flds)
    end do

    write(*, '("Final time = ", f8.4)') time
    write(*, '("Number of iterations = ", i0)') itime

    if (itime == 1) then
        do i = ions%n + 1, grid%ntotal
            grid%x0(i) = grid%xt(i)
        end do
    end if

    call write_profiles(grid, ions, flds, time, cfg%n0hot, els(1)%T, cfg%iline)
    call write_spectra(grid, ions, cfg%iline)
    call close_histories(diag)
    write(*, *) 'time', time

end program IonBoost
