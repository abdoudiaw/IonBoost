module mod_fields
    use mod_types, only: dp
    use mod_grid, only: grid_t
    use mod_ions, only: ion_species_t
    use mod_electrons, only: electron_species_t, e_density, e_ddensity_dphi, e_pressure
    use mod_linear_solver, only: linear_solver_t
    implicit none
    private
    public :: fields_t, init_fields, save_previous, solve_potential, &
              extend_sheath, compute_field, electron_work, &
              electrostatic_energy, total_electrons

    ! Electrostatic field state on the grid nodes (0:ntotal).
    ! Sign convention: phi is the normalized -e*Phi/Te, so E = +dphi/dx and
    ! Gauss's law reads dE/dx = ni - ne.
    type :: fields_t
        real(dp), allocatable :: phi(:), phiold(:)
        real(dp), allocatable :: E(:), Eold(:)
        real(dp), allocatable :: E_half(:)          ! 1:ntotal, at node midpoints
        real(dp), allocatable :: ne(:), rho(:)
        real(dp), allocatable :: nhot(:), nhotold(:)
        real(dp), allocatable :: gradnh(:), ghold(:)
        real(dp), allocatable :: phot(:), pe(:)
        real(dp), allocatable :: unode(:)           ! node velocity (advection in Whot)
        real(dp), allocatable :: dW(:), SW(:)
        real(dp) :: Whot1 = 0.0_dp, Whot2 = 0.0_dp  ! electron work: plasma / sheath
        real(dp) :: Whot1old = 0.0_dp, Whot2old = 0.0_dp
        real(dp) :: Whot = 0.0_dp
        ! tridiagonal assembly workspace
        real(dp), allocatable :: a(:), b(:), c(:), f(:)
    end type fields_t

contains

    subroutine init_fields(flds, grid, els)
        type(fields_t), intent(out) :: flds
        type(grid_t), intent(in) :: grid
        type(electron_species_t), intent(in) :: els(:)
        integer :: i, n, nt
        real(dp) :: alpha

        n = grid%ncell
        nt = grid%ntotal
        allocate(flds%phi(0:nt), flds%phiold(0:nt), flds%E(0:nt), flds%Eold(0:nt))
        allocate(flds%E_half(1:nt))
        allocate(flds%ne(0:nt), flds%rho(0:nt), flds%nhot(0:nt), flds%nhotold(0:nt))
        allocate(flds%gradnh(0:nt), flds%ghold(0:nt), flds%phot(0:nt), flds%pe(0:nt))
        allocate(flds%unode(0:nt), flds%dW(0:nt), flds%SW(0:nt))
        allocate(flds%a(1:nt), flds%b(0:nt), flds%c(0:nt - 1), flds%f(0:nt))

        flds%E = 0.0_dp
        flds%Eold = 0.0_dp
        flds%phiold = 0.0_dp
        flds%rho = 0.0_dp
        flds%gradnh = 0.0_dp
        flds%ghold = 0.0_dp
        flds%phot = 0.0_dp
        flds%pe = 0.0_dp
        flds%unode = 0.0_dp
        flds%dW = 0.0_dp
        flds%SW = 0.0_dp
        flds%E_half = 0.0_dp

        ! initial guess: exponential potential matching the sharp-edge field
        ! E(edge) = sqrt(2 pe(edge)) with phi(edge) = <T> (Eq. 8 of the paper)
        flds%phi(n) = e_pressure(els, 0.0_dp) / e_density(els, 0.0_dp)
        flds%E(n) = sqrt(2.0_dp * e_pressure(els, flds%phi(n)))
        alpha = flds%E(n) / flds%phi(n)
        do i = 0, n
            flds%phi(i) = flds%phi(n) * exp(alpha * (grid%x0(i) - grid%x0(n)))
        end do
        do i = n + 1, nt
            flds%phi(i) = flds%phi(n)
        end do
        flds%nhot = e_density(els, 0.0_dp)
        flds%ne = flds%nhot
        do i = n + 1, nt
            flds%nhot(i) = 0.0_dp
            flds%ne(i) = 0.0_dp
        end do
        flds%E(0) = 0.0_dp
    end subroutine init_fields

    subroutine save_previous(flds, grid)
        type(fields_t), intent(inout) :: flds
        type(grid_t), intent(inout) :: grid
        integer :: i
        do i = 0, grid%ntotal
            flds%nhotold(i) = flds%nhot(i)
            flds%phiold(i) = flds%phi(i)
            flds%Eold(i) = flds%E(i)
            flds%ghold(i) = flds%gradnh(i)
        end do
        do i = grid%ncell + 1, grid%ntotal
            grid%xtold(i) = grid%xt(i)
        end do
        do i = 1, grid%ntotal
            grid%dxtold(i) = grid%dxt(i)
        end do
        flds%Whot1old = flds%Whot1
        flds%Whot2old = flds%Whot2
    end subroutine save_previous

    ! Newton iteration on the nonlinear Poisson-Boltzmann equation
    !   d2phi/dx2 = ni - ne(phi)
    ! with dphi/dx = 0 at the wall and the linearized sharp-front condition
    ! E = sqrt(2 pe(phi)) at node ncell (Eq. 8 of the paper). Each iteration
    ! assembles the Jacobian system and hands it to the linear solver.
    subroutine solve_potential(flds, grid, ions, els, solver, niter)
        type(fields_t), intent(inout) :: flds
        type(grid_t), intent(in) :: grid
        type(ion_species_t), intent(in) :: ions
        type(electron_species_t), intent(in) :: els(:)
        class(linear_solver_t), intent(inout) :: solver
        integer, intent(in) :: niter
        integer :: i, it, n
        real(dp) :: neh, dne, peold, neold

        n = grid%ncell
        do it = 1, niter
            neh = e_density(els, flds%phi(0))
            dne = e_ddensity_dphi(els, flds%phi(0))
            flds%b(0) = -2.0_dp / max(grid%dxt(1) * grid%dxt(1), 1.0e-12_dp) + dne
            flds%c(0) = 2.0_dp / max(grid%dxt(1) * grid%dxt(1), 1.0e-12_dp)
            flds%f(0) = ions%dens(0) - neh + flds%phi(0) * dne

            do i = 1, n - 1
                neh = e_density(els, flds%phi(i))
                dne = e_ddensity_dphi(els, flds%phi(i))
                flds%a(i) = 2.0_dp / max(grid%dxt(i) * (grid%dxt(i) + grid%dxt(i + 1)), 1.0e-12_dp)
                flds%b(i) = -2.0_dp / max(grid%dxt(i) * grid%dxt(i + 1), 1.0e-12_dp) + dne
                flds%c(i) = 2.0_dp / max(grid%dxt(i + 1) * (grid%dxt(i) + grid%dxt(i + 1)), 1.0e-12_dp)
                flds%f(i) = ions%dens(i) - neh + flds%phi(i) * dne
            end do

            neold = e_density(els, flds%phi(n))
            dne = e_ddensity_dphi(els, flds%phi(n))
            peold = e_pressure(els, flds%phi(n))
            flds%a(n) = 2.0_dp / max(grid%dxt(n) * grid%dxt(n), 1.0e-12_dp)
            flds%b(n) = -flds%a(n) - sqrt(2.0_dp / max(peold, 1.0e-12_dp)) * neold / &
                        max(grid%dxt(n), 1.0e-12_dp) + dne
            flds%f(n) = ions%dens(n) - neold + flds%phi(n) * dne - &
                        sqrt(8.0_dp * max(peold, 1.0e-12_dp)) / max(grid%dxt(n), 1.0e-12_dp) * &
                        (1.0_dp + 0.5_dp * neold * flds%phi(n) / max(peold, 1.0e-12_dp))

            call solver%solve(n, flds%a(1:n), flds%b(0:n), flds%c(0:n - 1), &
                              flds%f(0:n), flds%phi(0:n))
        end do

        ! electron density and charge separation in the plasma
        do i = 0, n
            flds%nhot(i) = e_density(els, flds%phi(i))
            flds%ne(i) = flds%nhot(i)
            flds%rho(i) = ions%dens(i) - flds%ne(i)
        end do
        flds%phot(n) = flds%nhot(n) * els(1)%T
        flds%pe(n) = e_pressure(els, flds%phi(n))
    end subroutine solve_potential

    ! Continue phi beyond the ion front with the exact Boltzmann vacuum-sheath
    ! solution phi = phi_f + 2 T ln(1 + k x), k = sqrt(n0 e^{-phi_f/T} / 2T).
    ! The closed form is single-population; with several populations the
    ! hottest one is used, which is exact in the limit where the colder
    ! densities are negligible beyond the front (phi_front >> T_cold).
    ! The grid spacing follows the local Debye length.
    subroutine extend_sheath(flds, grid, els)
        type(fields_t), intent(inout) :: flds
        type(grid_t), intent(inout) :: grid
        type(electron_species_t), intent(in) :: els(:)
        integer :: i, n, ih
        real(dp) :: kvide, uphi, Th

        n = grid%ncell
        ih = maxloc(els%T, dim=1)
        Th = els(ih)%T
        kvide = sqrt(els(ih)%n0 / 2.0_dp / Th)
        uphi = 0.0_dp
        do i = n + 1, grid%ntotal
            if (i == n + 1) then
                grid%dxt(i) = 0.0_dp
            else
                grid%dxt(i) = 0.025_dp * sqrt(Th / max(flds%nhot(i - 1), 1.0e-12_dp))
            end if
            grid%xt(i) = grid%xt(i - 1) + grid%dxt(i)
            uphi = uphi + kvide * grid%dxt(i)
            flds%phi(i) = flds%phi(n) + 2.0_dp * Th * &
                          log(1.0_dp + uphi * exp(-0.5_dp * flds%phi(n) / Th))
            flds%nhot(i) = e_density(els, flds%phi(i))
            flds%ne(i) = flds%nhot(i)
            flds%rho(i) = -flds%ne(i)
            flds%phot(i) = flds%nhot(i) * Th
            flds%pe(i) = flds%phot(i)
            flds%E(i) = sqrt(2.0_dp * flds%pe(i))
        end do
    end subroutine extend_sheath

    ! E at the nodes from centered differences of phi; at the front node the
    ! sharp-front relation E = sqrt(2 pe) is exact and replaces the stencil
    subroutine compute_field(flds, grid)
        type(fields_t), intent(inout) :: flds
        type(grid_t), intent(in) :: grid
        integer :: i, n
        n = grid%ncell
        do i = 1, n
            flds%E_half(i) = (flds%phi(i) - flds%phi(i - 1)) / max(grid%dxt(i), 1.0e-12_dp)
        end do
        do i = 1, n - 1
            flds%E(i) = (grid%dxt(i + 1) * flds%E_half(i) + grid%dxt(i) * flds%E_half(i + 1)) / &
                        max(grid%dxt(i) + grid%dxt(i + 1), 1.0e-12_dp)
        end do
        flds%E(n) = sqrt(2.0_dp * flds%pe(n))
    end subroutine compute_field

    ! Energy handed by the (isothermal) electrons to the field and ions:
    ! dW = <phi> d(ne)|_comoving integrated over the grid, accumulated in
    ! Whot1 (plasma) and Whot2 (vacuum sheath). Requires unode and the
    ! electron-density gradients gradnh.
    subroutine electron_work(flds, grid, dt)
        type(fields_t), intent(inout) :: flds
        type(grid_t), intent(in) :: grid
        real(dp), intent(in) :: dt
        integer :: i, n, nt

        n = grid%ncell
        nt = grid%ntotal

        ! centered electron-density gradients, extrapolated at the edges
        do i = 1, nt - 1
            flds%gradnh(i) = (flds%nhot(i + 1) - flds%nhot(i - 1)) / &
                             max(grid%dxt(i) + grid%dxt(i + 1), 1.0e-12_dp)
        end do
        flds%gradnh(0) = 0.0_dp
        flds%gradnh(n) = flds%gradnh(n - 1) * flds%nhot(n) / max(flds%nhot(n - 1), 1.0e-12_dp)
        flds%gradnh(nt) = flds%gradnh(nt - 1) * flds%nhot(nt) / max(flds%nhot(nt - 1), 1.0e-12_dp)
        flds%gradnh(n + 1) = flds%gradnh(n)

        if (dt <= 0.0_dp) return

        flds%dW(0) = 0.5_dp * (flds%phiold(0) + flds%phi(0)) * &
                     (flds%nhot(0) - flds%nhotold(0) - flds%unode(0) * dt * &
                      0.5_dp * (flds%gradnh(0) + flds%ghold(0)))
        flds%Whot1 = flds%Whot1old + 0.25_dp * flds%dW(0) * (grid%dxtold(1) + grid%dxt(1))
        flds%SW(0) = 0.25_dp * flds%dW(0) * (grid%dxtold(1) + grid%dxt(1))
        do i = 1, n
            flds%dW(i) = 0.5_dp * (flds%phiold(i) + flds%phi(i)) * &
                         (flds%nhot(i) - flds%nhotold(i) - flds%unode(i) * dt * &
                          0.5_dp * (flds%gradnh(i) + flds%ghold(i)))
            flds%Whot1 = flds%Whot1 + 0.25_dp * flds%dW(i) * &
                         (grid%dxtold(i) + grid%dxt(i) + grid%dxtold(i + 1) + grid%dxt(i + 1))
            flds%SW(i) = flds%SW(i - 1) + 0.25_dp * flds%dW(i) * &
                         (grid%dxtold(i) + grid%dxt(i) + grid%dxtold(i + 1) + grid%dxt(i + 1))
        end do
        flds%Whot2 = flds%Whot2old
        do i = n + 1, nt - 1
            flds%dW(i) = 0.5_dp * (flds%phiold(i) + flds%phi(i)) * &
                         (flds%nhot(i) - flds%nhotold(i) - flds%unode(i) * dt * &
                          0.5_dp * (flds%gradnh(i) + flds%ghold(i)))
            flds%Whot2 = flds%Whot2 + 0.25_dp * flds%dW(i) * &
                         (grid%dxtold(i) + grid%dxt(i) + grid%dxtold(i + 1) + grid%dxt(i + 1))
            flds%SW(i) = flds%SW(i - 1) + 0.25_dp * flds%dW(i) * &
                         (grid%dxtold(i) + grid%dxt(i) + grid%dxtold(i + 1) + grid%dxt(i + 1))
        end do
        flds%dW(nt) = 0.5_dp * (flds%phiold(nt) + flds%phi(nt)) * &
                      (flds%nhot(nt) - flds%nhotold(nt) - flds%unode(nt) * dt * &
                       0.5_dp * (flds%gradnh(nt) + flds%ghold(nt)))
        flds%Whot2 = flds%Whot2 + 0.25_dp * flds%dW(nt) * (grid%dxtold(nt) + grid%dxt(nt))
        flds%SW(nt) = flds%SW(nt - 1) + 0.25_dp * flds%dW(nt) * (grid%dxtold(nt) + grid%dxt(nt))
        flds%Whot = flds%Whot1 + flds%Whot2
    end subroutine electron_work

    ! Field energy over the grid plus the analytic tail beyond the last node:
    ! for the Boltzmann sheath, int E^2/2 dx = T * E(ntotal)
    function electrostatic_energy(flds, grid, Th) result(En)
        type(fields_t), intent(in) :: flds
        type(grid_t), intent(in) :: grid
        real(dp), intent(in) :: Th
        real(dp) :: En
        integer :: i
        En = 0.0_dp
        do i = 1, grid%ntotal - 1
            En = En + 0.25_dp * flds%E(i)**2 * (grid%dxt(i) + grid%dxt(i + 1))
        end do
        En = En + 0.25_dp * flds%E(grid%ntotal)**2 * grid%dxt(grid%ntotal)
        En = En + flds%E(grid%ntotal) * Th
    end function electrostatic_energy

    ! Total electron number, with the analytic tail int ne dx = E(ntotal)
    function total_electrons(flds, grid) result(nte)
        type(fields_t), intent(in) :: flds
        type(grid_t), intent(in) :: grid
        real(dp) :: nte
        integer :: i
        nte = 0.0_dp
        do i = 1, grid%ntotal
            nte = nte + 0.5_dp * (flds%nhot(i - 1) + flds%nhot(i)) * grid%dxt(i)
        end do
        nte = nte + flds%E(grid%ntotal)
    end function total_electrons

end module mod_fields
