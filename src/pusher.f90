module mod_pusher
    use mod_types, only: dp
    use mod_grid, only: grid_t
    use mod_ions, only: ion_species_t
    use mod_fields, only: fields_t
    use mod_linear_solver, only: linear_solver_t
    implicit none
    private
    public :: first_half_kick, advance_velocity, sync_full_velocity, push_positions, sort_sheets

contains

    ! Leapfrog start: put vint half a step ahead of v
    subroutine first_half_kick(ions, flds, dt)
        type(ion_species_t), intent(inout) :: ions
        type(fields_t), intent(in) :: flds
        real(dp), intent(in) :: dt
        integer :: i
        do i = 1, ions%n
            ions%vint(i) = ions%v(i) + 0.5_dp * dt * flds%E(i) * ions%Z(i)
        end do
    end subroutine first_half_kick

    ! Advance the half-step velocity across a (possibly changing) step:
    ! vint_new = vint_old + 0.5(dtold+dt) Z E, solved implicitly so an optional
    ! viscosity nu couples neighboring sheets (nu = 0 reduces to the explicit
    ! kick). Sheet 0 is the fixed rear wall.
    subroutine advance_velocity(ions, grid, flds, solver, nu, dtold, dt)
        type(ion_species_t), intent(inout) :: ions
        type(grid_t), intent(in) :: grid
        type(fields_t), intent(inout) :: flds
        class(linear_solver_t), intent(inout) :: solver
        real(dp), intent(in) :: nu, dtold, dt
        integer :: i, n
        real(dp) :: delt

        n = ions%n
        delt = 0.5_dp * (dtold + dt)
        flds%b(0) = 1.0_dp / max(delt, 1.0e-12_dp)
        flds%c(0) = 0.0_dp
        flds%f(0) = 0.0_dp
        do i = 1, n - 1
            flds%a(i) = -2.0_dp * nu / max(grid%dxt(i) * (grid%dxt(i) + grid%dxt(i + 1)), 1.0e-12_dp)
            flds%b(i) = 1.0_dp / max(delt, 1.0e-12_dp) + &
                        2.0_dp * nu / max(grid%dxt(i) * grid%dxt(i + 1), 1.0e-12_dp)
            flds%c(i) = -2.0_dp * nu / max(grid%dxt(i + 1) * (grid%dxt(i) + grid%dxt(i + 1)), 1.0e-12_dp)
            flds%f(i) = flds%E(i) * ions%Z(i) + ions%vint(i) / max(delt, 1.0e-12_dp)
        end do
        flds%a(n) = 0.0_dp
        flds%b(n) = 1.0_dp / max(delt, 1.0e-12_dp)
        flds%f(n) = flds%E(n) * ions%Z(n) + ions%vint(n) / max(delt, 1.0e-12_dp)
        call solver%solve(n, flds%a(1:n), flds%b(0:n), flds%c(0:n - 1), &
                          flds%f(0:n), ions%vint(0:n))
    end subroutine advance_velocity

    ! Whole-step velocity from the half-step one: v(t) = vint(t-dt/2)
    ! + (dt/2) a(t-dt/4), with the field linearly interpolated in time
    subroutine sync_full_velocity(ions, flds, dt)
        type(ion_species_t), intent(inout) :: ions
        type(fields_t), intent(in) :: flds
        real(dp), intent(in) :: dt
        integer :: i
        do i = 1, ions%n
            ions%v(i) = ions%vint(i) + dt * (3.0_dp * flds%E(i) + flds%Eold(i)) * ions%Z(i) / 8.0_dp
        end do
    end subroutine sync_full_velocity

    ! Move the sheets; reflect any that cross the rear wall
    subroutine push_positions(ions, grid, dt)
        type(ion_species_t), intent(inout) :: ions
        type(grid_t), intent(inout) :: grid
        real(dp), intent(in) :: dt
        integer :: i
        do i = 1, ions%n
            grid%xt(i) = grid%xt(i) + dt * ions%vint(i)
            if (grid%xt(i) < grid%x0(0)) then
                grid%xt(i) = 2.0_dp * grid%x0(0) - grid%xt(i)
                ions%vint(i) = -ions%vint(i)
            end if
        end do
    end subroutine push_positions

    ! Keep the sheets ordered in x (crossing sheets exchange labels).
    ! Stable insertion sort on position; every per-sheet quantity, including
    ! the field values carried as the next Newton initial guess, follows the
    ! same permutation.
    subroutine sort_sheets(ions, grid, flds)
        type(ion_species_t), intent(inout) :: ions
        type(grid_t), intent(inout) :: grid
        type(fields_t), intent(inout) :: flds
        integer :: i, j, n
        real(dp) :: xxt, xx0, xxtold, vv, vvint, qq, nn, pp, EE, gg, dd, ZZ
        integer :: ii

        n = ions%n
        do j = 2, n
            xxt = grid%xt(j)
            xx0 = grid%x0(j)
            xxtold = grid%xtold(j)
            vv = ions%v(j)
            vvint = ions%vint(j)
            ii = ions%id0(j)
            qq = ions%q(j)
            nn = flds%nhot(j)
            pp = flds%phi(j)
            EE = flds%E(j)
            gg = flds%gradnh(j)
            dd = grid%dxt(j)
            ZZ = ions%Z(j)
            do i = j - 1, 1, -1
                if (grid%xt(i) <= xxt) exit
                grid%xt(i + 1) = grid%xt(i)
                grid%x0(i + 1) = grid%x0(i)
                grid%xtold(i + 1) = grid%xtold(i)
                ions%v(i + 1) = ions%v(i)
                ions%vint(i + 1) = ions%vint(i)
                ions%id0(i + 1) = ions%id0(i)
                ions%q(i + 1) = ions%q(i)
                flds%nhot(i + 1) = flds%nhot(i)
                flds%phi(i + 1) = flds%phi(i)
                flds%E(i + 1) = flds%E(i)
                flds%gradnh(i + 1) = flds%gradnh(i)
                grid%dxt(i + 1) = grid%dxt(i)
                ions%Z(i + 1) = ions%Z(i)
            end do
            grid%xt(i + 1) = xxt
            grid%x0(i + 1) = xx0
            grid%xtold(i + 1) = xxtold
            ions%v(i + 1) = vv
            ions%vint(i + 1) = vvint
            ions%id0(i + 1) = ii
            ions%q(i + 1) = qq
            flds%nhot(i + 1) = nn
            flds%phi(i + 1) = pp
            flds%E(i + 1) = EE
            flds%gradnh(i + 1) = gg
            grid%dxt(i + 1) = dd
            ions%Z(i + 1) = ZZ
        end do
    end subroutine sort_sheets

end module mod_pusher
