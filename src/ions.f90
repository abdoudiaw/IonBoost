module mod_ions
    use mod_types, only: dp
    use mod_grid, only: grid_t
    implicit none
    private
    public :: ion_species_t, init_ions, rebuild_density, ion_kinetic_energy

    ! Cold Lagrangian ion sheets. Sheet i carries a fixed charge q(i) between
    ! nodes; the density is rebuilt from the node spacing every step.
    ! A second ion species would be a second instance of this type with its own
    ! Z (and sheet charge), pushed by the same field; the Poisson assembly then
    ! sums the species densities.
    type :: ion_species_t
        integer :: n                          ! sheets indexed 0:n (n = ncell)
        real(dp), allocatable :: Z(:)         ! 0:n charge state per sheet
        real(dp), allocatable :: q(:)         ! 0:n fixed charge per sheet (qiSS)
        real(dp), allocatable :: n_init(:)    ! 0:n initial density (niSS)
        real(dp), allocatable :: v(:)         ! 0:n velocity at whole steps
        real(dp), allocatable :: vint(:)      ! 0:n velocity at half steps
        integer,  allocatable :: id0(:)       ! 0:n original sheet index
        real(dp), allocatable :: dens(:)      ! 0:n node density ni
        real(dp), allocatable :: dens_half(:) ! 1:n midpoint density ni1s2
    end type ion_species_t

contains

    subroutine init_ions(ions, grid, ni0, Zion, prog)
        type(ion_species_t), intent(out) :: ions
        type(grid_t), intent(in) :: grid
        real(dp), intent(in) :: ni0, Zion, prog
        integer :: i, n

        n = grid%ncell
        ions%n = n
        allocate(ions%Z(0:n), ions%q(0:n), ions%n_init(0:n))
        allocate(ions%v(0:n), ions%vint(0:n), ions%id0(0:n))
        allocate(ions%dens(0:n), ions%dens_half(1:n))

        ions%Z = Zion
        ions%n_init = ni0
        ions%v = 0.0_dp
        ions%vint = 0.0_dp
        do i = 0, n
            ions%id0(i) = i
        end do

        ions%q(0) = ions%n_init(0) * grid%dx0(1)
        do i = 1, n - 1
            ions%q(i) = ions%n_init(i) * (grid%dx0(i) + grid%dx0(i + 1)) / 2.0_dp
        end do
        ions%q(n) = ions%n_init(n) * grid%dx0(n) * (1.0_dp + prog) / 2.0_dp
        ions%dens = ni0
    end subroutine init_ions

    subroutine rebuild_density(ions, grid)
        type(ion_species_t), intent(inout) :: ions
        type(grid_t), intent(in) :: grid
        integer :: i, n

        n = ions%n
        ions%dens(0) = ions%q(0) / max(grid%dxt(1), 1.0e-12_dp)
        do i = 1, n - 1
            ions%dens(i) = 2.0_dp * ions%q(i) / max(grid%dxt(i) + grid%dxt(i + 1), 1.0e-12_dp)
        end do
        ! one-sided extrapolation at the front sheet
        ions%dens(n) = 2.0_dp * ions%q(n) / max(grid%dxt(n), 1.0e-12_dp) / &
                       max(1.0_dp + 2.0_dp * grid%dxt(n) / max(grid%dxt(n - 1), 1.0e-12_dp) - &
                           grid%dxt(n - 1) / max(grid%dxt(n - 2), 1.0e-12_dp), 1.0e-12_dp)
        do i = 1, n
            ions%dens_half(i) = 0.5_dp * (ions%dens(i - 1) + ions%dens(i))
        end do
    end subroutine rebuild_density

    ! Kinetic energy; the front sheet of an initially sharp ('step') profile
    ! carries half weight (it represents half a cell)
    function ion_kinetic_energy(ions) result(En)
        type(ion_species_t), intent(in) :: ions
        real(dp) :: En
        integer :: i
        En = 0.0_dp
        do i = 1, ions%n - 1
            En = En + ions%q(i) * ions%v(i)**2 / 2.0_dp / max(ions%Z(i), 1.0e-12_dp)
        end do
        En = En + ions%q(ions%n) * ions%v(ions%n)**2 / 4.0_dp / max(ions%Z(ions%n), 1.0e-12_dp)
    end function ion_kinetic_energy

end module mod_ions
