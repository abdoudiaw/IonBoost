module mod_grid
    use mod_types, only: dp
    use mod_config, only: run_config
    implicit none
    private
    public :: grid_t, build_grid, update_ion_spacing

    ! Node-centered 1D Lagrangian grid. Nodes 0:ncell ride with the ion sheets
    ! (node 0 is the rear wall); nodes ncell+1:ntotal discretize the ion-free
    ! vacuum sheath and are rebuilt from the field solution every step.
    ! Going multi-D means replacing this type (and the assembly loops that
    ! consume dxt) while the physics modules keep their interfaces.
    type :: grid_t
        integer :: ncell, nvacuum, ntotal
        real(dp), allocatable :: x0(:)             ! 0:ntotal initial positions
        real(dp), allocatable :: xt(:), xtold(:)   ! 0:ntotal current / previous
        real(dp), allocatable :: dx0(:)            ! 1:ncell initial spacing
        real(dp), allocatable :: dxt(:), dxtold(:) ! 1:ntotal current / previous spacing
    end type grid_t

contains

    subroutine build_grid(grid, cfg)
        type(grid_t), intent(out) :: grid
        type(run_config), intent(in) :: cfg
        integer :: i

        grid%ncell = cfg%ncell
        grid%nvacuum = cfg%nvacuum
        grid%ntotal = cfg%ntotal
        allocate(grid%x0(0:grid%ntotal), grid%xt(0:grid%ntotal), grid%xtold(0:grid%ntotal))
        allocate(grid%dx0(1:grid%ncell))
        allocate(grid%dxt(1:grid%ntotal), grid%dxtold(1:grid%ntotal))

        ! geometric progression: uniform density, cells shrinking towards the
        ! front by the factor prog so the front stays resolved as it accelerates
        if (cfg%prog == 1.0_dp) then
            grid%dx0(1) = cfg%length / cfg%ncell
        else
            grid%dx0(1) = cfg%length * (1.0_dp - cfg%prog) / (1.0_dp - cfg%prog**cfg%ncell)
        end if
        grid%x0(0) = -cfg%length
        grid%x0(1) = grid%x0(0) + grid%dx0(1)
        do i = 2, grid%ncell
            if (i == 2) then
                grid%dx0(2) = (1.0_dp + cfg%prog) * grid%dx0(1) - grid%dx0(1)
            else
                grid%dx0(i) = (grid%dx0(i - 2) + grid%dx0(i - 1)) * cfg%prog - grid%dx0(i - 1)
            end if
            grid%x0(i) = grid%x0(i - 1) + grid%dx0(i)
        end do

        do i = 0, grid%ncell
            grid%xt(i) = grid%x0(i)
            grid%xtold(i) = grid%x0(i)
        end do
        do i = 1, grid%ncell
            grid%dxt(i) = grid%dx0(i)
            grid%dxtold(i) = grid%dx0(i)
        end do
        ! vacuum nodes start collapsed onto the ion front
        do i = grid%ncell + 1, grid%ntotal
            grid%x0(i) = grid%xt(grid%ncell)
            grid%xt(i) = grid%xt(grid%ncell)
            grid%xtold(i) = grid%xt(i)
            grid%dxt(i) = 0.0_dp
            grid%dxtold(i) = 0.0_dp
        end do
    end subroutine build_grid

    subroutine update_ion_spacing(grid)
        type(grid_t), intent(inout) :: grid
        integer :: i
        do i = 1, grid%ncell
            grid%dxt(i) = grid%xt(i) - grid%xt(i - 1)
        end do
    end subroutine update_ion_spacing

end module mod_grid
