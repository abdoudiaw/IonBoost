module mod_linear_solver
    use mod_types, only: dp
    implicit none
    private
    public :: linear_solver_t, thomas_solver_t

    ! Seam between the physics (matrix assembly in mod_fields / mod_pusher) and
    ! the linear algebra. The 1D discretization is tridiagonal; a GPU or 3D
    ! backend (e.g. AmgX) implements this same interface with its own matrix
    ! storage, leaving the assembly and the Newton iteration untouched.
    type, abstract :: linear_solver_t
    contains
        procedure(solve_iface), deferred :: solve
    end type linear_solver_t

    abstract interface
        subroutine solve_iface(this, n, a, b, c, f, x)
            import :: dp, linear_solver_t
            class(linear_solver_t), intent(inout) :: this
            integer, intent(in) :: n
            real(dp), intent(in) :: a(1:n)      ! sub-diagonal
            real(dp), intent(in) :: b(0:n)      ! diagonal
            real(dp), intent(in) :: c(0:n-1)    ! super-diagonal
            real(dp), intent(in) :: f(0:n)      ! right-hand side
            real(dp), intent(inout) :: x(0:n)   ! solution
        end subroutine solve_iface
    end interface

    type, extends(linear_solver_t) :: thomas_solver_t
        real(dp), allocatable :: bx(:), fx(:)
    contains
        procedure :: solve => thomas_solve
    end type thomas_solver_t

contains

    subroutine thomas_solve(this, n, a, b, c, f, x)
        class(thomas_solver_t), intent(inout) :: this
        integer, intent(in) :: n
        real(dp), intent(in) :: a(1:n), b(0:n), c(0:n-1), f(0:n)
        real(dp), intent(inout) :: x(0:n)
        integer :: j

        if (.not. allocated(this%bx)) then
            allocate(this%bx(0:n), this%fx(0:n))
        else if (ubound(this%bx, 1) < n) then
            deallocate(this%bx, this%fx)
            allocate(this%bx(0:n), this%fx(0:n))
        end if

        this%bx(0) = b(0)
        this%fx(0) = f(0)
        if (abs(this%bx(0)) < 1.0e-12_dp) stop 'thomas_solve: zero pivot at j=0'
        do j = 1, n
            this%bx(j) = b(j) - a(j) * c(j - 1) / this%bx(j - 1)
            this%fx(j) = f(j) - a(j) * this%fx(j - 1) / this%bx(j - 1)
            if (abs(this%bx(j)) < 1.0e-12_dp) stop 'thomas_solve: zero pivot in forward sweep'
        end do

        x(n) = this%fx(n) / this%bx(n)
        do j = n - 1, 0, -1
            x(j) = (this%fx(j) - c(j) * x(j + 1)) / this%bx(j)
        end do
    end subroutine thomas_solve

end module mod_linear_solver
