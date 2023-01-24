module mod_solver
    implicit none
    private
    public :: gauss
    contains
      subroutine gauss(nmax,ncell,A,B,C,F,bx,fx,phi)
        ! Determine the electrostatic potential
        implicit none
        integer, intent(in) :: nmax, ncell
        real, intent(in) :: A(1:nmax), B(0:nmax), C(0:nmax-1), F(0:nmax)
        real, intent(inout) :: bx(0:nmax), fx(0:nmax), phi(0:nmax)
        integer :: j

        ! Initialize Bx and Fx
        bx(0) = B(0)
        fx(0) = F(0)
        
        ! Error handling
        if (bx(0) == 0) then
            print *, "Error: Division by zero in Bx(0)"
            stop
        end if

        ! Forward loop
        do j = 1,ncell
            bx(j) = B(j) - A(j)*C(j-1) /bx(j-1)
            fx(j) = F(j) - A(j)*fx(j-1)/bx(j-1)
            ! Error handling
            if (bx(j) == 0) then
                print *, "Error: Division by zero in Bx(" ,j,")"
                stop
            end if
        end do
        
        phi(ncell) = fx(ncell)/bx(ncell)
        
        ! Error handling
        if (bx(ncell) == 0) then
            print *, "Error: Division by zero in Bx(" ,ncell,")"
            stop
        end if
        ! Backward loop
        do j = ncell-1, 0, -1
            phi(j) = (fx(j) - C(j)*phi(j+1)) / bx(j)
             ! Error handling
            if (bx(j) == 0) then
                print *, "Error: Division by zero in Bx(" ,j,")"
                stop
            end if
        end do
      end subroutine gauss
end module mod_solver

