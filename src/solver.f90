module mod_solver
    use mod_types
    implicit none
    private
    public :: Thot_, Tcold_, gauss
    contains
        pure function Thot_(multiphase,trise,time) result(res)
        real(dp), intent (in):: time
        real(dp), intent (in):: trise
        real(dp) :: res
        logical, intent(in) :: multiphase
        !
        if(multiphase.and.time.le.trise) then
         res=0.01_dp+0.99_dp*time/trise
        else
        res=1._dp
        endif
        end function Thot_
      !
       pure function Tcold_(time) result(res)
        real(dp), intent (in):: time
        real(dp) :: res
        if(time.le.9.9_dp) then
          res=1._dp
        else
          res=1._dp
        endif
       end function Tcold_
      !
      pure subroutine gauss(nmax,ncell,A,B,C,F,bx,fx,phi)
        ! determine the electrostatic potential
        !
        real(dp), intent(in out) :: phi(0:nmax)
        real(dp), intent(in) :: A(1:nmax), C(0:nmax-1)
        real(dp), intent(in) :: B(0:nmax), F(0:nmax)
        real(dp), intent(in out) :: bx(0:nmax)
        real(dp), intent(in out) :: fx(0:nmax)
        integer(sp), intent(in) :: nmax
        integer(sp), intent(in) :: ncell
        integer(sp) :: j
        !
        Bx(0) = B(0)
        Fx(0) = F(0)
        do concurrent(j = 1:ncell)
            Bx(j) = B(j) - A(j)*C(j-1) /Bx(j-1)
            Fx(j) = F(j) - A(j)*Fx(j-1)/Bx(j-1)
        end do
            phi(ncell) = Fx(ncell)/Bx(ncell)
        do concurrent(j=ncell-1:0:-1)
            phi(j) = ( Fx(j)-C(j)*phi(j+1) )/Bx(j)
        end do
      end subroutine gauss
end module mod_solver
