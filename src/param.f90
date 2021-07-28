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
        integer :: iunit,ierr
        open(newunit=iunit,file='input.in',status='old',action='read',iostat=ierr)
        if( ierr == 0 ) then
            read(iunit,*)   ncell,nvacuum,prog
            read(iunit,*)   itmax,tmax
            read(iunit,*)   iter0,iter1,iter2,iter3
            read(iunit,*)   dti
            read(iunit,*)   n0cold,n0hot
            read(iunit,*)   Thmax,Tcmax
            read(iunit,*)   lfini,lmax
            read(iunit,*)   nb_cons,En_cons
            read(iunit,*)   T_MeV,LSS,nLSS
        else
            write(*,*) '*** Error reading the input file *** '
            write(*,*) 'Aborting...'
            stop
        endif
        close(iunit)
        if((.not.lfini).and.(nb_cons.or.En_cons))then
              nb_cons=.false.
              En_cons=.false.
              write(*,*)    'lfini=.false. therefore En_cons=.false. and nb_cons=.false.'
        endif
        if((.not.nb_cons).and.En_cons) then
          nb_cons=.true.
          write(*,*)    'En_cons=.true.  therefore  nb_cons=.true.'
        endif
        if(nb_cons.and.(n0cold.lt.1.e-04))then
          nb_cons=.false.
          write(*,*)  'n0cold.lt.1.d-04  therefore  nb_cons=.false.'
        endif
       ntotal=ncell+nvacuum
       if(ntotal.gt.nmax)then
         write(*,*) 'ntotal=ncell+nvacuum is bigger than nmax'
         stop
       endif
        nstep=(tmax/dti+0.5_dp)
        nsort=1
        if(nstep.gt.2000) nsort=nstep/2000
        laststep=.false.
        iline=1
        if(ncell.gt.3999) iline=ncell/2000
    !   Determine the length of the plasma
      if(.not.lfini) then
        time=0._dp
        dt=dti
        cs2=(n0cold+n0hot)/(n0cold/Tcmax+n0hot/Thmax)
        cs=sqrt(cs2)
        length=25._dp*cs
!
        if((time+dt).gt.(tmax-1.d-5)) dt=tmax-time
            cs2old=cs2
            time=time+dt
            cs2=(n0cold+n0hot)/(n0cold/Tcmax+n0hot/Thmax)
            cs=0.5_dp*(sqrt(cs2old)+sqrt(cs2))
            length=length+cs*dt
        endif
      if(lfini) length=lmax
        !   mesh
        rgauss=2._dp*length/sqrt(pi)
        ni0=n0cold+n0hot
              
        if(prog.eq.1._dp) then
             dx0(1)=length/ncell
        else
             dx0(1)=length*(1.-prog)/(1.-prog**ncell)
        endif
          x0(0)=-length
          x0(1)=x0(0)+dx0(1)
        !
        do i = 0,1
           niSS(i)=ni0
        end do
        do i = 2,ncell
           if(i.eq.2) then
                dx0(2)=(1._dp+prog)*dx0(1)*niSS(0)/niSS(1)-dx0(1)
           else
                dx0(i)=(dx0(i-2)+dx0(i-1))*prog*niSS(i-2)/niSS(i-1)-dx0(i-1)
           endif
           x0(i)=x0(i-1)+dx0(i)
           niSS(i)=ni0
        end do
        qiSS(0)=niSS(0)*dx0(1)
        do i = 1,ncell-1
           qiSS(i)=niSS(i)*(dx0(i)+dx0(i+1))/2._dp
        end do
        qiSS(ncell)=niSS(ncell)*dx0(ncell)*(1._dp+prog)/2.
    !  Saving the initial order of the ions to deal with breaking
        do i = 0,ncell
           idebut(i)=i
           irang(i)=i
        end do
    end subroutine read_input
end module mod_io
