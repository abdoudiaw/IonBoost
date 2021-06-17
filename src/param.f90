module mod_io
use iso_fortran_env, only: int32, real32, int64, real64
use mod_math, only: Thot_, Tcold_
use mod_share

implicit none

private
public :: read_input

contains

    subroutine read_input()
        ! Read the initial conditions from a data file
        character(len=*), intent(in) :: filename
        integer :: iunit,ierr
        open(newunit=iunit,file='IonBoost28.in',status='old',action='read',iostat=ierr)
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
            read(iunit,*)   profil,VTT
            read(iunit,*)   Ztest,multicouche,lay1,lay2
            read(iunit,*)   mix,charge2,multiphase,trise
        else
            write(*,*) '*** Error reading the input file *** '
            write(*,*) 'Aborting...'
            error stop
        endif  

        close(iunit)
    !    Check whether everything is fine in the file

        if(profil.ne.'sack'.and.profil.ne.'step'.and.profil.ne.'expo'.and.profil.ne.'gaus') then
              write (*,*) 'What is the initial density profile'
              stop
        endif
           
        if((.not.lfini).and.(nb_cons.or.En_cons))then
              nb_cons=.false.
              En_cons=.false.
              write(*,*)    'lfini=.false. therefore En_cons=.false. and nb_cons=.false.'
        endif
            
        if((.not.nb_cons).and.En_cons) then
          nb_cons=.true.
          write(*,*)    'En_cons=.true.  therefore  nb_cons=.true.'
        endif

        if(nb_cons.and.(n0cold.lt.1.d-04))then
          nb_cons=.false.
          write(*,*)  'n0cold.lt.1.d-04  therefore  nb_cons=.false.'
        endif
        
       ntotal=ncell+nvacuum
       
       if(ntotal.gt.nmax)then
         write(*,*) 'ntotal=ncell+nvacuum is bigger than nmax'
         stop
       endif

        nstep=(tmax/dti+0.5d0)
        nsort=1
        if(nstep.gt.2000) nsort=nstep/2000
        laststep=.false.

        iline=1
        if(ncell.gt.3999) iline=ncell/2000

    !   Determine the length of the plasma

      if(.not.lfini) then
         time=0.d0
         dt=dti
         Th=Thmax
         Tc=Tcmax
         cs2=(n0cold+n0hot)/(n0cold/Tc+n0hot/Th)
         cs=sqrt(cs2)

         if(profil.eq.'sack')length=10.d0*LSS+25.d0*cs
         if(profil.eq.'expo')length=LSS+25.d0*cs
         if(profil.eq.'step')length=25.d0*cs

         if(itmax.eq.1) goto 91

         Th=Thmax*Thot_(multiphase,trise,0.d0)
         Tc=Tcmax*Tcold_(0.d0)


        do 90 itime=1,itmax
           if(time.ge.(tmax-1.d-5)) goto 91
           if(itime.eq.itmax) goto 91
           if((time+dt).gt.(tmax-1.d-5)) dt=tmax-time

           cs2old=cs2
           time=time+dt
           Th=Thmax*Thot_(multiphase,trise,time)
           Tc=Tcmax*Tcold_(time)
           cs2=(n0cold+n0hot)/(n0cold/Tc+n0hot/Th)

           cs=0.5d0*(sqrt(cs2old)+sqrt(cs2))
           length=length+cs*dt

    90      continue
    91      continue
        endif

         if(lfini) length=lmax

    !   mesh 

        PI=4.D0*DATAN(1.D0)
        rgauss=2.d0*length/sqrt(PI)
        ni0=n0cold+n0hot
              
        if(prog.eq.1.d0)then
             dx0(1)=length/ncell
        else
             dx0(1)=length*(1.-prog)/(1.-prog**ncell)
        endif
          x0(0)=-length
          x0(1)=x0(0)+dx0(1)

        do i = 0,1
           niSS(i)=ni0
           if(profil.eq.'expo'.and.x0(i).gt.-LSS) niSS(i)=ni0*exp(-(x0(i)+LSS)/LSS)
           if(profil.eq.'sack') niSS(i)=ni0*(2.d0/pi)*DATAN(exp(-x0(i)/LSS))
           if(profil.eq.'gaus')  niSS(i)=ni0*exp(-((x0(i)+length)/rgauss)**2)
        end do

        do i = 2,ncell
           if(i.eq.2) then
                dx0(2)=(1.d0+prog)*dx0(1)*niSS(0)/niSS(1)-dx0(1)
           else
                dx0(i)=(dx0(i-2)+dx0(i-1))*prog*niSS(i-2)/niSS(i-1)-dx0(i-1)
           endif
           if((profil.eq.'expo'.or.profil.eq.'sack').and.dx0(i).gt.dx0(i-1).and.dx0(i).gt.LSS/5.d0) dx0(i)=LSS/5.d0
           if(profil.eq.'gaus' .and.dx0(i).gt.dx0(i-1).and.dx0(i).gt.lmax/5.d0) dx0(i)=lmax/5.d0

           x0(i)=x0(i-1)+dx0(i)
           niSS(i)=ni0
           if(profil.eq.'expo'.and.x0(i).gt.-LSS) niSS(i)=ni0*exp(-(x0(i)+LSS)/LSS)
           if(profil.eq.'sack') niSS(i)=ni0*(2.d0/pi)*DATAN(exp(-x0(i)/LSS))
           if(profil.eq.'gaus') niSS(i)=ni0*exp(-((x0(i)+length)/rgauss)**2)
        end do


        qiSS(0)=niSS(0)*dx0(1)
        do i = 1,ncell-1
           qiSS(i)=niSS(i)*(dx0(i)+dx0(i+1))/2.d0
        end do
        qiSS(ncell)=niSS(ncell)*dx0(ncell)*(1.d0+prog)/2.

    !     Saving the initial order of the ions to deal with breaking

        do i = 0,ncell
           idebut(i)=i
           irang(i)=i
        end do
    end subroutine read_input
    
end module mod_io
