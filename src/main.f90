program Mora
    use iso_fortran_env, only: int32, real32, int64, real64
    use mod_share
    use mod_solver, only: gauss
    use mod_io, only: read_input
    implicit none
    !
    call read_input()
    !
    Tnorm = T_MeV / 0.511
    e0 = 0.5 * Tnorm * (1 + 9 * Tnorm / 4 + 3 * Tnorm ** 2 / 4) / (1 + 3 * Tnorm / 2 + 3 * Tnorm ** 2 / 8)
    Tn0 = 2 * e0
    Th = Thmax
    !
    if (Th == 0) then
        write(*,*) 'Hot temperature should not be zero.'
        write(*,*) 'Aborting...'
        stop
    endif
    !
    phi(ncell) =  n0hot * Th / n0hot
    E(ncell) = n0hot * Th * exp(-phi(ncell) / Th)
    alpha = E(ncell) / phi(ncell)
    !
    do i = ncell + 1, ntotal
        v(i) = 0.
        ni(i) = 1.d-12
    enddo
    !
    E(0) = 0.d0
    grad_hot_density(0) = 0.d0
    time = 0.d0
    dt = dti
    nhm2 = n0hot
    nhm1 = n0hot
    Thm2 = Thmax
    Thm1 = Thm2
    Whot1 = 0.d0
    Whot2 = 0.d0
    !
    open(unit = 9, file = 'conservation.txt', status = 'unknown')
    open(unit = 10, file = 'historique.txt', status = 'unknown')
    write(9, *) '# time nti nthot nte n0hot &
                    En_ion Whot1 Whot2 Whot &
                    Th En_elec vmax vfinal'
    write(10, *) '# time xi v(ivmax) Energy_max E(ivmax) ne(ivmax) &
                     ni(ivmax) ni(0) nhot(0) lDebye lgrad_e_pressure &
                     lgrad_i_pressure ivmax idebut(ivmax)'
!
!  2 - time loop
    do 100 itime=1,itmax
        if(itime.le.3) then
           iter=iter0
        else
           iter=iter1
        endif
!
! 2.1.1 - Estimation of temperature in case where its dependence is fixed by the functions Thot
    Thold = Th
    if(.not.En_cons)then
        Th = Thmax
        idico = 0
        if(Th == 0)then
            write(*,*) 'Hot electron temperature should not be zero!'
            stop
        endif
        if((Th - Thold) / Thold > 0.05) then
            iter = iter0
        endif
        if((Th - Thold) / Thold < -0.05) then
            iter = iter0
        endif
    endif

! 2.1.2 - Estimation of the change in n0hot due to conservation of hot electrons number
    if(nb_cons)then
        nhm3 = nhm2
        nhm2 = nhm1
        nhm1 = n0hot
        n0hot = 3 * nhm1 - 3 * nhm2 + nhm3
    endif

! 2.1.3 - Estimation of the change of Th due to energy conservation
    if(En_cons)then
        Thm3 = Thm2
        Thm2 = Thm1
        Thm1 = Th
        Th = 3 * Thm1 - 3 * Thm2 + Thm3
    endif

    if(itime == (itime / nsort) * nsort) then
        print *, time
    endif
!
! 2.1.4 - Memory of old values
    if(itime > 1) then
        do i = 0, ntotal
            nhotold(i) = nhot(i)
            phiold(i) = phi(i)
            Eold(i) = E(i)
            ghold(i) = grad_hot_density(i)
        enddo
        do i = ncell + 1, ntotal
            xtold(i) = xt(i)
        enddo
        do i = 1, ntotal
            dxtold(i) = dxt(i)
        enddo
        Whot1old = Whot1
        Whot2old = Whot2
    endif
! 2.1.5 - Calculate new ionic density
    do i = 1, ncell
        dxt(i) = xt(i) - xt(i-1)
    end do
    !
    ni(0) = qiSS(0) / dxt(1)
    do i = 1, ncell-1
        ni(i) = 2 * qiSS(i) / (dxt(i) + dxt(i+1))
    end do
    !
    ni(ncell) = 2 * qiSS(ncell) / dxt(ncell) / (1 + 2 * dxt(ncell) / dxt(ncell-1) - dxt(ncell-1) / dxt(ncell-2))
    !
    do i = 1, ncell
        ni1s2(i) = (ni(i-1) + ni(i)) / 2
    end do
    !
! 2.2 - Iteration for the estimation of phi, Th
    do iphi = 1, iter
! 2.2.1 - Construction of arrays a, b, c, and f
        neh = n0hot * exp(-phi(0) / Th)
        b(0) = -2 / (dxt(1) * dxt(1)) - neh / Th
        c(0) = 2 / (dxt(1) * dxt(1))
        f(0) = ni(0) - neh * (1 + phi(0) / Th)
    !
    do i = 1, ncell-1
        neh = n0hot * exp(-phi(i) / Th)
        a(i) = 2 / dxt(i) / (dxt(i) + dxt(i+1))
        b(i) = -2 / (dxt(i) * dxt(i+1)) - neh / Th
        c(i) = 2 / dxt(i+1) / (dxt(i) + dxt(i+1))
        f(i) = ni(i) - neh * (1 + phi(i) / Th)
    enddo
    neh = n0hot * exp(-phi(ncell) / Th)
    neold = neh
    peold = neh * Th
    a(ncell) = 2 / (dxt(ncell) * dxt(ncell))
    b(ncell) = -a(ncell) - sqrt(2 / peold) * neold / dxt(ncell) - neh / Th
    f(ncell) = ni(ncell) - neh * (1 + phi(ncell) / Th) - sqrt(8 * peold) / dxt(ncell) * (1 + 0.5 * neold * phi(ncell) / peold)
! 2.2.2 - Inversion of the tridiagonal matrix
    call gauss(nmax,ncell,a,b,c,f,bx,fx,phi)
    if(iphi.ne.iter.and.itime.eq.1) cycle
    if(iphi.ne.iter.and.(((.not.nb_cons).and.(.not.En_cons)))) cycle


! 2.2.3 - Determine the electronic density
        do i = 0, ncell
            nhot(i) = n0hot * exp(-phi(i) / Th)
            ne(i) = nhot(i)
            rho(i) = ni(i) - ne(i)
        enddo
!
! 2.2.4 - Determine the electronic pressure and densities grad_i_pressureents
            phot(ncell) = nhot(ncell) * Th
            pe(ncell) = phot(ncell)
            grad_e_pressure = (log(ne(ncell)) - log(ne(ncell - 1))) / (xt(ncell) - xt(ncell - 1))
            grad_i_pressure = (log(ni(ncell)) - log(ni(ncell - 1))) / (xt(ncell) - xt(ncell - 1))
!
! 2.2.5 - Determine the electronic density in the vacuum sheath (empty of ions)
            kvide = sqrt(n0hot / 2.d0 / Th)
            uphi = 0.d0
            do i = ncell + 1, ntotal
                if (i == (ncell + 1)) then
                    dxt(ncell + 1) = 0.
                else
                    dxt(i) = 0.0250 * sqrt(Th / nhot(i - 1))
                endif
                xt(i) = xt(i - 1) + dxt(i)
                uphi = uphi + kvide * sqrt(1.d0 + pe(i - 1) / phot(i - 1)) * dxt(i)
                phi(i) = phi(ncell) + 2.d0 * Th * log(1.d0 + uphi * exp(-0.5d0 * phi(ncell) / Th))

                nhot(i) = n0hot * exp(-phi(i) / Th)
                ne(i) = nhot(i)
                rho(i) = -ne(i)
                phot(i) = nhot(i) * Th
                pe(i) = phot(i)
                E(i) = sqrt(2.d0 * pe(i))
           enddo
!
! Determine the total number of electrons and ions
            nti = 0.d0
            do i = 1, ncell
                nti = nti + ni1s2(i) * dxt(i)
            end do
            !
            nthot = 0.d0
            do i = 1, ntotal
                nthot = nthot + 0.5d0 * (nhot(i-1) + nhot(i)) * dxt(i)
            end do
            nthot = nthot + E(ntotal)
            nte = nthot
            !
! 2.2.7 - Adjust n0hot to conserve the total number of hot electrons
! Adjust n0hot to conserve the total number of hot electrons
            if (nb_cons) then
                if (itime == 1) then
                nthot0 = nthot
            else
                if (iphi > 1) n0hot = n0hot * (nthot0 / nthot)
                endif
            endif
            !
            if (itime > 1) then
                do i = ncell+1, ntotal
                    vint(i) = (xt(i) - xtold(i)) / dt
                end do
            endif
!
! Calculate the grad_i_pressureents of electron density
          do i=1,ntotal-1
             grad_hot_density(i)=(nhot(i+1)-nhot(i-1))/(dxt(i)+dxt(i+1))
          enddo
         !
         grad_hot_density(ncell)=grad_hot_density(ncell-1)*nhot(ncell)/nhot(ncell-1)
         grad_hot_density(ntotal)=grad_hot_density(ntotal-1)*nhot(ntotal)/nhot(ntotal-1)
         grad_hot_density(ncell+1)=grad_hot_density(ncell)
!
    !*     2.2.9 calcul de l'energie fournie par les electrons chauds [et
    !*           reajustement de Th (ancienne methode)]

          if(itime.gt.1.and.iphi.gt.1)then
             dWhot(0)=0.5d0*(phiold(0)+phi(0))*(nhot(0)-nhotold(0)-vint(0)*dt&
                *0.5d0*(grad_hot_density(0)+ghold(0)))
             Whot1=Whot1old+0.25d0*dWhot(0)*(dxtold(1)+dxt(1))
             SWhot(0)=0.25d0*dWhot(0)*(dxtold(1)+dxt(1))
             do i=1,ncell
                dWhot(i)=0.5d0*(phiold(i)+phi(i))*(nhot(i)-nhotold(i)&
                    -vint(i)*dt*0.5d0*(grad_hot_density(i)+ghold(i)))
                Whot1=Whot1+0.25d0*dWhot(i)*(dxtold(i)+dxt(i)+dxtold(i+1)+dxt(i+1))
                SWhot(i)=SWhot(i-1)+0.25d0*dWhot(i)*(dxtold(i)+dxt(i)+dxtold(i+1)+dxt(i+1))
             enddo
             Whot2=Whot2old
             do i=ncell+1,ntotal-1
                dWhot(i)=0.5d0*(phiold(i)+phi(i)) *(nhot(i)-nhotold(i)&
                    -vint(i)*dt*0.5d0*(grad_hot_density(i)+ghold(i)))
                Whot2=Whot2+0.25d0*dWhot(i)*(dxtold(i)+dxt(i)+dxtold(i+1)+dxt(i+1))
                SWhot(i)=SWhot(i-1)+0.25d0*dWhot(i)*(dxtold(i)+dxt(i)+dxtold(i+1)+dxt(i+1))
             enddo
             dWhot(ntotal)=0.5d0*(phiold(ntotal)+phi(ntotal))&
                    *(nhot(ntotal)-nhotold(ntotal)-vint(ntotal)&
                    *dt*0.5d0*(grad_hot_density(ntotal)+ghold(ntotal)))
             Whot2=Whot2+0.25d0*dWhot(ntotal) *(dxtold(ntotal)+dxt(ntotal))
             SWhot(ntotal)=SWhot(ntotal-1)+0.25d0*dWhot(ntotal)*(dxtold(ntotal)+dxt(ntotal))
             Whot=Whot1+Whot2
           endif
! * 2.2.11 calculate electric field throughout the relaxation
        do i = 1,ncell
            E1s2(i)=(phi(i)-phi(i-1))/dxt(i)
        end do
        do i = 1,ncell-1
            E(i)=(dxt(i+1)*E1s2(i)+dxt(i)*E1s2(i+1))/(dxt(i)+dxt(i+1))
        end do
        E(ncell)=sqrt(2.0*pe(ncell))
        ! 2.2.12 adjust velocity at 'whole' times
        if(itime.gt.1) then
            do i = 1,ncell
                v(i)=vint(i)+dt*(3.0*E(i)+Eold(i))charge(i)/8.
            end do
        endif
        ! 2.2.13 calculate kinetic energy of ions
        En_ion=0.
        do i = 1,ncell-1
            En_ion=En_ion+qiSS(i)*v(i)**2/2.d0/charge(i)
        end do
        if(profil.eq.'step') En_ion = En_ion+qiSS(ncell)*v(ncell)2/4.d0/charge(i)
        if(profil.ne.'step') En_ion = En_ion+qiSS(ncell)v(ncell)2/2.d0/charge(i)
        ! 2.2.14 calculate electrostatic energy
        En_elec=0.
        do i = 1,ntotal-1
            En_elec=En_elec+0.25d0(E(i)**2)(dxt(i)+dxt(i+1))
        end do
        En_elec=En_elec+0.25d0(E(ntotal)**2)dxt(ntotal)
        En_elec=En_elec+E(ntotal)Th
        ! 2.2.15 calculate total energy
        if(itime.eq.1)then
            En_totale=En_ion+En_elec+En_hot0
        endif
        !
        !* 2.2.16 Adjustment of Th
        !* (replaces previous adjustment)
        !
        if(En_cons.and.itime.gt.1.and.iphi.gt.1) then
            enew=(En_totale-En_ion-En_elec)Tn0/nthot/Thmax
            Thnew=2.Thmaxenew/Tn0
            !
            if(iphi.ne.iter) then
                Th=(iter2Th+Thnew)/(iter2+1)
            else
                Th=(iter3*Th+Thnew)/(iter3+1)
            endif
        endif
10 continue
        !* 2.3 Calculation of the theoretical electrostatic field
        !* 2.3.1 Normalization of the fields with self-similar fields
        if(itime.eq.1)then
            Enorm=0.d0
        else
            Ess=sqrt(Th)/time
            Enorm=E(ncell)/Ess
        endif
        !* 2.3.2 Comparison with the 'theoretical' field, supposed to be 'fit' at the boundary at all times
        Ebord0=sqrt(2.d0n0hot)/exp(0.5d0)
        Ebordth=2.d0Ebord0sqrt(Th)/sqrt(4.d0+(Ebord0time)2)
        ! 2.4 Stop test based on the number of iterations
        if(itmax.eq.1) exit
        ! 2.5 Time-dependent functions
        vmax=0.d0
        ivmax=ncell
        do i=1,ncell
            if(v(i)>vmax) then
                vmax=v(i)
                ivmax=i
            endif
        end do
        vfinal=v(ivmax)+E(ivmax)*charge(ivmax)*time
        lDebye=sqrt(Th/ne(ncell))
        lgrad_e_pressure=-1.0/grad_e_pressure
        if(itime.eq.1)then
            lgrad_i_pressure=0.
        else
            lgrad_i_pressure=-1.0/grad_i_pressure
        endif
        nstep=(tmax/dt)
        nsort=1
        if(nstep>2000)nsort=nstep/2000
        if(itime-1 != ((itime-1)/nsort)nsort .and. .not.laststep) go to 102
        !
        write(9,'(f8.2,13(f10.5))') time,nti,nthot,nte,n0hot,En_ion,&
        Whot1,Whot2,Whot,Th,En_elec,&
        En_elec+En_ion-Whot,vmax,vfinal
        !
        write(10,'(f8.2,f12.5,f7.4,f8.4,f7.5,2(f10.7),2(f8.6),3(f10.5),i6,i6)') time, &
        xt(ivmax),v(ivmax),0.5*v(ivmax)**2,E(ivmax),ne(ivmax),ni(ivmax),&
        ni(0),nhot(0),lDebye,lgrad_e_pressure,lgrad_i_pressure,ivmax,idebut(ivmax)
        !* 2.6 End of the time loop
        101 continue
        !
        !* 2.7 Write final results
        write(9,'(f8.2,13(f10.5))') time,nti,nthot,nte,n0hot,En_ion,&
        Whot1,Whot2,Whot,Th,En_elec,&
        En_elec+En_ion-Whot,vmax,vfinal
        write(10,'(f8.2,f12.5,f7.4,f8.4,f7.5,2(f10.7),2(f8.6),3(f10.5),i6,i6)') time, &
        xt(ivmax),v(ivmax),0.5*v(ivmax)**2,E(ivmax),ne(ivmax),ni(ivmax),&
        ni(0),nhot(0),lDebye,lgrad_e_pressure,lgrad_i_pressure,ivmax,idebut(ivmax)
        Rs=rgauss/ni(0)
        xi=(xt(ivmax)+lmax)/Rs
        tsR0=time/rgauss
        RsR0=Rs/rgauss
        Eq42=2.d0xi*exp(xi**2/2.d0)/rgauss/sqrt(RsR0)
        end
!*         2.6 Test for stopping the calculation
        if(time >= tmax-1e-5) exit
        if(itime == itmax) exit
        if( (itime-1) / nsort == (itime-1) / nsort) cycle
!
!*         2.7 Increment time
        dtold=dt
        dt=dti
        if((time+dt) > (tmax-1e-5)) then
           dt=tmax-time
           laststep=.true.
        endif
        time=time+dt
!
!*         2.8 modification of the velocity at intermediate times
!*             'half-integers'
        if(itime == 1) then

            do i = 1,ncell
                  vint(i) = v(i) + 0.5d0*dt*E(i)*charge(i)
               end do
        else
    !*            2.8.1 - construction of the arrays a, b, c, and f for
    !*              the implicit calculation of velocity with viscosity
           delt = 0.5 * (dtold + dt)
           b(1) = 1.d0/delt
           c(1) = 0.d0
           f(1) = 0.d0
           do i = 1, ncell - 1
              a(i) = -2.d0 * nu / dxt(i) / (dxt(i) + dxt(i + 1))
              b(i) = 1.d0 / delt + 2.d0 * nu / (dxt(i) * dxt(i + 1))
              c(i) = -2.d0 * nu / dxt(i + 1) / (dxt(i) + dxt(i + 1))
              f(i) = E(i) * charge(i) + vint(i) / delt
           end do
           a(ncell) = 0.d0
           b(ncell) = 1.d0 / delt
           f(ncell) = E(ncell) * charge(ncell) + vint(ncell) / delt
    !*          2.8.2 - inversion of the tridiagonal matrix
           call gauss(ncell, a, b, c, f, vint)
      endif
!*         2.9 modification of the position
    do i = 1,ncell
        xt(i) += dt * vint(i)
        if(xt(i) < x0(1)) then
            xt(i) = 2.d0 * x0(1) - xt(i)
            vint(i) = -vint(i)
        endif
    end do
!
!*         2.10 rearrangement of ions numbers
    do j = 2, ncell
        xxt = xt(j)
        xx0 = x0(j)
        xxtold = xtold(j)
        vv = v(j)
        vvint = vint(j)
        iidebut = idebut(j)
        qqiSS = qiSS(j)
        nnhot = nhot(j)
        pphi = phi(j)
        EE = E(j)
        ggrad_hot_density = grad_hot_density(j)
        ddxt = dxt(j)
        ccharge = charge(j)
        do i = j - 1, 1, -1
            if(xt(i) <= xxt) exit
            xt(i + 1) = xt(i)
            x0(i + 1) = x0(i)
            xtold(i + 1) = xtold(i)
            v(i + 1) = v(i)
            vint(i + 1) = vint(i)
            idebut(i + 1) = idebut(i)
            qiSS(i+1)=qqiSS
            nhot(i+1)=nnhot
            phi(i+1)=pphi
            E(i+1)=EE
            grad_hot_density(i+1)=ggrad_hot_density
            dxt(i+1)=ddxt
            charge(i+1)=ccharge
        end do
!* 2.11 incrementation du nombre d'iterations
        itime=itime+1
    end do
!* 3.2 final results display
    write(*, '("Final time = ", f8.4)') time
    write(*, '("Number of iterations = ", i0)') itime

!
!*    3 - Estimate self-similar density
      if(itime == 1) then
         do i=ncell+1,ntotal
            x0(i) = xt(i)
         enddo
      endif
    if(itmax == 1 .or. time <= 1.d-04) return
    do i=0,ntotal
       xi = xt(i)/sqrt(Th)/time
       if(xi < -1.d0) then
          nss(i) = n0hot
       else
          nss(i) = n0hot*exp(-(xi+1.d0))
       endif
    enddo
    do i=0,ntotal
       difth(i) = (ni(i) - nss(i))/nss(i)
    enddo
    print *, "Step 3 successfully completed!"
!
!*        4 - outputs
    open(unit=11,file='profil.txt',status='replace')
    write(11,*) '# x0 xt v phi E ni ne rho nhot nss difnor dWhot SWhot charge'
    !
    do i=0,ncell,iline
       write(11,'(f12.4,x,f12.4,x,f7.4,2(x,f9.6),5(x,f12.8))') &
            x0(i),xt(i),v(i),phi(i),&
            E(i),ni(i),ne(i),rho(i),&
            nhot(i),nss(i)
    enddo
    !
    do i=ncell+1,ntotal
        write(11,'(f12.4,x,f12.4,x,f7.4,2(x,f9.6),5(x,f12.8))') &
            x0(i),xt(i),v(i),phi(i),&
            E(i),ni(i),ne(i),rho(i),&
            nhot(i),nss(i)
    enddo
    close(11)

!
!* 4.2 Velocity and energy spectra
    open(12,file='spectres.txt', status='replace')
    write(12,*) '# vmoy Emoy dndv dndE'
    !
    do i=1,ncell
       vmoy(i)=0.5*(v(i-1)+v(i))
       Emoy(i)=0.25*(v(i-1)**2+v(i)**2)
       dndv(i)=(niSS(i-1)+niSS(i))*dx0(i)/(v(i)-v(i-1))/2.d0
       dndv(i)=abs(dndv(i))
       dndE(i)=(niSS(i-1)+niSS(i))*dx0(i)/(v(i)**2-v(i-1)**2)
       dndE(i)=abs(dndE(i))
    end do
    !
    do i=1,ncell,iline
       if(dndE(i)<1.d-3) then
          write(12,'(f7.4,a,f8.4,a,f10.6,a,f10.6)') vmoy(i),char(9),Emoy(i),char(9),dndv(i),char(9),dndE(i)
       endif
    end do
    close(12)
    !
    write(*,*) 'time', time
end

