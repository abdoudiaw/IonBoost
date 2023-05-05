program Ionboost

      ! implicit none
      parameter (nmax=150100,njmax=200)

      implicit double precision (a-h,o-z)
      implicit integer(kind=4) (i-n)
      
      integer(kind=4) idebut(0:nmax),irang(0:nmax),m,djout

      double precision length,ni0,n0cold,n0hot,neh,nec,neold,kvide
      double precision nti,nthot,ntcold,nte,lmax,nte1,nte2, nthot0,LSS,nombre
      double precision nhm3,nhm2,nhm1,nu,nnhot,nncold
      double precision lDebye,lgrade,lgradi,nLSS,tb,vmax0,nii,vv,qq,k1,k2
      double precision k3,k4,nhot1,vmax,Emax,C1,C2,D1,D2,a0,a1,a2,vfinal,vfinal2
      double precision borne2,borne1,Integ,a3,a4,xr,temp,ekmoy,Inv2,vintr,E1s2r
      double precision p1,p2,p3,p4,Hold,Hnew,venew,phi0,phiconv,slop,Et1,Et2,t1,t2
      double precision phimoy,path,xbar,Rbar,phibar,Integ3,cs3,cscin,csmoyen

      double precision erfc
    
      double precision   x0 (0:nmax), xt (0:nmax), xtold(0:nmax)
     double precision v  (0:nmax), vint (0:nmax), vemoy(0:nmax)
     double precision v2(0:nmax), vintold(0:nmax)
     double precision E  (0:nmax),  Eold (0:nmax), E2(0:nmax)
     double precision phi (0:nmax),  ni (0:nmax)
     double precision nilisse (0:nmax),  nelisse (0:nmax)
     double precision phiold (0:nmax), phi2old(0:nmax)
     double precision phi1m(0:nmax), phi2m(0:nmax), phi3m(0:nmax)
     double precision ne (0:nmax),  rho (0:nmax)
     double precision nhot (0:nmax),  ncold (0:nmax)
     double precision nhotold (0:nmax),  ncoldold (0:nmax)
     double precision gradnh (0:nmax),  gradnc (0:nmax)
     double precision ghold (0:nmax),  gcold (0:nmax)
     double precision pe (0:nmax)
     double precision nss (0:nmax)
     double precision dWhot (0:nmax), SWhot(0:nmax)
     double precision phi2(0:nmax), xt2(0:nmax), Enhot(0:nmax), csc(0:nmax)
     double precision En(0:nmax), Encold(0:nmax)
     
      double precision vmoy (1:nmax), Emoy (1:nmax)
     double precision dndv (1:nmax), dndE (1:nmax)

      double precision dx0  (1:nmax),  dxt (1:nmax), dxtold(1:nmax)
     double precision ni1s2 (1:nmax), E1s2 (1:nmax), E1s2old(1:nmax)
     double precision  niSS (0:nmax), qiSS(0:nmax), dxt2(1:nmax)
      
      double precision a (1:nmax), b(0:nmax), c(0:nmax-1), f(0:nmax)
     double precision bx(0:nmax), fx(0:nmax)

       double precision ve (0:njmax), fe0(0:njmax), fe(0:njmax) 
     double precision  Te(0:njmax), Te_0(0:njmax)
     double precision  veold(0:njmax), ue(0:njmax)
     double precision   ve1m(0:njmax), ve2m(0:njmax),ve3m(0:njmax)
     double precision   Inv0(0:njmax)

       double precision fe0_cold(0:njmax), fe_cold(0:njmax) 
     double precision  Te_cold(0:njmax)

       double precision fe0_hot(0:njmax), fe_hot(0:njmax)
     double precision  Te_hot(0:njmax)
     double precision  dve2s2_0(0:njmax)

       double precision dve, nephi, gnephi, pephi2, pephi, gpephi, H
     double precision period, Spe, Spe0, deltaE, poly, disc, Spe1  

       double precision nehot, necold
     
      logical lfini,nb_cons,En_cons,EOS,VTS,laststep,break,bitemp
      
      character*10 profil
     
      external Thot,Tcold,erfc

      namelist/par/ncell,nvide,nlisse,ve_max,Nbe,prog,itmax,tmax ,iter0,iter1,fiter2,fiter3,dti, &
       n0cold,n0hot,Thmax,Tcmax,lfini,lmax,nb_cons,&
       En_cons,EOS,T_MeV,LSS,nLSS,nu,profil,VTS,djout
      
! c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
      	
! *               1 - demarrage et conditions initiales
	
! *         1.1 lecture des donnees
      
      open (13,file='boost72.par',status='unknown',form='formatted')
      read (13,par)
      close(13)

      if(prog.lt.0.9d0)then
         prog=10.d0**(-5.d0/ncell)
      endif
      
      write(*,par)
      
      if(profil.ne.'sack'.and.profil.ne.'step'.and.profil.ne.'expo'.and.profil.ne.'gaus')then
        write (*,*) 'c''est quoi le profil?'
        stop
      endif
      
      if((.not.lfini).and.(nb_cons.or.En_cons))then
        nb_cons=.false.
        En_cons=.false.
        write(*,*)  'lfini=.false. et donc En_cons=.false. et nb_cons=.false.'
      endif
      
      if((.not.nb_cons).and.En_cons) then
        nb_cons=.true.
        write(*,*) 'En_cons=.true.  et donc  nb_cons=.true.'
      endif
      
      if(nb_cons.and.(n0cold.lt.1.d-04))then
        nb_cons=.false.
        write(*,*) 'n0cold.lt.1.d-04  et donc  nb_cons=.false.'
      endif
      
      if((.not.lfini).and.VTS) then
        VTS=.false.
        write(*,*) 'lfini=.false. et donc VTT=.false.'
      endif
      
      ntotal=ncell+nvide
      
      if(ntotal.gt.nmax)then
        write(*,*) 'ntotal=ncell+nvide est plus grand que nmax'
        stop
      endif

	omegadt=sqrt(n0cold+n0hot)*dti
	if(omegadt.gt.2.) then
	  write(*,*)'pas de temps reduit pour cause de stabilite'
	  dti=2./sqrt(n0cold+n0hot)
	  write(*,*)'dti=',dti
	endif

	bitemp=.false.
	if(Tcmax.GT.0.and.n0cold.gt.0.) bitemp=.true.
		
	nstep=(tmax/dti+0.5d0)
	nsort=1
	if(nstep.gt.2000)nsort=nstep/2000
	laststep=.false.
	
	iline=1
	iline_vide=1
	if(ncell.gt.3999)iline=ncell/2000
	if(nvide.gt.799)iline_vide=nvide/400

	
! *         1.2 construction du maillage
      
! *     1.2.1 determination de la longueur du plasma

      if(.not.lfini) then
        time=0.d0
        dt=dti
        Th=Thmax
        Tc=Tcmax
        cs2=(n0cold+n0hot)/(n0cold/Tc+n0hot/Th)
        cs=sqrt(cs2)
        if(profil.eq.'sack')length=(10.d0*LSS+25.d0)*cs
        if(profil.eq.'expo')length=(LSS+25.d0)*cs
        if(profil.eq.'step')length=25.d0*cs
         
        if(itmax.eq.1) goto 91
	
        Th=Thmax*Thot(dble(0.))
        Tc=Tcmax*Tcold(dble(0.))


        do 90 itime=1,itmax
	
          if(time.ge.(tmax-1.d-5)) goto 91
          if(itime.eq.itmax) goto 91
          if((time+dt).gt.(tmax-1.d-5)) dt=tmax-time
	
          cs2old=cs2
	
          time=time+dt
	
          Th=Thmax*Thot(time)
          Tc=Tcmax*Tcold(time)
          cs2=(n0cold+n0hot)/(n0cold/Tc+n0hot/Th)
	
          cs=0.5d0*(sqrt(cs2old)+sqrt(cs2))
          length=length+cs*dt

 90     continue
 91     continue
      endif
	
      if(lfini)length=lmax
      
! *     1.2.2 maillage spatial

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
	  if(profil.eq.'expo'.and.x0(i).gt.-LSS)   niSS(i)=ni0*exp(-(x0(i)+LSS)/LSS)
	  if(profil.eq.'sack')   niSS(i)=ni0*(2.d0/pi)*DATAN(exp(-x0(i)/LSS))
	  if(profil.eq.'gaus')   niSS(i)=ni0*exp(-((x0(i)+length)/rgauss)**2)
	end do

	do i = 2,ncell
	  if(i.eq.2)then
	    dx0(2)=(1.d0+prog)*dx0(1)*niSS(0)/niSS(1)-dx0(1)
	  else
	    dx0(i)=(dx0(i-2)+dx0(i-1))*prog*niSS(i-2)/niSS(i-1)-dx0(i-1)
	  endif
	  if((profil.eq.'expo'.or.profil.eq.'sack').and.dx0(i).gt.dx0(i-1).and.dx0(i).gt.LSS/5.d0) dx0(i)=LSS/5.d0
	  if(profil.eq.'gaus'.and.dx0(i).gt.dx0(i-1).and.dx0(i).gt.lmax/5.d0) dx0(i)=lmax/5.d0
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
      
! *     1.2.3 mise en memoire de l'ordre initial des ions

	do i = 0,ncell
	   idebut(i)=i
	   irang(i)=i
	end do
	
! *         1.3 conditions initiales et estimation du
! *             potentiel phi initial
	
        break=.false. !pas de deferlement au debut

	!Th=Thmax*Thot(0.)
	!Tc=Tcmax*Tcold(0.)

        Th = Thmax * Thot(dble(0.0d0))
        Tc = Tcmax * Tcold(dble(0.0d0))



      if(Th.eq.0.d0)then
         write(*,*) 'la temperature chaude ne doit pas s''annuler'
         stop
      endif

	phi(ncell)=(n0cold*Tc+n0hot*Th)/(n0cold+n0hot)
	E(ncell)=sqrt(2.*(n0cold*Tc*exp(-phi(ncell)/Tc)+n0hot*Th*exp(-phi(ncell)/Th)))
      alpha=E(ncell)/phi(ncell)

	do i = 0,ncell
	   xt(i)=x0(i)
	   v(i)=0.d0
           vemoy(i)=0.d0
	   phi(i)=phi(ncell)*exp(alpha*(x0(i)-x0(ncell)))
	end do
	
	do i =ncell+1,ntotal
	   v(i)=0.d0
           vemoy(i)=0.d0
           v2(i)=0.d0
	   ni(i)=1.d-12
	enddo
	
	E(0)=0.d0
	gradnh(0)=0.d0
	gradnc(0)=0.d0
	time=0.d0
	dt=dti
	nhm2=n0hot
	nhm1=n0hot
	Whot1=0.d0
	Whot2=0.d0
	Wcold=0.d0
	idico=0


! *     1.3.1 initialisation de la fonction de distribution des electrons
! *             et des temperatures locales

	xvide= ve_max**2/sqrt(8.d0)+0.3d0
        ve_max=ve_max*sqrt(Thmax)
	fe_min=n0hot/sqrt(Thmax)*exp(-ve_max**2/Thmax/2.d0)

	if(Thmax.eq.Tcmax.or.(.not.bitemp)) then
	   prog_v=1.d0
	   dve=ve_max/(Nbe-1)
	else
	   prog_v=(Thmax/Tcmax)**(0.5d0/(Nbe-2))
	   dve=ve_max*(prog_v-1.d0)/(prog_v**(Nbe-1)-1.d0)
        endif

	jcoldmax=Nbe
       do j=0,Nbe-1
	  if(j.eq.0)then
	    ve(0)=0.d0
	  else
	    ve(j)=ve(j-1)+dve
	    dve=dve*prog_v
	  endif
          fe0_hot(j)=n0hot/sqrt(Thmax)*exp(-(ve(j)**2)/Thmax/2.d0)
	  if(bitemp) then
	    fe0_cold(j)=n0cold/sqrt(Tcmax)*exp(-(ve(j)**2)/Tcmax/2.d0)
          else
            fe0_cold(j)=0.d0
	  endif
	  if(fe0_cold(j).lt.fe_min.and.j.lt.jcoldmax)then
	    jcoldmax=j
	  endif
          fe0(j)=fe0_hot(j)+fe0_cold(j)
	  fe_hot(j)=fe0_hot(j)
	  fe_cold(j)=fe0_cold(j)
	  fe(j)=fe0(j)
        enddo

        do j=0,Nbe-2
	  dve2s2=0.5d0*(ve(j+1)**2-ve(j)**2)
	  dve2s2_0(j)=dve2s2
          Te(j)=-dve2s2/(log(fe(j+1))-log(fe(j)))
          Te_0(j)=Te(j)
          Te_hot(j)= Thmax
          Te_cold(j)= Tcmax
        enddo 
        Te(Nbe-1)=Thmax
        Te_hot(Nbe-1)=Thmax
        Te_cold(Nbe-1)=Tcmax

        do j=0,Nbe-1
          write(*,*) ve(j)
          write(*,*) fe0(j), fe0_hot(j), fe0_cold(j)
          write(*,*) Te(j), Te_hot(j), Te_cold(j)
          write(*,*) 
        enddo 

        do j=0,Nbe-1
           ve2m(j)=ve(j)
           ve1m(j)=ve(j)
        enddo     
	
	! open(9,file='conservation',status='unknown')
	! write(9,'(37a)')'time',char(9), 'nti',char(9),'nthot', char(9), 'ntcold', char(9), 'nte',&
   !   *   , char(9), 'n0hot', char(9), 'n0cold'
   !   *   , char(9), 'En_ion', char(9), 'Whot1'
   !   *   , char(9), 'Whot2', char(9), 'Whot', char(9), 'Wcold'
   !   *   , char(9), 'temp', char(9), 'Tc', char(9), 'En_elec'
   !   *   , char(9), 'En_delta'
   !   *   , char(9), 'vmax', char(9), 'vfinal'

! 	open(10,file='historique',status='unknown')
! 	write(10,'(29a)')'time',char(9), 'xi',char(9)
!      *   ,'v(ivmax)', char(9), 'Energy_max', char(9)
!      *   ,'E(ivmax)'
!      *   , char(9), 'ne(ivmax)', char(9),'ni(ivmax)'
!      *   , char(9),'ni(0)'
!      *   , char(9),'nhot(0)',char(9),'ncold(0)'
!      *   , char(9),'lDebye'
! c     *	  , char(9),'lgrade', char(9),'lgradi'
! c     *   , char(9),'ivmax',char(9),'idebut(ivmax)'
	
	! open(14,file='histobis',status='unknown')
	! write(14,'(9a)')'time',char(9), 't/R0',char(9), 'xi',char(9)
   !   *   ,'R/R0', char(9), 'Eq42'

	! open(15,file='cs',status='unknown')
	! write(15,'(11a)')'time',char(9), 'cs(0)',char(9)
   !   *   ,'csmoyen', char(9), 'x_rar', char(9),'Te(x=0)', char(9)
   !   *   ,'vte(x=0)'

! *	open(800,file='timejHve',status='unknown')
! *	write(800,'(13a)')'time',char(9), 'j',char(9)
! *     *   ,'H', char(9), 've(j)', char(9)
! *     *   ,'period'
! *     *   , char(9), '1/phimoy', char(9),'deltaE'

!         open(100,file='Emax_v_slop_temp',status='unknown')
!         write(100,'(17a)') 'time',char(9),'Emax',char(9),'vmax'
!      *      ,char(9),'vfinal',char(9),'vfinal2',char(9),'slop'
!      *      ,char(9),'temp'
!      *      ,char(9),'temphot'
!      *      ,char(9),'tempcold'

!         open(250,file='fe(t)',status='unknown')
!       write(250,'(9a)') 'time',char(9), 've(j)'
!      *       ,char(9),'fe0(j)', char(9),'E(j)'
!      *       ,char(9),'Te(j)'

	
! *               2 - boucle temporelle
	
        do 100 itime=1,itmax
	
	if(itime.le.3) then
	   iter=iter0
	else
	   iter=iter1
	endif
                 
! *         2.1 determination ou estimation des temperatures et de n0hot
! *             et calcul de la nouvelle densite ionique
		
! *	2.1.2 - estimation de la modification de n0hot liee a la 
! *             conservation du nombre d'electrons chauds

	if(nb_cons)then
	   nhm3=nhm2
	   nhm2=nhm1
	   nhm1=n0hot
	   n0hot=3.d0*nhm1-3.d0*nhm2+nhm3
	endif
	
! *****************************
! c	if(itime.eq.(itime/nsort)*nsort) then
	   write(*,*)
	   write(*,*)time
! c	endif
! *	write(*,*)Th, Tc, n0hot
! *****************************

! *	2.1.4 - mise en memoire des anciennes valeurs de la densite
! *             electronique, du potentiel, du champ electrique,
! *             des gradients de densite electronique,
! *             et de la position dans la gaine vide d'ions

      if(itime.gt.1) then
         do i=0,ntotal
            nhotold(i)=nhot(i)
            ncoldold(i)=ncold(i)
            phiold(i)=phi(i)
            Eold(i)=E(i)
            E1s2old(i+1)=E1s2(i+1)
            ghold(i)=gradnh(i)
            gcold(i)=gradnc(i)
         enddo
         do i=ncell+1,ntotal
            xtold(i)=xt(i)
         enddo
         do i=1,ntotal
            dxtold(i)=dxt(i)
         enddo
         Whot1old=Whot1
         Whot2old=Whot2
         Wcoldold=Wcold

         do i=1,Nbe
            veold(i)=ve(i)
         enddo

         do i=0,ncell
            phi2old(i)=phiold(i)
         enddo
         do i=ncell+1,ntotal-1
            phi2old(i)=phiold(i+1)
         enddo   
   
      endif

! *	2.1.4 bis - estimation de la modification des vej 

        do j=0,Nbe-1
           ve3m(j)=ve2m(j)
           ve2m(j)=ve1m(j)           
           ve1m(j)=ve(j)
! c           ve(j)=3.d0*ve1m(j)-3.d0*ve2m(j)+ve3m(j)
           ve(j)=2.d0*ve1m(j)-1.d0*ve2m(j)
        enddo 
	
! *     2.1.5 calcul de la nouvelle densite ionique
	
	do i = 1,ncell
	   dxt(i)=xt(i)-xt(i-1)
! *	   ni1s2(i)=(qiSS(i-1)+qiSS(i))/dxt(i)/2.
	end do
	
	ni(0)=qiSS(0)/dxt(1)
	
! *	do i = 1,ncell-1
! *	   ni(i)=(dxt(i+1)*ni1s2(i)+dxt(i)*ni1s2(i+1))
! *     *         /(dxt(i)+dxt(i+1))
! *	end do
	
	do i = 1,ncell-1
	   ni(i)=2.d0*qiSS(i)/(dxt(i)+dxt(i+1))
	end do
	
! *	ni(ncell)=1.5d0*ni1s2(ncell)-0.5d0*ni1s2(ncell-1)
! * 	ni(ncell)=( (dxt(ncell-1)+2.d0*dxt(ncell))*ni1s2(ncell)
! *     *             -dxt(ncell)*ni1s2(ncell-1) )
! *     *		 /(dxt(ncell-1)+dxt(ncell))
! *      if(ni(ncell).lt.0.d0)ni(ncell)=1.d-10
	
 	ni(ncell)=2.d0*qiSS(ncell)/dxt(ncell)/(1.d0+2.d0*dxt(ncell)/dxt(ncell-1)-dxt(ncell-1)/dxt(ncell-2))

	do i = 1,ncell
	   ni1s2(i)=(ni(i-1)+ni(i))/2.d0
	end do
        
        
! *	   2.2 - iteration pour le calcul de phi, Th, Tc

	do 10 iphi=1,iter

! *-------------------------------------------------------------------
! *	2.2.1 - construction des tableaux a,b,c, et f
! *-------------------------------------------------------------------
! c----- initialisation des fj et des Tj c--------------------------------------

! c       reajustement de fe au premier pas de temps (remarque:
! c          phi(0) est tres proche de zero)


        if (itime.eq.1) then
          do j=0,Nbe-1
            fe_hot(j)=fe_hot(j)*exp(-phi(0)/Thmax)
	    if(bitemp)then
              fe_cold(j)=fe_cold(j)*exp(-phi(0)/Tcmax)
	    endif
            fe(j)=fe_hot(j)+fe_cold(j)
          enddo
	  do i=1,ncell
	    phi(i)=phi(i)-phi(0)
	  enddo
	  phi(0)=0.d0
          do j=0,Nbe-2
            Te_0(j)=-dve2s2_0(j)/(log(fe(j+1))-log(fe(j)))
          enddo 
        endif  

        do j=0,Nbe-2
	  dve2s2=0.5d0*(ve(j+1)**2-ve(j)**2)
          cooling=dve2s2/dve2s2_0(j)
          Te(j)=Te_0(j)*cooling
	  if(iphi.eq.iter) then
            Te_hot(j)=Thmax*cooling
	    if(bitemp)then
              Te_cold(j)=Tcmax*cooling
	    endif
	  endif
        enddo 

	if(itime.eq.1) phi0=phi(0)
         
! c---------------------------------------------------
! c------ Calcul i quelconque c---------------------------------------------------
	
	do i=0,ncell

        nephi=0.d0
        gnephi=0.d0
        pephi2=0.d0
        j=Nbe-2
        k=0

        do while ((ve(j)**2-2.d0*(phi(i)-phi(0))).ge.0.d0.and.(j.gt.0))

        ue(j+1)=sqrt(ve(j+1)**2-2*(phi(i)-phi(0)))
        ue(j)=sqrt(ve(j)**2-2*(phi(i)-phi(0)))
    
        nephi=nephi+fe(j)*exp((ve(j)**2/2+phi0-phi(i))/Te(j))*(-erfc(ue(j+1)/(sqrt(2.d0*Te(j)))) &
         +erfc(ue(j)/(sqrt(2.d0*Te(j))))) *sqrt(Te(j)) 

        gnephi=gnephi+fe(j)*exp((ve(j)**2/2+phi0-phi(i))/Te(j)) &
     *sqrt(Te(j))*((-erfc(ue(j+1)/(sqrt(2.d0*Te(j)))) &
     +erfc(ue(j)/(sqrt(2.d0*Te(j))))) *(-1.d0/Te(j)) &
     -(exp(-ue(j+1)**2/(2*Te(j)))/(ue(j+1)/sqrt(2.d0)) &
     -exp(-ue(j)**2/(2*Te(j)))/(ue(j)/sqrt(2.d0)))/sqrt(pi*Te(j)))
  
	if(i.eq.ncell)then
          pephi2=pephi2+fe(j)*exp((ve(j)**2/2+phi0)/Te(j)) &
    *sqrt(Te(j))*exp(-phi(ncell)/Te(j))*Te(j)*((-erfc(ue(j+1)/(sqrt(2.d0*Te(j)))) &
     +erfc(ue(j)/(sqrt(2.d0*Te(j))))) -(exp(-ue(j+1)**2/(2.d0*Te(j)))*ue(j+1) &
     -exp(-ue(j)**2/(2.d0*Te(j)))*ue(j))*sqrt(2.d0/(pi*Te(j))))
	endif
  
        j=j-1
        enddo
        
        k=j+1
        if (k.ge.1) then  
        nephi=nephi+fe(k)*exp((ve(k)**2/2+phi0-phi(i))/Te(k-1)) *(1.d0-erfc(ue(k)/(sqrt(2.d0*Te(k-1))))) *sqrt(Te(k-1))  
 
        gnephi=gnephi+fe(k)*exp((ve(k)**2/2+phi0-phi(i))/Te(k-1)) &
      *sqrt(Te(k-1))*( (-1.d0/Te(k-1))*(1.d0-erfc(ue(k)/sqrt(2.d0*Te(k-1)))) &
      -(exp(-ue(k)**2/(2*Te(k-1)))/(ue(k)/sqrt(2.d0)))   /sqrt(pi*Te(k-1)))

	if(i.eq.ncell)then
          pephi2=pephi2+fe(k)*exp((ve(k)**2/2+phi0)/Te(k-1)) &
        *sqrt(Te(k-1))*exp(-phi(ncell)/Te(k-1))*Te(k-1)*( &
         1.d0-erfc(ue(k)/sqrt(2.d0*Te(k-1))) &
         -exp(-ue(k)**2/(2.d0*Te(k-1)))*ue(k) &
        *sqrt(2.d0/(pi*Te(k-1)))) 
	endif

        endif

	j=Nbe-1

      nephi=nephi+fe(j)*exp((ve(j)**2/2+phi0-phi(i))/Te(j)) *erfc(ue(j)/sqrt(2.d0*Te(j)))*sqrt(Te(j)) 

        gnephi=gnephi+fe(j)*exp((ve(j)**2/2+phi0-phi(i))/Te(j)) &
     *sqrt(Te(j))*( &
       erfc(ue(j)/sqrt(2.d0*Te(j))) &
         *(-1.d0/Te(j)) & 
          -(0.d0 &
    -exp(-ue(j)**2/(2*Te(j)))/(ue(j)/sqrt(2.d0)))/sqrt(pi*Te(j)))
  
	if(i.eq.ncell)then
          pephi2=pephi2+fe(j)*exp((ve(j)**2/2+phi0)/Te(j))&
              *sqrt(Te(j))*exp(-phi(ncell)/Te(j))*Te(j)*(& 
                 (erfc(ue(j)/sqrt(2.d0*Te(j)))) &
                     -(0.d0 & 
                        -exp(-ue(j)**2/(2.d0*Te(j)))*ue(j))*sqrt(2.d0/(pi*Te(j))))  
          pephi=sqrt(2.d0*pephi2)
          gpephi=-nephi/pephi
	endif

	if(i.eq.0) then
	   b(0)=-2.d0/dxt(1)**2+gnephi
	   c(0)=2.d0/dxt(1)**2
	   f(0)=ni(0)-nephi+phi(0)*gnephi
	elseif(i.gt.0.and.i.lt.ncell) then
  	   a(i)=2.d0/dxt(i)/(dxt(i)+dxt(i+1))
	   b(i)=-2.d0/(dxt(i)*dxt(i+1))+gnephi
	   c(i)=2.d0/dxt(i+1)/(dxt(i)+dxt(i+1))
	   f(i)=ni(i)-nephi+phi(i)*gnephi  
	elseif(i.eq.ncell) then
	   a(ncell)=2.d0/dxt(ncell)**2
	   b(ncell)=-a(ncell)
   !   *         +gpephi*2.d0/dxt(ncell)
   !   *         +gnephi
           f(ncell)=ni(ncell)-(nephi-phi(ncell)*gnephi)
   !   *         -(pephi-phi(ncell)*gpephi)*2.d0/dxt(ncell)
	endif

	enddo

! c---------------------------------------------------------------------
       
! *	2.2.2 - inversion de la matrice tridiagonale
 
	call gauss(nmax,ncell,a,b,c,f,bx,fx,phi)

! c---------------------------------------------------------------------
	
      if(iphi.ne.iter.and.itime.eq.1) goto 10
      if(iphi.ne.iter.and.(.not.nb_cons).and.(.not.En_cons)) goto 10
     
! *     2.2.3 - calcul de la densite electronique dans le plasma

	do i=0,ncell

           nephi=0.d0
           nehot=0.d0
           necold=0.d0
           pephi=0.d0
           j=Nbe-2
           k=0

        do while ((ve(j)**2-2.d0*(phi(i)-phi(0))).gt.0.d0.and.(j.gt.0))
 
        ue(j+1)=sqrt(ve(j+1)**2-2.d0*(phi(i)-phi(0)))
        ue(j)=sqrt(ve(j)**2-2.d0*(phi(i)-phi(0)))
    
        nephi=nephi+fe(j)*exp((ve(j)**2/2+phi0-phi(i))/Te(j)) &
          *(-erfc(ue(j+1)/sqrt(2.d0*Te(j)))+erfc(ue(j)/sqrt(2.d0*Te(j)))) *sqrt(Te(j))

        nehot=nehot+fe_hot(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_hot(j))& 
         *(-erfc(ue(j+1)/sqrt(2.d0*Te_hot(j))) & 
          +erfc(ue(j)/sqrt(2.d0*Te_hot(j))))*sqrt(Te_hot(j))

	if(bitemp.and.j.lt.jcoldmax)then
        necold=necold +fe_cold(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_cold(j)) &
       *(-erfc(ue(j+1)/sqrt(2.d0*Te_cold(j))) +erfc(ue(j)/sqrt(2.d0*Te_cold(j))))*sqrt(Te_cold(j))
	endif

        j=j-1  
        enddo

        k=j+1

        if (k.ge.1) then   

        nephi=nephi+fe(k)*exp((ve(k)**2/2+phi0-phi(i))/Te(k-1)) &
     *(1.d0-erfc(ue(k)/sqrt(2.d0*Te(k-1)))) *sqrt(Te(k-1)) 

        nehot=nehot+fe_hot(k)*exp((ve(k)**2/2+phi0-phi(i))/Te_hot(k-1)) &
     *(1.d0-erfc(ue(k)/sqrt(2.d0*Te_hot(k-1)))) *sqrt(Te_hot(k-1)) 

	if(bitemp.and.k.lt.jcoldmax)then
        necold=necold +fe_cold(k)*exp((ve(k)**2/2+phi0-phi(i))/Te_cold(k-1)) &
     *(1.d0-erfc(ue(k)/sqrt(2.d0*Te_cold(k-1))))*sqrt(Te_cold(k-1))
	endif

        endif

	j=Nbe-1

        nephi=nephi+fe(j)*exp((ve(j)**2/2+phi0-phi(i))/Te(j)) &
      *erfc(ue(j)/sqrt(2.d0*Te(j))) *sqrt(Te(j))

        nehot=nehot+fe_hot(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_hot(j)) &
          *erfc(ue(j)/sqrt(2.d0*Te_hot(j))) *sqrt(Te_hot(j))

	if(bitemp.and.j.lt.jcoldmax)then
        necold=necold &
     +fe_cold(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_cold(j)) &
     *erfc(ue(j)/sqrt(2.d0*Te_cold(j))) *sqrt(Te_cold(j))
	endif

        ne(i)=nephi
        nhot(i)= nehot
        ncold(i)= necold
        rho(i)=ni(i)-ne(i)
	enddo
         
! *	2.2.4 - calcul de la pression electronique
! *		  et des gradients de densite
! *		  electronique et ionique au bord
	
	do i=ncell,ncell
           pephi=0.d0
           j=Nbe-2
           k=0

        do while ((ve(j)**2-2.d0*(phi(i)-phi(0))).ge.0.d0.and.(j.gt.0))
 
        ue(j+1)=sqrt(ve(j+1)**2-2.d0*(phi(i)-phi(0)))
        ue(j)=sqrt(ve(j)**2-2.d0*(phi(i)-phi(0)))


        pephi=pephi+fe(j)*exp((ve(j)**2/2+phi0)/Te(j)) &
          *sqrt(Te(j))*exp(-phi(ncell)/Te(j))*Te(j)*( &
            (-erfc(ue(j+1)/(sqrt(2.d0*Te(j)))) &
               +erfc(ue(j)/(sqrt(2.d0*Te(j))))) & 
                -(exp(-ue(j+1)**2/(2.d0*Te(j)))*ue(j+1) & 
                 -exp(-ue(j)**2/(2.d0*Te(j)))*ue(j))*sqrt(2.d0/(pi*Te(j))))
             
        j=j-1  
           enddo

        k=j+1

        if (k.ge.1) then 

        pephi=pephi+fe(k)*exp((ve(k)**2/2+phi0)/Te(k-1)) & 
         *sqrt(Te(k-1))*exp(-phi(ncell)/Te(k-1))*Te(k-1)* &
           (1.d0-erfc(ue(k)/sqrt(2.d0*Te(k-1))) &
             -exp(-ue(k)**2/(2.d0*Te(k-1)))*ue(k) &
               *sqrt(2.d0/(pi*Te(k-1)))) 
  
        endif

	j=Nbe-1

        pephi=pephi+fe(j)*exp((ve(j)**2/2+phi0)/Te(j)) &
          *sqrt(Te(j))*exp(-phi(ncell)/Te(j))*Te(j)* &
            (erfc(ue(j)/sqrt(2.d0*Te(j))) &
              -(0.d0 -exp(-ue(j)**2/(2.d0*Te(j)))*ue(j))*sqrt(2.d0/(pi*Te(j))))
             
        pe(i)=pephi
      
	enddo
   	
! *     2.2.5 calcul de la densite electronique dans la gaine vide d'ions

! c--------- calcul du potentiel dans la gaine ----------------------------------
    
      kvide=1.d0/sqrt(2.d0)
      uphi=0.d0
 
      do i=ncell+1,ntotal
        
         Spe=0.d0
         Spe0=0.d0
         Spe1=0.d0

         if(i.eq.(ncell+1)) then
            dxt(ncell+1)=0.
         else
! c            dxt(i)=0.0018*sqrt(1.d0/nhot(i-1))
            dxt(i)=xvide*sqrt(Thmax/ne(i-1))/nvide
         endif

         xt(i)=xt(i-1)+dxt(i)

       j=Nbe-2
       k=0

       do while ((ve(j)**2-2.d0*(phi(i-1)-phi(0))).gt.0.d0.and.(j.gt.0))
 
        ue(j+1)=sqrt(ve(j+1)**2-2.d0*(phi(i-1)-phi(0)))
        ue(j)=sqrt(ve(j)**2-2.d0*(phi(i-1)-phi(0)))

        Spe=Spe+fe(j)*exp((ve(j)**2/2+phi0)/Te(j)) &
          *sqrt(Te(j))*exp(-phi(i-1)/Te(j))*Te(j)*(& 
           (-erfc(ue(j+1)/(sqrt(2.d0*Te(j)))) & 
             +erfc(ue(j)/(sqrt(2.d0*Te(j))))) &
               -(exp(-ue(j+1)**2/(2.d0*Te(j)))*ue(j+1) &
                 -exp(-ue(j)**2/(2.d0*Te(j)))*ue(j))*sqrt(2.d0/(pi*Te(j))))

        j=j-1
  
        enddo

        k=j+1

        if (k.ge.1) then 

        Spe=Spe+fe(k)*exp((ve(k)**2/2+phi0)/Te(k-1)) &  
        *sqrt(Te(k-1))*exp(-phi(i-1)/Te(k-1))*Te(k-1)* & 
         (1.d0-erfc(ue(k)/sqrt(2.d0*Te(k-1))) &
           -exp(-ue(k)**2/(2.d0*Te(k-1)))*ue(k) &
             *sqrt(2.d0/(pi*Te(k-1))))

        endif

	j=Nbe-1

        Spe=Spe+fe(j)*exp((ve(j)**2/2+phi0)/Te(j)) & 
         *sqrt(Te(j))*exp(-phi(i-1)/Te(j))*Te(j)* & 
          (erfc(ue(j)/sqrt(2.d0*Te(j))) & 
           -(0.d0 -exp(-ue(j)**2/(2.d0*Te(j)))*ue(j))*sqrt(2.d0/(pi*Te(j))))

        Spe=Spe*exp(phi(i-1))    !+Spe0

        uphi=uphi+kvide*sqrt(Spe)*dxt(i)   
         phi(i)=phi(ncell)
   !   *   +2.d0*log(1.d0+uphi*exp(-0.5d0*phi(ncell)))

! c---------calcul de la densite et du champ dans la gaine

        nephi=0.d0
        nehot=0.d0
        necold=0.d0
        pephi=0.d0

        j=Nbe-2
        k=0

        do while ((ve(j)**2-2.d0*(phi(i)-phi(0))).gt.0.d0.and.(j.gt.0))
 
        ue(j+1)=sqrt(ve(j+1)**2-2.d0*(phi(i)-phi(0)))
        ue(j)=sqrt(ve(j)**2-2.d0*(phi(i)-phi(0)))
    
        nephi=nephi+fe(j)*exp((ve(j)**2/2+phi0-phi(i))/Te(j)) & 
         *(-erfc(ue(j+1)/sqrt(2.d0*Te(j)))+erfc(ue(j)/sqrt(2.d0*Te(j)))) *sqrt(Te(j))

        nehot=nehot+fe_hot(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_hot(j)) &
       *(-erfc(ue(j+1)/sqrt(2.d0*Te_hot(j))) &
      +erfc(ue(j)/sqrt(2.d0*Te_hot(j))))*sqrt(Te_hot(j))

	if(bitemp.and.j.lt.jcoldmax)then
        necold=necold  & 
     +fe_cold(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_cold(j))  & 
       *(-erfc(ue(j+1)/sqrt(2.d0* Te_cold(j)))  & 
       +erfc(ue(j)/sqrt(2.d0* Te_cold(j))))*sqrt(Te_cold(j))
	endif

        pephi=pephi+fe(j)*exp((ve(j)**2/2+phi0)/Te(j)) &
       *sqrt(Te(j))*exp(-phi(i)/Te(j))*Te(j)*( & 
       (-erfc(ue(j+1)/sqrt(2.d0*Te(j))) & 
        +erfc(ue(j)/sqrt(2.d0*Te(j)))) & 
       -(exp(-ue(j+1)**2/(2.d0*Te(j)))*ue(j+1) & 
       -exp(-ue(j)**2/(2.d0*Te(j)))*ue(j))*sqrt(2.d0/(pi*Te(j))))

        j=j-1

        enddo

        k=j+1

        if (k.ge.1) then   

        nephi=nephi+fe(k)*exp((ve(k)**2/2+phi0-phi(i))/Te(k-1)) & 
      *(1.d0-erfc(ue(k)/sqrt(2.d0*Te(k-1))))*sqrt(Te(k-1)) 

        nehot=nehot+fe_hot(k)*exp((ve(k)**2/2+phi0-phi(i))/Te_hot(k-1))  & 
      *(1.d0-erfc(ue(k)/sqrt(2.d0*Te_hot(k-1))))*sqrt(Te_hot(k-1)) 

	if(bitemp.and.k.lt.jcoldmax)then
        necold=necold  & 
       +fe_cold(k)*exp((ve(k)**2/2+phi0-phi(i))/Te_cold(k-1))  & 
      *(1.d0-erfc(ue(k)/sqrt(2.d0*Te_cold(k-1))))  & 
       *sqrt(Te_cold(k-1))
	endif

        pephi=pephi+fe(k)*exp((ve(k)**2/2+phi0)/Te(k-1)) & 
       *sqrt(Te(k-1))*exp(-phi(i)/Te(k-1))*Te(k-1)* & 
       (1.d0-erfc(ue(k)/sqrt(2.d0*Te(k-1))) & 
       -exp(-ue(k)**2/(2.d0*Te(k-1)))*ue(k) & 
       *sqrt(2.d0/(pi*Te(k-1)))) 
      
        endif

	j=Nbe-1

        nephi=nephi+fe(j)*exp((ve(j)**2/2+phi0-phi(i))/Te(j)) & 
      *erfc(ue(j)/sqrt(2.d0*Te(j))) *sqrt(Te(j))

        nehot=nehot+fe_hot(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_hot(j))  & 
     *erfc(ue(j)/sqrt(2.d0*Te_hot(j)))*sqrt(Te_hot(j))

	if(bitemp.and.j.lt.jcoldmax)then
        necold=necold +fe_cold(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_cold(j)) &
          *erfc(ue(j)/sqrt(2.d0* Te_cold(j)))*sqrt(Te_cold(j))
	endif

        pephi=pephi+fe(j)*exp((ve(j)**2/2+phi0)/Te(j)) &
      *sqrt(Te(j))*exp(-phi(i)/Te(j))*Te(j)* &
       (erfc(ue(j)/sqrt(2.d0*Te(j))) &
       -(0.d0 &
       -exp(-ue(j)**2/(2.d0*Te(j)))*ue(j))*sqrt(2.d0/(pi*Te(j))))

        ne(i)=nephi
        nhot(i)= nehot
        ncold(i)= necold
        rho(i)=-ne(i)
        pe(i)=pephi
        E(i)=sqrt(2.d0*pe(i))
        
        enddo

! *     2.2.6 calcul du nombre total d'electrons et d'ions

      nti=0.d0
      do i=1,ncell
         nti=nti+ni1s2(i)*dxt(i)
      enddo

      nte=0.d0
      nthot=0.d0
      ntcold=0.d0
      do i=1,ntotal
          nte = nte +0.5d0*(ne(i-1)+ne(i))*dxt(i)
          nthot = nthot +0.5d0*(nhot(i-1)+ nhot(i))*dxt(i)
          ntcold = ntcold +0.5d0*(ncold(i-1)+ ncold(i))*dxt(i)
      enddo
      nte = nte +E(ntotal)
      if(Thmax.eq.Tcmax)then
        nthot = nthot +E(ntotal)*n0hot/(n0hot+n0cold)
        ntcold = ntcold +E(ntotal)*n0cold/(n0hot+n0cold)
      else
        nthot = nthot +E(ntotal)
      endif

      if (itime.eq.1) then
         nombre= nte
      endif   
      
      nte1=0.d0
      do i=1,ncell
         nte1=nte1+0.5d0*(ne(i-1)+ne(i))*dxt(i)
      enddo

      nte2=0.d0
      do i=ncell+1,ntotal 
         nte2=nte2+0.5d0*(ne(i-1)+ne(i))*dxt(i)
      enddo
      nte2=nte2+E(ntotal)

	if(iphi.eq.iter)then

      do i=0,ntotal
       ekmoy=0.d0
       ekhot=0.d0
       ekcold=0.d0
 
      j=Nbe-2

       do while ((ve(j)**2-2.d0*(phi(i)-phi(0))).gt.0.d0.and.(j.gt.0))

        ue(j+1)=sqrt(ve(j+1)**2-2.d0*(phi(i)-phi(0)))
        ue(j)=sqrt(ve(j)**2-2.d0*(phi(i)-phi(0)))

       ekmoy=ekmoy+fe(j)*exp((ve(j)**2/2+phi0-phi(i))/Te(j))*Te(j) & 
       *(sqrt(pi/2.d0)*sqrt(Te(j))* &
        (-erfc(ue(j+1)/(sqrt(2.d0*Te(j)))) &
          +erfc(ue(j)/(sqrt(2.d0*Te(j)))))- &
           (ue(j+1)*exp(-(ue(j+1)**2/(2*Te(j)))) &
             -ue(j)*exp(-(ue(j)**2/(2*Te(j))))))

       ekhot=ekhot & 
       +fe_hot(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_hot(j))*Te_hot(j) &
        *(sqrt(pi/2.d0)*sqrt(Te_hot(j))* &
         (-erfc(ue(j+1)/(sqrt(2.d0*Te_hot(j)))) &
           +erfc(ue(j)/(sqrt(2.d0*Te_hot(j)))))- &
            (ue(j+1)*exp(-(ue(j+1)**2/(2*Te_hot(j)))) &
              -ue(j)*exp(-(ue(j)**2/(2*Te_hot(j))))))

	if(bitemp.and.j.lt.jcoldmax)then
       ekcold=ekcold & 
       +fe_cold(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_cold(j))*Te_cold(j) &
        *(sqrt(pi/2.d0)*sqrt(Te_cold(j))* &
         (-erfc(ue(j+1)/(sqrt(2.d0* Te_cold(j)))) &
           +erfc(ue(j)/(sqrt(2.d0* Te_cold(j)))))- &
            (ue(j+1)*exp(-(ue(j+1)**2/(2* Te_cold(j)))) &
              -ue(j)*exp(-(ue(j)**2/(2* Te_cold(j))))))
	endif

       j=j-1

       enddo 

       k=j+1
       if (k.ge.1) then 

         ekmoy=ekmoy+fe(j)*exp((ve(j)**2/2+phi0-phi(i))/Te(j))*Te(j) &
          *(sqrt(pi/2.d0)*sqrt(Te(j))* &
          erfc(ue(j)/sqrt(2.d0*Te(j)))- &
          (0.d0-ue(j)*exp(-(ue(j)**2/(2*Te(j))))))
    
           ekhot=ekhot &
         +fe_hot(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_hot(j))*Te_hot(j) &
         *(sqrt(pi/2.d0)*sqrt(Te_hot(j))* &
         erfc(ue(j)/sqrt(2.d0*Te_hot(j)))- &
          (0.d0-ue(j)*exp(-(ue(j)**2/(2*Te_hot(j))))))
    
       if(bitemp.and.j.lt.jcoldmax)then
           ekcold=ekcold  &
         +fe_cold(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_cold(j))*Te_cold(j)  &
         *(sqrt(pi/2.d0)*sqrt(Te_cold(j))*  &
         erfc(ue(j)/sqrt(2.d0* Te_cold(j)))-  &
         (0.d0-ue(j)*exp(-(ue(j)**2/(2* Te_cold(j))))))
       endif
    
       endif

       j=Nbe-1

       ekmoy = ekmoy + fe(j)*exp((ve(j)**2/2+phi0-phi(i))/Te(j))*Te(j) &
       * (sqrt(pi/2.d0)*sqrt(Te(j)) &
       * erfc(ue(j)/sqrt(2.d0*Te(j))) - &
       (0.d0 - ue(j)*exp(-(ue(j)**2/(2*Te(j))))))
  
  ekhot = ekhot &
       + fe_hot(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_hot(j))*Te_hot(j) &
       * (sqrt(pi/2.d0)*sqrt(Te_hot(j)) &
       * erfc(ue(j)/sqrt(2.d0*Te_hot(j))) - &
       (0.d0 - ue(j)*exp(-(ue(j)**2/(2*Te_hot(j))))))
  
  if(bitemp.and.j.lt.jcoldmax) then
      ekcold = ekcold &
          + fe_cold(j)*exp((ve(j)**2/2+phi0-phi(i))/Te_cold(j))*Te_cold(j) &
          * (sqrt(pi/2.d0)*sqrt(Te_cold(j)) &
          * erfc(ue(j)/sqrt(2.d0* Te_cold(j))) - &
          (0.d0 - ue(j)*exp(-(ue(j)**2/(2* Te_cold(j)))))) 
  endif
  


        En(i)=ekmoy/sqrt(2.d0*pi)
        Enhot(i)= ekhot/sqrt(2.d0*pi)
        Encold(i)= ekcold/sqrt(2.d0*pi)
  
      enddo

      En_0=0.d0
      En_hot0=0.d0
      En_cold0=0.d0

      do i=1,ntotal
         En_0=En_0+0.5d0*(En(i-1)+En(i))*dxt(i)
         En_hot0=En_hot0+0.5d0*(Enhot(i-1)+Enhot(i))*dxt(i)
         En_cold0=En_cold0+0.5d0*(Encold(i-1)+Encold(i))*dxt(i)
      enddo   

      En_0=En_0+0.5d0*E(ntotal)*Thmax
      if(Thmax.eq.Tcmax)then
        En_hot0=En_hot0+0.5d0*E(ntotal)*n0hot/(n0hot+n0cold)
        En_cold0=En_cold0+0.5d0*E(ntotal)*n0cold/(n0hot+n0cold)
      else
        En_hot0=En_hot0+0.5d0*E(ntotal)
      endif

      temp=2.d0*En_0/nte
      temphot=2.d0*En_hot0/nthot
	if(bitemp)then
      tempcold=2.d0*En_cold0/ntcold
	else
      tempcold=0.d0
	endif

      write(*,*) 'Nb_ions =',nti,'nte1 =',nte1
      write(*,*) 'nte2 =',nte2

! *     2.2.6.1 calcul de la vitesse acoustique ionique cinetique

      do i = 0, ntotal
         cs3 = 0.d0
         j = Nbe - 2
     
         do while ((ve(j)**2 - 2.d0*(phi(i) - phi(0))) > 0.d0 .and. (j > 0))
             ue(j+1) = sqrt(ve(j+1)**2 - 2.d0*(phi(i) - phi(0)))
             ue(j) = sqrt(ve(j)**2 - 2.d0*(phi(i) - phi(0)))
             cs3 = cs3 + fe(j)*exp((ve(j)**2/2 + phi0 - phi(i))/Te(j)) &
                 *(-erfc(ue(j+1)/sqrt(2.d0*Te(j))) + erfc(ue(j)/sqrt(2.d0*Te(j)))) &
                 /sqrt(Te(j))
             j = j - 1
         end do 
     
         k = j + 1
         if (k >= 1) then 
             cs3 = cs3 + fe(k)*exp((ve(k)**2/2 + phi0 - phi(i))/Te(k-1)) &
                 *(1.d0 - erfc(ue(k)/sqrt(2.d0*Te(k-1)))) &
                 /sqrt(Te(k-1))
         end if
     
         j = Nbe - 1
         cs3 = cs3 + fe(j)*exp((ve(j)**2/2 + phi0 - phi(i))/Te(j)) &
             *erfc(ue(j)/sqrt(2.d0*Te(j))) &
             /sqrt(Te(j))
         csc(i) = cs3
     end do 
     
     cscin = sqrt(ni(0)/csc(0))

     

! *     2.2.7 ajustement de n0hot pour conserver le nombre
! *           total d'electrons chauds (alors le nombre d'electrons
! *           froids est automatiquement conserve, car le potentiel
! *           en 0 s'ajuste en consequence)

! c      if(nb_cons)then
! c         if(itime.eq.1)then
! c            nthot0=nthot
! c         else
! c            if(iphi.gt.1) n0hot=n0hot*(nthot0/nthot)
! c         endif
! c      endif
      
	if(En_cons.and.itime.eq.1)then
! c	   En_cold0=En_cold
	   write(*,*)
	   write(*,*)'En_hot0=', En_hot0
	   write(*,*)'En_cold0=', En_cold0
	   write(*,*)
	endif
	
	endif


! *     2.2.8 calcul de la vitesse d'interface des cellules vides d'ions

      if(itime.gt.1)then
         do i=ncell+1,ntotal
            vint(i)=(xt(i)-xtold(i))/dt
            vintold(i)=vint(i)
         enddo
      endif
            
! *     2.2.11 calcul du champ electrique dans toute la detente

        do i = 1,ncell
	   E1s2(i)=(phi(i)-phi(i-1))/dxt(i)
	end do
		
        do i=ncell+2,ntotal
	   E1s2(i-1)=(phi(i)-phi(i-1))/dxt(i)
	end do  

	do i = 1,ncell-1
	   E(i)=(dxt(i+1)*E1s2(i)+dxt(i)*E1s2(i+1))/(dxt(i)+dxt(i+1))
	end do
	
 	Enou=E(ncell-1)+dxt(ncell)*(ni(ncell)-ne(ncell))
        E(ncell)=sqrt(2.d0*pe(ncell))

! *     2.2.12 ajustement de la vitesse aux temps 'entiers'
	
	if(itime.gt.1) then
	   do i = 1,ncell
	      v(i)=vint(i)+dt*(3.d0*E(i)+Eold(i))/8.d0
	   end do	
	endif

! *     2.2.12.1 calcul de la vitesse moyenne electronique

	if(itime.gt.1) then

           do i=1,ncell
           vemoy(i)=(v(i)+(E(i)-Eold(i))/(ne(i)*dt))
           enddo

           do i=ncell+1,ntotal
           v2(i)=1.5d0*vint(i)-0.5d0*vintold(i)
           vemoy(i)=(v2(i)+(E(i)-Eold(i))/(ne(i)*dt))
           enddo

        endif   

! *     2.2.13 calcul de l'energie cinetique des ions

      En_ion=0.
	do i = 1,ncell-1
	   En_ion=En_ion+qiSS(i)*v(i)**2/2.d0
	end do	
	if(profil.eq.'step') En_ion=En_ion+qiSS(ncell)*v(ncell)**2/4.d0
	if(profil.ne.'step') En_ion=En_ion+qiSS(ncell)*v(ncell)**2/2.d0

! *     2.2.14 calcul de l'energie electrostatique

   En_elec=0.
	do i = 1,ntotal-1
	   En_elec=En_elec+0.25d0*(E(i)**2)*(dxt(i)+dxt(i+1))
	end do
	En_elec=En_elec+0.25d0*(E(ntotal)**2)*dxt(ntotal)
	En_elec=En_elec+E(ntotal)*Th
			
! *     2.2.15 calcul de l'energie totale

	if(itime.eq.1)then
	   En_totale=En_ion+En_elec+En_0
	endif

! *-------------------------------------------------------------------
! *  calcul de la variation d'energie des electrons dans le potentiel
! *--------------------------------------------------------------------

       if(itime.gt.1.and.iphi.gt.1) then

       xt2(0)=-lmax

      do i=0,ncell
         phi2(i)=phi(i)
         dxt2(i+1)=dxt(i+1)
         xt2(i+1)=xt2(i)+dxt2(i+1)
      enddo  

      do i=ncell+1,ntotal-1
         phi2(i)=phi(i+1)
         dxt2(i)=dxt(i+1)
         xt2(i)=xt(i+1)
      enddo 
      
      do j=1,Nbe-1

      i=0
      m=0 
      period=0.d0
      deltaE=0.d0
      Inv2=0.d0
      Hold=veold(j)**2/2.d0+phi2old(0)
      Hnew=ve(j)**2/2.d0+phi2(0)

      phimoy=0.d0
      path=0.d0

      if (iphi.eq.2) then
         H=Hold
      else
         H=Hnew
      endif 
      
! c--------------- boucle d'integration --------------------

      if(phi2(1).lt.H) then

      do while (phi2(i+1).lt.H)
         m=i
         if (i.eq.0) then

      Integ=dxt2(i+1)/sqrt(H-(phi2(i+1)+phi2(i))/2.d0)

       period=period+Integ

      deltaE=deltaE+(phi2(0)-phi2old(0) - vint(0)*0.5d0*(E1s2(1)+E1s2old(1))*dt)*Integ

      Inv2=Inv2+sqrt(H-phi2(i+1))*dxt2(i+1) 

         else

	if(i.ge.(ntotal-1)) then
	   write(*,*) 'attention i egal ou sup a ntotal-1 dans la'
	   write(*,*) 'boucle d''integration pour le calcul de la'
	   write(*,*) 'variation d''energie des electrons'
	   write(*,*) 'j=',j,'   ve(j)=',ve(j)
	   close(9)
	   close(10)
	   close(14)
! *	   close(800)
	   close(100)
	   close(250)
	   stop
	endif

        C1=xt2(i-1)-xt2(i)-(xt2(i-1)**2-xt2(i)**2)/(xt2(i)+xt2(i+1))
        D1=phi2(i)-phi2(i-1)+(phi2(i)-phi2(i+1)) & 
                *(xt2(i-1)**2-xt2(i)**2)/(xt2(i)**2-xt2(i+1)**2) 
        C2=xt2(i+2)-xt2(i+1)-(xt2(i+2)**2-xt2(i+1)**2)/(xt2(i)+xt2(i+1))
        D2=phi2(i+1)-phi2(i+2)+(phi2(i)-phi2(i+1)) &
                 *(xt2(i+2)**2-xt2(i+1)**2)/(xt2(i)**2-xt2(i+1)**2) 

            a1=-(C1*D1+C2*D2)/(C1**2+C2**2)
            a2=(phi2(i)-phi2(i+1)-a1*(xt2(i)-xt2(i+1)))/(xt2(i)**2-xt2(i+1)**2)
            a0=phi2(i)-a1*xt2(i)-a2*xt2(i)**2

         xbar=(xt2(i)+xt2(i+1))/2.d0
         phibar=a0+a1*xbar+a2*xbar**2
         Rbar=H-phibar


          if (abs(a2)*dxt2(i+1)**2.lt.Rbar*1.d-10) then
            if (abs(a1+2*a2*xbar)*dxt2(i+1).lt.Rbar*1.d-5) then
              Integ=dxt2(i+1)/sqrt(Rbar)
            else
              Integ=2*(sqrt(H-phi2(i))-sqrt(H-phi2(i+1)))/(a1 +2.d0*a2*xbar)
            endif        
         else 
            
            disc=-4.d0*(H-a0)*a2-a1**2

            if ((disc.gt.0.d0).and.(a2.lt.0.d0)) then
               
            borne2=(2.d0*a2*xt2(i+1)+a1)/sqrt(disc)
            borne1=(2.d0*a2*xt2(i)+a1)/sqrt(disc)

            Integ=abs((asinh(borne2)-asinh(borne1))/sqrt(-a2))

            elseif ((disc.lt.0.d0).and.(a2.gt.0.d0)) then 

            borne2=(2.d0*a2*xt2(i+1)+a1)/sqrt(-disc)
            borne1=(2.d0*a2*xt2(i)+a1)/sqrt(-disc)

            Integ=abs((asin(borne2)-asin(borne1))/sqrt(a2))

            elseif ((disc.lt.0.d0).and.(a2.lt.0.d0)) then 

         borne2=abs(2.d0*sqrt((H-(a0+a1*xt2(i+1)+a2*xt2(i+1)**2))*(-a2)) + 2.d0*(-a2)*xt2(i+1)-a1)
         borne1=abs(2.d0*sqrt((H-(a0+a1*xt2(i)+a2*xt2(i)**2))*(-a2)) + 2.d0*(-a2)*xt2(i)-a1)

            Integ=abs((log(borne2)-log(borne1))/sqrt(-a2))

            endif 

      endif

       period=period+Integ
 
       deltaE=deltaE+(phi2(i)-phi2old(i) -vint(i)*0.5d0*(E1s2(i+1)+E1s2old(i+1))*dt)*Integ

      phimoy=phimoy+((phi2(i)-phi2old(i))/dt-vint(i)*0.5d0*(E1s2(i)+E1s2old(i)))*dxt2(i)/phi2old(i)

      path=path+dxt2(i)
    
      endif

      i=i+1

      enddo

! c--------- correction de l'integration --------------------

      if (phi2(m+1).ne.H) then
         a3=(phi2(m+2)-phi2(m+1))/(xt2(m+2)-xt2(m+1))
         a4=phi2(m+1)-a3*xt2(m+1)
         xr=(H-a4)/a3
                  
        C1=xt2(m)-xt2(m+1)-(xt2(m)**2-xt2(m+1)**2)/(xt2(m+1)+xr)
        D1=phi2(m+1)-phi2(m)+(phi2(m+1)-H) & 
                *(xt2(m)**2-xt2(m+1)**2)/(xt2(m+1)**2-xr**2) 
        C2=xt2(i+2)-xr-(xt2(i+2)**2-xr**2)/(xt2(m+1)+xr)
        D2=H-phi2(i+2)+(phi2(m+1)-H) &
                 *(xt2(i+2)**2-xr**2)/(xt2(m+1)**2-xr**2) 


            a1=-(C1*D1+C2*D2)/(C1**2+C2**2)
            a2=(phi2(m+1)-H-a1*(xt2(m+1)-xr)) /(xt2(m+1)**2-xr**2)
            a0=phi2(m+1)-a1*xt2(m+1)-a2*xt2(m+1)**2

            disc=-4.d0*(H-a0)*a2-a1**2

            if ((disc.gt.0.d0).and.(a2.lt.0.d0)) then
               
        borne1=(2.d0*a2*xt2(m+1)+a1)/sqrt(disc)
        borne2=(2.d0*a2*xr+a1)/sqrt(disc)

        Integ=abs((asinh(borne2)-asinh(borne1))/sqrt(-a2))

        elseif ((disc.lt.0.d0).and.(a2.gt.0.d0)) then             

       borne1=(2.d0*a2*xt2(m+1)+a1)/sqrt(-disc)
       borne2=(2.d0*a2*xr+a1)/sqrt(-disc)

      if(borne2.ge.1d0) then
	borne2=1.d0
      endif
  
       Integ=abs((asin(borne2)-asin(borne1))/sqrt(a2))

            elseif ((disc.lt.0.d0).and.(a2.lt.0.d0)) then 

         borne2=abs(2.d0*(-a2)*xr-a1)
         borne1=abs(2.d0*sqrt((H-(a0+a1*xt2(m+1)+a2*xt2(m+1)**2))*(-a2)) + 2.d0*(-a2)*xt2(m+1)-a1)

            Integ=abs((log(borne2)-log(borne1))/sqrt(-a2))

            endif

       period=period+Integ

      deltaE=deltaE+(phi2(m+1)-phi2old(m+1) -vint(m+1)*0.5d0*(E1s2(m+2)+E1s2old(m+2))*dt)*Integ

      phimoy=phimoy+((phi2(m+1)-phi2old(m+1))/dt- vint(m+1)*0.5d0*(E1s2(m+2)+E1s2old(m+2)))*dxt2(m+1) &
      /phi2old(m+1)

      path=path+dxt2(m+1)
               
      Inv2=Inv2+sqrt(H-phi2(m+1))*(xr-xt2(m+1))
    
       endif
! c ---------------- calcul des invariants   

       if (itime.eq.2) then
         Inv0(j)=Inv2
       endif


! c --------------- calcul des nouveaux vej 

       deltaE=(deltaE/period)
        
       if ((2.d0*(deltaE-phi2(0)+phi2old(0))+veold(j)**2).ge.0.d0) then  
       venew=sqrt(2.d0*(deltaE-phi2(0)+phi2old(0))+veold(j)**2)
       endif

         if(iphi.ne.iter) then
            ve(j)=(fiter2*ve(j)+venew)/(fiter2+1.d0)	     
         else
            ve(j)=(fiter3*ve(j)+venew)/(fiter3+1.d0)	     
         endif

         endif 

	if(phi2(1).lt.H) then         
         phimoy=abs(phimoy)/abs(path)
	else
	 phimoy=(phi2(0)+phi2(1))/2.d0
	endif

! c ------------------ SORTIES ---------------------------------------
                
! *        if (iphi.eq.iter) then
! *
! *	write(800,'(f8.2,a,i2,5(a,e12.5))')
! *     *    time,char(9), j,char(9)
! *     *   ,H, char(9), ve(j), char(9)
! *     *   ,period
! *     *   , char(9), 1/phimoy, char(9),deltaE
! *
! *       endif

      enddo   
      endif

	
 10    continue

! *         2.3 champs 'theoriques'

! *	2.3.1 normalisation au champ theorique du modele self-similaire

        if(itime.eq.1)then
 	   Enorm=0.d0
 	else
 	   Ess=1.d0*sqrt(Th)/time
 	   Enorm=E(ncell)/Ess
 	endif

! *	2.3.2 comparaison au champ 'theorique' suppose 'fite' le
! *	      champ au bord a tout instant

 	Ebord0=sqrt(2.d0*n0hot)/exp(0.5d0)
 	Ebordth=2.d0*Ebord0*sqrt(Th)/sqrt(4.d0+(Ebord0*time)**2)

	
! *         2.4 test d'arret sur le nombre d'iterations
	
	if(itmax.eq.1) goto 101


! *         2.5 sorties de fonctions dependant du temps

	vmax=0.d0
	ivmax=ncell
	do i=1,ncell
	   if(v(i).gt.vmax) then
	     vmax=v(i)
	     ivmax=i
	   endif
	enddo


	Emax=0.d0
	do i=1,ncell
	   if(E(i).gt.Emax) then
	     Emax=E(i)
             iEmax=i
	   endif
	enddo

        slop=0.d0
        if (time.gt.lmax) then
           Et1=Eold(ivmax)
           Et2=E(ivmax)
           slop=log(Et2/Et1)/log(1-dt/time)
	vfinal=v(ivmax)+E(ivmax)*time/(slop-1)
	vfinal2=v(ivmax)+E(ivmax)*time/(0.3)
        else
        vfinal=v(ivmax)+E(ivmax)*time
        vfinal2=v(ivmax)+E(ivmax)*time        
        endif   
           
   !      write(100,'(f8.2,8(a,e12.5))') 
   !   *    time,char(9),Emax,char(9),vmax
   !   *    ,char(9),vfinal,char(9),vfinal2
   !   *    ,char(9),slop,char(9),temp
   !   *    ,char(9),temphot,char(9),tempcold

        rhomin=0.d0
        xmin=0.d0
        do i=1,ncell
           if (xt(i).gt.0.d0) then
           if (abs(rho(i)).gt.rhomin) then
              rhomin=abs(rho(i))
              xmin=xt(i)
           endif
           endif
        enddo

	if(itime.eq.1) then
	   x_rar=0.d0
	   csmoyen=cscin
	else
	   x_rar=x_rar+(csold+cscin)*dt/2.d0
	   csmoyen=x_rar/time
	endif

	lDebye=sqrt(Th/ne(ncell))

	nstep=(tmax/dt)
	nsort=1
	if(nstep.gt.2000)nsort=nstep/2000

	nsortf=1
	if(nstep.gt.100)nsortf=nstep/100

	! if(itime-1.ne.((itime-1)/nsort)*nsort.and..not.laststep)
!      *   go to 102
! 	write(9,
!      *'(f8.2,17(a,f10.5))')
!      *   time,char(9),nti,char(9),nthot,char(9),ntcold
!      *   ,char(9),nte,char(9),n0hot,char(9),n0cold
!      *   ,char(9),En_ion,char(9),Whot1
!      *   ,char(9),Whot2,char(9),Whot,char(9),Wcold,char(9),temp
!      *   ,char(9),Tc,char(9),En_elec,char(9),En_elec+En_ion-Whot
!      *   ,char(9),vmax,char(9),vfinal

! 	write(10,
!      *'(f8.2,a,f12.5,a,f7.4,a,f8.4,a,f7.5
!      *  ,a,f10.7,a,f10.7,3(a,f8.6),a,f10.5,a,f10.5,a,f10.5
!      *  ,a,i6,a,i6)')
!      *   time,char(9),xt(ivmax),char(9),v(ivmax),char(9),.5*v(ivmax)**2
!      *   ,char(9),E(ivmax)
!      *   ,char(9),ne(ivmax),char(9),ni(ivmax)
!      *   ,char(9),ni(0)
!      *   ,char(9),nhot(0),char(9),ncold(0)
!      *   ,char(9),lDebye
! c     *   ,char(9),lgrade,char(9),lgradi
! c     *   ,char(9),ivmax,char(9),idebut(ivmax)
	
	Rs=rgauss/ni(0)
	xi=(xt(ivmax)+lmax)/Rs
	tsR0=time/rgauss
	RsR0=Rs/rgauss
	Eq42=2.d0*xi*exp(xi**2/2.d0)/rgauss/sqrt(RsR0)

	! write(14,
   !   *'(f8.2,a,f7.3,a,f7.4,a,f8.3,a,f8.4,a,f7.4)')
   !   *   time,char(9),tsR0,char(9),xi,char(9),RsR0,char(9),Eq42

   !     write(15,'(f8.2,a,f7.4,a,f7.4,a,f8.3,a,f7.5,a,f7.5)')
   !   &  time,char(9), cscin,char(9), csmoyen,char(9), x_rar,char(9), 
   !   &  2.d0*En(0)/ne(0),char(9), sqrt(2.d0*En(0)/ne(0))
	
	! if(itime-1.ne.((itime-1)/nsortf)*nsortf.and..not.laststep)
   !   *   go to 102
      do j=0,Nbe-1,djout
        write(250,'(f8.2,a,f6.4,a,e12.5,a,f7.4,a,f7.4)')  time,char(9),ve(j),char(9),fe0(j),char(9),ve(j)**2/2  ,char(9),Te(j)  
      enddo 

102	continue

! *         2.6 test d'arret sur le temps (ici s'arrete normalement
! *		  un calcul standard avec itmax.gt.1 et tmax.gt.0 apres
! *		  un nombre fini d'iterations)
	
	if(time.ge.(tmax-1.d-5)) goto 101
	if(itime.eq.itmax) goto 101


! *         2.7 incrementation du temps
	! 
	csold=cscin

	dtold=dt
	dt=dti
	if(VTS) dt=dti/sqrt(ni(0)/ni0)
	if(idico.eq.1)dt=dti/100.d0
	if((time+dt).gt.(tmax-1.d-5)) then
	   dt=tmax-time
	   laststep=.true.
	endif
	time=time+dt

! *         2.8 modification de la vitesse aux temps intermediaires
! *             'demi-entiers'
	
	if(itime.eq.1) then

	   do i = 1,ncell
	      vint(i)=v(i)+0.5d0*dt*E(i)
	   end do

	else

! *	   do i = 1,ncell
! *	      vint(i)=vint(i)+0.5*(dtold+dt)*E(i)
! *	   end do

! *            2.8.1 - construction des tableaux a,b,c, et f pour 
! *	       le calcul implicite de la vitesse avec viscosite

	   delt=0.5*(dtold+dt)
	   b(0)=1.d0/delt
	   c(0)=0.d0
	   f(0)=0.d0
	
	   do i=1,ncell-1
	      a(i)=-2.d0*nu/dxt(i)/(dxt(i)+dxt(i+1))
	      b(i)=1.d0/delt+2.d0*nu/(dxt(i)*dxt(i+1))
	      c(i)=-2.d0*nu/dxt(i+1)/(dxt(i)+dxt(i+1))
	      f(i)=E(i)+vint(i)/delt
	   enddo
	
	   a(ncell)=0.d0
	   b(ncell)=1.d0/delt
           f(ncell)=E(ncell)+vint(ncell)/delt
     
! *	      2.8.2 - inversion de la matrice tridiagonale

	   call gauss(nmax,ncell,a,b,c,f,bx,fx,vint)
	
	endif
	

! *         2.9 modification de la position
	
	do i = 1,ncell
	   xt(i)=xt(i)+dt*vint(i)
	end do
	

! *         2.10 rearrangement des numeros des ions

! c------- temps de deferlement
       
	do i=1,ncell-1
        if ((xt(i).gt.xt(i+1)).and.(.not.break)) then
           tb=time-dt
           break=.true.
        endif
        enddo

	if(.not.break) goto 100
! c-------	
         do j=2,ncell
	   xxt=xt(j)
	   xx0=x0(j)
	   xxtold=xtold(j)
	   vv=v(j)
	   vvint=vint(j)
	   iidebut=idebut(j)
	   qqiSS=qiSS(j)
	   nnhot=nhot(j)
	   nncold=ncold(j)
	   pphi=phi(j)
	   EE=E(j)
	   ggradnh=gradnh(j)
	   ggradnc=gradnc(j)
	   ddxt=dxt(j)
	   do i=j-1,1,-1
	      if(xt(i).le.xxt) goto 95
	      xt(i+1)=xt(i)
	      x0(i+1)=x0(i)
	      xtold(i+1)=xtold(i)
	      v(i+1)=v(i)
	      vint(i+1)=vint(i)
	      idebut(i+1)=idebut(i)
		qiSS(i+1)=qiSS(i)
		nhot(i+1)=nhot(i)
		ncold(i+1)=ncold(i)
		phi(i+1)=phi(i)
		E(i+1)=E(i)
		gradnh(i+1)=gradnh(i)
		gradnc(i+1)=gradnc(i)
		dxt(i+1)=dxt(i)
	   enddo
	   i=0
 95	   xt(i+1)=xxt
	   x0(i+1)=xx0
	   xtold(i+1)=xxtold
	   v(i+1)=vv
	   vint(i+1)=vvint
	   idebut(i+1)=iidebut
	   qiSS(i+1)=qqiSS
	   nhot(i+1)=nnhot
	   ncold(i+1)=nncold
	   phi(i+1)=pphi
	   E(i+1)=EE
	   gradnh(i+1)=ggradnh
	   gradnc(i+1)=ggradnc
	   dxt(i+1)=ddxt
	enddo
	
	do j=1,ncell
	   irang(idebut(j))=j
	enddo

	
100     continue 
101	continue
	close(9)
	close(10)
	close(14)
        close(15)
! *        close(800)
	close(100)
	close(250)

	write(*,*) 'OK etape2'

                	
! *               3 - calcul de la densite theorique (modele
! *		    self-similaire)
	
      if(itime.eq.1)then
         do i=ncell+1,ntotal
            x0(i)=xt(i)
         enddo
      endif
	
	if(itmax.eq.1.or.time.le.1.d-04) goto 201

	do i=0,ntotal
	   xi=xt(i)/sqrt(Th)/time
	   if(xi.lt.-1.d0) then
	      nss(i)=n0hot
	   elseif(xi.lt.40.d0) then
	      nss(i)=n0hot*exp(-(xi+1.d0))
	   else
	      nss(i)=0.
	   endif
	enddo
	
201	continue

	write(*,*) 'OK etape3'

! c-------------------------------------------------------
        write(*,*) 'vmax=',vmax
        write(*,*) 'vfinal=',vfinal,vfinal2
        write(*,*) 'Lss=',LSS
        write(*,*) 'length=',length
        write(*,*) 'Tdef=',tb
        write(*,*) 'Encell=',E(ncell), E(ncell-1), Enou
        write(*,*) 'phi(0)=', phi(0), phi0
! c--------------------------------------------------------- 
! *		4.2 spectres en vitesse et en energie
	
	do i=1,ncell
	   vmoy(i)=0.5*(v(i-1)+v(i))
	   Emoy(i)=0.25*(v(i-1)**2+v(i)**2)
	   dndv(i)=(niSS(i-1)+niSS(i))*dx0(i)/(v(i)-v(i-1))/2.d0
	   dndv(i)=abs(dndv(i))
	   dndE(i)=(niSS(i-1)+niSS(i))*dx0(i)/(v(i)**2-v(i-1)**2)
	   dndE(i)=abs(dndE(i))
	enddo
	
	
	write(*,*) 
	write(*,*) time
	write(*,*) 
	stop
	end      

      Subroutine Gauss(nmax,ncell,A,B,C,F,bx,fx,phi)



      double precision phi (0:nmax)
	double precision a (1:nmax),  b(0:nmax), c(0:nmax-1), f(0:nmax), bx(0:nmax), fx(0:nmax)


      j=0
      Bx(j) = B(j)
      Fx(j) = F(j)
	
      do j=1,ncell
         Bx(j) = B(j) - A(j)*C(j-1) /Bx(j-1) 
         Fx(j) = F(j) - A(j)*Fx(j-1)/Bx(j-1)   
      end do
	
! C... descending recursion for solution

      j = ncell
      phi(j) = Fx(j)/Bx(j)
	
      do j=ncell-1,0,-1
         phi(j) = ( Fx(j)-C(j)*phi(j+1) )/Bx(j)
      end do

      Return
      End
	
! c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
! c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**

      function Thot(time)
      
      implicit double precision (a-h,o-z)
      if(time.le.9.9d0) then
	   Thot=1.0d0
	else
	   Thot=1.0d0
	endif
      
      Return
      End
	
! c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
! c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**

      function Tcold(time)
      
      implicit double precision (a-h,o-z)
      if(time.le.9.9d0) then
	   Tcold=1.d0
	else
	   Tcold=1.d0
	endif
      
      Return
      End
	
! c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
! c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**

      function asinh(x)

      double precision asinh,x
      asinh=log(x+sqrt(x**2+1))
      end
     
! c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
! c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**

      function erfc(x)

!
!! ERFC evaluates the error function ERFC(X).
!
!  Modified:
!
!    09 May 2003 - 31 March 2008 (by P M)
!
!  Author: (modified by Patrick M. to output ERFC instead of ERF)
!
!    W J Cody,
!    Mathematics and Computer Science Division,
!    Argonne National Laboratory,
!    Argonne, Illinois, 60439.
!
!  Reference:
!
!    W J Cody,
!    "Rational Chebyshev approximations for the error function",
!    Mathematics of Computation, 
!    1969, pages 631-638.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of ERFC.
!
!    Output, real ( kind = 8 ) ERFC, the value of ERFC(X).
!
      implicit none

      real(kind=8), save, dimension(5) :: a = (/ &
      3.16112374387056560D+00, &
      1.13864154151050156D+02, &
      3.77485237685302021D+02, &
      3.20937758913846947D+03, &
      1.85777706184603153D-01 /)
  real(kind=8), save, dimension(4) :: b = (/ &
      2.36012909523441209D+01, &
      2.44024637934444173D+02, &
      1.28261652607737228D+03, &
      2.84423683343917062D+03 /)
  real(kind=8), save, dimension(9) :: c = (/ &
      5.64188496988670089D-01, &
      8.88314979438837594D+00, &
      6.61191906371416295D+01, &
      2.98635138197400131D+02, &
      8.81952221241769090D+02, &
      1.71204761263407058D+03, &
      2.05107837782607147D+03, &
      1.23033935479799725D+03, &
      2.15311535474403846D-08 /)
  real(kind=8), save, dimension(8) :: d = (/ &
      1.57449261107098347D+01, &
      1.17693950891312499D+02, &
      5.37181101862009858D+02, &
      1.62138957456669019D+03, &
      3.29079923573345963D+03, &
      4.36261909014324716D+03, &
      3.43936767414372164D+03, &
      1.23033935480374942D+03 /)
  real(kind=8), save, dimension(6) :: p = (/ &
      3.05326634961232344D-01, &
      3.60344899949804439D-01, &
      1.25781726111229246D-01, &
      1.60837851487422766D-02, &
      6.58749161529837803D-04, &
      1.63153871373020978D-02 /)
  real(kind=8), save, dimension(5) :: q = (/ &
      2.56852019228982242D+00, &
      1.87295284992346047D+00, &
      5.27905102951428412D-01, &
      6.05183413124413191D-02, &
      2.33520497626869185D-03 /)

      real(kind=8), parameter :: SQRPI = 0.56418958354775628695D+00
      real(kind=8), parameter :: THRESH = 0.46875D+00
      real(kind=8), parameter :: XBIG = 26.543D+00
      real(kind=8) del
      real(kind=8) erf
      real(kind=8) erfc

      ! real ( kind = 8 ), parameter :: SQRPI = 0.56418958354775628695D+00
      ! real ( kind = 8 ), parameter :: THRESH = 0.46875D+00
      real ( kind = 8 ) x
      real ( kind = 8 ) xabs
      ! real ( kind = 8 ), parameter :: XBIG = 26.543D+00
      real ( kind = 8 ) xden
      real ( kind = 8 ) xnum
      real ( kind = 8 ) xsq
      integer :: i 

      xabs = abs ( x )
!
!  Evaluate ERF(X) for |X| <= 0.46875.
!
      if ( xabs <= THRESH ) then

         if ( epsilon ( xabs ) < xabs ) then
            xsq = xabs * xabs
         else
            xsq = 0.0D+00
         end if

         xnum = a(5) * xsq
         xden = xsq
         do i = 1, 3
            xnum = ( xnum + a(i) ) * xsq
            xden = ( xden + b(i) ) * xsq
         end do

         erf = x * ( xnum + a(4) ) / ( xden + b(4) )

	erfc=1.d0-erf
!
!  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
!
      else if ( xabs <= 4.0D+00 ) then

         xnum = c(9) * xabs
         xden = xabs
         do i = 1, 7
            xnum = ( xnum + c(i) ) * xabs
            xden = ( xden + d(i) ) * xabs
         end do

         erf = ( xnum + c(8) ) / ( xden + d(8) )
         xsq = aint ( xabs * 16.0D+00 ) / 16.0D+00
         del = ( xabs - xsq ) * ( xabs + xsq )
         erfc = exp ( - xsq * xsq ) * exp ( - del ) * erf

! *         erf = ( 0.5D+00 - erfc ) + 0.5D+00
	
         if ( x < 0.0D+00 ) then
            erfc = 2.d0 - erfc
         end if
!
!  Evaluate ERFC(X) for 4.0 < |X|.
!
      else

         if ( XBIG <= xabs ) then

            if ( 0.0D+00 < x ) then
               erfc = 0.0D+00
            else
               erfc = 2.0D+00
            end if

         else

            xsq = 1.0D+00 / ( xabs * xabs )

            xnum = p(6) * xsq
            xden = xsq
            do i = 1, 4
               xnum = ( xnum + p(i) ) * xsq
               xden = ( xden + q(i) ) * xsq
            end do

            erf = xsq * ( xnum + p(5) ) / ( xden + q(5) )
            erf = ( SQRPI -  erf ) / xabs
            xsq = aint ( xabs * 16.0D+00 ) / 16.0D+00
            del = ( xabs - xsq ) * ( xabs + xsq )
            erfc = exp ( - xsq * xsq ) * exp ( - del ) * erf

! *            erf = ( 0.5D+00 - erfc ) + 0.5D+00
            if ( x < 0.0D+00 ) then
            erfc = 2.d0 - erfc
            end if

         end if

      end if

      return
      end

! c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
! c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
