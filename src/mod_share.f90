module mod_share
use mod_types
implicit none
public
!
! parameters
!
integer, parameter :: nmax=80000
integer  :: ncell, nvacuum,Nbe, itmax, jcoldmax
integer  :: ntotal, nsort, iline, nstep, iline_vide, itime,ivmax
integer  :: idebut(0:nmax),irang(0:nmax),itest(0:nmax),m,djout
integer  :: i, j, idico, iter, iphi, iidebut
real(dp) :: length,ni0,n0cold,n0hot,neh,nec,neold,kvide
real(dp) :: nti,nthot,ntcold,nte,lmax
real(dp) :: nthot0nhm3,nhm2,nhm1,LSS,nu,nnhot,nncold, nhm3, Tcm3, Thm3
real(dp) :: lDebye,lgrade,lgradi,lay1,lay2,mix,mixp
real(dp) :: x0 (0:nmax), xt (0:nmax), xtold(0:nmax)
real(dp) :: v  (0:nmax),  vint (0:nmax)
real(dp) :: E  (0:nmax),  Eold (0:nmax)
real(dp) :: phi (0:nmax),  ni (0:nmax)
real(dp) :: nilisse (0:nmax),  nelisse (0:nmax)
real(dp) :: phiold (0:nmax)
real(dp) :: ne (0:nmax),  rho (0:nmax)
real(dp) :: nhot (0:nmax),  ncold (0:nmax)
real(dp) :: nhotold (0:nmax),  ncoldold (0:nmax)
real(dp) :: gradnh (0:nmax),  gradnc (0:nmax)
real(dp) :: ghold (0:nmax),  gcold (0:nmax)
real(dp) :: phot (0:nmax),  pcold (0:nmax)
real(dp) :: pe (0:nmax)
real(dp) :: nss (0:nmax),  difth (0:nmax)
real(dp) :: dWhot (0:nmax), SWhot(0:nmax)
real(dp) :: vmoy (1:nmax), Emoy (1:nmax), dndv (1:nmax), dndE (1:nmax)
real(dp) :: dx0  (1:nmax),  dxt (1:nmax), dxtold(1:nmax)
real(dp) :: ni1s2 (1:nmax), E1s2 (1:nmax), niSS (0:nmax), qiSS(0:nmax), charge(0:nmax)
real(dp) :: a (1:nmax), b(0:nmax), c(0:nmax-1), f(0:nmax), bx(0:nmax), fx(0:nmax)
real(dp) :: ve_max, prog, iter0, iter1, iter2, iter3, dti,tmax
real(dp) :: Tcmax,Thmax, T_MeV, nLSS, charge2, p_negatif
real(dp) :: omegadt, Th, Tc, dt, cs, cs2, time, cs2old, rgauss
real(dp) :: e0, Tnorm, Tn0, alpha, whot
real(dp) :: Thm1, Thm2, Whot2, Whot1, Tcm1,Tcm2,Wcold, Tcold, Thot, Thold,Tcoold, En_ion, En_totale, En_elec
real(dp) :: enew, thnew
real(dp) :: whot2old, whot1old, wcoldold, peold, grade, gradi, uphi, en_hot0, En_cold, En_cold0,nthot0,trise
real(dp) :: xxtold, xxt, xx0, vvint, vv, qqiss, pphi, ggradnh, ggradnc, EE, ccharge, xi
real(dp) :: Ess,Enorm, Ebord0, Ebordth, vfinal, vmax, dtold, delt, ddxt
real(dp) :: Eq42, RsR0, tsR0, Rs,Ztest
real(dp) :: x0test (0:nmax), xttest (0:nmax), xttestold(0:nmax)
real(dp) :: vtest (0:nmax), vtestint (0:nmax)
real(dp) :: Etest  (0:nmax),  Etestold (0:nmax)
integer  :: fileunit
logical  :: lfini,nb_cons,En_cons,EOS,VTT,laststep,multicouche,multiphase
character(len=10) :: profil
!
end module mod_share

