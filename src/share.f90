!module mod_share
!    use mod_types
!    implicit none
!    private
!    public :: nmax, ncell, nvacuum, Nbe, itmax,  &
!              ntotal, nsort, iline, nstep, iline_vide, itime, ivmax, &
!              idebut, irang, itest, m, djout, i, j, idico, iter, iphi, &
!              x0, xt, xtold, v, vint, E, Eold, phi, ni, nilisse, &
!              nelisse, phiold, ne, rho, nhot,  nhotold, &
!              gradnh, ghold,  phot, pe, nss, difth, &
!              dWhot, SWhot, niSS, qiSS, charge, b, f, bx, fx, vtest, &
!              vtestint, Etest, Etestold, x0test, xttest, xttestold, vmoy, &
!              Emoy, dndv, dndE, dx0, dxt, dxtold, ni1s2, E1s2, a, c, &
!              length, ni0, neh, neold, k
!    real(dp), parameter :: n0hot=0.0
!    logical  :: lfini,nb_cons,En_cons,EOS,VTT,laststep
!    character(len=10) :: profil
!end module mod_share
!

module mod_share
use mod_types
implicit none
public
!
! parameters
!
integer, parameter :: nmax=80000
integer  :: ncell, nvacuum,Nbe, itmax, jcoldmax, &
            ntotal, nsort, iline, nstep, iline_vide, itime,ivmax, &
            idebut(0:nmax),irang(0:nmax),itest(0:nmax),m,djout, &
            i, j, idico, iter, iphi, iidebut
real(8), dimension(0:nmax) :: x0, xt, xtold, v, vint, E, Eold, &
                               phi, ni, nilisse, nelisse, phiold, &
                               ne, rho, nhot,ncold, nhotold, ncoldold, &
                               gradnh, gradnc, ghold, gcold, phot, &
                               pcold, pe, nss, difth, dWhot, SWhot, &
                               niSS, qiSS, charge,b, f, bx, fx, &
                               vtest, vtestint, Etest, Etestold, x0test, xttest, xttestold
real(8), dimension(1:nmax) :: vmoy, Emoy, dndv, dndE, dx0, dxt, dxtold, ni1s2, E1s2, a
real(8), dimension(0:nmax-1) :: c
real(8) :: length,ni0,n0cold,n0hot,neh,nec,neold,kvide, &
            nti,nthot,ntcold,nte,lmax, &
            nthot0nhm3,nhm2,nhm1,LSS,nu,nnhot,nncold, nhm3, Tcm3, Thm3, &
            lDebye,lgrade,lgradi,lay1,lay2,mix,mixp, &
            ve_max, prog, iter0, iter1, iter2, iter3, dti,tmax, &
            Tcmax,Thmax, T_MeV, nLSS, charge2, p_negatif, &
            omegadt, Th, Tc, dt, cs, cs2, time, cs2old, rgauss, &
            e0, Tnorm, Tn0, alpha, whot, Thm1, Thm2, Whot2, Whot1, &
            Tcm1,Tcm2,Wcold, Tcold, Thot, Thold,Tcoold, En_ion, En_totale, En_elec, &
            whot2old, whot1old, peold, enew, thnew, &
            grade, gradi, uphi, en_hot0, En_cold, &
            En_cold0,nthot0,trise, xxtold, xxt, xx0, &
            vvint, vv, qqiss, pphi, ggradnh, ggradnc, &
            EE, ccharge, xi, Ess,Enorm, Ebord0, Ebordth, &
            vfinal, vmax, dtold, delt, ddxt, Eq42, RsR0, tsR0, Rs,Ztest
integer  :: fileunit
logical  :: lfini,nb_cons,En_cons,EOS,VTT,laststep,multicouche,multiphase
character(len=10) :: profil
end module mod_share
