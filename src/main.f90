program IonBoost
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mod_share
    use mod_io, only: read_input
    implicit none

    call read_input()

    profil = 'step'
    nu = 0.d0
    VTT = .true.
    En_hot0 = 0.d0
    E = 0.d0
    Eold = 0.d0
    phi = 0.d0
    phiold = 0.d0
    v = 0.d0
    vint = 0.d0
    nhot = 0.d0
    nhotold = 0.d0
    ne = 0.d0
    ni = 0.d0
    rho = 0.d0
    gradnh = 0.d0
    ghold = 0.d0
    phot = 0.d0
    pe = 0.d0
    dWhot = 0.d0
    SWhot = 0.d0
    charge = 1.d0
    idebut = 0

    do i = 0, ncell
        xt(i) = x0(i)
        xtold(i) = x0(i)
        nhot(i) = n0hot
        ne(i) = n0hot
        niSS(i) = ni0
        idebut(i) = i
        irang(i) = i
    end do

    qiSS(0) = niSS(0) * dx0(1)
    do i = 1, ncell - 1
        qiSS(i) = niSS(i) * (dx0(i) + dx0(i + 1)) / 2.d0
    end do
    qiSS(ncell) = niSS(ncell) * dx0(ncell) * (1.d0 + prog) / 2.d0

    do i = 1, ncell
        dxt(i) = dx0(i)
        dxtold(i) = dx0(i)
    end do

    Tnorm = T_MeV / 0.511d0
    e0 = 0.5d0 * Tnorm * (1.d0 + 9.d0 * Tnorm / 4.d0 + 3.d0 * Tnorm**2 / 4.d0) / &
         (1.d0 + 3.d0 * Tnorm / 2.d0 + 3.d0 * Tnorm**2 / 8.d0)
    Tn0 = 2.d0 * e0
    Th = Thmax

    if (Th == 0.d0) then
        write(*, *) 'Hot temperature should not be zero.'
        stop
    end if

    phi(ncell) = Th
    E(ncell) = sqrt(2.d0 * n0hot * Th * exp(-phi(ncell) / Th))
    alpha = E(ncell) / phi(ncell)

    do i = 0, ncell
        phi(i) = phi(ncell) * exp(alpha * (x0(i) - x0(ncell)))
    end do

    do i = ncell + 1, ntotal
        xt(i) = xt(ncell)
        xtold(i) = xt(i)
        v(i) = 0.d0
        vint(i) = 0.d0
        ni(i) = 1.d-12
        ne(i) = 0.d0
        nhot(i) = 0.d0
        dxt(i) = 0.d0
        dxtold(i) = 0.d0
    end do

    E(0) = 0.d0
    gradnh(0) = 0.d0
    time = 0.d0
    dt = dti
    dtold = dt
    nhm2 = n0hot
    nhm1 = n0hot
    Thm2 = Thmax
    Thm1 = Thm2
    Whot1 = 0.d0
    Whot2 = 0.d0

    open(unit = 9, file = 'conservation.txt', status = 'unknown')
    open(unit = 10, file = 'historique.txt', status = 'unknown')
    write(9, '(a)') '# time nti nthot nte n0hot En_ion Whot1 Whot2 Whot Th En_elec En_balance vmax vfinal'
    write(10, '(a)') '# time xt(ivmax) v(ivmax) Energy_max E(ivmax) ne(ivmax) ni(ivmax) ni(0) nhot(0) lDebye lgrad_e_pressure lgrad_i_pressure ivmax idebut(ivmax)'

    do itime = 1, itmax
        if (itime <= 3) then
            iter = int(iter0)
        else
            iter = int(iter1)
        end if

        Thold = Th
        if (.not. En_cons) then
            Th = Thmax
            idico = 0
            if (Th == 0.d0) then
                write(*, *) 'Hot electron temperature should not be zero!'
                stop
            end if
            if (abs((Th - Thold) / max(Thold, 1.d-12)) > 0.05d0) iter = int(iter0)
        end if

        if (nb_cons) then
            nhm3 = nhm2
            nhm2 = nhm1
            nhm1 = n0hot
            n0hot = 2.d0 * nhm1 - nhm2
            n0hot = max(n0hot, 1.d-12)
        end if

        if (En_cons) then
            Thm3 = Thm2
            Thm2 = Thm1
            Thm1 = Th
            Th = 2.d0 * Thm1 - Thm2
            Th = max(Th, 1.d-12)
        end if

        if (mod(itime, max(nsort, 1)) == 0) print *, time

        if (itime > 1) then
            do i = 0, ntotal
                nhotold(i) = nhot(i)
                phiold(i) = phi(i)
                Eold(i) = E(i)
                ghold(i) = gradnh(i)
            end do
            do i = ncell + 1, ntotal
                xtold(i) = xt(i)
            end do
            do i = 1, ntotal
                dxtold(i) = dxt(i)
            end do
            Whot1old = Whot1
            Whot2old = Whot2
        end if

        do i = 1, ncell
            dxt(i) = xt(i) - xt(i - 1)
        end do

        ni(0) = qiSS(0) / max(dxt(1), 1.d-12)
        do i = 1, ncell - 1
            ni(i) = 2.d0 * qiSS(i) / max(dxt(i) + dxt(i + 1), 1.d-12)
        end do
        ni(ncell) = 2.d0 * qiSS(ncell) / max(dxt(ncell), 1.d-12) / &
                    max(1.d0 + 2.d0 * dxt(ncell) / max(dxt(ncell - 1), 1.d-12) - &
                        dxt(ncell - 1) / max(dxt(ncell - 2), 1.d-12), 1.d-12)

        do i = 1, ncell
            ni1s2(i) = 0.5d0 * (ni(i - 1) + ni(i))
        end do

        do iphi = 1, iter
            neh = n0hot * exp(-phi(0) / Th)
            b(0) = -2.d0 / max(dxt(1) * dxt(1), 1.d-12) - neh / Th
            c(0) = 2.d0 / max(dxt(1) * dxt(1), 1.d-12)
            f(0) = ni(0) - neh * (1.d0 + phi(0) / Th)

            do i = 1, ncell - 1
                neh = n0hot * exp(-phi(i) / Th)
                a(i) = 2.d0 / max(dxt(i) * (dxt(i) + dxt(i + 1)), 1.d-12)
                b(i) = -2.d0 / max(dxt(i) * dxt(i + 1), 1.d-12) - neh / Th
                c(i) = 2.d0 / max(dxt(i + 1) * (dxt(i) + dxt(i + 1)), 1.d-12)
                f(i) = ni(i) - neh * (1.d0 + phi(i) / Th)
            end do

            neh = n0hot * exp(-phi(ncell) / Th)
            neold = neh
            peold = neh * Th
            a(ncell) = 2.d0 / max(dxt(ncell) * dxt(ncell), 1.d-12)
            b(ncell) = -a(ncell) - sqrt(2.d0 / max(peold, 1.d-12)) * neold / max(dxt(ncell), 1.d-12) - neh / Th
            f(ncell) = ni(ncell) - neh * (1.d0 + phi(ncell) / Th) - &
                       sqrt(8.d0 * max(peold, 1.d-12)) / max(dxt(ncell), 1.d-12) * &
                       (1.d0 + 0.5d0 * neold * phi(ncell) / max(peold, 1.d-12))

            call gauss_solve(nmax, ncell, a, b, c, f, bx, fx, phi)
            if (iphi /= iter .and. itime == 1) cycle
            if (iphi /= iter .and. (.not. nb_cons) .and. (.not. En_cons)) cycle

            do i = 0, ncell
                nhot(i) = n0hot * exp(-phi(i) / Th)
                ne(i) = nhot(i)
                rho(i) = ni(i) - ne(i)
            end do

            phot(ncell) = nhot(ncell) * Th
            pe(ncell) = phot(ncell)
            grade = (log(max(ne(ncell), 1.d-12)) - log(max(ne(ncell - 1), 1.d-12))) / &
                    max(xt(ncell) - xt(ncell - 1), 1.d-12)
            gradi = (log(max(ni(ncell), 1.d-12)) - log(max(ni(ncell - 1), 1.d-12))) / &
                    max(xt(ncell) - xt(ncell - 1), 1.d-12)

            kvide = sqrt(n0hot / 2.d0 / Th)
            uphi = 0.d0
            do i = ncell + 1, ntotal
                if (i == ncell + 1) then
                    dxt(ncell + 1) = 0.d0
                else
                    dxt(i) = 0.025d0 * sqrt(Th / max(nhot(i - 1), 1.d-12))
                end if
                xt(i) = xt(i - 1) + dxt(i)
                ! Boltzmann electrons: pe*exp(phi/Th) = n0hot*Th is constant, so the
                ! exact sheath solution phi = phi(ncell) + 2*Th*log(1 + kvide*x*exp(-phi(ncell)/(2*Th)))
                ! needs uphi = kvide*x with kvide = sqrt(n0hot/(2*Th))
                uphi = uphi + kvide * dxt(i)
                phi(i) = phi(ncell) + 2.d0 * Th * log(1.d0 + uphi * exp(-0.5d0 * phi(ncell) / Th))
                nhot(i) = n0hot * exp(-phi(i) / Th)
                ne(i) = nhot(i)
                rho(i) = -ne(i)
                phot(i) = nhot(i) * Th
                pe(i) = phot(i)
                E(i) = sqrt(2.d0 * pe(i))
            end do

            nti = 0.d0
            do i = 1, ncell
                nti = nti + ni1s2(i) * dxt(i)
            end do

            nthot = 0.d0
            do i = 1, ntotal
                nthot = nthot + 0.5d0 * (nhot(i - 1) + nhot(i)) * dxt(i)
            end do
            nthot = nthot + E(ntotal)
            nte = nthot

            if (itime == 1) En_hot0 = nthot * (e0 / Tn0) * Th

            if (nb_cons) then
                if (itime == 1) then
                    nthot0 = nthot
                else if (iphi > 1) then
                    n0hot = n0hot * (nthot0 / max(nthot, 1.d-12))
                    n0hot = max(n0hot, 1.d-12)
                end if
            end if

            if (itime > 1) then
                do i = ncell + 1, ntotal
                    vint(i) = (xt(i) - xtold(i)) / max(dt, 1.d-12)
                end do
            end if

            do i = 1, ntotal - 1
                gradnh(i) = (nhot(i + 1) - nhot(i - 1)) / max(dxt(i) + dxt(i + 1), 1.d-12)
            end do
            gradnh(ncell) = gradnh(ncell - 1) * nhot(ncell) / max(nhot(ncell - 1), 1.d-12)
            gradnh(ntotal) = gradnh(ntotal - 1) * nhot(ntotal) / max(nhot(ntotal - 1), 1.d-12)
            gradnh(ncell + 1) = gradnh(ncell)

            if (itime > 1 .and. iphi > 1) then
                dWhot(0) = 0.5d0 * (phiold(0) + phi(0)) * &
                           (nhot(0) - nhotold(0) - vint(0) * dt * 0.5d0 * (gradnh(0) + ghold(0)))
                Whot1 = Whot1old + 0.25d0 * dWhot(0) * (dxtold(1) + dxt(1))
                SWhot(0) = 0.25d0 * dWhot(0) * (dxtold(1) + dxt(1))
                do i = 1, ncell
                    dWhot(i) = 0.5d0 * (phiold(i) + phi(i)) * &
                               (nhot(i) - nhotold(i) - vint(i) * dt * 0.5d0 * (gradnh(i) + ghold(i)))
                    Whot1 = Whot1 + 0.25d0 * dWhot(i) * (dxtold(i) + dxt(i) + dxtold(i + 1) + dxt(i + 1))
                    SWhot(i) = SWhot(i - 1) + 0.25d0 * dWhot(i) * (dxtold(i) + dxt(i) + dxtold(i + 1) + dxt(i + 1))
                end do
                Whot2 = Whot2old
                do i = ncell + 1, ntotal - 1
                    dWhot(i) = 0.5d0 * (phiold(i) + phi(i)) * &
                               (nhot(i) - nhotold(i) - vint(i) * dt * 0.5d0 * (gradnh(i) + ghold(i)))
                    Whot2 = Whot2 + 0.25d0 * dWhot(i) * (dxtold(i) + dxt(i) + dxtold(i + 1) + dxt(i + 1))
                    SWhot(i) = SWhot(i - 1) + 0.25d0 * dWhot(i) * (dxtold(i) + dxt(i) + dxtold(i + 1) + dxt(i + 1))
                end do
                dWhot(ntotal) = 0.5d0 * (phiold(ntotal) + phi(ntotal)) * &
                                (nhot(ntotal) - nhotold(ntotal) - vint(ntotal) * dt * 0.5d0 * (gradnh(ntotal) + ghold(ntotal)))
                Whot2 = Whot2 + 0.25d0 * dWhot(ntotal) * (dxtold(ntotal) + dxt(ntotal))
                SWhot(ntotal) = SWhot(ntotal - 1) + 0.25d0 * dWhot(ntotal) * (dxtold(ntotal) + dxt(ntotal))
                Whot = Whot1 + Whot2
            end if

            do i = 1, ncell
                E1s2(i) = (phi(i) - phi(i - 1)) / max(dxt(i), 1.d-12)
            end do
            do i = 1, ncell - 1
                E(i) = (dxt(i + 1) * E1s2(i) + dxt(i) * E1s2(i + 1)) / max(dxt(i) + dxt(i + 1), 1.d-12)
            end do
            E(ncell) = sqrt(2.d0 * pe(ncell))

            if (itime > 1) then
                do i = 1, ncell
                    v(i) = vint(i) + dt * (3.d0 * E(i) + Eold(i)) * charge(i) / 8.d0
                end do
            end if

            En_ion = 0.d0
            do i = 1, ncell - 1
                En_ion = En_ion + qiSS(i) * v(i)**2 / 2.d0 / max(charge(i), 1.d-12)
            end do
            if (profil == 'step') En_ion = En_ion + qiSS(ncell) * v(ncell)**2 / 4.d0 / max(charge(ncell), 1.d-12)
            if (profil /= 'step') En_ion = En_ion + qiSS(ncell) * v(ncell)**2 / 2.d0 / max(charge(ncell), 1.d-12)

            En_elec = 0.d0
            do i = 1, ntotal - 1
                En_elec = En_elec + 0.25d0 * E(i)**2 * (dxt(i) + dxt(i + 1))
            end do
            En_elec = En_elec + 0.25d0 * E(ntotal)**2 * dxt(ntotal)
            En_elec = En_elec + E(ntotal) * Th

            if (itime == 1) En_totale = En_ion + En_elec + En_hot0

            if (En_cons .and. itime > 1 .and. iphi > 1) then
                enew = (En_totale - En_ion - En_elec) * Tn0 / max(nthot * Thmax, 1.d-12)
                Thnew = 2.d0 * Thmax * enew / max(Tn0, 1.d-12)
                if (iphi /= iter) then
                    Th = (iter2 * Th + Thnew) / (iter2 + 1.d0)
                else
                    Th = (iter3 * Th + Thnew) / (iter3 + 1.d0)
                end if
                Th = max(Th, 1.d-12)
            end if

        end do

        if (itime == 1) then
            Enorm = 0.d0
        else
            Ess = sqrt(Th) / max(time, 1.d-12)
            Enorm = E(ncell) / max(Ess, 1.d-12)
        end if
        Ebord0 = sqrt(2.d0 * n0hot) / exp(0.5d0)
        Ebordth = 2.d0 * Ebord0 * sqrt(Th) / sqrt(4.d0 + (Ebord0 * time)**2)

        if (itmax == 1) exit

        vmax = 0.d0
        ivmax = ncell
        do i = 1, ncell
            if (v(i) > vmax) then
                vmax = v(i)
                ivmax = i
            end if
        end do
        vfinal = v(ivmax) + E(ivmax) * charge(ivmax) * time
        lDebye = sqrt(Th / max(ne(ncell), 1.d-12))
        lgrade = -1.d0 / sign(max(abs(grade), 1.d-12), grade)
        if (itime == 1) then
            lgradi = 0.d0
        else
            lgradi = -1.d0 / sign(max(abs(gradi), 1.d-12), gradi)
        end if

        nstep = int(tmax / max(dt, 1.d-12))
        nsort = 1
        if (nstep > 2000) nsort = max(nstep / 2000, 1)
        if (mod(itime - 1, nsort) == 0 .or. laststep) then
            write(9, '(*(g0,:,","))') time, nti, nthot, nte, n0hot, En_ion, &
                Whot1, Whot2, Whot, Th, En_elec, En_elec + En_ion - Whot, vmax, vfinal
            write(10, '(*(g0,:,","))') time, xt(ivmax), v(ivmax), 0.5d0 * v(ivmax)**2, E(ivmax), ne(ivmax), ni(ivmax), &
                ni(0), nhot(0), lDebye, lgrade, lgradi, ivmax, idebut(ivmax)
        end if

        Rs = rgauss / max(ni(0), 1.d-12)
        xi = (xt(ivmax) + lmax) / max(Rs, 1.d-12)
        tsR0 = time / max(rgauss, 1.d-12)
        RsR0 = Rs / max(rgauss, 1.d-12)
        Eq42 = 2.d0 * xi * exp(xi**2 / 2.d0) / max(rgauss * sqrt(max(RsR0, 1.d-12)), 1.d-12)

        if (time >= tmax - 1.d-5) exit
        if (itime == itmax) exit

        dtold = dt
        dt = dti
        if (VTT) dt = dti / sqrt(max(ni(0) / ni0, 1.d-12))
        if ((time + dt) > (tmax - 1.d-5)) then
            dt = tmax - time
            laststep = .true.
        end if
        time = time + dt

        if (itime == 1) then
            do i = 1, ncell
                vint(i) = v(i) + 0.5d0 * dt * E(i) * charge(i)
            end do
        else
            delt = 0.5d0 * (dtold + dt)
            b(0) = 1.d0 / max(delt, 1.d-12)
            c(0) = 0.d0
            f(0) = 0.d0
            do i = 1, ncell - 1
                a(i) = -2.d0 * nu / max(dxt(i) * (dxt(i) + dxt(i + 1)), 1.d-12)
                b(i) = 1.d0 / max(delt, 1.d-12) + 2.d0 * nu / max(dxt(i) * dxt(i + 1), 1.d-12)
                c(i) = -2.d0 * nu / max(dxt(i + 1) * (dxt(i) + dxt(i + 1)), 1.d-12)
                f(i) = E(i) * charge(i) + vint(i) / max(delt, 1.d-12)
            end do
            a(ncell) = 0.d0
            b(ncell) = 1.d0 / max(delt, 1.d-12)
            f(ncell) = E(ncell) * charge(ncell) + vint(ncell) / max(delt, 1.d-12)
            call gauss_solve(nmax, ncell, a, b, c, f, bx, fx, vint)
        end if

        do i = 1, ncell
            xt(i) = xt(i) + dt * vint(i)
            if (xt(i) < x0(0)) then
                xt(i) = 2.d0 * x0(0) - xt(i)
                vint(i) = -vint(i)
            end if
        end do

        do j = 2, ncell
            xxt = xt(j)
            xx0 = x0(j)
            xxtold = xtold(j)
            vv = v(j)
            vvint = vint(j)
            iidebut = idebut(j)
            qqiss = qiSS(j)
            nnhot = nhot(j)
            pphi = phi(j)
            EE = E(j)
            ggradnh = gradnh(j)
            ddxt = dxt(j)
            ccharge = charge(j)

            do i = j - 1, 1, -1
                if (xt(i) <= xxt) exit
                xt(i + 1) = xt(i)
                x0(i + 1) = x0(i)
                xtold(i + 1) = xtold(i)
                v(i + 1) = v(i)
                vint(i + 1) = vint(i)
                idebut(i + 1) = idebut(i)
                qiSS(i + 1) = qiSS(i)
                nhot(i + 1) = nhot(i)
                phi(i + 1) = phi(i)
                E(i + 1) = E(i)
                gradnh(i + 1) = gradnh(i)
                dxt(i + 1) = dxt(i)
                charge(i + 1) = charge(i)
            end do

            xt(i + 1) = xxt
            x0(i + 1) = xx0
            xtold(i + 1) = xxtold
            v(i + 1) = vv
            vint(i + 1) = vvint
            idebut(i + 1) = iidebut
            qiSS(i + 1) = qqiss
            nhot(i + 1) = nnhot
            phi(i + 1) = pphi
            E(i + 1) = EE
            gradnh(i + 1) = ggradnh
            dxt(i + 1) = ddxt
            charge(i + 1) = ccharge
        end do
    end do

    write(*, '("Final time = ", f8.4)') time
    write(*, '("Number of iterations = ", i0)') itime

    if (itime == 1) then
        do i = ncell + 1, ntotal
            x0(i) = xt(i)
        end do
    end if

    if (.not. (itmax == 1 .or. time <= 1.d-04)) then
        do i = 0, ntotal
            xi = xt(i) / max(sqrt(Th) * time, 1.d-12)
            if (xi < -1.d0) then
                nss(i) = n0hot
            else
                nss(i) = n0hot * exp(-(xi + 1.d0))
            end if
        end do
        do i = 0, ntotal
            difth(i) = (ni(i) - nss(i)) / max(nss(i), 1.d-12)
        end do
    end if

    open(unit = 11, file = 'profil.txt', status = 'replace')
    write(11, '(a)') '# x0 xt v phi E ni ne rho nhot nss'
    do i = 0, ncell, iline
        write(11, '(*(g0,:,","))') &
            x0(i), xt(i), v(i), phi(i), E(i), ni(i), ne(i), rho(i), nhot(i), nss(i)
    end do
    do i = ncell + 1, ntotal
        write(11, '(*(g0,:,","))') &
            x0(i), xt(i), v(i), phi(i), E(i), ni(i), ne(i), rho(i), nhot(i), nss(i)
    end do
    close(11)

    open(unit = 12, file = 'spectres.txt', status = 'replace')
    write(12, '(a)') '# vmoy Emoy dndv dndE'
    do i = 1, ncell
        vmoy(i) = 0.5d0 * (v(i - 1) + v(i))
        Emoy(i) = 0.25d0 * (v(i - 1)**2 + v(i)**2)
        dndv(i) = (niSS(i - 1) + niSS(i)) * dx0(i) / max(v(i) - v(i - 1), 1.d-12) / 2.d0
        dndv(i) = abs(dndv(i))
        dndE(i) = (niSS(i - 1) + niSS(i)) * dx0(i) / max(v(i)**2 - v(i - 1)**2, 1.d-12)
        dndE(i) = abs(dndE(i))
    end do
    do i = 1, ncell, iline
        if (ieee_is_finite(dndE(i)) .and. ieee_is_finite(dndv(i))) then
            write(12, '(*(g0,:,","))') vmoy(i), Emoy(i), dndv(i), dndE(i)
        end if
    end do
    close(12)

    close(9)
    close(10)
    write(*, *) 'time', time

contains

    subroutine gauss_solve(nmax_local, ncell_local, a_local, b_local, c_local, f_local, bx_local, fx_local, phi_local)
        implicit none
        integer, intent(in) :: nmax_local, ncell_local
        real(8), intent(in) :: a_local(1:nmax_local), b_local(0:nmax_local), c_local(0:nmax_local-1), f_local(0:nmax_local)
        real(8), intent(inout) :: bx_local(0:nmax_local), fx_local(0:nmax_local), phi_local(0:nmax_local)
        integer :: j_local

        bx_local(0) = b_local(0)
        fx_local(0) = f_local(0)
        if (abs(bx_local(0)) < 1.d-12) stop 'Division by zero in gauss_solve at j=0'

        do j_local = 1, ncell_local
            bx_local(j_local) = b_local(j_local) - a_local(j_local) * c_local(j_local - 1) / bx_local(j_local - 1)
            fx_local(j_local) = f_local(j_local) - a_local(j_local) * fx_local(j_local - 1) / bx_local(j_local - 1)
            if (abs(bx_local(j_local)) < 1.d-12) stop 'Division by zero in gauss_solve forward sweep'
        end do

        phi_local(ncell_local) = fx_local(ncell_local) / bx_local(ncell_local)
        do j_local = ncell_local - 1, 0, -1
            phi_local(j_local) = (fx_local(j_local) - c_local(j_local) * phi_local(j_local + 1)) / bx_local(j_local)
        end do
    end subroutine gauss_solve
end program IonBoost
