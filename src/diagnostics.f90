module mod_diagnostics
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mod_types, only: dp
    use mod_grid, only: grid_t
    use mod_ions, only: ion_species_t
    use mod_fields, only: fields_t
    implicit none
    private
    public :: diagnostics_t, open_histories, write_history, close_histories, &
              write_profiles, write_spectra

    type :: diagnostics_t
        integer :: ucons = -1, uhist = -1
    end type diagnostics_t

contains

    subroutine open_histories(diag)
        type(diagnostics_t), intent(inout) :: diag
        open(newunit=diag%ucons, file='conservation.txt', status='unknown')
        open(newunit=diag%uhist, file='historique.txt', status='unknown')
        write(diag%ucons, '(a)') '# time nti nthot nte n0hot En_ion Whot1 Whot2 Whot Th En_elec En_balance vmax vfinal'
        write(diag%uhist, '(a)') '# time xt(ivmax) v(ivmax) Energy_max E(ivmax) ne(ivmax) ni(ivmax) ' // &
            'ni(0) nhot(0) lDebye lgrad_e_pressure lgrad_i_pressure ivmax idebut(ivmax)'
    end subroutine open_histories

    subroutine write_history(diag, grid, ions, flds, time, n0hot, Th, &
                             nti, nthot, En_ion, En_elec, first_step)
        type(diagnostics_t), intent(in) :: diag
        type(grid_t), intent(in) :: grid
        type(ion_species_t), intent(in) :: ions
        type(fields_t), intent(in) :: flds
        real(dp), intent(in) :: time, n0hot, Th, nti, nthot, En_ion, En_elec
        logical, intent(in) :: first_step
        integer :: i, ivmax, n
        real(dp) :: vmax, vfinal, lDebye, grade, gradi, lgrade, lgradi

        n = grid%ncell
        vmax = 0.0_dp
        ivmax = n
        do i = 1, n
            if (ions%v(i) > vmax) then
                vmax = ions%v(i)
                ivmax = i
            end if
        end do
        vfinal = ions%v(ivmax) + flds%E(ivmax) * ions%Z(ivmax) * time
        lDebye = sqrt(Th / max(flds%ne(n), 1.0e-12_dp))

        ! edge log-gradient scale lengths of ne and ni (paper Eqs. 18 and 20)
        grade = (log(max(flds%ne(n), 1.0e-12_dp)) - log(max(flds%ne(n - 1), 1.0e-12_dp))) / &
                max(grid%xt(n) - grid%xt(n - 1), 1.0e-12_dp)
        lgrade = -1.0_dp / sign(max(abs(grade), 1.0e-12_dp), grade)
        if (first_step) then
            lgradi = 0.0_dp
        else
            gradi = (log(max(ions%dens(n), 1.0e-12_dp)) - log(max(ions%dens(n - 1), 1.0e-12_dp))) / &
                    max(grid%xt(n) - grid%xt(n - 1), 1.0e-12_dp)
            lgradi = -1.0_dp / sign(max(abs(gradi), 1.0e-12_dp), gradi)
        end if

        write(diag%ucons, '(*(g0,:,","))') time, nti, nthot, nthot, n0hot, En_ion, &
            flds%Whot1, flds%Whot2, flds%Whot, Th, En_elec, &
            En_elec + En_ion - flds%Whot, vmax, vfinal
        write(diag%uhist, '(*(g0,:,","))') time, grid%xt(ivmax), ions%v(ivmax), &
            0.5_dp * ions%v(ivmax)**2, flds%E(ivmax), flds%ne(ivmax), ions%dens(ivmax), &
            ions%dens(0), flds%nhot(0), lDebye, lgrade, lgradi, ivmax, ions%id0(ivmax)
    end subroutine write_history

    subroutine close_histories(diag)
        type(diagnostics_t), intent(inout) :: diag
        close(diag%ucons)
        close(diag%uhist)
    end subroutine close_histories

    subroutine write_profiles(grid, ions, flds, time, n0hot, Th, iline)
        type(grid_t), intent(in) :: grid
        type(ion_species_t), intent(in) :: ions
        type(fields_t), intent(in) :: flds
        real(dp), intent(in) :: time, n0hot, Th
        integer, intent(in) :: iline
        integer :: i, iunit
        real(dp) :: xi
        real(dp), allocatable :: nss(:), dens_full(:), v_full(:)

        ! self-similar reference density n0 exp(-(x/cs t + 1))
        allocate(nss(0:grid%ntotal), dens_full(0:grid%ntotal), v_full(0:grid%ntotal))
        nss = 0.0_dp
        if (time > 1.0e-4_dp) then
            do i = 0, grid%ntotal
                xi = grid%xt(i) / max(sqrt(Th) * time, 1.0e-12_dp)
                if (xi < -1.0_dp) then
                    nss(i) = n0hot
                else
                    nss(i) = n0hot * exp(-(xi + 1.0_dp))
                end if
            end do
        end if
        dens_full(0:ions%n) = ions%dens
        dens_full(ions%n + 1:grid%ntotal) = 1.0e-12_dp
        v_full(0:ions%n) = ions%v
        v_full(ions%n + 1:grid%ntotal) = 0.0_dp

        open(newunit=iunit, file='profil.txt', status='replace')
        write(iunit, '(a)') '# x0 xt v phi E ni ne rho nhot nss'
        do i = 0, ions%n, iline
            write(iunit, '(*(g0,:,","))') grid%x0(i), grid%xt(i), v_full(i), flds%phi(i), &
                flds%E(i), dens_full(i), flds%ne(i), flds%rho(i), flds%nhot(i), nss(i)
        end do
        do i = ions%n + 1, grid%ntotal
            write(iunit, '(*(g0,:,","))') grid%x0(i), grid%xt(i), v_full(i), flds%phi(i), &
                flds%E(i), dens_full(i), flds%ne(i), flds%rho(i), flds%nhot(i), nss(i)
        end do
        close(iunit)
    end subroutine write_profiles

    subroutine write_spectra(grid, ions, iline)
        type(grid_t), intent(in) :: grid
        type(ion_species_t), intent(in) :: ions
        integer, intent(in) :: iline
        integer :: i, iunit
        real(dp) :: vmoy, Emoy, dndv, dndE

        open(newunit=iunit, file='spectres.txt', status='replace')
        write(iunit, '(a)') '# vmoy Emoy dndv dndE'
        do i = 1, ions%n, iline
            vmoy = 0.5_dp * (ions%v(i - 1) + ions%v(i))
            Emoy = 0.25_dp * (ions%v(i - 1)**2 + ions%v(i)**2)
            dndv = abs((ions%n_init(i - 1) + ions%n_init(i)) * grid%dx0(i) / &
                       max(ions%v(i) - ions%v(i - 1), 1.0e-12_dp) / 2.0_dp)
            dndE = abs((ions%n_init(i - 1) + ions%n_init(i)) * grid%dx0(i) / &
                       max(ions%v(i)**2 - ions%v(i - 1)**2, 1.0e-12_dp))
            if (ieee_is_finite(dndv) .and. ieee_is_finite(dndE)) then
                write(iunit, '(*(g0,:,","))') vmoy, Emoy, dndv, dndE
            end if
        end do
        close(iunit)
    end subroutine write_spectra

end module mod_diagnostics
