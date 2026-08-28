module mod_electrons
    use mod_types, only: dp
    implicit none
    private
    public :: electron_species_t, e_density, e_ddensity_dphi, e_pressure

    ! Isothermal Boltzmann electron population: n(phi) = n0 exp(-phi/T).
    ! The solver only ever sees the summed closure functions below, so adding a
    ! population (e.g. cold electrons for a two-temperature expansion) means
    ! appending one entry to the species array in the driver. The vacuum-sheath
    ! extension in mod_fields currently assumes a single population (exact
    ! solution); generalize it before adding a second one.
    type :: electron_species_t
        real(dp) :: n0
        real(dp) :: T
    end type electron_species_t

contains

    pure function e_density(els, phi) result(ne)
        type(electron_species_t), intent(in) :: els(:)
        real(dp), intent(in) :: phi
        real(dp) :: ne
        integer :: s
        ne = 0.0_dp
        do s = 1, size(els)
            ne = ne + els(s)%n0 * exp(-phi / els(s)%T)
        end do
    end function e_density

    pure function e_ddensity_dphi(els, phi) result(dne)
        type(electron_species_t), intent(in) :: els(:)
        real(dp), intent(in) :: phi
        real(dp) :: dne
        integer :: s
        dne = 0.0_dp
        do s = 1, size(els)
            dne = dne - els(s)%n0 / els(s)%T * exp(-phi / els(s)%T)
        end do
    end function e_ddensity_dphi

    pure function e_pressure(els, phi) result(pe)
        type(electron_species_t), intent(in) :: els(:)
        real(dp), intent(in) :: phi
        real(dp) :: pe
        integer :: s
        pe = 0.0_dp
        do s = 1, size(els)
            pe = pe + els(s)%n0 * els(s)%T * exp(-phi / els(s)%T)
        end do
    end function e_pressure

end module mod_electrons
