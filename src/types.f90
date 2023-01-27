module mod_types
    implicit none
    private
    public :: rp, dp
    integer, parameter :: rp = selected_real_kind(p=24)
    integer, parameter :: dp = selected_real_kind(p=53)
end module mod_types
