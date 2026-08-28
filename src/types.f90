module mod_types
    implicit none
    private
    public :: dp
    integer, parameter :: dp = selected_real_kind(15, 307)   ! IEEE double
end module mod_types
