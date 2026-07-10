module reverse_shock_state
    implicit none
    private
    public :: wait_phase, precross_phase, postcross_phase
    public :: rs_db3, rs_tcross, rs_rcross, rs_e3cross, rs_gam20, rs_u3cross
    public :: rs_v3cross, rs_m3cross, rs_gammcross, rs_b3ordered, rs_nstate

    integer, parameter :: wait_phase = -1, precross_phase = 1
    integer, parameter :: postcross_phase = 2
    integer, parameter :: rs_db3 = 1, rs_tcross = 2, rs_rcross = 3, rs_e3cross = 4
    integer, parameter :: rs_gam20 = 5, rs_u3cross = 6, rs_v3cross = 7
    integer, parameter :: rs_m3cross = 8, rs_gammcross = 9, rs_b3ordered = 10
    integer, parameter :: rs_nstate = 10
end module reverse_shock_state
