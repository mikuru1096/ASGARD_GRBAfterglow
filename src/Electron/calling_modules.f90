module get_Y
  use electron_radiation_kernel, only: besselk, get_syn, get_syn_simpson, get_syn_selected, get_nu_a, &
                                        get_nu_a_nonuniform, get_syn_adaptive
  use electron_cooling_kernel, only: get_forward_cooling, get_SSA_numerical, get_IC_numerical, &
                                      get_Y_Nakar, get_Y_Fan
  implicit none
  private

  public :: besselk, get_syn, get_syn_simpson, get_syn_selected, get_forward_cooling, get_nu_a
  public :: get_nu_a_nonuniform, get_syn_adaptive
  public :: get_SSA_numerical, get_IC_numerical, get_Y_Nakar, get_Y_Fan

end module get_Y
