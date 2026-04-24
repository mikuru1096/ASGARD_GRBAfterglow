!f2py: skip
module hadronic_radiation_kernel
    use constants
    use hadronic_common
    implicit none

contains

real(8) function hadronic_syn_kernel_ultrarel(V_nu,gam_p,B_field_g)
    implicit real(8)(A-H,O-Z)
    real(8), intent(in) :: V_nu,gam_p,B_field_g
    real(8), parameter :: B_crit=4.41d13
    real(8) :: mass_p_gev,mass_e_gev,energy_photon_gev,energy_proton_gev
    real(8) :: b_dimless,mass_ratio,xbar,x_arg,y,poly

    mass_p_gev=Para_m_p_E*Para_erg2eV*1d-9
    mass_e_gev=Para_m_energy*Para_erg2eV*1d-9
    energy_photon_gev=Para_h_GeV*V_nu
    energy_proton_gev=gam_p*mass_p_gev
    b_dimless=B_field_g/B_crit
    mass_ratio=mass_e_gev/mass_p_gev
    xbar=energy_photon_gev*mass_p_gev/(3d0*energy_proton_gev*energy_proton_gev*b_dimless*mass_ratio*mass_ratio)
    x_arg=two*xbar

    if (x_arg < 1d-2) then
        hadronic_syn_kernel_ultrarel=1.80842d0*xbar**(1d0/3d0)*two**(-2d0/3d0)
    else if (x_arg < one) then
        y=dlog10(x_arg)
        poly=-0.35775237d0 &
            -0.83695385d0*y &
            -1.1449608d0*y*y &
            -0.68137283d0*y**3 &
            -0.22754737d0*y**4 &
            -0.031967334d0*y**5
        hadronic_syn_kernel_ultrarel=ten**poly/two
    else if (x_arg < ten) then
        y=dlog10(x_arg)
        poly=-0.35842494d0 &
            -0.79652041d0*y &
            -1.6113032d0*y*y &
            +0.26055213d0*y**3 &
            -1.6979017d0*y**4 &
            +0.032955035d0*y**5
        hadronic_syn_kernel_ultrarel=ten**poly/two
    else if (x_arg < 1d2) then
        hadronic_syn_kernel_ultrarel=(pi/4d0)*dexp(-x_arg)*(one-99d0/(162d0*x_arg))
    else
        hadronic_syn_kernel_ultrarel=zero
    end if
end function hadronic_syn_kernel_ultrarel

subroutine hadronic_get_proton_syn_state(R_loc,B_field_g,Num_gam_p,Num_nu,gam_p,dN_gam_p,V_seed,P_had_syn,Seed_had_syn)
    implicit real(8)(A-H,O-Z)
    integer, intent(in) :: Num_gam_p,Num_nu
    real(8), intent(in) :: R_loc,B_field_g,gam_p(Num_gam_p),dN_gam_p(Num_gam_p),V_seed(Num_nu)
    real(8), intent(out) :: P_had_syn(Num_nu),Seed_had_syn(Num_nu)
    integer :: I_nu,I_gam
    real(8) :: temp_syn,Rariv2,V_cal,dInteg,Fx
    real(8) :: dN_mid(Num_gam_p),dln_gam(Num_gam_p),temp_para,gam_mid(Num_gam_p)

    temp_syn=dsqrt(3d0)*Para_e*Para_e*Para_e/Para_m_p_E
    Rariv2=R_loc*R_loc
    P_had_syn=zero
    Seed_had_syn=zero

    do I_gam=1,Num_gam_p-1
        gam_mid(I_gam)=dsqrt(gam_p(I_gam)*gam_p(I_gam+1))
        dN_mid(I_gam)=(dN_gam_p(I_gam)+dN_gam_p(I_gam+1))/two
        dln_gam(I_gam)=dlog(gam_p(I_gam+1)/gam_p(I_gam))
    end do

    do I_nu=1,Num_nu
        V_cal=V_seed(I_nu)
        dInteg=zero
        do I_gam=1,Num_gam_p-1
            Fx=hadronic_syn_kernel_ultrarel(V_cal,gam_mid(I_gam),B_field_g)
            dInteg=dInteg+dN_mid(I_gam)*Fx*gam_mid(I_gam)*dln_gam(I_gam)
        end do
        P_had_syn(I_nu)=temp_syn*B_field_g*dInteg
        Seed_had_syn(I_nu)=P_had_syn(I_nu)/(Rariv2*max(V_cal,1d-60))
    end do

    temp_para=4d0*pi*Para_c*Para_h
    Seed_had_syn=Seed_had_syn/temp_para
end subroutine hadronic_get_proton_syn_state

end module hadronic_radiation_kernel
