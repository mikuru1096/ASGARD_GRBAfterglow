!f2py: skip
module hadronic_secondary_radiation_kernel
    use constants
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    real(8), parameter :: am3_c_cgs = Para_c
    real(8), parameter :: am3_sigma_t_cgs = Para_sigmaT
    real(8), parameter :: am3_mass_electron_gev = Para_m_energy*Para_erg2eV*1d-9
    real(8), parameter :: am3_mass_pion_charged_gev = 1.396d-1
    real(8), parameter :: am3_mass_muon_gev = 1.0566d-1
    real(8), parameter :: am3_b_crit_gauss = 4.41d13
    real(8), parameter :: am3_erg_to_gev = 6.24d2

    public :: hadronic_secondary_initialize_synchrotron_kernel
    public :: hadronic_secondary_initialize_inverse_compton_kernel
    public :: hadronic_secondary_pion_synchrotron_rate
    public :: hadronic_secondary_muon_synchrotron_rate
    public :: hadronic_secondary_pion_inverse_compton_rate
    public :: hadronic_secondary_muon_inverse_compton_rate
    public :: hadronic_secondary_radiation_operator

contains

subroutine hadronic_secondary_radiation_operator(num_had,hadron_energy_gev,num_ph,photon_energy_gev, &
                                                 pion_plus_per_gev,pion_minus_per_gev, &
                                                 muon_minus_left_per_gev,muon_minus_right_per_gev, &
                                                 muon_plus_left_per_gev,muon_plus_right_per_gev, &
                                                 photons_on_had_grid_per_gev, &
                                                 ind_min_energy_pho_hadgrid,magnetic_field_g, &
                                                 pion_synch_rate_per_gev,muon_synch_rate_per_gev, &
                                                 pion_ic_rate_per_gev,muon_ic_rate_per_gev, &
                                                 synch_dln_energy,kernel_pion,kernel_muon, &
                                                 ic_dln_energy,delta_e_pi,jmax_pi,delta_e_mu,jmax_mu)
    integer, intent(in) :: num_had,num_ph,ind_min_energy_pho_hadgrid
    real(8), intent(in) :: hadron_energy_gev(num_had),photon_energy_gev(num_ph)
    real(8), intent(in) :: pion_plus_per_gev(num_had),pion_minus_per_gev(num_had)
    real(8), intent(in) :: muon_minus_left_per_gev(num_had),muon_minus_right_per_gev(num_had)
    real(8), intent(in) :: muon_plus_left_per_gev(num_had),muon_plus_right_per_gev(num_had)
    real(8), intent(in) :: photons_on_had_grid_per_gev(num_ph),magnetic_field_g
    real(8), intent(out) :: pion_synch_rate_per_gev(num_ph),muon_synch_rate_per_gev(num_ph)
    real(8), intent(out) :: pion_ic_rate_per_gev(num_ph),muon_ic_rate_per_gev(num_ph)
    real(8), intent(out) :: synch_dln_energy,ic_dln_energy
    real(8), intent(out) :: kernel_pion(num_ph,num_had),kernel_muon(num_ph,num_had)
    integer, intent(out) :: delta_e_pi(num_had),jmax_pi(num_had),delta_e_mu(num_had),jmax_mu(num_had)

    call hadronic_secondary_initialize_synchrotron_kernel(num_had,hadron_energy_gev,num_ph,photon_energy_gev, &
                                                          magnetic_field_g,synch_dln_energy,kernel_pion,kernel_muon)
    call hadronic_secondary_initialize_inverse_compton_kernel(num_had,hadron_energy_gev,num_ph,photon_energy_gev, &
                                                              ind_min_energy_pho_hadgrid,ic_dln_energy, &
                                                              delta_e_pi,jmax_pi,delta_e_mu,jmax_mu)

    call hadronic_secondary_apply_radiation_kernels(num_had,pion_plus_per_gev,pion_minus_per_gev, &
                                                    muon_minus_left_per_gev,muon_minus_right_per_gev, &
                                                    muon_plus_left_per_gev,muon_plus_right_per_gev, &
                                                    num_ph,photons_on_had_grid_per_gev, &
                                                    synch_dln_energy,kernel_pion,kernel_muon, &
                                                    ic_dln_energy,delta_e_pi,jmax_pi,delta_e_mu,jmax_mu, &
                                                    pion_synch_rate_per_gev,muon_synch_rate_per_gev, &
                                                    pion_ic_rate_per_gev,muon_ic_rate_per_gev)
end subroutine hadronic_secondary_radiation_operator

subroutine hadronic_secondary_apply_radiation_kernels(num_had,pion_plus_per_gev,pion_minus_per_gev, &
                                                      muon_minus_left_per_gev,muon_minus_right_per_gev, &
                                                      muon_plus_left_per_gev,muon_plus_right_per_gev, &
                                                      num_ph,photons_on_had_grid_per_gev, &
                                                      synch_dln_energy,kernel_pion,kernel_muon, &
                                                      ic_dln_energy,delta_e_pi,jmax_pi,delta_e_mu,jmax_mu, &
                                                      pion_synch_rate_per_gev,muon_synch_rate_per_gev, &
                                                      pion_ic_rate_per_gev,muon_ic_rate_per_gev)
    integer, intent(in) :: num_had,num_ph
    integer, intent(in) :: delta_e_pi(num_had),jmax_pi(num_had),delta_e_mu(num_had),jmax_mu(num_had)
    real(8), intent(in) :: pion_plus_per_gev(num_had),pion_minus_per_gev(num_had)
    real(8), intent(in) :: muon_minus_left_per_gev(num_had),muon_minus_right_per_gev(num_had)
    real(8), intent(in) :: muon_plus_left_per_gev(num_had),muon_plus_right_per_gev(num_had)
    real(8), intent(in) :: photons_on_had_grid_per_gev(num_ph)
    real(8), intent(in) :: synch_dln_energy,ic_dln_energy
    real(8), intent(in) :: kernel_pion(num_ph,num_had),kernel_muon(num_ph,num_had)
    real(8), intent(out) :: pion_synch_rate_per_gev(num_ph),muon_synch_rate_per_gev(num_ph)
    real(8), intent(out) :: pion_ic_rate_per_gev(num_ph),muon_ic_rate_per_gev(num_ph)
    real(8) :: pion_total(num_had),muon_total(num_had),coeff_pi,coeff_mu

    call hadronic_secondary_validate_density(num_had,pion_plus_per_gev,"pion_plus_per_gev")
    call hadronic_secondary_validate_density(num_had,pion_minus_per_gev,"pion_minus_per_gev")
    pion_total = pion_plus_per_gev + pion_minus_per_gev
    pion_synch_rate_per_gev = synch_dln_energy*matmul(kernel_pion,pion_total)

    call hadronic_secondary_validate_density(num_had,muon_minus_left_per_gev,"muon_minus_left_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_minus_right_per_gev,"muon_minus_right_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_plus_left_per_gev,"muon_plus_left_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_plus_right_per_gev,"muon_plus_right_per_gev")
    muon_total = muon_minus_left_per_gev + muon_minus_right_per_gev + muon_plus_left_per_gev + &
                 muon_plus_right_per_gev
    muon_synch_rate_per_gev = synch_dln_energy*matmul(kernel_muon,muon_total)

    call hadronic_secondary_validate_density(num_ph,photons_on_had_grid_per_gev, &
                                             "photons_on_had_grid_per_gev")
    call hadronic_secondary_validate_density(num_had,pion_plus_per_gev,"pion_plus_per_gev")
    call hadronic_secondary_validate_density(num_had,pion_minus_per_gev,"pion_minus_per_gev")
    coeff_pi = hadronic_secondary_ic_coeff(am3_mass_pion_charged_gev)
    call hadronic_secondary_compute_ic_channel(num_ph,photons_on_had_grid_per_gev,num_had,pion_total, &
                                               delta_e_pi,jmax_pi,ic_dln_energy,coeff_pi,pion_ic_rate_per_gev)

    call hadronic_secondary_validate_density(num_ph,photons_on_had_grid_per_gev, &
                                             "photons_on_had_grid_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_minus_left_per_gev,"muon_minus_left_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_minus_right_per_gev,"muon_minus_right_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_plus_left_per_gev,"muon_plus_left_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_plus_right_per_gev,"muon_plus_right_per_gev")
    coeff_mu = hadronic_secondary_ic_coeff(am3_mass_muon_gev)
    call hadronic_secondary_compute_ic_channel(num_ph,photons_on_had_grid_per_gev,num_had,muon_total, &
                                               delta_e_mu,jmax_mu,ic_dln_energy,coeff_mu,muon_ic_rate_per_gev)
end subroutine hadronic_secondary_apply_radiation_kernels

subroutine hadronic_secondary_initialize_synchrotron_kernel(num_had,hadron_energy_gev,num_ph,photon_energy_gev, &
                                                            magnetic_field_g,dln_energy,kernel_pion,kernel_muon)
    integer, intent(in) :: num_had,num_ph
    real(8), intent(in) :: hadron_energy_gev(num_had),photon_energy_gev(num_ph),magnetic_field_g
    real(8), intent(out) :: dln_energy,kernel_pion(num_ph,num_had),kernel_muon(num_ph,num_had)
    integer :: i,j

    call hadronic_secondary_validate_positive_log_grid(num_had,hadron_energy_gev, &
                                                       "hadron_energy_gev",dln_energy)
    call hadronic_secondary_validate_positive_log_grid(num_ph,photon_energy_gev, &
                                                       "photon_energy_gev")
    if (magnetic_field_g <= zero) then
        error stop "secondary synchrotron requires magnetic_field_g > 0."
    end if

    do i=1,num_ph
        do j=1,num_had
            kernel_pion(i,j) = hadronic_secondary_syn_kernel_ultrarel(photon_energy_gev(i), &
                                                                      hadron_energy_gev(j), &
                                                                      am3_mass_pion_charged_gev, &
                                                                      magnetic_field_g)
            kernel_muon(i,j) = hadronic_secondary_syn_kernel_ultrarel(photon_energy_gev(i), &
                                                                      hadron_energy_gev(j), &
                                                                      am3_mass_muon_gev, &
                                                                      magnetic_field_g)
        end do
    end do
end subroutine hadronic_secondary_initialize_synchrotron_kernel

subroutine hadronic_secondary_initialize_inverse_compton_kernel(num_had,hadron_energy_gev, &
                                                                num_ph,photon_energy_gev, &
                                                                ind_min_energy_pho_hadgrid,dln_energy, &
                                                                delta_e_pi,jmax_pi,delta_e_mu,jmax_mu)
    integer, intent(in) :: num_had,num_ph,ind_min_energy_pho_hadgrid
    real(8), intent(in) :: hadron_energy_gev(num_had),photon_energy_gev(num_ph)
    real(8), intent(out) :: dln_energy
    integer, intent(out) :: delta_e_pi(num_had),jmax_pi(num_had),delta_e_mu(num_had),jmax_mu(num_had)
    real(8) :: dln_had,dln_ph

    call hadronic_secondary_validate_positive_log_grid(num_had,hadron_energy_gev, &
                                                       "hadron_energy_gev",dln_had)
    call hadronic_secondary_validate_positive_log_grid(num_ph,photon_energy_gev, &
                                                       "photon_energy_gev",dln_ph)
    if (dabs(dln_had-dln_ph) > dmax1(1d-12,1d-10*dabs(dln_had))) then
        error stop "secondary IC requires hadron/photon grids with one shared logarithmic spacing."
    end if

    dln_energy = dln_had
    call hadronic_secondary_build_ic_species_kernel(num_had,hadron_energy_gev,dln_energy,num_ph, &
                                                    ind_min_energy_pho_hadgrid,am3_mass_pion_charged_gev, &
                                                    delta_e_pi,jmax_pi)
    call hadronic_secondary_build_ic_species_kernel(num_had,hadron_energy_gev,dln_energy,num_ph, &
                                                    ind_min_energy_pho_hadgrid,am3_mass_muon_gev, &
                                                    delta_e_mu,jmax_mu)
end subroutine hadronic_secondary_initialize_inverse_compton_kernel

subroutine hadronic_secondary_pion_synchrotron_rate(num_had,pion_plus_per_gev,pion_minus_per_gev, &
                                                    num_ph,dln_energy,kernel_pion,pion_synch_rate_per_gev)
    integer, intent(in) :: num_had,num_ph
    real(8), intent(in) :: pion_plus_per_gev(num_had),pion_minus_per_gev(num_had)
    real(8), intent(in) :: dln_energy,kernel_pion(num_ph,num_had)
    real(8), intent(out) :: pion_synch_rate_per_gev(num_ph)
    real(8) :: pion_total(num_had)

    call hadronic_secondary_validate_density(num_had,pion_plus_per_gev,"pion_plus_per_gev")
    call hadronic_secondary_validate_density(num_had,pion_minus_per_gev,"pion_minus_per_gev")
    pion_total = pion_plus_per_gev + pion_minus_per_gev
    pion_synch_rate_per_gev = dln_energy*matmul(kernel_pion,pion_total)
end subroutine hadronic_secondary_pion_synchrotron_rate

subroutine hadronic_secondary_muon_synchrotron_rate(num_had,muon_minus_left_per_gev, &
                                                    muon_minus_right_per_gev,muon_plus_left_per_gev, &
                                                    muon_plus_right_per_gev,num_ph,dln_energy,kernel_muon, &
                                                    muon_synch_rate_per_gev)
    integer, intent(in) :: num_had,num_ph
    real(8), intent(in) :: muon_minus_left_per_gev(num_had),muon_minus_right_per_gev(num_had)
    real(8), intent(in) :: muon_plus_left_per_gev(num_had),muon_plus_right_per_gev(num_had)
    real(8), intent(in) :: dln_energy,kernel_muon(num_ph,num_had)
    real(8), intent(out) :: muon_synch_rate_per_gev(num_ph)
    real(8) :: muon_total(num_had)

    call hadronic_secondary_validate_density(num_had,muon_minus_left_per_gev,"muon_minus_left_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_minus_right_per_gev,"muon_minus_right_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_plus_left_per_gev,"muon_plus_left_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_plus_right_per_gev,"muon_plus_right_per_gev")
    muon_total = muon_minus_left_per_gev + muon_minus_right_per_gev + muon_plus_left_per_gev + &
                 muon_plus_right_per_gev
    muon_synch_rate_per_gev = dln_energy*matmul(kernel_muon,muon_total)
end subroutine hadronic_secondary_muon_synchrotron_rate

subroutine hadronic_secondary_pion_inverse_compton_rate(num_ph,photons_on_had_grid_per_gev,num_had, &
                                                        pion_plus_per_gev,pion_minus_per_gev,dln_energy, &
                                                        delta_e_pi,jmax_pi,pion_ic_rate_per_gev)
    integer, intent(in) :: num_ph,num_had,delta_e_pi(num_had),jmax_pi(num_had)
    real(8), intent(in) :: photons_on_had_grid_per_gev(num_ph),pion_plus_per_gev(num_had)
    real(8), intent(in) :: pion_minus_per_gev(num_had),dln_energy
    real(8), intent(out) :: pion_ic_rate_per_gev(num_ph)
    real(8) :: pion_total(num_had),coeff_pi

    call hadronic_secondary_validate_density(num_ph,photons_on_had_grid_per_gev, &
                                             "photons_on_had_grid_per_gev")
    call hadronic_secondary_validate_density(num_had,pion_plus_per_gev,"pion_plus_per_gev")
    call hadronic_secondary_validate_density(num_had,pion_minus_per_gev,"pion_minus_per_gev")
    pion_total = pion_plus_per_gev + pion_minus_per_gev
    coeff_pi = hadronic_secondary_ic_coeff(am3_mass_pion_charged_gev)
    call hadronic_secondary_compute_ic_channel(num_ph,photons_on_had_grid_per_gev,num_had,pion_total, &
                                               delta_e_pi,jmax_pi,dln_energy,coeff_pi,pion_ic_rate_per_gev)
end subroutine hadronic_secondary_pion_inverse_compton_rate

subroutine hadronic_secondary_muon_inverse_compton_rate(num_ph,photons_on_had_grid_per_gev,num_had, &
                                                        muon_minus_left_per_gev,muon_minus_right_per_gev, &
                                                        muon_plus_left_per_gev,muon_plus_right_per_gev, &
                                                        dln_energy,delta_e_mu,jmax_mu,muon_ic_rate_per_gev)
    integer, intent(in) :: num_ph,num_had,delta_e_mu(num_had),jmax_mu(num_had)
    real(8), intent(in) :: photons_on_had_grid_per_gev(num_ph),dln_energy
    real(8), intent(in) :: muon_minus_left_per_gev(num_had),muon_minus_right_per_gev(num_had)
    real(8), intent(in) :: muon_plus_left_per_gev(num_had),muon_plus_right_per_gev(num_had)
    real(8), intent(out) :: muon_ic_rate_per_gev(num_ph)
    real(8) :: muon_total(num_had),coeff_mu

    call hadronic_secondary_validate_density(num_ph,photons_on_had_grid_per_gev, &
                                             "photons_on_had_grid_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_minus_left_per_gev,"muon_minus_left_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_minus_right_per_gev,"muon_minus_right_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_plus_left_per_gev,"muon_plus_left_per_gev")
    call hadronic_secondary_validate_density(num_had,muon_plus_right_per_gev,"muon_plus_right_per_gev")
    muon_total = muon_minus_left_per_gev + muon_minus_right_per_gev + muon_plus_left_per_gev + &
                 muon_plus_right_per_gev
    coeff_mu = hadronic_secondary_ic_coeff(am3_mass_muon_gev)
    call hadronic_secondary_compute_ic_channel(num_ph,photons_on_had_grid_per_gev,num_had,muon_total, &
                                               delta_e_mu,jmax_mu,dln_energy,coeff_mu,muon_ic_rate_per_gev)
end subroutine hadronic_secondary_muon_inverse_compton_rate

real(8) function hadronic_secondary_syn_kernel_ultrarel(photon_energy_gev,particle_energy_gev, &
                                                        particle_mass_gev,magnetic_field_g)
    real(8), intent(in) :: photon_energy_gev,particle_energy_gev,particle_mass_gev,magnetic_field_g
    real(8) :: norm,b_dimless,mass_ratio,xbar,two_xbar,y,poly

    norm = (3d0*dsqrt(3d0)/pi) * am3_sigma_t_cgs * am3_c_cgs * magnetic_field_g * &
           am3_b_crit_gauss/(8d0*pi) * am3_erg_to_gev/particle_mass_gev
    b_dimless = magnetic_field_g/am3_b_crit_gauss
    mass_ratio = am3_mass_electron_gev/particle_mass_gev
    xbar = photon_energy_gev*particle_mass_gev / &
           (3d0*particle_energy_gev*particle_energy_gev*b_dimless*mass_ratio*mass_ratio)
    two_xbar = 2d0*xbar

    if (two_xbar < 1d-2) then
        hadronic_secondary_syn_kernel_ultrarel = norm*1.80842d0*xbar**(1d0/3d0)*2d0**(-2d0/3d0)
    else if (two_xbar < one) then
        y = dlog10(two_xbar)
        poly = -0.35775237d0 - 0.83695385d0*y - 1.1449608d0*y*y - 0.68137283d0*y**3 - &
               0.22754737d0*y**4 - 0.031967334d0*y**5
        hadronic_secondary_syn_kernel_ultrarel = norm*(ten**poly)/2d0
    else if (two_xbar < 1d1) then
        y = dlog10(two_xbar)
        poly = -0.35842494d0 - 0.79652041d0*y - 1.6113032d0*y*y + 0.26055213d0*y**3 - &
               1.6979017d0*y**4 + 0.032955035d0*y**5
        hadronic_secondary_syn_kernel_ultrarel = norm*(ten**poly)/2d0
    else if (two_xbar < 1d2) then
        hadronic_secondary_syn_kernel_ultrarel = norm*(pi/4d0)*dexp(-two_xbar)*(one-99d0/(162d0*two_xbar))
    else
        hadronic_secondary_syn_kernel_ultrarel = zero
    end if
end function hadronic_secondary_syn_kernel_ultrarel

subroutine hadronic_secondary_build_ic_species_kernel(num_had,hadron_energy_gev,dln_energy,num_ph, &
                                                      ind_min_energy_pho_hadgrid,mass_gev,delta_e,jmax)
    integer, intent(in) :: num_had,num_ph,ind_min_energy_pho_hadgrid
    real(8), intent(in) :: hadron_energy_gev(num_had),dln_energy,mass_gev
    integer, intent(out) :: delta_e(num_had),jmax(num_had)
    integer :: i,candidate,jmax1
    real(8) :: gamma

    do i=1,num_had
        gamma = hadron_energy_gev(i)/mass_gev
        delta_e(i) = int(dlog(gamma*gamma)/dln_energy)
        jmax1 = int(dlog(0.5d0*mass_gev/gamma)/dln_energy) + ind_min_energy_pho_hadgrid
        candidate = jmax1 + delta_e(i)
        if (candidate > num_ph) then
            jmax(i) = num_ph
        else
            jmax(i) = candidate
        end if
    end do
end subroutine hadronic_secondary_build_ic_species_kernel

subroutine hadronic_secondary_compute_ic_channel(num_ph,photons_on_had_grid_per_gev,num_had, &
                                                 hadron_density_per_gev,delta_e,jmax,dln_energy, &
                                                 coeff_cgs,rate_per_gev)
    integer, intent(in) :: num_ph,num_had,delta_e(num_had),jmax(num_had)
    real(8), intent(in) :: photons_on_had_grid_per_gev(num_ph),hadron_density_per_gev(num_had)
    real(8), intent(in) :: dln_energy,coeff_cgs
    real(8), intent(out) :: rate_per_gev(num_ph)
    integer :: i,j,j0,src_idx
    real(8) :: z

    rate_per_gev = zero
    do j=1,num_ph
        j0 = j - 1
        z = zero
        do i=1,num_had
            if (j0 < delta_e(i) .or. j0 > jmax(i)) cycle
            src_idx = j - delta_e(i)
            if (src_idx < 1 .or. src_idx > num_ph) then
                error stop "secondary IC kernel maps to an out-of-grid photon source index."
            end if
            z = z + photons_on_had_grid_per_gev(src_idx)*hadron_density_per_gev(i)
        end do
        rate_per_gev(j) = z*dln_energy*coeff_cgs
    end do
end subroutine hadronic_secondary_compute_ic_channel

real(8) function hadronic_secondary_ic_coeff(mass_gev)
    real(8), intent(in) :: mass_gev
    real(8) :: mass_ratio

    mass_ratio = mass_gev/am3_mass_electron_gev
    hadronic_secondary_ic_coeff = am3_c_cgs*am3_sigma_t_cgs/(mass_ratio*mass_ratio)
end function hadronic_secondary_ic_coeff

subroutine hadronic_secondary_validate_positive_log_grid(num_grid,energy_grid,name,dln_energy)
    integer, intent(in) :: num_grid
    character(*), intent(in) :: name
    real(8), intent(in) :: energy_grid(num_grid)
    real(8), intent(out), optional :: dln_energy
    integer :: i
    real(8) :: dln_local,dln_i

    if (num_grid < 2) then
        error stop "secondary radiation requires at least two grid points."
    end if
    if (.not. ieee_is_finite(energy_grid(1)) .or. energy_grid(1) <= zero) then
        error stop "secondary radiation requires positive finite grid energies."
    end if

    if (.not. ieee_is_finite(energy_grid(2)) .or. energy_grid(2) <= zero) then
        error stop "secondary radiation requires positive finite grid energies."
    end if
    if (energy_grid(2) <= energy_grid(1)) then
        error stop "secondary radiation requires strictly increasing grids."
    end if

    dln_local = dlog(energy_grid(2)/energy_grid(1))
    do i=2,num_grid
        if (.not. ieee_is_finite(energy_grid(i)) .or. energy_grid(i) <= zero) then
            error stop "secondary radiation requires positive finite grid energies."
        end if
        if (energy_grid(i) <= energy_grid(i-1)) then
            error stop "secondary radiation requires strictly increasing grids."
        end if
        if (i >= 3) then
            dln_i = dlog(energy_grid(i)/energy_grid(i-1))
            if (dabs(dln_i-dln_local) > dmax1(1d-12,1d-6*dabs(dln_local))) then
                error stop "secondary radiation requires logarithmically uniform grids."
            end if
        end if
    end do

    if (present(dln_energy)) dln_energy = dln_local
end subroutine hadronic_secondary_validate_positive_log_grid

subroutine hadronic_secondary_validate_density(num_grid,values,name)
    integer, intent(in) :: num_grid
    character(*), intent(in) :: name
    real(8), intent(in) :: values(num_grid)
    integer :: i

    do i=1,num_grid
        if (.not. ieee_is_finite(values(i))) then
            error stop "secondary radiation requires finite densities."
        end if
    end do
end subroutine hadronic_secondary_validate_density

end module hadronic_secondary_radiation_kernel
