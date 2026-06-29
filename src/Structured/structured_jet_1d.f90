subroutine structured_jet_flux_1d(Boundary,E_iso_grid,Gamma0_grid,active_grid,V_seed,V_obs,Tobs, &
                                  n,Num_theta_patch,Num_phi_patch,Num_nu,Num_nu_obs,Num_Tobs,Num_R,Num_gam_e, &
                                  Num_theta_sed,Num_phi_sed,index_dyn,index_Y,index_syn_intger, &
                                  include_reverse_sync,include_forward_ssc,include_hadronic,include_proton_synch, &
                                  include_pg,include_neutrino,num_gam_p,num_nu_nu, &
                                  reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e,reverse_sigma, &
                                  hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc,axisymmetric, &
                                  n_threads_outer,n_threads_inner,adaptive_substeps,substep_rtol, &
                                  substep_min,substep_max,thermal_electrons,electron_solver_id, &
                                  fwd_sync_obs,fwd_ssc_obs,fwd_hadronic_obs,rev_sync_obs,total_obs, &
                                  track_tobs,track_gamma,track_radius,track_mass,track_bfield,track_nu_m,track_nu_c,track_nu_a)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_external_density_profile, density_jump_max
    implicit none
    integer, intent(in) :: n,Num_theta_patch,Num_phi_patch,Num_nu,Num_nu_obs,Num_Tobs,Num_R,Num_gam_e
    integer, intent(in) :: Num_theta_sed,Num_phi_sed,index_dyn,index_Y,index_syn_intger
    integer, intent(in) :: include_reverse_sync,include_forward_ssc,include_hadronic,include_proton_synch
    integer, intent(in) :: include_pg,include_neutrino,num_gam_p,num_nu_nu
    integer, intent(in) :: axisymmetric,n_threads_outer,n_threads_inner,adaptive_substeps,substep_min,substep_max
    integer, intent(in) :: thermal_electrons,electron_solver_id,active_grid(Num_theta_patch,Num_phi_patch)
    real(8), intent(in) :: Boundary(n),E_iso_grid(Num_theta_patch,Num_phi_patch),Gamma0_grid(Num_theta_patch,Num_phi_patch)
    real(8), intent(in) :: V_seed(Num_nu),V_obs(Num_nu_obs),Tobs(Num_Tobs),substep_rtol
    real(8), intent(in) :: reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e,reverse_sigma
    real(8), intent(in) :: hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc
    real(8), intent(out) :: fwd_sync_obs(Num_nu_obs,Num_Tobs),fwd_ssc_obs(Num_nu_obs,Num_Tobs)
    real(8), intent(out) :: fwd_hadronic_obs(Num_nu_obs,Num_Tobs),rev_sync_obs(Num_nu_obs,Num_Tobs),total_obs(Num_nu_obs,Num_Tobs)
    real(8), intent(out) :: track_tobs(Num_R),track_gamma(Num_R),track_radius(Num_R),track_mass(Num_R),track_bfield(Num_R)
    real(8), intent(out) :: track_nu_m(Num_R),track_nu_c(Num_R),track_nu_a(Num_R)

    real(8) :: Boundary_sed(n)
    integer :: track_set,n_threads_projection

    fwd_sync_obs=zero; fwd_ssc_obs=zero; fwd_hadronic_obs=zero; rev_sync_obs=zero; total_obs=zero
    track_tobs=zero; track_gamma=zero; track_radius=zero; track_mass=zero; track_bfield=zero
    track_nu_m=zero; track_nu_c=zero; track_nu_a=zero; track_set=0
    n_threads_projection=max(n_threads_outer,n_threads_inner)
    Boundary_sed=Boundary

    !$ call omp_set_dynamic(.false.)
    if (n_threads_outer > 1 .and. n_threads_inner > 1) then
        !$ call omp_set_max_active_levels(2)
    else
        !$ call omp_set_max_active_levels(1)
    end if

    if (axisymmetric /= 0) then
        call run_axisymmetric()
    else
        call run_nonaxisymmetric()
    end if

    total_obs=fwd_sync_obs+fwd_ssc_obs+fwd_hadronic_obs+rev_sync_obs

contains

    subroutine run_axisymmetric()
        implicit none
        real(8), allocatable :: rt_axis(:,:),rg_axis(:,:),rr_axis(:,:),sync_axis(:,:,:),ssc_axis(:,:,:)
        real(8), allocatable :: had_axis(:,:,:),rev_axis(:,:,:)

        allocate(rt_axis(Num_R,Num_theta_patch),rg_axis(Num_R,Num_theta_patch),rr_axis(Num_R,Num_theta_patch), &
                 sync_axis(Num_nu,Num_R,Num_theta_patch),ssc_axis(Num_nu,Num_R,Num_theta_patch), &
                 had_axis(Num_nu,Num_R,Num_theta_patch),rev_axis(Num_nu,Num_R,Num_theta_patch))
        rt_axis=zero; rg_axis=one; rr_axis=zero; sync_axis=zero; ssc_axis=zero; had_axis=zero; rev_axis=zero
        call structured_solve_axisymmetric(Boundary,E_iso_grid,Gamma0_grid,active_grid,V_seed,n,Num_theta_patch, &
                                           Num_phi_patch,Num_nu,Num_R,Num_gam_e,index_dyn,index_Y,index_syn_intger, &
                                           include_reverse_sync,include_forward_ssc,include_hadronic,include_proton_synch, &
                                           include_pg,include_neutrino,num_gam_p,num_nu_nu, &
                                           reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e, &
                                           reverse_sigma, &
                                           hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc, &
                                           n_threads_outer,n_threads_inner,adaptive_substeps,substep_rtol,substep_min, &
                                           substep_max,thermal_electrons,electron_solver_id,rt_axis,rg_axis,rr_axis, &
                                           sync_axis,ssc_axis,had_axis,rev_axis, &
                                           track_tobs,track_gamma,track_radius,track_mass,track_bfield,track_nu_m,track_nu_c, &
                                           track_nu_a,track_set)
        call sed_interpolation_structured(Boundary_sed,zero,rt_axis,rg_axis,rr_axis,sync_axis,V_seed,V_obs,Tobs, &
                                          n,Num_nu,Num_nu_obs,Num_Tobs,Num_theta_patch,Num_R,Num_phi_sed, &
                                          n_threads_projection,fwd_sync_obs)
        if (include_forward_ssc /= 0) call sed_interpolation_structured(Boundary_sed,zero,rt_axis,rg_axis,rr_axis, &
                                          ssc_axis,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs,Num_Tobs, &
                                          Num_theta_patch,Num_R,Num_phi_sed,n_threads_projection,fwd_ssc_obs)
        if (include_hadronic /= 0 .and. (include_proton_synch /= 0 .or. include_pg /= 0)) &
            call sed_interpolation_structured(Boundary_sed,zero, &
                                          rt_axis,rg_axis,rr_axis,had_axis,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs, &
                                          Num_Tobs,Num_theta_patch,Num_R,Num_phi_sed,n_threads_projection,fwd_hadronic_obs)
        if (include_reverse_sync /= 0) call sed_interpolation_structured(Boundary_sed,zero,rt_axis,rg_axis,rr_axis, &
                                          rev_axis,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs,Num_Tobs, &
                                          Num_theta_patch,Num_R,Num_phi_sed,n_threads_projection,rev_sync_obs)
    end subroutine run_axisymmetric

    subroutine run_nonaxisymmetric()
        implicit none
        real(8), allocatable :: rt_phi(:,:,:),rg_phi(:,:,:),rr_phi(:,:,:),sync_phi(:,:,:,:),ssc_phi(:,:,:,:)
        real(8), allocatable :: had_phi(:,:,:,:),rev_phi(:,:,:,:)

        allocate(rt_phi(Num_R,Num_theta_patch,Num_phi_patch),rg_phi(Num_R,Num_theta_patch,Num_phi_patch), &
                 rr_phi(Num_R,Num_theta_patch,Num_phi_patch),sync_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch), &
                 ssc_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch),had_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch), &
                 rev_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch))
        rt_phi=zero; rg_phi=one; rr_phi=zero; sync_phi=zero; ssc_phi=zero; had_phi=zero; rev_phi=zero
        call structured_solve_nonaxisymmetric(Boundary,E_iso_grid,Gamma0_grid,active_grid,V_seed,n,Num_theta_patch, &
                                              Num_phi_patch,Num_nu,Num_R,Num_gam_e,index_dyn,index_Y,index_syn_intger, &
                                              include_reverse_sync,include_forward_ssc,include_hadronic,include_proton_synch, &
                                              include_pg,include_neutrino,num_gam_p,num_nu_nu, &
                                              reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e, &
                                              reverse_sigma, &
                                              hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc, &
                                              n_threads_outer,n_threads_inner,adaptive_substeps,substep_rtol,substep_min, &
                                              substep_max,thermal_electrons,electron_solver_id,rt_phi,rg_phi,rr_phi, &
                                              sync_phi,ssc_phi,had_phi,rev_phi, &
                                              track_tobs,track_gamma,track_radius,track_mass,track_bfield,track_nu_m, &
                                              track_nu_c,track_nu_a,track_set)
        call sed_interpolation_structured_phi(Boundary_sed,rt_phi,rg_phi,rr_phi,sync_phi,V_seed,V_obs,Tobs, &
                                              n,Num_nu,Num_nu_obs,Num_Tobs,Num_theta_patch,Num_phi_patch,Num_R, &
                                              n_threads_projection,fwd_sync_obs)
        if (include_forward_ssc /= 0) call sed_interpolation_structured_phi(Boundary_sed,rt_phi,rg_phi,rr_phi, &
                                              ssc_phi,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs,Num_Tobs, &
                                              Num_theta_patch,Num_phi_patch,Num_R,n_threads_projection,fwd_ssc_obs)
        if (include_hadronic /= 0 .and. (include_proton_synch /= 0 .or. include_pg /= 0)) &
            call sed_interpolation_structured_phi(Boundary_sed, &
                                              rt_phi,rg_phi,rr_phi,had_phi,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs, &
                                              Num_Tobs,Num_theta_patch,Num_phi_patch,Num_R,n_threads_projection,fwd_hadronic_obs)
        if (include_reverse_sync /= 0) call sed_interpolation_structured_phi(Boundary_sed,rt_phi,rg_phi,rr_phi, &
                                              rev_phi,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs,Num_Tobs, &
                                              Num_theta_patch,Num_phi_patch,Num_R,n_threads_projection,rev_sync_obs)
    end subroutine run_nonaxisymmetric
end subroutine structured_jet_flux_1d

subroutine structured_solve_axisymmetric(Boundary,E_iso_grid,Gamma0_grid,active_grid,V_seed,n,Num_theta_patch, &
                                         Num_phi_patch,Num_nu,Num_R,Num_gam_e,index_dyn,index_Y,index_syn_intger, &
                                         include_reverse_sync,include_forward_ssc,include_hadronic,include_proton_synch, &
                                         include_pg,include_neutrino,num_gam_p,num_nu_nu, &
                                         reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e, &
                                         reverse_sigma, &
                                         hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc,n_threads_outer, &
                                         n_threads_inner,adaptive_substeps,substep_rtol,substep_min, &
                                         substep_max,thermal_electrons,electron_solver_id,rt_axis,rg_axis,rr_axis, &
                                         sync_axis,ssc_axis,had_axis,rev_axis, &
                                         track_tobs,track_gamma,track_radius,track_mass,track_bfield,track_nu_m,track_nu_c, &
                                         track_nu_a,track_set)
    !$ use omp_lib
    use constants
    implicit none
    integer, intent(in) :: n,Num_theta_patch,Num_phi_patch,Num_nu,Num_R,Num_gam_e,index_dyn,index_Y,index_syn_intger
    integer, intent(in) :: include_reverse_sync,include_forward_ssc,include_hadronic,include_proton_synch
    integer, intent(in) :: include_pg,include_neutrino,num_gam_p,num_nu_nu
    integer, intent(in) :: n_threads_outer,n_threads_inner,adaptive_substeps,substep_min,substep_max,thermal_electrons
    integer, intent(in) :: electron_solver_id
    integer, intent(in) :: active_grid(Num_theta_patch,Num_phi_patch)
    integer, intent(inout) :: track_set
    real(8), intent(in) :: Boundary(n),E_iso_grid(Num_theta_patch,Num_phi_patch),Gamma0_grid(Num_theta_patch,Num_phi_patch)
    real(8), intent(in) :: V_seed(Num_nu),substep_rtol,hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc
    real(8), intent(in) :: reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e,reverse_sigma
    real(8), intent(out) :: rt_axis(Num_R,Num_theta_patch),rg_axis(Num_R,Num_theta_patch),rr_axis(Num_R,Num_theta_patch)
    real(8), intent(out) :: sync_axis(Num_nu,Num_R,Num_theta_patch),ssc_axis(Num_nu,Num_R,Num_theta_patch)
    real(8), intent(out) :: had_axis(Num_nu,Num_R,Num_theta_patch),rev_axis(Num_nu,Num_R,Num_theta_patch)
    real(8), intent(out) :: track_tobs(Num_R),track_gamma(Num_R),track_radius(Num_R),track_mass(Num_R),track_bfield(Num_R)
    real(8), intent(out) :: track_nu_m(Num_R),track_nu_c(Num_R),track_nu_a(Num_R)
    integer, allocatable :: solve_index(:), solve_reps(:)
    integer :: it,iu,rep_idx,unique_count

    allocate(solve_index(Num_theta_patch),solve_reps(Num_theta_patch))
    solve_index=0; solve_reps=0; unique_count=0
    do it=1,Num_theta_patch
        if (active_grid(it,1) /= 0) call register_axis_patch(it)
    end do

    !$OMP PARALLEL DO num_threads(n_threads_outer) schedule(dynamic,1)
    do iu=1,unique_count
        call solve_axis_patch(solve_reps(iu))
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO num_threads(n_threads_outer) schedule(static) private(rep_idx)
    do it=1,Num_theta_patch
        rep_idx=solve_index(it)
        if (rep_idx /= 0 .and. rep_idx /= it) call copy_axis_patch(rep_idx,it)
    end do
    !$OMP END PARALLEL DO

    deallocate(solve_index,solve_reps)

contains

    subroutine register_axis_patch(it)
        implicit none
        integer, intent(in) :: it
        integer :: jr,rep

        do jr=1,unique_count
            rep=solve_reps(jr)
            if (E_iso_grid(it,1) == E_iso_grid(rep,1) .and. Gamma0_grid(it,1) == Gamma0_grid(rep,1)) then
                solve_index(it)=rep
                return
            end if
        end do
        unique_count=unique_count+1
        solve_reps(unique_count)=it
        solve_index(it)=it
    end subroutine register_axis_patch

    subroutine solve_axis_patch(it)
        implicit none
        integer, intent(in) :: it

        call structured_solve_element(Boundary,E_iso_grid(it,1),Gamma0_grid(it,1),V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                      index_dyn,index_Y,index_syn_intger,include_forward_ssc,include_hadronic, &
                                      include_reverse_sync,include_proton_synch,include_pg,include_neutrino, &
                                      num_gam_p,num_nu_nu,n_threads_inner, &
                                      reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e, &
                                      reverse_sigma,hadronic_p_p, &
                                      hadronic_epsilon_p,hadronic_eta_acc,adaptive_substeps, &
                                      substep_rtol,substep_min,substep_max,thermal_electrons,electron_solver_id, &
                                      rt_axis(:,it),rg_axis(:,it),rr_axis(:,it),sync_axis(:,:,it),ssc_axis(:,:,it), &
                                      had_axis(:,:,it),rev_axis(:,:,it), &
                                      track_tobs,track_gamma,track_radius,track_mass,track_bfield,track_nu_m,track_nu_c, &
                                      track_nu_a,track_set)
    end subroutine solve_axis_patch

    subroutine copy_axis_patch(src,dst)
        implicit none
        integer, intent(in) :: src,dst

        rt_axis(:,dst)=rt_axis(:,src); rg_axis(:,dst)=rg_axis(:,src); rr_axis(:,dst)=rr_axis(:,src)
        sync_axis(:,:,dst)=sync_axis(:,:,src); ssc_axis(:,:,dst)=ssc_axis(:,:,src)
        had_axis(:,:,dst)=had_axis(:,:,src); rev_axis(:,:,dst)=rev_axis(:,:,src)
    end subroutine copy_axis_patch
end subroutine structured_solve_axisymmetric

subroutine structured_solve_nonaxisymmetric(Boundary,E_iso_grid,Gamma0_grid,active_grid,V_seed,n,Num_theta_patch, &
                                            Num_phi_patch,Num_nu,Num_R,Num_gam_e,index_dyn,index_Y,index_syn_intger, &
                                            include_reverse_sync,include_forward_ssc,include_hadronic,include_proton_synch, &
                                            include_pg,include_neutrino,num_gam_p,num_nu_nu, &
                                            reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e, &
                                            reverse_sigma, &
                                            hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc,n_threads_outer, &
                                            n_threads_inner,adaptive_substeps,substep_rtol,substep_min, &
                                            substep_max,thermal_electrons,electron_solver_id,rt_phi,rg_phi,rr_phi, &
                                            sync_phi,ssc_phi,had_phi,rev_phi, &
                                            track_tobs,track_gamma,track_radius,track_mass,track_bfield,track_nu_m,track_nu_c, &
                                            track_nu_a,track_set)
    !$ use omp_lib
    use constants
    implicit none
    integer, intent(in) :: n,Num_theta_patch,Num_phi_patch,Num_nu,Num_R,Num_gam_e,index_dyn,index_Y,index_syn_intger
    integer, intent(in) :: include_reverse_sync,include_forward_ssc,include_hadronic,include_proton_synch
    integer, intent(in) :: include_pg,include_neutrino,num_gam_p,num_nu_nu
    integer, intent(in) :: n_threads_outer,n_threads_inner,adaptive_substeps,substep_min,substep_max,thermal_electrons
    integer, intent(in) :: electron_solver_id
    integer, intent(in) :: active_grid(Num_theta_patch,Num_phi_patch)
    integer, intent(inout) :: track_set
    real(8), intent(in) :: Boundary(n),E_iso_grid(Num_theta_patch,Num_phi_patch),Gamma0_grid(Num_theta_patch,Num_phi_patch)
    real(8), intent(in) :: V_seed(Num_nu),substep_rtol,hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc
    real(8), intent(in) :: reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e,reverse_sigma
    real(8), intent(out) :: rt_phi(Num_R,Num_theta_patch,Num_phi_patch),rg_phi(Num_R,Num_theta_patch,Num_phi_patch)
    real(8), intent(out) :: rr_phi(Num_R,Num_theta_patch,Num_phi_patch),sync_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch)
    real(8), intent(out) :: ssc_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch),had_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch)
    real(8), intent(out) :: rev_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch)
    real(8), intent(out) :: track_tobs(Num_R),track_gamma(Num_R),track_radius(Num_R),track_mass(Num_R),track_bfield(Num_R)
    real(8), intent(out) :: track_nu_m(Num_R),track_nu_c(Num_R),track_nu_a(Num_R)
    integer, allocatable :: solve_index(:), solve_reps(:)
    integer :: it,ip,flat,iu,rep_flat,rep_it,rep_ip,unique_count

    allocate(solve_index(Num_theta_patch*Num_phi_patch),solve_reps(Num_theta_patch*Num_phi_patch))
    solve_index=0; solve_reps=0; unique_count=0
    do ip=1,Num_phi_patch
        do it=1,Num_theta_patch
            flat=flatten_patch(it,ip)
            if (active_grid(it,ip) /= 0) call register_phi_patch(it,ip,flat)
        end do
    end do

    !$OMP PARALLEL DO num_threads(n_threads_outer) schedule(dynamic,1) private(rep_flat,rep_it,rep_ip)
    do iu=1,unique_count
        rep_flat=solve_reps(iu)
        call unflatten_patch(rep_flat,rep_it,rep_ip)
        call solve_phi_patch(rep_it,rep_ip)
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO num_threads(n_threads_outer) collapse(2) schedule(static) private(flat,rep_flat,rep_it,rep_ip)
    do ip=1,Num_phi_patch
        do it=1,Num_theta_patch
            flat=flatten_patch(it,ip)
            rep_flat=solve_index(flat)
            if (rep_flat /= 0 .and. rep_flat /= flat) then
                call unflatten_patch(rep_flat,rep_it,rep_ip)
                call copy_phi_patch(rep_it,rep_ip,it,ip)
            end if
        end do
    end do
    !$OMP END PARALLEL DO

    deallocate(solve_index,solve_reps)

contains

    integer function flatten_patch(it,ip)
        implicit none
        integer, intent(in) :: it,ip

        flatten_patch=(ip-1)*Num_theta_patch+it
    end function flatten_patch

    subroutine unflatten_patch(flat,it,ip)
        implicit none
        integer, intent(in) :: flat
        integer, intent(out) :: it,ip

        it=mod(flat-1,Num_theta_patch)+1
        ip=(flat-1)/Num_theta_patch+1
    end subroutine unflatten_patch

    subroutine register_phi_patch(it,ip,flat)
        implicit none
        integer, intent(in) :: it,ip,flat
        integer :: jr,rep,rep_it,rep_ip

        do jr=1,unique_count
            rep=solve_reps(jr)
            call unflatten_patch(rep,rep_it,rep_ip)
            if (E_iso_grid(it,ip) == E_iso_grid(rep_it,rep_ip) .and. &
                Gamma0_grid(it,ip) == Gamma0_grid(rep_it,rep_ip)) then
                solve_index(flat)=rep
                return
            end if
        end do
        unique_count=unique_count+1
        solve_reps(unique_count)=flat
        solve_index(flat)=flat
    end subroutine register_phi_patch

    subroutine solve_phi_patch(it,ip)
        implicit none
        integer, intent(in) :: it,ip

        call structured_solve_element(Boundary,E_iso_grid(it,ip),Gamma0_grid(it,ip),V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                      index_dyn,index_Y,index_syn_intger,include_forward_ssc,include_hadronic, &
                                      include_reverse_sync,include_proton_synch,include_pg,include_neutrino, &
                                      num_gam_p,num_nu_nu,n_threads_inner, &
                                      reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e, &
                                      reverse_sigma,hadronic_p_p, &
                                      hadronic_epsilon_p,hadronic_eta_acc,adaptive_substeps, &
                                      substep_rtol,substep_min,substep_max,thermal_electrons,electron_solver_id, &
                                      rt_phi(:,it,ip),rg_phi(:,it,ip),rr_phi(:,it,ip),sync_phi(:,:,it,ip), &
                                      ssc_phi(:,:,it,ip),had_phi(:,:,it,ip),rev_phi(:,:,it,ip), &
                                      track_tobs,track_gamma,track_radius,track_mass,track_bfield, &
                                      track_nu_m,track_nu_c,track_nu_a,track_set)
    end subroutine solve_phi_patch

    subroutine copy_phi_patch(src_it,src_ip,dst_it,dst_ip)
        implicit none
        integer, intent(in) :: src_it,src_ip,dst_it,dst_ip

        rt_phi(:,dst_it,dst_ip)=rt_phi(:,src_it,src_ip)
        rg_phi(:,dst_it,dst_ip)=rg_phi(:,src_it,src_ip)
        rr_phi(:,dst_it,dst_ip)=rr_phi(:,src_it,src_ip)
        sync_phi(:,:,dst_it,dst_ip)=sync_phi(:,:,src_it,src_ip)
        ssc_phi(:,:,dst_it,dst_ip)=ssc_phi(:,:,src_it,src_ip)
        had_phi(:,:,dst_it,dst_ip)=had_phi(:,:,src_it,src_ip)
        rev_phi(:,:,dst_it,dst_ip)=rev_phi(:,:,src_it,src_ip)
    end subroutine copy_phi_patch
end subroutine structured_solve_nonaxisymmetric

subroutine structured_solve_element(Boundary,E_iso,Gamma0,V_seed,n,Num_nu,Num_R,Num_gam_e,index_dyn,index_Y,index_syn_intger, &
                                    include_forward_ssc,include_hadronic,include_reverse_sync,include_proton_synch, &
                                    include_pg,include_neutrino,num_gam_p,num_nu_nu,n_threads, &
                                    reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e, &
                                    reverse_sigma,hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc,adaptive_substeps, &
                                    substep_rtol, &
                                    substep_min,substep_max,thermal_electrons,electron_solver_id,R_Tobs,R_Gamma,R, &
                                    sync_abs,ssc_abs,had_abs,rev_abs,track_tobs,track_gamma,track_radius,track_mass,track_bfield, &
                                    track_nu_m,track_nu_c,track_nu_a,track_set)
    !$ use omp_lib
    use constants
    use dynamics_common, only: dynamics_external_density_profile, density_jump_max, rs_shell_matter_fraction
    use electron_shell_transport_common, only: electron_solver_dg_1d, electron_solver_fullhide_1d
    use electron_radiation_kernel, only: get_syn_selected_state
    use electron_reverse_kernel, only: electron_reverse_evolve
    implicit none
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_dyn,index_Y,index_syn_intger,include_forward_ssc
    integer, intent(in) :: include_hadronic,include_proton_synch,include_pg,include_neutrino
    integer, intent(in) :: num_gam_p,num_nu_nu,n_threads,adaptive_substeps
    integer, intent(in) :: include_reverse_sync,substep_min,substep_max,thermal_electrons,electron_solver_id
    integer, intent(inout) :: track_set
    real(8), intent(in) :: Boundary(n),E_iso,Gamma0,V_seed(Num_nu),substep_rtol
    real(8), intent(in) :: reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e,reverse_sigma
    real(8), intent(in) :: hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc
    real(8), intent(out) :: R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R)
    real(8), intent(out) :: sync_abs(Num_nu,Num_R),ssc_abs(Num_nu,Num_R),had_abs(Num_nu,Num_R),rev_abs(Num_nu,Num_R)
    real(8), intent(out) :: track_tobs(Num_R),track_gamma(Num_R),track_radius(Num_R),track_mass(Num_R),track_bfield(Num_R)
    real(8), intent(out) :: track_nu_m(Num_R),track_nu_c(Num_R),track_nu_a(Num_R)
    real(8) :: B_local(n)
    real(8) :: R_mass(Num_R),nu_m_local(Num_R),nu_c_local(Num_R),nu_a_local(Num_R)
    real(8) :: M3(Num_R),B3(Num_R),U3(Num_R),V3(Num_R),Gamma34(Num_R)
    real(8) :: Secondary_M3(density_jump_max,Num_R),Secondary_U3(density_jump_max,Num_R)
    real(8) :: Secondary_V3(density_jump_max,Num_R),Secondary_B3(density_jump_max,Num_R)
    real(8) :: Secondary_M3_total(Num_R),Secondary_U3_total(Num_R),Secondary_V3_total(Num_R),Secondary_B3_total(Num_R)
    real(8) :: Secondary_pressure_total(Num_R),Secondary_enthalpy_density_total(Num_R)
    real(8) :: Secondary_gamma_contact(Num_R),Secondary_pressure_3(Num_R),Secondary_gamma_43(Num_R)
    real(8) :: Secondary_beta_rs(Num_R),Secondary_u_diss(Num_R),Secondary_dissipated_energy(Num_R)
    real(8) :: Secondary_electron_injected_energy(Num_R),Secondary_branch_gamma_m(density_jump_max,Num_R)
    real(8) :: Secondary_branch_gamma_contact(density_jump_max,Num_R),Secondary_branch_gamma_43(density_jump_max,Num_R)
    real(8) :: Secondary_branch_compression(density_jump_max,Num_R)
    real(8) :: Secondary_branch_beta_rs(density_jump_max,Num_R),Secondary_branch_u_diss(density_jump_max,Num_R)
    real(8) :: Secondary_nu_m(Num_R),Secondary_nu_c(Num_R)
    real(8) :: Secondary_start_radius(density_jump_max),Secondary_end_radius(density_jump_max)
    real(8) :: Secondary_start_tobs_axis(density_jump_max),Secondary_end_tobs_axis(density_jump_max)
    real(8) :: T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gam_m_cross,B3_ordered_cross
    real(8), allocatable :: gam_e(:),dN_gam_e(:,:),P_syn(:,:),Seed_syn(:,:),P_ssc(:,:),Seed_ssc(:,:),B_field(:),shell_energy(:)
    real(8) :: dNe,prev_radius,shell_volume,P_emit_tmp(Num_nu),Tau_syn_tmp(Num_nu)
    logical :: Secondary_event_active(density_jump_max)
    integer :: ir

    B_local=Boundary; B_local(1)=Gamma0; B_local(14)=E_iso
    allocate(gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R),P_syn(Num_nu,Num_R),Seed_syn(Num_nu,Num_R), &
             P_ssc(Num_nu,Num_R),Seed_ssc(Num_nu,Num_R),B_field(Num_R),shell_energy(Num_R))
    if (include_reverse_sync /= 0) then
        call dynamics_reverse(reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e, &
                              reverse_sigma,B_local,n,Num_R,T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross, &
                              M3_cross,gam_m_cross,B3_ordered_cross,R_Tobs,R_Gamma,R,R_mass,M3,B3,U3,V3,Gamma34, &
                              Secondary_M3,Secondary_U3,Secondary_V3,Secondary_B3, &
                              Secondary_M3_total,Secondary_U3_total,Secondary_V3_total,Secondary_B3_total, &
                              Secondary_pressure_total,Secondary_enthalpy_density_total, &
                              Secondary_gamma_contact,Secondary_pressure_3,Secondary_gamma_43,Secondary_beta_rs, &
                              Secondary_u_diss,Secondary_dissipated_energy,Secondary_electron_injected_energy, &
                              Secondary_branch_gamma_m,Secondary_branch_gamma_contact,Secondary_branch_gamma_43, &
                              Secondary_branch_compression,Secondary_branch_beta_rs,Secondary_branch_u_diss, &
                              Secondary_nu_m,Secondary_nu_c, &
                              Secondary_event_active,Secondary_start_radius,Secondary_end_radius, &
                              Secondary_start_tobs_axis,Secondary_end_tobs_axis)
    else
        call dynamics_forward(B_local,n,Num_R,index_dyn,R_Tobs,R_Gamma,R,R_mass)
        M3=zero; B3=zero; U3=zero; V3=zero; Gamma34=zero
        T_cross=zero; R_cross=zero; U3_cross=zero; M3_cross=zero
    end if
    select case(electron_solver_id)
    case(electron_solver_fullhide_1d)
        call fs_electron_fullhide_1d(B_local,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                     index_Y,index_syn_intger,n_threads,adaptive_substeps,substep_rtol, &
                                     substep_min,substep_max,thermal_electrons,gam_e,dN_gam_e,P_syn,Seed_syn, &
                                     nu_m_local,nu_c_local,nu_a_local)
    case(electron_solver_dg_1d)
        if (thermal_electrons /= 0) error stop 'structured_solve_element: dg_1d does not support thermal electrons.'
        call fs_electron_dg_1d(B_local,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y, &
                               index_syn_intger,n_threads,gam_e,dN_gam_e,P_syn,Seed_syn,nu_m_local,nu_c_local,nu_a_local)
    case default
        error stop 'structured_solve_element: unsupported electron solver id.'
    end select
    P_ssc=zero; Seed_ssc=zero; had_abs=zero; rev_abs=zero
    if (include_forward_ssc /= 0) call ssc_spec(R,gam_e,dN_gam_e,V_seed,Seed_syn,Num_nu,Num_R,Num_gam_e,n_threads,P_ssc,Seed_ssc)

    do ir=1,Num_R
        call dynamics_external_density_profile(B_local(12),B_local(11),R(ir),B_local(n),1,B_local(21),B_local(22),B_local(23),dNe)
        B_field(ir)=0.39d0*dsqrt(B_local(6)*dNe*(R_Gamma(ir)*(R_Gamma(ir)-one)))
        prev_radius=zero
        if (ir > 1) prev_radius=R(ir-1)
        shell_volume=(4d0/3d0)*pi*(R(ir)**3-prev_radius**3)
        shell_energy(ir)=zero
        if (include_hadronic /= 0) shell_energy(ir)=hadronic_epsilon_p*(R_Gamma(ir)-one)*shell_volume*dNe*Para_m_p*Para_c**2
    end do

    if (include_hadronic /= 0 .and. (include_proton_synch /= 0 .or. include_pg /= 0)) then
        block
            real(8), allocatable :: gam_p(:),dN_gam_p(:,:),P_had_syn(:,:),Seed_had_syn(:,:),P_pg(:,:),V_nu(:),P_nu(:,:)

        allocate(gam_p(num_gam_p),dN_gam_p(num_gam_p,Num_R),P_had_syn(Num_nu,Num_R),Seed_had_syn(Num_nu,Num_R), &
                 P_pg(Num_nu,Num_R),V_nu(num_nu_nu),P_nu(num_nu_nu,Num_R))
        call fs_hadronic_1d(R_Tobs,R_Gamma,R,shell_energy,B_field,V_seed,Seed_syn,hadronic_p_p, &
                            hadronic_epsilon_p,hadronic_eta_acc, &
                            include_proton_synch,include_pg,include_neutrino,Num_nu,Num_R,num_gam_p,num_nu_nu,n_threads, &
                            gam_p,dN_gam_p,P_had_syn,Seed_had_syn,P_pg,V_nu,P_nu)
        Seed_syn=Seed_syn+Seed_had_syn
        had_abs=P_had_syn+P_pg
        end block
    end if

    if (include_reverse_sync /= 0) then
        block
            real(8), allocatable :: gam_e_rev(:),dN_gam_e_rev(:,:),P_rev_syn(:,:),Seed_rev_syn(:,:)
            real(8) :: para_m_ej,reverse_total,reverse_target
            integer :: ig

        allocate(gam_e_rev(Num_gam_e),dN_gam_e_rev(Num_gam_e,Num_R),P_rev_syn(Num_nu,Num_R),Seed_rev_syn(Num_nu,Num_R))
        para_m_ej=E_iso*rs_shell_matter_fraction(reverse_sigma)/Gamma0/Para_c**2
        call electron_reverse_evolve(reverse_delta_t_s*Para_c,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e, &
                                     Gamma0,B_local(5),B_local(6),B_local(8),B_local(12),B_local(11),para_m_ej, &
                                     B_local(21),B_local(22),B_local(23),B_local(n), &
                                     T_cross,R_cross,U3_cross,V3_cross,M3_cross,R_Tobs,R_Gamma,R,B3,M3,U3,V3,V_seed, &
                                     Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads,gam_e_rev,dN_gam_e_rev, &
                                     electron_solver_id)
        do ir=1,Num_R
            reverse_total=zero
            do ig=2,Num_gam_e
                reverse_total=reverse_total+0.5d0*(dN_gam_e_rev(ig,ir)+dN_gam_e_rev(ig-1,ir)) &
                              *(gam_e_rev(ig)-gam_e_rev(ig-1))
            end do
            reverse_target=M3(ir)/Para_m_p*reverse_f_e
            if (reverse_total > zero .and. reverse_target > zero) &
                dN_gam_e_rev(:,ir)=dN_gam_e_rev(:,ir)*reverse_target/reverse_total
            call get_syn_selected_state(index_syn_intger,R(ir),B3(ir),Num_gam_e,Num_nu,n_threads,gam_e_rev, &
                                        dN_gam_e_rev(:,ir),V_seed,P_emit_tmp,P_rev_syn(:,ir), &
                                        Seed_rev_syn(:,ir),Tau_syn_tmp)
        end do
        rev_abs=P_rev_syn
        end block
    end if

    call structured_apply_observer_absorption(Boundary,R_Gamma,R,V_seed,Seed_syn,Seed_ssc,P_syn,P_ssc, &
                                              Num_nu,Num_R,n_threads,include_reverse_sync,include_hadronic, &
                                              include_proton_synch,include_pg,sync_abs,ssc_abs,rev_abs,had_abs)

    !$OMP CRITICAL(structured_track_copy)
    if (track_set == 0) then
        track_tobs=R_Tobs; track_gamma=R_Gamma; track_radius=R; track_mass=R_mass; track_bfield=B_field
        track_nu_m=nu_m_local; track_nu_c=nu_c_local; track_nu_a=nu_a_local; track_set=1
    end if
    !$OMP END CRITICAL(structured_track_copy)
    deallocate(gam_e,dN_gam_e,P_syn,Seed_syn,P_ssc,Seed_ssc,B_field,shell_energy)
end subroutine structured_solve_element

subroutine structured_apply_observer_absorption(Boundary,R_Gamma,R,V_seed,Seed_syn,Seed_ssc,P_syn,P_ssc, &
                                                Num_nu,Num_R,n_threads,include_reverse_sync,include_hadronic, &
                                                include_proton_synch,include_pg,sync_abs,ssc_abs,rev_abs,had_abs)
    use constants
    implicit none
    integer, intent(in) :: Num_nu,Num_R,n_threads,include_reverse_sync,include_hadronic,include_proton_synch,include_pg
    real(8), intent(in) :: Boundary(*),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(in) :: Seed_syn(Num_nu,Num_R),Seed_ssc(Num_nu,Num_R),P_syn(Num_nu,Num_R),P_ssc(Num_nu,Num_R)
    real(8), intent(inout) :: rev_abs(Num_nu,Num_R),had_abs(Num_nu,Num_R)
    real(8), intent(out) :: sync_abs(Num_nu,Num_R),ssc_abs(Num_nu,Num_R)
    real(8) :: prefactor(Num_nu,Num_R),tau_extra(Num_nu,Num_R),absorption(Num_nu,Num_R)

    tau_extra=zero
    call annihilation(R_Gamma,R,V_seed,Seed_syn,Seed_ssc,tau_extra,Num_nu,Num_R,n_threads,absorption)
    associate(luminosity_distance_cm => Boundary(13), redshift_factor => one+Boundary(8))
        prefactor=absorption/(4d0*pi*luminosity_distance_cm**2)*redshift_factor
    end associate
    sync_abs=P_syn*prefactor
    ssc_abs=P_ssc*prefactor
    if (include_reverse_sync /= 0) rev_abs=rev_abs*prefactor
    if (include_hadronic /= 0 .and. (include_proton_synch /= 0 .or. include_pg /= 0)) had_abs=had_abs*prefactor
end subroutine structured_apply_observer_absorption
