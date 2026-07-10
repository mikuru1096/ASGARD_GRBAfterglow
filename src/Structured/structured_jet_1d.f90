subroutine jet_flux_1d(Boundary,E_iso_grid,Gamma0_grid,active_grid,V_seed,V_obs,Tobs, &
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
    use dynamics_density_profile, only: density_profile, jump_max
    implicit none
    integer, intent(in) :: n,Num_theta_patch,Num_phi_patch,Num_nu,Num_nu_obs,Num_Tobs,Num_R,Num_gam_e
    integer, intent(in) :: Num_theta_sed,Num_phi_sed,index_dyn,index_Y,index_syn_intger
    integer, intent(in) :: include_reverse_sync,include_forward_ssc,include_hadronic,include_proton_synch
    integer, intent(in) :: include_pg,include_neutrino,num_gam_p,num_nu_nu
    integer, intent(in) :: axisymmetric,n_threads_outer,n_threads_inner,adaptive_substeps,substep_min,substep_max
    integer, intent(in), dimension(Num_theta_patch,Num_phi_patch) :: active_grid
    integer, intent(in) :: thermal_electrons,electron_solver_id
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_theta_patch,Num_phi_patch) :: E_iso_grid,Gamma0_grid
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu_obs) :: V_obs
    real(8), intent(in), dimension(Num_Tobs) :: Tobs
    real(8), intent(in) :: substep_rtol
    real(8), intent(in) :: reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e,reverse_sigma
    real(8), intent(in) :: hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc
    real(8), intent(out), dimension(Num_nu_obs,Num_Tobs) :: fwd_sync_obs,fwd_ssc_obs
    real(8), intent(out), dimension(Num_nu_obs,Num_Tobs) :: fwd_hadronic_obs,rev_sync_obs,total_obs
    real(8), intent(out), dimension(Num_R) :: track_tobs,track_gamma,track_radius,track_mass,track_bfield
    real(8), intent(out), dimension(Num_R) :: track_nu_m,track_nu_c,track_nu_a

    real(8), dimension(n) :: Boundary_sed
    integer :: track_set,n_threads_projection

    fwd_sync_obs=0d0; fwd_ssc_obs=0d0; fwd_hadronic_obs=0d0; rev_sync_obs=0d0; total_obs=0d0
    track_tobs=0d0; track_gamma=0d0; track_radius=0d0; track_mass=0d0; track_bfield=0d0
    track_nu_m=0d0; track_nu_c=0d0; track_nu_a=0d0; track_set=0
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
        real(8), allocatable, dimension(:,:) :: rt_axis,rg_axis,rr_axis
        real(8), allocatable, dimension(:,:,:) :: sync_axis,ssc_axis
        real(8), allocatable, dimension(:,:,:) :: had_axis,rev_axis

        allocate(rt_axis(Num_R,Num_theta_patch),rg_axis(Num_R,Num_theta_patch),rr_axis(Num_R,Num_theta_patch), &
                 sync_axis(Num_nu,Num_R,Num_theta_patch),ssc_axis(Num_nu,Num_R,Num_theta_patch), &
                 had_axis(Num_nu,Num_R,Num_theta_patch),rev_axis(Num_nu,Num_R,Num_theta_patch))
        rt_axis=0d0; rg_axis=1d0; rr_axis=0d0; sync_axis=0d0; ssc_axis=0d0; had_axis=0d0; rev_axis=0d0
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
        call sed_interpolation_structured(Boundary_sed,0d0,rt_axis,rg_axis,rr_axis,sync_axis,V_seed,V_obs,Tobs, &
                                          n,Num_nu,Num_nu_obs,Num_Tobs,Num_theta_patch,Num_R,Num_phi_sed, &
                                          n_threads_projection,fwd_sync_obs)
        if (include_forward_ssc /= 0) call sed_interpolation_structured(Boundary_sed,0d0,rt_axis,rg_axis,rr_axis, &
                                          ssc_axis,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs,Num_Tobs, &
                                          Num_theta_patch,Num_R,Num_phi_sed,n_threads_projection,fwd_ssc_obs)
        if (include_hadronic /= 0 .and. (include_proton_synch /= 0 .or. include_pg /= 0)) &
            call sed_interpolation_structured(Boundary_sed,0d0, &
                                          rt_axis,rg_axis,rr_axis,had_axis,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs, &
                                          Num_Tobs,Num_theta_patch,Num_R,Num_phi_sed,n_threads_projection,fwd_hadronic_obs)
        if (include_reverse_sync /= 0) call sed_interpolation_structured(Boundary_sed,0d0,rt_axis,rg_axis,rr_axis, &
                                          rev_axis,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs,Num_Tobs, &
                                          Num_theta_patch,Num_R,Num_phi_sed,n_threads_projection,rev_sync_obs)
    end subroutine run_axisymmetric

    subroutine run_nonaxisymmetric()
        implicit none
        real(8), allocatable, dimension(:,:,:) :: rt_phi,rg_phi,rr_phi
        real(8), allocatable, dimension(:,:,:,:) :: sync_phi,ssc_phi
        real(8), allocatable, dimension(:,:,:,:) :: had_phi,rev_phi

        allocate(rt_phi(Num_R,Num_theta_patch,Num_phi_patch),rg_phi(Num_R,Num_theta_patch,Num_phi_patch), &
                 rr_phi(Num_R,Num_theta_patch,Num_phi_patch),sync_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch), &
                 ssc_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch),had_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch), &
                 rev_phi(Num_nu,Num_R,Num_theta_patch,Num_phi_patch))
        rt_phi=0d0; rg_phi=1d0; rr_phi=0d0; sync_phi=0d0; ssc_phi=0d0; had_phi=0d0; rev_phi=0d0
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
        call sed_structured_phi(Boundary_sed,rt_phi,rg_phi,rr_phi,sync_phi,V_seed,V_obs,Tobs, &
                                              n,Num_nu,Num_nu_obs,Num_Tobs,Num_theta_patch,Num_phi_patch,Num_R, &
                                              n_threads_projection,fwd_sync_obs)
        if (include_forward_ssc /= 0) call sed_structured_phi(Boundary_sed,rt_phi,rg_phi,rr_phi, &
                                              ssc_phi,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs,Num_Tobs, &
                                              Num_theta_patch,Num_phi_patch,Num_R,n_threads_projection,fwd_ssc_obs)
        if (include_hadronic /= 0 .and. (include_proton_synch /= 0 .or. include_pg /= 0)) &
            call sed_structured_phi(Boundary_sed, &
                                              rt_phi,rg_phi,rr_phi,had_phi,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs, &
                                              Num_Tobs,Num_theta_patch,Num_phi_patch,Num_R,n_threads_projection,fwd_hadronic_obs)
        if (include_reverse_sync /= 0) call sed_structured_phi(Boundary_sed,rt_phi,rg_phi,rr_phi, &
                                              rev_phi,V_seed,V_obs,Tobs,n,Num_nu,Num_nu_obs,Num_Tobs, &
                                              Num_theta_patch,Num_phi_patch,Num_R,n_threads_projection,rev_sync_obs)
    end subroutine run_nonaxisymmetric
end subroutine jet_flux_1d

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
    use electron_cooling_kernel, only: cooling_reset
    implicit none
    integer, intent(in) :: n,Num_theta_patch,Num_phi_patch,Num_nu,Num_R,Num_gam_e,index_dyn,index_Y,index_syn_intger
    integer, intent(in) :: include_reverse_sync,include_forward_ssc,include_hadronic,include_proton_synch
    integer, intent(in) :: include_pg,include_neutrino,num_gam_p,num_nu_nu
    integer, intent(in) :: n_threads_outer,n_threads_inner,adaptive_substeps,substep_min,substep_max,thermal_electrons
    integer, intent(in) :: electron_solver_id
    integer, intent(in), dimension(Num_theta_patch,Num_phi_patch) :: active_grid
    integer, intent(inout) :: track_set
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_theta_patch,Num_phi_patch) :: E_iso_grid,Gamma0_grid
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in) :: substep_rtol,hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc
    real(8), intent(in) :: reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e,reverse_sigma
    real(8), intent(out), dimension(Num_R,Num_theta_patch) :: rt_axis,rg_axis,rr_axis
    real(8), intent(out), dimension(Num_nu,Num_R,Num_theta_patch) :: sync_axis,ssc_axis
    real(8), intent(out), dimension(Num_nu,Num_R,Num_theta_patch) :: had_axis,rev_axis
    real(8), intent(out), dimension(Num_R) :: track_tobs,track_gamma,track_radius,track_mass,track_bfield
    real(8), intent(out), dimension(Num_R) :: track_nu_m,track_nu_c,track_nu_a
    integer, allocatable, dimension(:) :: solve_index,solve_reps
    integer :: it,iu,rep_idx,unique_count

    allocate(solve_index(Num_theta_patch),solve_reps(Num_theta_patch))
    solve_index=0; solve_reps=0; unique_count=0
    do it=1,Num_theta_patch
        if (active_grid(it,1) /= 0) call register_axis_patch(it)
    end do

    !$OMP PARALLEL num_threads(n_threads_outer)
    !$ if (n_threads_outer > 1) call cooling_reset()
    !$OMP DO schedule(dynamic,1)
    do iu=1,unique_count
        call solve_axis_patch(solve_reps(iu))
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP PARALLEL DO num_threads(n_threads_outer) schedule(static) private(rep_idx)
    do it=1,Num_theta_patch
        rep_idx=solve_index(it)
        if (rep_idx /= 0 .and. rep_idx /= it) then
            rt_axis(:,it)=rt_axis(:,rep_idx)
            rg_axis(:,it)=rg_axis(:,rep_idx)
            rr_axis(:,it)=rr_axis(:,rep_idx)
            sync_axis(:,:,it)=sync_axis(:,:,rep_idx)
            ssc_axis(:,:,it)=ssc_axis(:,:,rep_idx)
            had_axis(:,:,it)=had_axis(:,:,rep_idx)
            rev_axis(:,:,it)=rev_axis(:,:,rep_idx)
        end if
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
    use electron_cooling_kernel, only: cooling_reset
    implicit none
    integer, intent(in) :: n,Num_theta_patch,Num_phi_patch,Num_nu,Num_R,Num_gam_e,index_dyn,index_Y,index_syn_intger
    integer, intent(in) :: include_reverse_sync,include_forward_ssc,include_hadronic,include_proton_synch
    integer, intent(in) :: include_pg,include_neutrino,num_gam_p,num_nu_nu
    integer, intent(in) :: n_threads_outer,n_threads_inner,adaptive_substeps,substep_min,substep_max,thermal_electrons
    integer, intent(in) :: electron_solver_id
    integer, intent(in), dimension(Num_theta_patch,Num_phi_patch) :: active_grid
    integer, intent(inout) :: track_set
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_theta_patch,Num_phi_patch) :: E_iso_grid,Gamma0_grid
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in) :: substep_rtol,hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc
    real(8), intent(in) :: reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e,reverse_sigma
    real(8), intent(out), dimension(Num_R,Num_theta_patch,Num_phi_patch) :: rt_phi,rg_phi
    real(8), intent(out), dimension(Num_R,Num_theta_patch,Num_phi_patch) :: rr_phi
    real(8), intent(out), dimension(Num_nu,Num_R,Num_theta_patch,Num_phi_patch) :: sync_phi
    real(8), intent(out), dimension(Num_nu,Num_R,Num_theta_patch,Num_phi_patch) :: ssc_phi,had_phi
    real(8), intent(out), dimension(Num_nu,Num_R,Num_theta_patch,Num_phi_patch) :: rev_phi
    real(8), intent(out), dimension(Num_R) :: track_tobs,track_gamma,track_radius,track_mass,track_bfield
    real(8), intent(out), dimension(Num_R) :: track_nu_m,track_nu_c,track_nu_a
    integer, allocatable, dimension(:) :: solve_index,solve_reps
    integer :: it,ip,flat,iu,rep_flat,rep_it,rep_ip,unique_count

    allocate(solve_index(Num_theta_patch*Num_phi_patch),solve_reps(Num_theta_patch*Num_phi_patch))
    solve_index=0; solve_reps=0; unique_count=0
    do ip=1,Num_phi_patch
        do it=1,Num_theta_patch
            flat=(ip-1)*Num_theta_patch+it
            if (active_grid(it,ip) /= 0) call register_phi_patch(it,ip,flat)
        end do
    end do

    !$OMP PARALLEL num_threads(n_threads_outer) private(rep_flat,rep_it,rep_ip)
    !$ if (n_threads_outer > 1) call cooling_reset()
    !$OMP DO schedule(dynamic,1)
    do iu=1,unique_count
        rep_flat=solve_reps(iu)
        rep_it=mod(rep_flat-1,Num_theta_patch)+1
        rep_ip=(rep_flat-1)/Num_theta_patch+1
        call solve_phi_patch(rep_it,rep_ip)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP PARALLEL DO num_threads(n_threads_outer) collapse(2) schedule(static) private(flat,rep_flat,rep_it,rep_ip)
    do ip=1,Num_phi_patch
        do it=1,Num_theta_patch
            flat=(ip-1)*Num_theta_patch+it
            rep_flat=solve_index(flat)
            if (rep_flat /= 0 .and. rep_flat /= flat) then
                rep_it=mod(rep_flat-1,Num_theta_patch)+1
                rep_ip=(rep_flat-1)/Num_theta_patch+1
                call copy_phi_patch(rep_it,rep_ip,it,ip)
            end if
        end do
    end do
    !$OMP END PARALLEL DO

    deallocate(solve_index,solve_reps)

contains

    subroutine register_phi_patch(it,ip,flat)
        implicit none
        integer, intent(in) :: it,ip,flat
        integer :: jr,rep,rep_it,rep_ip

        do jr=1,unique_count
            rep=solve_reps(jr)
            rep_it=mod(rep-1,Num_theta_patch)+1
            rep_ip=(rep-1)/Num_theta_patch+1
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
    use dynamics_density_profile, only: density_profile, jump_max
    use electron_shell_transport, only: solver_dg, solver_fullhide
    use electron_radiation_kernel, only: syn_state
    use electron_reverse_kernel, only: electron_reverse_evolve
    implicit none
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_dyn,index_Y,index_syn_intger,include_forward_ssc
    integer, intent(in) :: include_hadronic,include_proton_synch,include_pg,include_neutrino
    integer, intent(in) :: num_gam_p,num_nu_nu,n_threads,adaptive_substeps
    integer, intent(in) :: include_reverse_sync,substep_min,substep_max,thermal_electrons,electron_solver_id
    integer, intent(inout) :: track_set
    real(8), intent(in), dimension(n) :: Boundary
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in) :: E_iso,Gamma0,substep_rtol
    real(8), intent(in) :: reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e,reverse_sigma
    real(8), intent(in) :: hadronic_p_p,hadronic_epsilon_p,hadronic_eta_acc
    real(8), intent(out), dimension(Num_R) :: R_Tobs,R_Gamma,R
    real(8), intent(out), dimension(Num_nu,Num_R) :: sync_abs,ssc_abs,had_abs,rev_abs
    real(8), intent(out), dimension(Num_R) :: track_tobs,track_gamma,track_radius,track_mass,track_bfield
    real(8), intent(out), dimension(Num_R) :: track_nu_m,track_nu_c,track_nu_a
    real(8), dimension(n) :: B_local
    real(8), dimension(Num_R) :: R_mass,nu_m_local,nu_c_local,nu_a_local
    real(8), dimension(Num_R) :: M3,B3,U3,V3,Gamma34
    real(8), dimension(jump_max,Num_R) :: branch_m3,branch_u3
    real(8), dimension(jump_max,Num_R) :: branch_v3,branch_b3
    real(8), dimension(Num_R) :: total_m3,total_u3,total_v3,total_b3
    real(8), dimension(Num_R) :: total_p,total_w
    real(8), dimension(Num_R) :: avg_gc,avg_p3,avg_g43
    real(8), dimension(Num_R) :: avg_brs,total_ud,diss_e
    real(8), dimension(Num_R) :: inj_e
    real(8), dimension(jump_max,Num_R) :: branch_gm
    real(8), dimension(jump_max,Num_R) :: branch_gc,branch_g43
    real(8), dimension(jump_max,Num_R) :: branch_comp
    real(8), dimension(jump_max,Num_R) :: branch_brs,branch_ud
    real(8), dimension(Num_R) :: nu_m_rs,nu_c_rs
    real(8), dimension(jump_max) :: start_r,end_r
    real(8), dimension(jump_max) :: start_t,end_t
    real(8) :: T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross,M3_cross,gmcross,b3ordx
    real(8), allocatable, dimension(:) :: gam_e,B_field,shell_energy
    real(8), allocatable, dimension(:,:) :: dN_gam_e,P_syn,Seed_syn,P_ssc,Seed_ssc
    real(8), dimension(Num_nu) :: P_emit_tmp,Tau_syn_tmp
    real(8) :: dNe,prev_radius,shell_volume
    logical, dimension(jump_max) :: event_on
    integer :: ir

    B_local=Boundary; B_local(1)=Gamma0; B_local(14)=E_iso
    allocate(gam_e(Num_gam_e),dN_gam_e(Num_gam_e,Num_R),P_syn(Num_nu,Num_R),Seed_syn(Num_nu,Num_R), &
             P_ssc(Num_nu,Num_R),Seed_ssc(Num_nu,Num_R),B_field(Num_R),shell_energy(Num_R))
    if (include_reverse_sync /= 0) then
        call dynamics_reverse(reverse_delta_t_s,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e, &
                              reverse_sigma,B_local,n,Num_R,T_cross,R_cross,e3_cross,gam20,U3_cross,V3_cross, &
                              M3_cross,gmcross,b3ordx,R_Tobs,R_Gamma,R,R_mass,M3,B3,U3,V3,Gamma34, &
                              branch_m3,branch_u3,branch_v3,branch_b3, &
                              total_m3,total_u3,total_v3,total_b3, &
                              total_p,total_w, &
                              avg_gc,avg_p3,avg_g43,avg_brs, &
                              total_ud,diss_e,inj_e, &
                              branch_gm,branch_gc,branch_g43, &
                              branch_comp,branch_brs,branch_ud, &
                              nu_m_rs,nu_c_rs, &
                              event_on,start_r,end_r, &
                              start_t,end_t)
    else
        call dynamics_forward(B_local,n,Num_R,index_dyn,R_Tobs,R_Gamma,R,R_mass)
        M3=0d0; B3=0d0; U3=0d0; V3=0d0; Gamma34=0d0
        T_cross=0d0; R_cross=0d0; U3_cross=0d0; M3_cross=0d0
    end if
    select case(electron_solver_id)
    case(solver_fullhide)
        call fs_fullhide_1d(B_local,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e, &
                                     index_Y,index_syn_intger,n_threads,adaptive_substeps,substep_rtol, &
                                     substep_min,substep_max,thermal_electrons,gam_e,dN_gam_e,P_syn,Seed_syn, &
                                     nu_m_local,nu_c_local,nu_a_local)
    case(solver_dg)
        call fs_dg_1d(B_local,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y, &
                               index_syn_intger,n_threads,thermal_electrons,gam_e,dN_gam_e, &
                               P_syn,Seed_syn,nu_m_local,nu_c_local,nu_a_local)
    case default
        error stop 'structured_solve_element: unsupported electron solver id.'
    end select
    P_ssc=0d0; Seed_ssc=0d0; had_abs=0d0; rev_abs=0d0
    if (include_forward_ssc /= 0) call ssc_spec(R,gam_e,dN_gam_e,V_seed,Seed_syn,Num_nu,Num_R,Num_gam_e,n_threads,P_ssc,Seed_ssc)

    do ir=1,Num_R
        call density_profile(B_local(12),B_local(11),R(ir),B_local(n),1,B_local(21),B_local(22),B_local(23),dNe)
        B_field(ir)=0.39d0*dsqrt(B_local(6)*dNe*(R_Gamma(ir)*(R_Gamma(ir)-1d0)))
        prev_radius=0d0
        if (ir > 1) prev_radius=R(ir-1)
        shell_volume=(4d0/3d0)*pi*(R(ir)**3-prev_radius**3)
        shell_energy(ir)=0d0
        if (include_hadronic /= 0) shell_energy(ir)=hadronic_epsilon_p*(R_Gamma(ir)-1d0)*shell_volume*dNe*Para_m_p*Para_c**2
    end do

    if (include_hadronic /= 0 .and. (include_proton_synch /= 0 .or. include_pg /= 0)) then
        block
            real(8), allocatable, dimension(:) :: gam_p,V_nu
            real(8), allocatable, dimension(:,:) :: dN_gam_p,P_had_syn,Seed_had_syn,P_pg,P_nu

        allocate(gam_p(num_gam_p),dN_gam_p(num_gam_p,Num_R),P_had_syn(Num_nu,Num_R),Seed_had_syn(Num_nu,Num_R), &
                 P_pg(Num_nu,Num_R),V_nu(num_nu_nu),P_nu(num_nu_nu,Num_R))
        call hadronic_1d(R_Tobs,R_Gamma,R,shell_energy,B_field,V_seed,Seed_syn,hadronic_p_p, &
                            hadronic_epsilon_p,hadronic_eta_acc, &
                            include_proton_synch,include_pg,include_neutrino,Num_nu,Num_R,num_gam_p,num_nu_nu,n_threads, &
                            gam_p,dN_gam_p,P_had_syn,Seed_had_syn,P_pg,V_nu,P_nu)
        Seed_syn=Seed_syn+Seed_had_syn
        had_abs=P_had_syn+P_pg
        end block
    end if

    if (include_reverse_sync /= 0) then
        block
            real(8), allocatable, dimension(:) :: gam_e_rev
            real(8), allocatable, dimension(:,:) :: dN_gam_e_rev,P_rev_syn,Seed_rev_syn
            real(8) :: para_m_ej,reverse_total,reverse_target
            integer :: ig

        allocate(gam_e_rev(Num_gam_e),dN_gam_e_rev(Num_gam_e,Num_R),P_rev_syn(Num_nu,Num_R),Seed_rev_syn(Num_nu,Num_R))
        if (reverse_sigma <= 0d0) then
            para_m_ej=E_iso/Gamma0/Para_c**2
        else
            para_m_ej=E_iso/(1d0+reverse_sigma)/Gamma0/Para_c**2
        end if
        call electron_reverse_evolve(reverse_delta_t_s*Para_c,reverse_epsilon_e,reverse_epsilon_b,reverse_p,reverse_f_e, &
                                     Gamma0,B_local(5),B_local(6),B_local(8),B_local(12),B_local(11),para_m_ej, &
                                     B_local(21),B_local(22),B_local(23),B_local(n), &
                                     T_cross,R_cross,U3_cross,V3_cross,M3_cross,R_Tobs,R_Gamma,R,B3,M3,U3,V3,V_seed, &
                                     Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads,gam_e_rev,dN_gam_e_rev, &
                                     electron_solver_id)
        do ir=1,Num_R
            reverse_total=0d0
            do ig=2,Num_gam_e
                reverse_total=reverse_total+0.5d0*(dN_gam_e_rev(ig,ir)+dN_gam_e_rev(ig-1,ir)) &
                              *(gam_e_rev(ig)-gam_e_rev(ig-1))
            end do
            reverse_target=M3(ir)/Para_m_p*reverse_f_e
            if (reverse_total > 0d0 .and. reverse_target > 0d0) &
                dN_gam_e_rev(:,ir)=dN_gam_e_rev(:,ir)*reverse_target/reverse_total
            call syn_state(index_syn_intger,R(ir),B3(ir),Num_gam_e,Num_nu,n_threads,gam_e_rev, &
                                        dN_gam_e_rev(:,ir),V_seed,P_emit_tmp,P_rev_syn(:,ir), &
                                        Seed_rev_syn(:,ir),Tau_syn_tmp)
        end do
        rev_abs=P_rev_syn
        end block
    end if

    call apply_absorption(Boundary,R_Gamma,R,V_seed,Seed_syn,Seed_ssc,P_syn,P_ssc, &
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

subroutine apply_absorption(Boundary,R_Gamma,R,V_seed,Seed_syn,Seed_ssc,P_syn,P_ssc, &
                                                Num_nu,Num_R,n_threads,include_reverse_sync,include_hadronic, &
                                                include_proton_synch,include_pg,sync_abs,ssc_abs,rev_abs,had_abs)
    use constants
    implicit none
    integer, intent(in) :: Num_nu,Num_R,n_threads,include_reverse_sync,include_hadronic,include_proton_synch,include_pg
    real(8), intent(in), dimension(*) :: Boundary
    real(8), intent(in), dimension(Num_R) :: R_Gamma,R
    real(8), intent(in), dimension(Num_nu) :: V_seed
    real(8), intent(in), dimension(Num_nu,Num_R) :: Seed_syn,Seed_ssc,P_syn,P_ssc
    real(8), intent(inout), dimension(Num_nu,Num_R) :: rev_abs,had_abs
    real(8), intent(out), dimension(Num_nu,Num_R) :: sync_abs,ssc_abs
    real(8), dimension(Num_nu,Num_R) :: prefactor,tau_extra,absorption

    tau_extra=0d0
    call annihilation(R_Gamma,R,V_seed,Seed_syn,Seed_ssc,tau_extra,Num_nu,Num_R,n_threads,absorption)
    associate(luminosity_distance_cm => Boundary(13), redshift_factor => 1d0+Boundary(8))
        prefactor=absorption/(4d0*pi*luminosity_distance_cm**2)*redshift_factor
    end associate
    sync_abs=P_syn*prefactor
    ssc_abs=P_ssc*prefactor
    if (include_reverse_sync /= 0) rev_abs=rev_abs*prefactor
    if (include_hadronic /= 0 .and. (include_proton_synch /= 0 .or. include_pg /= 0)) had_abs=had_abs*prefactor
end subroutine apply_absorption
