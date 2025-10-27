! Created by  on 2021/1/31.
subroutine annihilation(R_gamma,R,V_seed,seed_syn,seed_ssc,Num_nu,Num_R,n_threads, absorption)
    !$ use omp_lib
    use constants
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: Num_R,Num_nu,n_threads
    real(8), intent(in) :: R_gamma(Num_R),R(Num_R),V_seed(Num_nu),seed_syn(Num_nu,Num_R),seed_ssc(Num_nu,Num_R)
    real(8), intent(out) :: absorption(Num_nu,Num_R)

    allocatable :: seed_tot(:,:),seed_tot_mean(:,:),ep1(:,:),ep2(:,:),ep2ep1(:,:),dVloc(:), &
                dRariv_Sigma(:),beta(:),seed_tot_mean_dVloc(:,:)
                
    allocate (seed_tot(Num_nu,Num_R),dRariv_Sigma(Num_R),ep1(1,Num_nu),ep2(Num_nu-1,1),ep2ep1(Num_nu-1,Num_nu), &
            dVloc(Num_nu-1),seed_tot_mean(Num_nu,Num_R),beta(Num_R),seed_tot_mean_dVloc(Num_nu,Num_R))

    absorption=zero

    Num_cos=50
    dcos_bin=two/Num_cos
    Cross_Area=3.0d0/16.0d0*Para_sigmaT
    absorption=zero
    seed_tot=zero

    seed_tot=seed_syn+seed_ssc
    
    para_hEme=Para_h/para_m_energy
    dRariv_Sigma=Cross_Area*R/(12.0d0*R_gamma)
    beta=dsqrt(one-R_gamma**(-2))
    
    ep1(1,:)=para_hEme*V_seed
    ep2(:,1)=para_hEme*(V_seed(1:Num_nu-1)+V_seed(2:Num_nu))/two
    dVloc=V_seed(2:Num_nu)-V_seed(1:Num_nu-1)
    seed_tot_mean=(seed_tot(1:Num_nu-1,:)+seed_tot(2:Num_nu,:))/two
    ep2ep1=matmul(ep2,ep1)
    
!    call system_clock(int1)

    !$ call omp_set_dynamic(.true.)
    !$OMP PARALLEL num_threads(n_threads)
    !$OMP DO
    do i_R=1,Num_R
        seed_tot_mean_dVloc(:,i_R)=seed_tot_mean(:,i_R)*dVloc
    end do
    !$OMP END DO
    
    !$OMP DO SIMD
    do I_R=1,Num_R
        do Nu_s1=1,Num_nu
            Tau=zero
            do i_cos=1,Num_cos+1
                d_cos=-one+dcos_bin*(i_cos-one)
                z=(one-d_cos)/two
                Tau1=zero
                do Nu_s2=1,Num_nu-1
                    Temp_s0=ep2ep1(Nu_s2,Nu_s1)*z
                    if (Temp_s0 <= one .or. Temp_s0 >=1.0d12) cycle
                    Temp_b02=one-one/Temp_s0
                    Temp_b0=dsqrt(Temp_b02)
                    Temp_log=dlog((one+Temp_b0)/(one-Temp_b0))
                    Sigma_Gamma=(one-Temp_b02)*((3.0d0-Temp_b02*Temp_b02)*Temp_log-two*Temp_b0*(two-Temp_b02))
                    Tau1=Tau1+Sigma_Gamma*seed_tot_mean_dVloc(Nu_s2,I_R)
                end do
                Tau=Tau+Tau1*dcos_bin*(one-beta(I_R)*d_cos)
            end do
            Tau=Tau*dRariv_Sigma(I_R)/two
            if ((Tau-1.0d-4) < 1.0d-5) Tau=1.0d-4
            temp_abs=(one-dexp(-Tau))/Tau
            absorption(Nu_s1,I_R)=absorption(Nu_s1,I_R)+temp_abs
        end do
    end do
    !$OMP END DO SIMD
    !$OMP END PARALLEL

!    call system_clock(int2)
!    print*, 'time=', (int2-int1)/1000.0

    deallocate(seed_tot,dRariv_Sigma,ep1,ep2,ep2ep1,dVloc,seed_tot_mean,beta,seed_tot_mean_dVloc)

    return
end subroutine annihilation
