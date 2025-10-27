subroutine ssc_spec(R,gam_e,dN_gam_e,V_seed,seed,Num_nu,Num_R,Num_gam_e,n_threads, P_SSC_spec,seed_SSC)
    use constants
    !$ use omp_lib
    IMPLICIT REAL(8)(A-H,O-Z)
    !***********************************************************
    integer, intent(in) :: Num_nu,Num_R,Num_gam_e,n_threads
    real(8), intent(in) :: R(Num_R),gam_e(Num_gam_e)
    real(8), intent(in) :: dN_gam_e(Num_gam_e,Num_R),V_seed(Num_nu),seed(Num_nu,Num_R)
    real(8), intent(out) :: P_SSC_spec(Num_nu,Num_R),seed_SSC(Num_nu,Num_R)

    allocatable :: simpson_weights(:), V_weights(:)

    allocate (simpson_weights(Num_gam_e), V_weights(Num_nu))

    P_SSC_spec=zero
    seed_SSC=zero
    
    para_hEme = Para_h/para_m_energy

    h_nu = log(V_seed(2))-log(V_seed(1))
    h_gam = log(gam_e(2))-log(gam_e(1))

    call compute_simpson_weights(simpson_weights, Num_gam_e)
    call compute_simpson_weights(V_weights, Num_nu)

    P_SSC_spec=zero
    seed_SSC=zero
    
    !$ call omp_set_dynamic(.true.)
    !$OMP PARALLEL num_threads(n_threads)
    !$OMP DO SIMD
    do I_R=1,Num_R
        do I_nu=1,Num_nu
            Vloc=V_seed(I_nu)
            Ephoton2eV=para_hEme*Vloc
            II=1
            do i=1,Num_gam_e
                if (gam_e(i)<=Ephoton2eV) then
                    II=II+1
                else
                    exit
                end if
            end do
            if (II==Num_gam_e) cycle
            
            dInteg=zero
            simpson_sum_nu = zero
            do Nu_s=1,Num_nu
               if (Vloc <= V_seed(Nu_s)) then
                  emission_int2 = zero
                  simpson_sum_gam = zero
                  do i_game=1,Num_gam_e
                     gam_val = gam_e(i_game)
                     if (Vloc > 0.25d0*V_seed(Nu_s)/gam_val/gam_val) then
                        fssc = Vloc/V_seed(Nu_s) - 0.25d0/gam_val/gam_val
                     else
                        fssc = zero
                     end if
                     val = dN_gam_e(i_game, I_R) * fssc / gam_val
                     weight = simpson_weights(i_game)
                     simpson_sum_gam = simpson_sum_gam + val * weight
                  end do
                  emission_int2 = (h_gam/3.0d0) * simpson_sum_gam
               else
                  emission_int2 = 0.0d0
                  simpson_sum_gam = 0.0d0
                  do i_game = II,Num_gam_e-1
                     gam_val = gam_e(i_game)
                     temp = gam_val - Ephoton2eV
                     if (temp <= 0) cycle
                     q = Vloc / (4.0d0 * gam_val * V_seed(Nu_s) * temp)
                     if (q >= one) cycle
                     q_gamma = Ephoton2eV / temp
                     fssc = two*q*(log(q)-q)+one+q+q_gamma*q_gamma/(two*(one+q_gamma))*(one-q)
                     val = dN_gam_e(i_game, I_R) * fssc / gam_val
                     weight = simpson_weights(i_game)
                     simpson_sum_gam = simpson_sum_gam + val * weight
                  end do
                  emission_int2 = (h_gam/3.0d0) * simpson_sum_gam
               end if
               val = seed(Nu_s, I_R) * emission_int2
               weight = V_weights(Nu_s)
               simpson_sum_nu = simpson_sum_nu + val * weight
            end do
            dInteg = (h_nu/3.0d0) * simpson_sum_nu
            
            P_v=dInteg*Vloc
            P_SSC_spec(I_nu,I_R)=P_SSC_spec(I_nu,I_R)+P_v
            F1=dInteg/R(I_R)/R(I_R)
            seed_SSC(I_nu,I_R)=seed_SSC(I_nu,I_R)+F1
        end do
    end do
    !$OMP END DO SIMD
    !$OMP END PARALLEL

    Temp_para=0.75d0*Para_c*Para_h*Para_SigmaT
    P_SSC_spec=P_SSC_spec*Temp_para
    
    Temp_para2=4.0d0*pi*Para_c*Para_h
    seed_SSC=seed_SSC/Temp_para2*Temp_para
    
    deallocate(simpson_weights, V_weights)

    return

contains

    subroutine compute_simpson_weights(weights, n)
        integer, intent(in) :: n
        real(8), intent(out) :: weights(n)
        integer :: i
        
        weights = 1.0d0
        if (n >= 3) then
            do i = 2, n-1
                if (mod(i,2) == 0) then
                    weights(i) = 4.0d0
                else
                    weights(i) = 2.0d0
                endif
            end do
        endif
    end subroutine compute_simpson_weights

end subroutine ssc_spec
