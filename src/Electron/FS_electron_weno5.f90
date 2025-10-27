module get_Y
  use constants
  private

  public :: besselk, get_syn, get_SSA_numerical, get_IC_numerical, get_Y_Nakar, get_Y_Fan

contains

!****************************************************************************************
!**************************** Besselk function interpolation ****************************
!****************************************************************************************
function besselk(var)
    REAL(8) :: val,var,besselk
    integer :: i
    REAL(8),dimension(100),parameter :: minus_theta= &
                 (/1.0000000e-04,1.1233240e-04,1.2618569e-04,1.4174742e-04,1.5922828e-04,1.7886495e-04,2.0092330e-04,& 
                 2.2570197e-04,2.5353645e-04,2.8480359e-04,3.1992671e-04,3.5938137e-04,4.0370173e-04,4.5348785e-04,&
                 5.0941380e-04,5.7223677e-04,6.4280731e-04,7.2208090e-04,8.1113083e-04,9.1116276e-04,1.0235310e-03,&
                 1.1497570e-03,1.2915497e-03,1.4508288e-03,1.6297508e-03,1.8307383e-03,2.0565123e-03,2.3101297e-03,&
                 2.5950242e-03,2.9150531e-03,3.2745492e-03,3.6783798e-03,4.1320124e-03,4.6415888e-03,5.2140083e-03,&
                 5.8570208e-03,6.5793322e-03,7.3907220e-03,8.3021757e-03,9.3260335e-03,1.0476158e-02,1.1768120e-02,&
                 1.3219411e-02,1.4849683e-02,1.6681005e-02,1.8738174e-02,2.1049041e-02,2.3644894e-02,2.6560878e-02,&
                 2.9836472e-02,3.3516027e-02,3.7649358e-02,4.2292429e-02,4.7508102e-02,5.3366992e-02,5.9948425e-02,&
                 6.7341507e-02,7.5646333e-02,8.4975344e-02,9.5454846e-02,1.0722672e-01,1.2045035e-01,1.3530478e-01,&
                 1.5199111e-01,1.7073526e-01,1.9179103e-01,2.1544347e-01,2.4201283e-01,2.7185882e-01,3.0538555e-01,&
                 3.4304693e-01,3.8535286e-01,4.3287613e-01,4.8626016e-01,5.4622772e-01,6.1359073e-01,6.8926121e-01,&
                 7.7426368e-01,8.6974900e-01,9.7700996e-01,1.0974988e+00,1.2328467e+00,1.3848864e+00,1.5556761e+00,&
                 1.7475284e+00,1.9630407e+00,2.2051307e+00,2.4770764e+00,2.7825594e+00,3.1257158e+00,3.5111917e+00,&
                 3.9442061e+00,4.4306215e+00,4.9770236e+00,5.5908102e+00,6.2802914e+00,7.0548023e+00,7.9248290e+00,&
                 8.9021509e+00,1.0000000e+01/)
    REAL(8),dimension(100),parameter :: besselk2= &
                 (/2.0000000e+08,1.5849658e+08,1.2560583e+08,9.9540471e+07,7.8884121e+07,6.2514316e+07,4.9541527e+07,&
                 3.9260813e+07,3.1113522e+07,2.4656934e+07,1.9540199e+07,1.5485273e+07,1.2271814e+07,9.7252027e+06,&
                 7.7070567e+06,6.1077105e+06,4.8402560e+06,3.8358200e+06,3.0398217e+06,2.4090066e+06,1.9090964e+06,&
                 1.5129262e+06,1.1989680e+06,9.5016153e+05,7.5298666e+05,5.9672895e+05,4.7289738e+05,3.7476298e+05,&
                 2.9699315e+05,2.3536189e+05,1.8652017e+05,1.4781394e+05,1.1713992e+05,9.2831277e+04,7.3567095e+04,&
                 5.8300561e+04,4.6202094e+04,3.6614266e+04,2.9016076e+04,2.2994640e+04,1.8222755e+04,1.4441118e+04,&
                 1.1444235e+04,9.0692572e+03,7.1871275e+03,5.6955719e+03,4.5135397e+03,3.5767994e+03,2.8344487e+03,&
                 2.2461486e+03,1.7799308e+03,1.4104612e+03,1.1176629e+03,8.8562540e+02,7.0173970e+02,5.5601353e+02,&
                 4.4052817e+02,3.4900815e+02,2.7648028e+02,2.1900342e+02,1.7345426e+02,1.3735766e+02,1.0875212e+02,&
                 8.6083185e+01,6.8119011e+01,5.3883384e+01,4.2602693e+01,3.3663886e+01,2.6581150e+01,2.0969514e+01,&
                 1.6523922e+01,1.3002650e+01,1.0214166e+01,8.0067134e+00,6.2600579e+00,4.8789400e+00,3.7878895e+00,&
                 2.9271104e+00,2.2492185e+00,1.7166521e+00,1.2996171e+00,9.7445621e-01,7.2235470e-01,5.2831390e-01,&
                 3.8033947e-01,2.6880181e-01,1.8593568e-01,1.2545279e-01,8.2245762e-02,5.2165420e-02,3.1855005e-02,&
                 1.8626561e-02,1.0365680e-02,5.4526367e-03,2.6905368e-03,1.2347515e-03,5.2199962e-04,2.0112029e-04,&
                 6.9778918e-05,2.1509817e-05/)

   besselk=0.0d0
   
   do i=1,99
       if (var > minus_theta(i) .and. var <= minus_theta(i+1)) then
           val = (var-minus_theta(i))/(minus_theta(i+1)-minus_theta(i))
           besselk = besselk2(i)+val*(besselk2(i+1)-besselk2(i))
           exit
       end if
   end do

   return
end function besselk


!****************************************************************************************
!************************** get syn power and number density ****************************
!****************************************************************************************
subroutine get_syn(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e,V_seed, &
                       P_syn,Seed_syn)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: R_loc,DB,gam_e(Num_gam_e),dN_gam_e(Num_gam_e),V_seed(Num_nu)
real(8), intent(out) ::P_syn(Num_nu),Seed_syn(Num_nu)

real(8),allocatable,dimension (:) :: d_nu,dN1,ddN
allocate (dN1(Num_gam_e),ddN(Num_gam_e-1),d_nu(Num_nu-1))


    factor=(3.62d0/pi)**2
    Temp_syn=dsqrt(3d0)*para_e*para_e*para_e/Para_m_energy
    Rariv2=R_loc*R_loc
    d_nu=V_seed(2:Num_nu)-V_seed(1:Num_nu-1)
    dN1=dN_gam_e/(gam_e*gam_e)
    ddN=dN1(1:Num_gam_e-1)-dN1(2:Num_gam_e)
    
    !$ call omp_set_dynamic(.true.)
    !$OMP PARALLEL num_threads(n_threads)
    !$OMP DO SIMD
    do I_nu=1,Num_nu
       V_cal=V_seed(I_nu)
       dInteg=zero
       Tau=zero
       do I_gam_e=1,Num_gam_e-1
          gam_e_mean2=(gam_e(I_gam_e)+gam_e(I_gam_e+1))**2/4d0
          Vc=(4.2d6)*gam_e_mean2*DB !Which is $\nu_c$
          x=V_cal/Vc !Which is ($\nu/\nu_c$)
          Fx=1.81d0*dexp(-x)/dsqrt(x**(-2d0/3d0)+factor) !Approximate function of synchrotron radiation spectrum
!         Fx=2.149d0*x**(one/3.0d0)*dexp(-x) !!Another approximate function
          dN=(dN_gam_e(I_gam_e)+dN_gam_e(I_gam_e+1))/two
          dgam_e=gam_e(I_gam_e+1)-gam_e(I_gam_e)
          dInteg=dInteg+dN*Fx*dgam_e
          !====================  [SSA]  ======================
          Tau=Tau+gam_e_mean2*ddN(I_gam_e)*Fx
          !Synchrotron self absorption effect
          !===================================================
       end do
       P_v=Temp_syn*DB*dInteg ! with units in erg/Hz/s
       Tau=1.046d4*Tau*DB/(4d0*pi*Rariv2*V_cal*V_cal) !Synchrotron self absorption effect
       if ((Tau-1d-4) < 1d-5) Tau=1d-4
       P_syn(I_nu)=P_v*(one-dexp(-Tau))/Tau !Radiation transfer equation for the emission-absorption plasma
       Seed_syn(I_nu)=P_syn(I_nu)/(Rariv2*V_cal)
    end do
    !$OMP END DO SIMD
    !$OMP END PARALLEL

    temp_para=4d0*pi*Para_c*Para_h
    Seed_syn=Seed_syn/temp_para


deallocate (d_nu,dN1,ddN)
end subroutine get_syn

!****************************************************************************************
!**************** get SSA pile-up effect with Ghisellini & Svensson 1991 ****************
!****************************************************************************************
subroutine get_SSA_numerical(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn, dot_gam_e)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dot_gam_e(Num_gam_e)

real(8),allocatable,dimension (:) :: dF_nu,d_nu,sigma_SSA
allocate (dF_nu(Num_nu),d_nu(Num_nu-1),sigma_SSA(Num_nu-1))

    B_cr=4.4d13
    Temp1=2.5042d-22*B_cr/DB
    Temp2=7.787d-22*B_cr/DB
    Cyclotron_nu=para_e*DB/(two*pi*para_m_e*para_c)
    
    dot_gam_e=zero
    d_nu=V_seed(2:Num_nu)-V_seed(1:Num_nu-1)
    dF_nu=Seed_syn*para_h*V_seed*para_c
    
    do I_gam_e=1,Num_gam_e
       gam=gam_e(I_gam_e)
       gam2=gam*gam
       gam3=gam2*gam
       V_lowlim=Cyclotron_nu/gam
       V_uplim=1.5d0*gam2*Cyclotron_nu
       do I_nu=1,Num_nu-1
          if (V_lowlim < V_seed(I_nu) .and. V_seed(I_nu) <= V_uplim) then
             sigma_SSA(I_nu)=Temp1*(3d0*V_lowlim/V_seed(I_nu))**(5d0/3d0)
          else if (V_seed(I_nu) > V_uplim) then
             sigma_SSA(I_nu)=Temp2/gam3*(Cyclotron_nu/V_seed(I_nu))*exp(-V_seed(I_nu)/V_uplim)
          else
             sigma_SSA(I_nu)=zero
          end if
       end do
       dot_gam_e(I_gam_e)=sum(sigma_SSA*0.5d0*(dF_nu(2:Num_nu)+dF_nu(1:Num_nu-1))*d_nu)
    end do


deallocate (dF_nu,d_nu,sigma_SSA)
end subroutine get_SSA_numerical


!****************************************************************************************
!******************* get Compton cooling effiency fully numerical ***********************
!****************************************************************************************
subroutine get_IC_numerical(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn, dot_gam_e)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),Seed_syn(Num_nu)
real(8), intent(out) :: dot_gam_e(Num_gam_e)

real(8),allocatable,dimension (:) :: d_nu,photon_number
real(8),allocatable,dimension (:,:) :: gam_e_mean,E_seed1,gam_e_mean_E_seed
allocate (d_nu(Num_nu-1),photon_number(Num_nu-1),gam_e_mean(1,Num_gam_e),E_seed1(Num_nu,1), &
          gam_e_mean_E_seed(Num_nu-1,Num_gam_e-1))

    para_hEme = Para_h/para_m_energy

    d_nu=V_seed(2:Num_nu)-V_seed(1:Num_nu-1)

    gam_e_mean(1,:)=(gam_e(1:Num_gam_e-1)+gam_e(2:Num_gam_e))/two
    E_seed1(:,1)=V_seed*para_hEme
    gam_e_mean_E_seed=4d0*matmul(E_seed1,gam_e_mean)
    
    dot_gam_e=zero

    photon_number=(Seed_syn(1:Num_nu-1)+Seed_syn(2:Num_nu))/two
    
    !$ call omp_set_dynamic(.true.)
    !$OMP PARALLEL num_threads(n_threads)
    !$OMP DO SIMD
    do i_gam_e=1,Num_gam_e
       dInteg1=zero
       game=gam_e_mean(1,i_gam_e)
       game_pow=game*game
       var=0.25d0/game_pow
       do I_nu=1,Num_nu-1      !frequency circulation for seed photons
          dInteg2=zero
          V_t=V_seed(I_nu)
          E_t2eV=para_hEme*V_t
          do Nu_s=1,Num_nu-1     !frequency circulation for SSC photons
             fssc=zero
             Vloc=V_seed(Nu_s)
             E_Vloc2eV=para_hEme*Vloc
             if (Vloc > var*V_t .and. Vloc <= V_t) then
                fssc=Vloc/V_t-var
             else
                uplim=(4d0*game_pow*E_t2eV)/(one+gam_e_mean_E_seed(I_nu,i_gam_e))
                if (E_Vloc2eV > uplim) cycle
                temp=game-E_Vloc2eV
                q=E_Vloc2eV/(gam_e_mean_E_seed(I_nu,i_gam_e)*temp)
                fssc=two*q*(log(q)-q)+one+q+ &
                0.5d0*(one-q)*(4d0*game*E_t2eV*q)**2/(1+4d0*game*q*E_t2eV)
             end if
             dInteg2=dInteg2+Vloc*fssc*d_nu(Nu_s)
          end do
          dot_gam_e(i_gam_e)=dot_gam_e(i_gam_e)+photon_number(I_nu)/V_t*d_nu(I_nu)*dInteg2
       end do
    end do
    !$OMP END DO SIMD
    !$OMP END PARALLEL

    dot_gam_e=dot_gam_e/gam_e/gam_e*para_h*Para_h*Para_SigmaT/para_m_energy
    dot_gam_e(Num_gam_e)=0.99*dot_gam_e(Num_gam_e-1)

deallocate (d_nu,gam_e_mean,photon_number,E_seed1,gam_e_mean_E_seed)
end subroutine get_IC_numerical


!****************************************************************************************
!********************* get Compton Y with Nakar & Piran 2009 ****************************
!****************************************************************************************
subroutine get_Y_Nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn, Compton)
!$ use omp_lib
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e,Num_nu,n_threads
real(8), intent(in) :: gam_e(Num_gam_e),V_seed(Num_nu),P_syn(Num_nu)
real(8), intent(out) :: Compton(Num_gam_e)

real(8),allocatable,dimension (:) :: hat_nu,d_nu
allocate (hat_nu(Num_gam_e),d_nu(Num_nu-1))

    Compton=zero
    var_Compensation=zero
    hat_nu=Para_m_energy/Para_h/gam_e
    d_nu=V_seed(2:Num_nu)-V_seed(1:Num_nu-1)
    
    do I_Compton=1,Num_gam_e
       do I_nu=1,Num_nu
          if (hat_nu(I_Compton) < V_seed(I_nu)) then
             temp=(hat_nu(I_Compton)-V_seed(I_nu-1))/(V_seed(I_nu)-V_seed(I_nu-1))
             var_Compensation=temp*(P_syn(I_nu)-P_syn(I_nu-1))*(hat_nu(I_Compton)-V_seed(I_nu-1))
             Compton(I_Compton)=sum(0.5*(P_syn(1:I_nu-2)+P_syn(2:I_nu-1))*d_nu(1:I_nu-2))+var_Compensation
          exit
          end if
       end do
    end do
    

deallocate (hat_nu,d_nu)
end subroutine get_Y_Nakar

!****************************************************************************************
!************************ get Compton Y with Fan et al. 2008 ****************************
!****************************************************************************************
subroutine get_Y_Fan(Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e, Compton)
implicit REAL(8)(A-H,O-Z)
integer, intent(in) :: Num_gam_e
real(8), intent(in) :: Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,gam_e(Num_gam_e)
real(8), intent(out) :: Compton(Num_gam_e)


    eta=(Gam_e_m/Gam_e_c)**(p-two)
    if (eta-one > 0.001) eta=one

        do i_gam_e=1,Num_gam_e
            if (Num_gam_e > i_gam_e) then
!*****************************[The general inverse Compton effect]******************************
                hat_gam=5.4246D6/sqrt(DB*gam_e(i_gam_e+1))
                if (Gam_e_m > Gam_e_c) then
                    if (hat_gam < Gam_e_c) then
                        eta_NK=zero
                    else
                        if (hat_gam < Gam_e_m) then
                            if (p>2) then
                                Step1=(p-1)/(p-2)*Gam_e_m-Gam_e_c
                                eta_NK=(hat_gam-Gam_e_c)/Step1
                            else
                                Step1=Gam_e_m**(p-1)*Gam_e_max**(2-p)-(p-1)*Gam_e_m-(2-p)*Gam_e_c
                                eta_NK=(2-p)*(hat_gam-Gam_e_c)/Step1
                            end if
                        else
                            if (p>2) then
                                Step2=Gam_e_m**(p-1)*hat_gam**(2-p)
                                Step3=(p-1)*Gam_e_m-(p-2)*Gam_e_c
                                eta_NK=1-Step2/Step3
                            else
                                Step2=Gam_e_m**(p-1)*Gam_e_max**(2-p)-(p-1)*Gam_e_m-(2-p)*Gam_e_c
                                Step3=Gam_e_m**(p-1)*(Gam_e_max**(2-p)-hat_gam**(2-p))
                                eta_NK=1-Step2/Step3
                            end if
                        end if
                    end if
                else
                    if (hat_gam < Gam_e_m) then
                        eta_NK=zero
                    else
                        if (hat_gam < Gam_e_c) then
                            if (p>2) then
                                Step4=Gam_e_c**(3-p)/(p-2.0)-Gam_e_m**(3-p)
                                eta_NK=(hat_gam**(3-p)-Gam_e_m**(3-p))/Step4
                            else
                                Step4=(3-p)*Gam_e_c*Gam_e_max**(2-p)-Gam_e_c**(3-p)-(2-p)*Gam_e_m**(3-p)
                                eta_NK=(2-p)*(hat_gam**(3-p)-Gam_e_m**(3-p))/Step4
                            end if
                        else
                            if (p>2) then
                                Step5=(3-p)*Gam_e_c*hat_gam**(2-p)
                                Step6=Gam_e_c**(3.0-p)-(p-2)*Gam_e_m**(3.0-p)
                                eta_NK=1-Step5/Step6
                            else
                                Step5=(3-p)*Gam_e_c*(Gam_e_max**(2-p)-hat_gam**(2-p))
                                Step6=(3-p)*Gam_e_c*Gam_e_max**(2-p)-Gam_e_c**(3-p)-(2-p)*Gam_e_m**(3-p)
                                eta_NK=1-Step5/Step6
                            end if
                        end if
                    end if
                end if
                Compton(i_gam_e)=0.5d0*(-1.0+sqrt(1.0+4.0*eta*eta_NK*Epsilon_e/Epsilon_b))
            else
                Compton(i_gam_e)=0.99*Compton(i_gam_e-1)
            end if
        end do


end subroutine get_Y_Fan

end module get_Y

subroutine fs_electron_weno5(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,n_threads, &
                                gam_e,dN_gam_e,P_syn,Seed_syn)
    !$ use omp_lib
    use constants
    use get_Y
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_R,Num_gam_e
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R),Seed_syn(Num_nu,Num_R)

    real(8),allocatable,dimension (:) :: f_r_times_gam_e,dEl,para_minus_gam_e_p,dEl1,x, &
            dN_x,dF1,fp,flux,Compton,Compton1,dot_gam_e,dot_gam_e_SSA
    real(8),allocatable,dimension (:,:) :: temp_store,temp_store_extended
    real(8),allocatable,dimension (:) :: dN_x_extended,fp_extended,flux_extended,dEl1_extended
    allocate (f_r_times_gam_e(Num_gam_e),dEl(Num_gam_e),para_minus_gam_e_p(Num_gam_e), &
            dEl1(Num_gam_e),x(Num_gam_e),dN_x(Num_gam_e),dF1(Num_gam_e),fp(Num_gam_e), &
            flux(0:Num_gam_e), temp_store(3,Num_gam_e),Compton(Num_gam_e), &
              dot_gam_e(Num_gam_e),dot_gam_e_SSA(Num_gam_e),Compton1(Num_gam_e))
    allocate(dN_x_extended(1-2:Num_gam_e+3),temp_store_extended(3, 1-2:Num_gam_e+3),&
             fp_extended(1-2:Num_gam_e+3),flux_extended(1-2:Num_gam_e+3),dEl1_extended(1-2:Num_gam_e+3)) !ghost cells

    !***********************[Parameter Initial]**********************
    Eta_0=Boundary(1)
    R_ini=Boundary(4)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    E_iso=Boundary(14)
    f_e=Boundary(16)
    R0=Boundary(n)
    
    P_syn=zero
    Seed_syn=zero
    
    !*****************Part 1: given the boundary consition [Using the analytical approximation]*********************
    if (A_star > zero) then
        dNe_wind=A_star*3.0d35/R(1)**2
        Para_N_e_ini=4d0*pi*R_ini*A_star*3.0d35
        if (dNe_wind <= dNe_ISM/4d0) then
            dNe=dNe_ISM
        else
            dNe=dNe_wind
        end if
    else
        dNe=dNe_ISM
        Para_N_e_ini=4d0/3d0*pi*R_ini**3*dNe_ISM
    end if
    
    if (R(1)<R0) then
        dNe=A_star*3.0d35/R0**2*4
        Para_N_e_ini=4d0/3d0*pi*R_ini**3*dNe_ISM
    end if

    DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(1)*(R_Gamma(1)-one)))
    Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
    DB_min=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma(Num_R)*(R_Gamma(Num_R)-one)))
    Gam_e_max_max=3d0*Para_m_energy/dsqrt(8d0*DB_min*Para_e**3)
    Gam_e_m=(p-two)/(p-one)*Epsilon_e/f_e*1836d0*(R_Gamma(1)-one)+one
    if (p<2.05 .and. p>=2.0) then
        Gam_e_m=0.05d0/1.05d0*Epsilon_e/f_e*1836d0*(R_Gamma(1)-one)+one
    else if (p<2 .and. p>1) then
        Gam_e_m=((two-p)/(p-one)*Epsilon_e/f_e*1836d0*(R_Gamma(1)-one)*Gam_e_max**(p-two))**(one/(p-one))+one
    end if
    Gam_e_c=7.7d8/(one+dsqrt(Epsilon_e/Epsilon_b))/R_Gamma(1)/DB**2/(R_Tobs(1)/two)
    do I_gam_e=1,Num_gam_e
        Gam_e(I_gam_e)=3d0*ten**(dlog10(Gam_e_max_max)*(I_gam_e-1)/(Num_gam_e-1))
        if (Gam_e_m > Gam_e_c) then
            if (Gam_e_c > Gam_e(I_gam_e) .or. Gam_e_max < Gam_e(I_gam_e)) then
                dN_gam_e(I_gam_e,1)=zero
            else
                Q1=Para_N_e_ini*Gam_e_c
                if (Gam_e_m > Gam_e(I_gam_e)) then
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e(I_gam_e)**(-2)
                else
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e_m**(p-one)*Gam_e(I_gam_e)**(-(p+one))
                end if
            end if
        else
            if (Gam_e_m > Gam_e(I_gam_e) .or. Gam_e_max < Gam_e(I_gam_e)) then
                dN_gam_e(I_gam_e,1)=zero
            else
                Q1=Para_N_e_ini*Gam_e_m**(p-one)
                if (Gam_e_c > Gam_e(I_gam_e)) then
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e(I_gam_e)**(-p)
                else
                    dN_gam_e(I_gam_e,1)=Q1*Gam_e_c*Gam_e(I_gam_e)**(-(p+one))
                end if
            end if
        end if
    end do
    !*******************Part 1 is completed [has been checked and there is no bug]**********************************
    !*******************Part 2: To calculate the electron distribution**********************************************
    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
    d_x=dlog10(gam_e(2)/gam_e(1))
    factor_adv=Para_sigmaT/(6.0d0*pi*Para_m_energy)
    para_minus_gam_e_p=one/(gam_e-one)**p*gam_e*dlog(ten)
    
    do I_tobs=2,Num_R
        R_loc=R(I_tobs-1)
        R_Gamma_loc=(R_Gamma(I_tobs)+R_Gamma(I_tobs-1))/two
        if (A_star > zero) then
            dNe_wind=A_star*3.0d35/R_loc**2
            if (dNe_wind <= dNe_ISM/4d0) then
                dNe=dNe_ISM
            else
                dNe=dNe_wind
            end if
        else
            dNe=dNe_ISM
        end if
        
        if (R_loc<R0) then
            dNe=A_star*3.0d35/R0**2
        end if

        DB=0.39d0*dsqrt(Epsilon_b*dNe*(R_Gamma_loc*(R_Gamma_loc-one)))
        Gam_e_max=3d0*Para_m_energy/dsqrt(8d0*DB*Para_e**3)
        Gam_e_m=(p-two)/(p-one)*Epsilon_e*1836d0*(R_Gamma_loc-one)/f_e+one
        if (p<2.05 .and. p>=2.0) then
            Gam_e_m=0.05d0/1.05d0*Epsilon_e*1836d0*(R_Gamma_loc-one)/f_e+one
        else if (p<2 .and. p>1) then
            Gam_e_m=((two-p)/(p-one)*Epsilon_e/f_e*1836d0*(R_Gamma_loc-one)*Gam_e_max**(p-two))**(one/(p-one))+one
        end if
        Gam_e_m_p=(p-one)*(Gam_e_m-one)**(p-one)
        Gam_e_c=7.7d8*(one+z)/R_Gamma_loc/DB**2/R_Tobs(I_tobs)

        beta_Gam=dsqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1/(f_r*Gam_e_max+1.333/(R(I_tobs)+R(I_tobs-1)))
        !***********************[Here we have presented the choice on Delta_r]******************************************
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=Int(dDD/dDR)
        L1=L1+10
        dDR=dDD/L1
        CFL=dDR/d_x
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)
        
        call get_syn(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e(:,I_tobs-1),V_seed, &
                     P_syn(:,I_tobs),Seed_syn(:,I_tobs))
        
        call get_SSA_numerical(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn(:,I_tobs), dot_gam_e_SSA)
        
        select case(index_Y)
        
        case(1)
        call get_IC_numerical(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn(:,I_tobs), &
                       dot_gam_e)
        
        dEl=(f_r+(dot_gam_e-dot_gam_e_SSA)/beta_Gam/R_Gamma_loc/para_c)*gam_e
        
        case(2)
        call get_Y_Nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn(:,I_tobs), &
                         Compton)

        Q=4d0*pi*R_loc*R_loc*para_c
        Compton=one+Compton/Q/(4d0*R_Gamma_loc*R_Gamma_loc*dNe*Para_m_p_E)
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA/beta_Gam/R_Gamma_loc/para_c)*gam_e
        
        case(3)
        call get_Y_Fan(Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e, &
                       Compton)
        Compton=one+Compton
        Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
        dEl=(f_r*Compton-dot_gam_e_SSA/beta_Gam/R_Gamma_loc/para_c)*gam_e
        
        case default
         
        print*, 'invalid Compton case, check your chosen model!'
        stop
         
        end select
        
        dEl(1:Num_gam_e-1)=(dEl(1:Num_gam_e-1)+dEl(2:Num_gam_e))*0.5d0
        dEl(Num_gam_e)=dEl(Num_gam_e-1)*0.5d0
        
        do L=1,L1
            R_loc=R_loc+dDR
            dEl1=(dEl+one/R_loc)/dlog(ten)
            if (A_star >= zero) then
                dNe_wind=A_star*3.0d35/R_loc**2
                if (dNe_wind <= dNe_ISM/4d0) then
                    dNe=dNe_ISM
                else
                    dNe=dNe_wind
                end if
            else
                dNe=dNe_ISM
            end if
            
            R_loc=R_loc+dDR
            Q=4d0/3d0*pi*(3d0*R_loc**2+dDR*(3d0*R_loc+dDR))*dNe*f_e*Gam_e_m_p  !here Q is Q_0*\gamma_m**p
            dF1=zero
            where(gam_e<Gam_e_max .and. gam_e>Gam_e_m) dF1=Q*para_minus_gam_e_p
!            dF1=dF1+Q*para_maxwell/Gam_e_m_p*(one-f_e)
            
            
            dN_x_extended(1-2:0) = dN_x(1)
            dN_x_extended(1:Num_gam_e) = dN_x
            dN_x_extended(Num_gam_e+1:Num_gam_e+3) = dN_x(Num_gam_e)
            temp_store_extended(1,:) = dN_x_extended
            
            dEl1_extended(1-2:0) = dEl1(1)
            dEl1_extended(1:Num_gam_e) = dEl1
            dEl1_extended(Num_gam_e+1:Num_gam_e+3) = dEl1(Num_gam_e)

            do j=1,3

              call update_ghost_cells(dN_x_extended, Num_gam_e)
        
              fp_extended = dEl1_extended * dN_x_extended
        
              do i_gam_e = 1, Num_gam_e
                 if (dEl1_extended(i_gam_e) <= 0.0d0) then
                    flux_extended(i_gam_e) = fpx(fp_extended(i_gam_e-2:i_gam_e+2))
                 else
                    flux_extended(i_gam_e) = fmx(fp_extended(i_gam_e-1:i_gam_e+3))
                 end if
              end do
        
              if(j==1) then
                do i = 1, Num_gam_e
                   dN_x_extended(i) = temp_store_extended(1,i) + CFL*(flux_extended(i)-flux_extended(i-1))
                end do
                temp_store_extended(2,:) = dN_x_extended
              else if(j==2) then
                do i = 1, Num_gam_e
                   dN_x_extended(i) = 0.75d0*temp_store_extended(1,i) + 0.25d0*(temp_store_extended(2,i) + &
                                  CFL*(flux_extended(i)-flux_extended(i-1)))
                end do
                temp_store_extended(3,:) = dN_x_extended
              else if(j==3) then
                do i = 1, Num_gam_e
                   dN_x_extended(i) = (temp_store_extended(1,i) + 2.0d0*(temp_store_extended(3,i) + &
                                  CFL*(flux_extended(i)-flux_extended(i-1))))/3.0d0
                end do
              end if
           end do
    
    ! 更新内部点（不包括ghost cells）
           dN_x = dN_x_extended(1:Num_gam_e)
           dN_x = dN_x + dF1*dDR
            
           where(dN_x < 0.0d0) dN_x = 0.0d0
        
           if (L1 == L) then
              dN_gam_e(:,I_tobs)=dN_x/gam_e/dlog(ten)
           end if
        end do
    end do
 
    deallocate (f_r_times_gam_e,dEl,para_minus_gam_e_p,dEl1,x,dN_x,dF1,fp,flux,temp_store,Compton,Compton1,dot_gam_e, &
                dN_x_extended, temp_store_extended, fp_extended, flux_extended)
    
    return
end subroutine

subroutine update_ghost_cells(arr, n)
        integer, intent(in) :: n
        real(8), intent(inout) :: arr(1-2:n+3)
        
        arr(1-2:0) = arr(1)
        arr(n+1:n+3) = arr(n)
        
        ! 可选：周期性边界条件
        ! arr(1-2) = arr(n-1)
        ! arr(1-1) = arr(n)
        ! arr(n+1) = arr(1)
        ! arr(n+2) = arr(2)
    end subroutine update_ghost_cells

function fpx(fps)
    real(8) :: fps(-2:2), fpx
    real(8) :: omega(3), fu(3), beta(3)
    real(8) :: tao5, totalpha, alpha(3), fomega(3), eps

    omega(1) = 0.1d0;   omega(2) = 0.6d0;   omega(3) = 0.3d0
    eps=1d-30
    
    if(any(isnan(fps))) then
        fpx = 0.0d0
        return
    end if
    
    fu(1) =  1.0d0/3.0d0*fps(-2) - 7.0d0/6.0d0*fps(-1) + 11.0d0/6.0d0*fps(0)
    fu(2) = -1.0d0/6.0d0*fps(-1) + 5.0d0/6.0d0*fps(0)  + 1.0d0/3.0d0*fps(1)
    fu(3) =  1.0d0/3.0d0*fps(0)  + 5.0d0/6.0d0*fps(1)  - 1.0d0/6.0d0*fps(2)
    
    beta(1) = 13.0d0/12.0d0*( fps(-2) - 2.0d0*fps(-1) + fps(0) )**2 &
            + 0.25d0*( fps(-2) - 4.0d0*fps(-1) + 3.0d0*fps(0) )**2
    beta(2) = 13.0d0/12.0d0*( fps(-1) - 2.0d0*fps(0) + fps(1) )**2 &
            + 0.25d0*( fps(1) - fps(-1) )**2
    beta(3) = 13.0d0/12.0d0*( fps(0) - 2.0d0*fps(1) + fps(2) )**2 &
            + 0.25d0*( 3.0d0*fps(0) - 4.0d0*fps(1) + fps(2) )**2
    
    tao5 = abs(beta(1) - beta(3))
    
    alpha(:) = omega(:)*( 1.0d0 + tao5/(beta(:)+eps) )
    totalpha = alpha(1) + alpha(2) + alpha(3)
    
    if(totalpha < eps) then
        fpx = fu(2)
    else
        fomega(:) = alpha(:)/totalpha
        fpx = fu(1)*fomega(1) + fu(2)*fomega(2) + fu(3)*fomega(3)
    end if
    
    return
end function fpx

function fmx(fms)
    real(8) :: fms(-1:3), fmx
    real(8) :: omega(3), fu(3), beta(3)
    real(8) :: tao5, totalpha, alpha(3), fomega(3), eps
    
    omega(1) = 0.1d0;   omega(2) = 0.6d0;   omega(3) = 0.3d0
    eps=1d-30
    
    if(any(isnan(fms))) then
        fmx = 0.0d0
        return
    end if
    
    fu(1) =  1.0d0/3.0d0*fms(3) - 7.0d0/6.0d0*fms(2) + 11.0d0/6.0d0*fms(1)
    fu(2) = -1.0d0/6.0d0*fms(2) + 5.0d0/6.0d0*fms(1) + 1.0d0/3.0d0*fms(0)
    fu(3) =  1.0d0/3.0d0*fms(1) + 5.0d0/6.0d0*fms(0) - 1.0d0/6.0d0*fms(-1)
    
    beta(1) = 13.0d0/12.0d0*( fms(3) - 2.0d0*fms(2) + fms(1) )**2 &
            + 0.25d0*( fms(3) - 4.0d0*fms(2) + 3.0d0*fms(1) )**2
    beta(2) = 13.0d0/12.0d0*( fms(2) - 2.0d0*fms(1) + fms(0) )**2 &
            + 0.25d0*( fms(2) - fms(0) )**2
    beta(3) = 13.0d0/12.0d0*( fms(1) - 2.0d0*fms(0) + fms(-1) )**2 &
            + 0.25d0*( 3.0d0*fms(1) - 4.0d0*fms(0) + fms(-1) )**2
    
    tao5 = abs(beta(1) - beta(3))
    
    alpha(:) = omega(:)*( 1.0d0 + tao5/(beta(:)+eps) )
    totalpha = alpha(1) + alpha(2) + alpha(3)
    
    if(totalpha < eps) then
        fmx = fu(2)
    else
        fomega(:) = alpha(:)/totalpha
        fmx = fu(1)*fomega(1) + fu(2)*fomega(2) + fu(3)*fomega(3)
    end if
    
    return
end function fmx
