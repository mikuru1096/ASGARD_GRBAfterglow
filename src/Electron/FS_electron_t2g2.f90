!Calculate the electron distributions of forward shock.
!Second-order spatial accuracy with upwind scheme (Beam-Warming)
!Time discretization: first-order implicit Euler
!****************************************************************************************
!******************************* main program *******************************************
!****************************************************************************************
subroutine fs_electron_t2g2(Boundary,R_Tobs,R_Gamma,R,V_seed,n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger,n_threads, &
                                   gam_e,dN_gam_e,P_syn,Seed_syn,V_m,V_c,V_a)
    !$ use omp_lib
    use constants
    use get_Y
    IMPLICIT REAL(8)(A-H,O-Z)
    integer, intent(in) :: n,Num_nu,Num_R,Num_gam_e,index_Y,index_syn_intger
    real(8), intent(in) :: Boundary(n),R_Tobs(Num_R),R_Gamma(Num_R),R(Num_R),V_seed(Num_nu)
    real(8), intent(out) :: dN_gam_e(Num_gam_e,Num_R),gam_e(Num_gam_e),P_syn(Num_nu,Num_R), &
                            Seed_syn(Num_nu,Num_R), V_m(Num_R), V_c(Num_R), V_a(Num_R)
    
    ! Main work arrays
    real(8),allocatable,dimension (:) :: dEl,dF1,para_minus_gam_e_p,dot_gam_e_SSA, &
                                         dN_x,Compton,dot_gam_e
    real(8),allocatable,dimension (:) :: v_eff, dEL_mean
    real(8),allocatable,dimension (:) :: a_diag, b_diag, c_diag, rhs, sol
    
    ! Allocate work arrays
    allocate (dEl(Num_gam_e), dF1(Num_gam_e), para_minus_gam_e_p(Num_gam_e), &
              dot_gam_e_SSA(Num_gam_e), dN_x(Num_gam_e), Compton(Num_gam_e), &
              dot_gam_e(Num_gam_e), v_eff(Num_gam_e), dEL_mean(Num_gam_e))
    
    Eta_0=Boundary(1)
    R_ini=Boundary(4)
    Epsilon_e=Boundary(5)
    Epsilon_b=Boundary(6)
    p=Boundary(7)
    z=Boundary(8)
    dNe_ISM=Boundary(11)
    A_star=Boundary(12)
    E_iso=Boundary(14)
    T_log10_duration=Boundary(15)
    f_e=Boundary(16)
    R_tr=Boundary(21)
    f_jump=Boundary(22)
    f_wide=Boundary(23)
    R0=Boundary(n)
    
    P_syn=zero
    Seed_syn=zero
    V_m=zero
    V_c=zero
    V_a=zero

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
    if (p<2.01 .and. p>=2.0) then
        Gam_e_m=0.01d0/1.01d0*Epsilon_e/f_e*1836d0*(R_Gamma(1)-one)+one
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
    
    dN_x=dN_gam_e(:,1)*gam_e*dlog(ten)
    d_x=dlog10(gam_e(2)/gam_e(1))
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

        beta_Gam=sqrt(one-one/R_Gamma_loc**2)
        f_r=(1.35d-19)/beta_Gam/R_Gamma_loc*DB**2/pi
        dDR=0.1/(f_r*Gam_e_max+1.333/(R(I_tobs)+R(I_tobs-1)))
        
        dDD=R(I_tobs)-R(I_tobs-1)
        L1=max(100,min(1000,Int(dDD/dDR)))
        dDR=dDD/L1
        CFL=dDR/d_x
        
        dN_x=dN_gam_e(:,I_tobs-1)*gam_e*dlog(ten)
        
        V_m(I_tobs-1)=4.2d6*DB*Gam_e_m*Gam_e_m/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        V_c(I_tobs-1)=4.2d6*DB*Gam_e_c*Gam_e_c/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))
        
        call get_nu_a(R_loc,DB,Num_gam_e,gam_e,dN_gam_e(:,I_tobs-1), temp_val)
        V_a(I_tobs-1)=temp_val/(R_Gamma_loc*(1d0-beta_Gam)*(one+z))

        select case(index_syn_intger)
        case(1)
            call get_syn(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e(:,I_tobs-1),V_seed, &
                         P_syn(:,I_tobs),Seed_syn(:,I_tobs))
        case(2)
            call get_syn_simpson(R_loc,DB,Num_gam_e,Num_nu,n_threads,gam_e,dN_gam_e(:,I_tobs-1),V_seed, &
                         P_syn(:,I_tobs),Seed_syn(:,I_tobs))
        end select
        
        call get_SSA_numerical(DB,Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn(:,I_tobs), dot_gam_e_SSA)
        
        select case(index_Y)
        case(1)
            call get_IC_numerical(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,Seed_syn(:,I_tobs), dot_gam_e)
            dEl=(f_r+(dot_gam_e-dot_gam_e_SSA)/beta_Gam/R_Gamma_loc/para_c)*gam_e
        case(2)
            call get_Y_Nakar(Num_gam_e,Num_nu,n_threads,gam_e,V_seed,P_syn(:,I_tobs), Compton)
            Q=4d0*pi*R_loc*R_loc*para_c
            Compton=one+Compton/Q/(4d0*R_Gamma_loc*R_Gamma_loc*dNe*Para_m_p_E)
            Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
            dEl=(f_r*Compton-dot_gam_e_SSA/beta_Gam/R_Gamma_loc/para_c)*gam_e
        case(3)
            call get_Y_Fan(Epsilon_e,Epsilon_b,p,DB,Gam_e_m,Gam_e_c,Gam_e_max,Num_gam_e,gam_e, Compton)
            Compton=one+Compton
            Gam_e_max=Gam_e_max/sqrt(Compton(Num_gam_e))
            dEl=(f_r*Compton-dot_gam_e_SSA/beta_Gam/R_Gamma_loc/para_c)*gam_e
        case default
            print*, 'invalid Compton case, check your chosen model!'
            stop
        end select
        
        ! Effective cooling velocity v = dEl / gamma_e
        v_eff = dEl/gam_e
        
        do L=1,L1
            R_loc=R_loc+dDR
            
            if (A_star > zero) then
                dNe_wind=A_star*3.0d35/R_loc**2
                if (dNe_wind <= dNe_ISM/4d0) then
                    dNe=dNe_ISM
                else
                    dNe=dNe_wind
                end if
            else
                dNe=dNe_ISM*(1.0+(f_jump-1d0)*exp(-(log10(R_loc)-log10(R_tr))**2/(2*f_wide*f_wide)))
            end if
        
            if (R_loc<R0) then
                dNe=A_star*3.0d35/R0**2
            end if
            
            Q=4d0/3d0*pi*(3d0*R_loc**2+dDR*(3d0*R_loc+dDR))*dNe*f_e*Gam_e_m_p
            dF1=zero
            where(gam_e<Gam_e_max .and. gam_e>Gam_e_m) dF1=Q*para_minus_gam_e_p

            ! ============== Beam-Warming upwind solver ==============
            allocate(a_diag(Num_gam_e), b_diag(Num_gam_e), c_diag(Num_gam_e), &
                     rhs(Num_gam_e), sol(Num_gam_e))
            
            a_diag = zero
            b_diag = zero
            c_diag = zero
            rhs = zero
            sol = zero
            
            ! ---------- Boundary conditions ----------
            ! i=1: Dirichlet boundary, n_1 = 0
            b_diag(1) = 1.0d0
            rhs(1) = 0.0d0
            
            ! i=2: first-order upwind because i=0 is unavailable
            ! Equation: n_2 + dr/dx*v_2*n_2 - dr/dx*v_1*n_1 + dr/R*n_2 = n_2_old + dr*S_2
            ! Rearranged: (1 + CFL*v_2 + dDR/R) * n_2 - CFL*v_1*n_1 = n_2_old + dr*S_2
            a_diag(2) = -CFL * v_eff(1)  ! n_1 coefficient
            b_diag(2) = 1.0d0 + CFL * (v_eff(2) + one/R_loc/log(ten))  ! n_2 coefficient
            rhs(2) = dN_x(2) + dDR * dF1(2)
            
            ! ---------- Interior points: second-order upwind ----------
            ! For the cooling case (v < 0), use the Beam-Warming stencil.
            ! d(v*n)/dx ~= (3*v_i*n_i - 4*v_{i-1}*n_{i-1} + v_{i-2}*n_{i-2}) / (2*dx)
            ! Treat the n_{i-2} term explicitly with values from the previous step.
            do i = 3, Num_gam_e-1
                a_diag(i) = -(2.0d0 * CFL) * v_eff(i-1)  ! n_{i-1} coefficient
                b_diag(i) = 1.0d0 + (1.5d0 * CFL) * (v_eff(i) + one/R_loc/log(ten))  ! n_i coefficient
                rhs(i) = dN_x(i) + dDR * dF1(i) - (0.5d0 * CFL) * v_eff(i-2) * dN_x(i-2)
            end do
            
            ! ---------- High-energy boundary ----------
            ! i=Num_gam_e: first-order upwind at high-energy boundary
            a_diag(Num_gam_e) = -CFL * v_eff(Num_gam_e-1)
            b_diag(Num_gam_e) = 1.0d0 + CFL * (v_eff(Num_gam_e) + one/R_loc/log(ten))
            rhs(Num_gam_e) = dN_x(Num_gam_e) + dDR * dF1(Num_gam_e)
            
            ! ---------- Solve the tridiagonal system ----------
            ! solve tridiagonal system
            call thomas_solve_tridiag(Num_gam_e, a_diag, b_diag, c_diag, rhs, sol)
            
            ! Enforce non-negative solution
            where(sol < zero) sol = zero
            
            dN_x = sol
            
            ! Release local tridiagonal buffers
            deallocate(a_diag, b_diag, c_diag, rhs, sol)

            if (L1 == L) then
                dN_gam_e(:,I_tobs)=dN_x/gam_e/dlog(ten)
            end if
        end do
    end do

    deallocate (dEl, dF1, para_minus_gam_e_p, dot_gam_e_SSA, dN_x, &
                Compton, dot_gam_e, v_eff, dEL_mean)

    return
    
contains

    ! Thomas algorithm for a tridiagonal linear system
    subroutine thomas_solve_tridiag(n, a, b, c, d, x)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: a(n), b(n), c(n), d(n)
        real(8), intent(out) :: x(n)
        real(8), allocatable :: cp(:), dp(:)
        integer :: i
        
        allocate(cp(n), dp(n))
        
        ! Forward elimination
        cp(1) = c(1) / b(1)
        dp(1) = d(1) / b(1)
        
        do i = 2, n
            cp(i) = c(i) / (b(i) - a(i) * cp(i-1))
            dp(i) = (d(i) - a(i) * dp(i-1)) / (b(i) - a(i) * cp(i-1))
        end do
        
        ! Back substitution
        x(n) = dp(n)
        do i = n-1, 1, -1
            x(i) = dp(i) - cp(i) * x(i+1)
        end do
        
        deallocate(cp, dp)
    end subroutine thomas_solve_tridiag

end subroutine fs_electron_t2g2
