! Reverse-shock hydro jump conditions used by the main RS dynamics and secondary density-jump branches.
module reverse_jump_conditions
    implicit none
    private
    public :: secondary_reverse_contact_rh
contains

subroutine secondary_reverse_contact_rh(gamma4,n1,n4,e4,p4,gamma_c,p3,gamma43,n3_over_n4,beta_rs)
    use constants
    implicit none
    integer :: I
    real(8), intent(in) :: gamma4,n1,n4,e4,p4
    real(8), intent(out) :: gamma_c,p3,gamma43,n3_over_n4,beta_rs
    real(8) :: lo,hi,mid,f_lo,f_hi,f_mid,beta_c,beta4

    if (gamma4 <= one .or. n1 <= zero .or. n4 <= zero .or. e4 <= zero .or. p4 <= zero) &
        error stop 'secondary_reverse_contact_rh requires positive hot upstream state'
    lo=one+1d-10; hi=gamma4*(one-1d-10)
    call pressure_difference(lo,f_lo)
    call pressure_difference(hi,f_hi)
    if (f_lo*f_hi > zero) error stop 'secondary_reverse_contact_rh has no physical bracket'
    do I=1,80
        mid=0.5d0*(lo+hi)
        call pressure_difference(mid,f_mid)
        if (f_lo*f_mid <= zero) then
            hi=mid; f_hi=f_mid
        else
            lo=mid; f_lo=f_mid
        end if
    end do
    gamma_c=0.5d0*(lo+hi)
    beta_c=dsqrt(one-gamma_c**(-2)); beta4=dsqrt(one-gamma4**(-2))
    gamma43=gamma_c*gamma4*(one-beta_c*beta4)
    p3=p4+(4d0/3d0)*(gamma43*gamma43-one)*(e4+p4)
    call solve_reverse_shock_compression(n3_over_n4,beta_rs)

contains

    subroutine pressure_difference(gamma_trial, diff)
    implicit none
    real(8), intent(in) :: gamma_trial
    real(8), intent(out) :: diff
    real(8) :: beta_trial,beta_up,gamma_rel,p2_trial,p3_trial

        beta_trial=dsqrt(one-gamma_trial**(-2)); beta_up=dsqrt(one-gamma4**(-2))
        gamma_rel=gamma_trial*gamma4*(one-beta_trial*beta_up)
        p2_trial=(4d0/3d0)*(gamma_trial*gamma_trial-one)*n1*Para_m_p*Para_c**2
        p3_trial=p4+(4d0/3d0)*(gamma_rel*gamma_rel-one)*(e4+p4)
        diff=p2_trial-p3_trial
    end subroutine pressure_difference

    subroutine solve_reverse_shock_compression(comp,beta_rs_out)
    implicit none
    real(8), intent(out) :: comp,beta_rs_out
    integer :: K
    real(8) :: beta_s_lo,beta_s_hi,beta_s_mid,g_lo,g_hi,g_mid,eps_beta,beta_comp_one

        if (beta4 <= beta_c) error stop 'secondary_reverse_contact_rh reverse shock speed bracket is empty'
        beta_rs_out=zero
        eps_beta=1d-14*(beta4-beta_c)
        beta_s_lo=beta_c+eps_beta
        beta_s_hi=beta4-eps_beta
        call solve_compression_unity_speed(beta_s_lo,beta_s_hi,beta_comp_one)
        beta_s_hi=beta_comp_one-eps_beta
        call shock_momentum_difference(beta_s_lo,g_lo,comp)
        call shock_momentum_difference(beta_s_hi,g_hi,comp)
        if (g_lo*g_hi > zero) then
            comp=-one
            return
        end if
        do K=1,80
            beta_s_mid=0.5d0*(beta_s_lo+beta_s_hi)
            call shock_momentum_difference(beta_s_mid,g_mid,comp)
            if (g_lo*g_mid <= zero) then
                beta_s_hi=beta_s_mid; g_hi=g_mid
            else
                beta_s_lo=beta_s_mid; g_lo=g_mid
            end if
        end do
        beta_rs_out=0.5d0*(beta_s_lo+beta_s_hi)
        call shock_momentum_difference(beta_rs_out,g_mid,comp)
    end subroutine solve_reverse_shock_compression

    subroutine solve_compression_unity_speed(beta_lo,beta_hi,beta_one)
    implicit none
    real(8), intent(in) :: beta_lo,beta_hi
    real(8), intent(out) :: beta_one
    integer :: K
    real(8) :: lo_c,hi_c,mid_c,c_lo,c_mid,comp_tmp,diff_tmp

        lo_c=beta_lo; hi_c=beta_hi
        call shock_momentum_difference(lo_c,diff_tmp,comp_tmp)
        c_lo=comp_tmp-one
        do K=1,80
            mid_c=0.5d0*(lo_c+hi_c)
            call shock_momentum_difference(mid_c,diff_tmp,comp_tmp)
            c_mid=comp_tmp-one
            if (c_lo*c_mid <= zero) then
                hi_c=mid_c
            else
                lo_c=mid_c; c_lo=c_mid
            end if
        end do
        beta_one=0.5d0*(lo_c+hi_c)
    end subroutine solve_compression_unity_speed

    subroutine shock_momentum_difference(beta_s,diff,comp)
    implicit none
    real(8), intent(in) :: beta_s
    real(8), intent(out) :: diff,comp
    real(8) :: beta4_s,beta3_s,gamma4_s,gamma3_s,u4_s,u3_s,w4,w3,e3_trial

        beta4_s=(beta4-beta_s)/(one-beta4*beta_s)
        beta3_s=(beta_c-beta_s)/(one-beta_c*beta_s)
        gamma4_s=one/dsqrt(one-beta4_s*beta4_s)
        gamma3_s=one/dsqrt(one-beta3_s*beta3_s)
        u4_s=gamma4_s*dabs(beta4_s)
        u3_s=gamma3_s*dabs(beta3_s)
        comp=u4_s/u3_s
        e3_trial=3d0*p3
        w4=n4*Para_m_p*Para_c**2+e4+p4
        w3=comp*n4*Para_m_p*Para_c**2+e3_trial+p3
        diff=(w4*u4_s*u4_s+p4)-(w3*u3_s*u3_s+p3)
    end subroutine shock_momentum_difference
end subroutine secondary_reverse_contact_rh

end module reverse_jump_conditions
