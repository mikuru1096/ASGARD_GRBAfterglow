! Reverse-shock hydro jump conditions used by the main RS dynamics and secondary density-jump branches.
module reverse_jump_conditions
    implicit none
    private
    public :: reverse_contact
contains

subroutine reverse_contact(gamma4,n1,n4,e4,p4,gamma_c,p3,gamma43,comp,beta_rs)
    use constants
    implicit none
    integer :: I
    real(8), intent(in) :: gamma4,n1,n4,e4,p4
    real(8), intent(out) :: gamma_c,p3,gamma43,comp,beta_rs
    real(8) :: lo,hi,mid,f_lo,f_hi,f_mid,beta_c,beta4

    if (gamma4 <= 1d0 .or. n1 <= 0d0 .or. n4 <= 0d0 .or. e4 <= 0d0 .or. p4 <= 0d0) &
        error stop 'reverse_contact requires positive hot upstream state'
    lo=1d0+1d-10; hi=gamma4*(1d0-1d-10)
    call pressure_difference(lo,f_lo)
    call pressure_difference(hi,f_hi)
    if (f_lo*f_hi > 0d0) error stop 'reverse_contact has no physical bracket'
    do I=1,80
        mid=0.5d0*(lo+hi)
        call pressure_difference(mid,f_mid)
        if (f_lo*f_mid <= 0d0) then
            hi=mid; f_hi=f_mid
        else
            lo=mid; f_lo=f_mid
        end if
    end do
    gamma_c=0.5d0*(lo+hi)
    beta_c=dsqrt(1d0-gamma_c**(-2)); beta4=dsqrt(1d0-gamma4**(-2))
    gamma43=gamma_c*gamma4*(1d0-beta_c*beta4)
    p3=p4+(4d0/3d0)*(gamma43*gamma43-1d0)*(e4+p4)
    call solve_comp(comp,beta_rs)

contains

    subroutine pressure_difference(gamma_trial, diff)
    implicit none
    real(8), intent(in) :: gamma_trial
    real(8), intent(out) :: diff
    real(8) :: beta_trial,beta_up,gamma_rel,p2_trial,p3_trial

        beta_trial=dsqrt(1d0-gamma_trial**(-2)); beta_up=dsqrt(1d0-gamma4**(-2))
        gamma_rel=gamma_trial*gamma4*(1d0-beta_trial*beta_up)
        p2_trial=(4d0/3d0)*(gamma_trial*gamma_trial-1d0)*n1*Para_m_p*Para_c**2
        p3_trial=p4+(4d0/3d0)*(gamma_rel*gamma_rel-1d0)*(e4+p4)
        diff=p2_trial-p3_trial
    end subroutine pressure_difference

    subroutine solve_comp(comp,brs_out)
    implicit none
    real(8), intent(out) :: comp,brs_out
    integer :: K
    real(8) :: bs_lo,bs_hi,bs_mid,g_lo,g_hi,g_mid,eps_beta,comp1

        if (beta4 <= beta_c) error stop 'reverse_contact reverse shock speed bracket is empty'
        brs_out=0d0
        eps_beta=1d-14*(beta4-beta_c)
        bs_lo=beta_c+eps_beta
        bs_hi=beta4-eps_beta
        call unity_comp(bs_lo,bs_hi,comp1)
        bs_hi=comp1-eps_beta
        call shock_diff(bs_lo,g_lo,comp)
        call shock_diff(bs_hi,g_hi,comp)
        if (g_lo*g_hi > 0d0) then
            comp=-1d0
            return
        end if
        do K=1,80
            bs_mid=0.5d0*(bs_lo+bs_hi)
            call shock_diff(bs_mid,g_mid,comp)
            if (g_lo*g_mid <= 0d0) then
                bs_hi=bs_mid; g_hi=g_mid
            else
                bs_lo=bs_mid; g_lo=g_mid
            end if
        end do
        brs_out=0.5d0*(bs_lo+bs_hi)
        call shock_diff(brs_out,g_mid,comp)
    end subroutine solve_comp

    subroutine unity_comp(beta_lo,beta_hi,beta_one)
    implicit none
    real(8), intent(in) :: beta_lo,beta_hi
    real(8), intent(out) :: beta_one
    integer :: K
    real(8) :: lo_c,hi_c,mid_c,c_lo,c_mid,comp_tmp,diff_tmp

        lo_c=beta_lo; hi_c=beta_hi
        call shock_diff(lo_c,diff_tmp,comp_tmp)
        c_lo=comp_tmp-1d0
        do K=1,80
            mid_c=0.5d0*(lo_c+hi_c)
            call shock_diff(mid_c,diff_tmp,comp_tmp)
            c_mid=comp_tmp-1d0
            if (c_lo*c_mid <= 0d0) then
                hi_c=mid_c
            else
                lo_c=mid_c; c_lo=c_mid
            end if
        end do
        beta_one=0.5d0*(lo_c+hi_c)
    end subroutine unity_comp

    subroutine shock_diff(beta_s,diff,comp)
    implicit none
    real(8), intent(in) :: beta_s
    real(8), intent(out) :: diff,comp
    real(8) :: beta4_s,beta3_s,gamma4_s,gamma3_s,u4_s,u3_s,w4,w3,e3_trial

        beta4_s=(beta4-beta_s)/(1d0-beta4*beta_s)
        beta3_s=(beta_c-beta_s)/(1d0-beta_c*beta_s)
        gamma4_s=1d0/dsqrt(1d0-beta4_s*beta4_s)
        gamma3_s=1d0/dsqrt(1d0-beta3_s*beta3_s)
        u4_s=gamma4_s*dabs(beta4_s)
        u3_s=gamma3_s*dabs(beta3_s)
        comp=u4_s/u3_s
        e3_trial=3d0*p3
        w4=n4*Para_m_p*Para_c**2+e4+p4
        w3=comp*n4*Para_m_p*Para_c**2+e3_trial+p3
        diff=(w4*u4_s*u4_s+p4)-(w3*u3_s*u3_s+p3)
    end subroutine shock_diff
end subroutine reverse_contact

end module reverse_jump_conditions
