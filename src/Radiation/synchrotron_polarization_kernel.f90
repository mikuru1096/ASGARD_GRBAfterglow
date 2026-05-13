!f2py: skip
module synchrotron_polarization_kernel
    use constants
    implicit none
    private

    public :: synchrotron_fg_kernel, synchrotron_polarized_components

contains

! 计算同步辐射F(x)、G(x)核：小/大x用解析渐近式，中间区间用u变量Gauss积分。
subroutine synchrotron_fg_kernel(x,Fx,Gx)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: Fx,Gx
    real(8), parameter :: gamma_two_thirds=1.3541179394264004169d0
    real(8) :: xp,common,large_corr_f,large_corr_g

    if (x <= zero) error stop "synchrotron_fg_kernel requires x > 0."

    if (x < 1d-4) then
        xp=x**(one/3d0)
        Fx=two**(2d0/3d0)*gamma_two_thirds*xp
        Gx=two**(-one/3d0)*gamma_two_thirds*xp
    else if (x > 50d0) then
        common=dsqrt(pi*x/two)*dexp(-x)
        large_corr_f=one+55d0/(72d0*x)
        large_corr_g=one+7d0/(72d0*x)
        Fx=common*large_corr_f
        Gx=common*large_corr_g
    else
        call synchrotron_fg_integral(x,Fx,Gx)
    end if
end subroutine synchrotron_fg_kernel

! 计算单粒子两个线偏振分量：(F+G)/2 和 (F-G)/2。
subroutine synchrotron_polarized_components(x,P_perp_kernel,P_parallel_kernel)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: P_perp_kernel,P_parallel_kernel
    real(8) :: Fx,Gx

    call synchrotron_fg_kernel(x,Fx,Gx)
    P_perp_kernel=0.5d0*(Fx+Gx)
    P_parallel_kernel=0.5d0*(Fx-Gx)
end subroutine synchrotron_polarized_components

! 在有限u区间上做复合4点Gauss积分；上限由exp(-x cosh u)的尾部指数确定。
subroutine synchrotron_fg_integral(x,Fx,Gx)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: Fx,Gx
    integer, parameter :: n_cell=96
    real(8), parameter :: q1=0.3399810435848562648d0,q2=0.8611363115940525752d0
    real(8), parameter :: w1=0.6521451548625461426d0,w2=0.3478548451374538574d0
    integer :: i
    real(8) :: u_max,du,u_mid,u_half,sum_f,sum_g

    u_max=dlog((80d0/x)+dsqrt((80d0/x)*(80d0/x)+one))
    du=u_max/dble(n_cell)
    u_half=0.5d0*du
    sum_f=zero
    sum_g=zero

    do i=1,n_cell
        u_mid=(dble(i)-0.5d0)*du
        call add_fg_node(x,u_mid-u_half*q2,w2*u_half,sum_f,sum_g)
        call add_fg_node(x,u_mid-u_half*q1,w1*u_half,sum_f,sum_g)
        call add_fg_node(x,u_mid+u_half*q1,w1*u_half,sum_f,sum_g)
        call add_fg_node(x,u_mid+u_half*q2,w2*u_half,sum_f,sum_g)
    end do

    Fx=x*sum_f
    Gx=x*sum_g
end subroutine synchrotron_fg_integral

! 累加一个Gauss节点上的F/G被积函数值。
subroutine add_fg_node(x,u,weight,sum_f,sum_g)
    implicit none
    real(8), intent(in) :: x,u,weight
    real(8), intent(inout) :: sum_f,sum_g
    real(8) :: exp_arg,exp_val,cosh_u

    cosh_u=dcosh(u)
    exp_arg=-x*cosh_u
    exp_val=dexp(exp_arg)
    sum_f=sum_f+weight*exp_val*dcosh(5d0*u/3d0)/cosh_u
    sum_g=sum_g+weight*exp_val*dcosh(2d0*u/3d0)
end subroutine add_fg_node

end module synchrotron_polarization_kernel
