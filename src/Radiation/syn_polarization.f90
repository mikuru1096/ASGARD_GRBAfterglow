!f2py: skip
module syn_polarization
    use constants
    implicit none
    private

    public :: synchrotron_fg_kernel, synchrotron_polarized_components

contains

! 计算同步辐射 F(x)、G(x) 核。
! Compute synchrotron F(x) and G(x) kernels.
!
! 小/大 x 用解析渐近式，中间区间用 u 变量 Gauss 积分。
! Use analytic asymptotes at small/large x and Gauss integration over u in between.
subroutine synchrotron_fg_kernel(x,Fx,Gx)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: Fx,Gx
    real(8), parameter :: gam23=1.3541179394264004169d0
    real(8) :: xp,kernel,corrf,corrg

    if (x <= 0d0) error stop "synchrotron_fg_kernel requires x > 0."

    if (x < 1d-4) then
        xp=x**(1d0/3d0)
        Fx=2d0**(2d0/3d0)*gam23*xp
        Gx=2d0**(-1d0/3d0)*gam23*xp
    else if (x > 50d0) then
        kernel=dsqrt(pi*x/2d0)*dexp(-x)
        corrf=1d0+55d0/(72d0*x)
        corrg=1d0+7d0/(72d0*x)
        Fx=kernel*corrf
        Gx=kernel*corrg
    else
        call fg_integral(x,Fx,Gx)
    end if
end subroutine synchrotron_fg_kernel

! 计算单粒子的两个线偏振分量: (F+G)/2 和 (F-G)/2。
! Compute the 2 linear polarization components for a single particle.
subroutine synchrotron_polarized_components(x,P_perp_kernel,P_parallel_kernel)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: P_perp_kernel,P_parallel_kernel
    real(8) :: Fx,Gx

    call synchrotron_fg_kernel(x,Fx,Gx)
    P_perp_kernel=0.5d0*(Fx+Gx)
    P_parallel_kernel=0.5d0*(Fx-Gx)
end subroutine synchrotron_polarized_components

! 在有限 u 区间上做复合 4 点 Gauss 积分。
! Composite 4-point Gauss integration over a finite u interval.
subroutine fg_integral(x,Fx,Gx)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: Fx,Gx
    integer, parameter :: n_cell=96
    real(8), parameter :: q1=0.3399810435848562648d0,q2=0.8611363115940525752d0
    real(8), parameter :: w1=0.6521451548625461426d0,w2=0.3478548451374538574d0
    integer :: i
    real(8) :: umax,du,umid,uhalf,sumf,sumg

    umax=dlog((80d0/x)+dsqrt((80d0/x)*(80d0/x)+1d0))
    du=umax/dble(n_cell)
    uhalf=0.5d0*du
    sumf=0d0
    sumg=0d0

    do i=1,n_cell
        umid=(dble(i)-0.5d0)*du
        call add_node(x,umid-uhalf*q2,w2*uhalf,sumf,sumg)
        call add_node(x,umid-uhalf*q1,w1*uhalf,sumf,sumg)
        call add_node(x,umid+uhalf*q1,w1*uhalf,sumf,sumg)
        call add_node(x,umid+uhalf*q2,w2*uhalf,sumf,sumg)
    end do

    Fx=x*sumf
    Gx=x*sumg
end subroutine fg_integral

! 累加一个 Gauss 节点上的 F/G 被积函数值。
! Accumulate F/G integrands at a Gauss node.
subroutine add_node(x,u,weight,sumf,sumg)
    implicit none
    real(8), intent(in) :: x,u,weight
    real(8), intent(inout) :: sumf,sumg
    real(8) :: arg,eval,coshu

    coshu=dcosh(u)
    arg=-x*coshu
    eval=dexp(arg)
    sumf=sumf+weight*eval*dcosh(5d0*u/3d0)/coshu
    sumg=sumg+weight*eval*dcosh(2d0*u/3d0)
end subroutine add_node

end module syn_polarization
