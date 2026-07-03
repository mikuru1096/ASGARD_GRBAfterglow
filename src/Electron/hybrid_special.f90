module hybrid_special
    implicit none
    private

    public :: gamma_uic, bessel_k0, bessel_k1

    ! SLATEC 外部函数接口：只在这里保留原始符号名。
    ! SLATEC external interfaces: keep the original symbol names isolated here.
    interface
        double precision function dgamic(a, x)
            implicit none
            real(8), intent(in) :: a, x
        end function dgamic

        double precision function dbesk0(x)
            implicit none
            real(8), intent(in) :: x
        end function dbesk0

        double precision function dbesk1(x)
            implicit none
            real(8), intent(in) :: x
        end function dbesk1
    end interface

contains

    ! 上不完全 gamma 函数 Γ(s,x)，转发到 SLATEC dgamic。
    ! Upper incomplete gamma Γ(s,x), forwarded to SLATEC dgamic.
    function gamma_uic(s, x) result(val)
        implicit none
        real(8), intent(in) :: s, x
        real(8) :: val

        val = dgamic(s, x)
    end function gamma_uic

    ! 修正 Bessel K0，转发到 SLATEC dbesk0。
    ! Modified Bessel K0, forwarded to SLATEC dbesk0.
    function bessel_k0(z) result(val)
        implicit none
        real(8), intent(in) :: z
        real(8) :: val

        val = dbesk0(z)
    end function bessel_k0

    ! 修正 Bessel K1，转发到 SLATEC dbesk1。
    ! Modified Bessel K1, forwarded to SLATEC dbesk1.
    function bessel_k1(z) result(val)
        implicit none
        real(8), intent(in) :: z
        real(8) :: val

        val = dbesk1(z)
    end function bessel_k1

end module hybrid_special
