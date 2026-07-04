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

        double precision function dbsk0e(x)
            implicit none
            real(8), intent(in) :: x
        end function dbsk0e

        double precision function dbsk1e(x)
            implicit none
            real(8), intent(in) :: x
        end function dbsk1e
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

    ! 修正 Bessel K0，用 SLATEC scaled K0 避免大 z 下错误栈中止。
    ! Modified Bessel K0, evaluated through SLATEC scaled K0.
    function bessel_k0(z) result(val)
        implicit none
        real(8), intent(in) :: z
        real(8) :: val

        val = exp(-z)*dbsk0e(z)
    end function bessel_k0

    ! 修正 Bessel K1，用 SLATEC scaled K1 避免大 z 下错误栈中止。
    ! Modified Bessel K1, evaluated through SLATEC scaled K1.
    function bessel_k1(z) result(val)
        implicit none
        real(8), intent(in) :: z
        real(8) :: val

        val = exp(-z)*dbsk1e(z)
    end function bessel_k1

end module hybrid_special
