module hybrid_special_functions
    implicit none
    private

    public :: gammauic, besselk0, besselk1

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

    function gammauic(s, x) result(val)
        implicit none
        real(8), intent(in) :: s, x
        real(8) :: val

        val = dgamic(s, x)
    end function gammauic

    function besselk0(z) result(val)
        implicit none
        real(8), intent(in) :: z
        real(8) :: val

        val = dbesk0(z)
    end function besselk0

    function besselk1(z) result(val)
        implicit none
        real(8), intent(in) :: z
        real(8) :: val

        val = dbesk1(z)
    end function besselk1

end module hybrid_special_functions
