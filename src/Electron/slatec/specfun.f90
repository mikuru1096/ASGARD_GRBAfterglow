module specfun
    implicit none
    private

    public :: gammauic, besselk, besselk0, besselk1

    real(8), private, parameter :: euler_gamma = 0.5772156649015328

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

    function besselk(nu, z, n) result(val)
        implicit none
        real(8), intent(in) :: nu, z
        integer, intent(in) :: n
        real(8) :: val(n)

        integer :: nzero
        integer, parameter :: KODE = 1

        call dbesk(z, nu, KODE, n, &
            val, nzero)
        
        if (nzero .gt. 0) then
            print *, "Encountered underflow (nu, z, n) = (", nu, z, n, &
                "), some outputs are zero."
        end if
    end function besselk

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

end module specfun