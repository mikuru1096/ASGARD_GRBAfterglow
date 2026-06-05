! f2py: skip
! public: integral_thermal1, integral_thermal12, integral_cp, solve_theta, normalized_hybrid_spec, normalized_hybrid_spec_lg

module hybrid_spectrum_kernel_fast
   use specfun, only: gammauic, besselk0, besselk1
   implicit none
   private
   public :: integral_thermal1, integral_thermal12, &
      integral_cpl, solve_theta, &
      normalized_hybrid_spec, normalized_hybrid_spec_lg
   
   real(8), private, parameter :: m700 = -7.0d2
   real(8), private :: global_gmin, global_lnc

   ! coefficients for integral_thermal
   real(8), private, parameter :: coeff0 = -1.0d0
   real(8), private, parameter :: coeff1 = 1.0d0/2.0d0
   real(8), private, parameter :: coeff2 = 1.0d0/8.0d0
   real(8), private, parameter :: coeff3 = 1.0d0/1.6d1
   real(8), private, parameter :: coeff4 = 5.0d0/1.28d2
   real(8), private, parameter :: coeff5 = 7.0d0/2.56d2
   real(8), private, parameter :: coeff6 = 2.1d1/1.024d3
   real(8), private, parameter :: coeff7 = 3.3d1/2.048d3

   ! constants for integral_thermal
   real(8), private, parameter :: d2 = 1.0d0/2.0d0
   real(8), private, parameter :: d3 = 1.0d0/3.0d0
   real(8), private, parameter :: d4 = 1.0d0/4.0d0
   real(8), private, parameter :: d5 = 1.0d0/5.0d0
   real(8), private, parameter :: d6 = 1.0d0/6.0d0
   real(8), private, parameter :: d7 = 1.0d0/7.0d0
   real(8), private, parameter :: d8 = 1.0d0/8.0d0
   real(8), private, parameter :: d9 = 1.0d0/9.0d0
   real(8), private, parameter :: d10 = 1.0d0/1.0d1
   real(8), private, parameter :: d11 = 1.0d0/1.1d1
   real(8), private, parameter :: d2x3 = 1.0d0/6.0d0
   real(8), private, parameter :: d4x5 = 1.0d0/2.0d1
   real(8), private, parameter :: d6x7 = 1.0d0/4.2d1
   real(8), private, parameter :: d8x9 = 1.0d0/7.2d1
   real(8), private, parameter :: d10x11 = 1.0d0/1.1d2

   ! accuracy for integral_thermal
   real(8), private, parameter :: tol = 1.0d-15
   real(8), private, parameter :: floor_eps = 1.0d-30

   real(8), private, parameter :: ln_10 = log(1.0d1)
   
   contains

   subroutine integral_thermal1(gamma_min, theta, &
      it1)
      implicit none
      real(8), intent(in) :: gamma_min
      real(8), intent(in) :: theta
      real(8), intent(out) :: it1

      real(8) :: inv_theta, inv_theta2, theta2
      real(8) :: gdt, inv_gdt, inv_gdt2, exp_mgdt
      real(8) :: theta_term, gdt_term, uig_term, cor_term
      real(8) :: besk0, besk1

      if (gamma_min < 1.0d0) then
         error stop "Error in integral_thermal1: gamma_min must be >= 1"
      end if
      if (theta <= 0.0d0) then
         error stop "Error in integral_thermal1: theta must be > 0"
      end if

      theta2 = theta*theta
      inv_theta = 1.0d0/theta
      inv_theta2 = inv_theta*inv_theta
      gdt = gamma_min*inv_theta
      exp_mgdt = exp(-gdt)
      inv_gdt = 1.0d0/gdt
      inv_gdt2 = inv_gdt*inv_gdt

      besk0 = besselk0(inv_theta)
      besk1 = besselk1(inv_theta)
      it1 = 2.0d0*theta*besk1 + besk0

      ! 0 order correction
      uig_term = (2.0d0 + 2.0d0*gdt + gdt*gdt)*exp_mgdt ! Gamma(3;gdt,\infty)
      cor_term = coeff0*theta2*uig_term ! (n=1, k=0) -1 * theta^{2} * Gamma(3;gdt,\infty)
      it1 = it1 + cor_term
      
      ! 1 order correction
      uig_term = exp_mgdt ! Gamma(1;gdt,\infty)
      cor_term = coeff1*uig_term ! (n=1, k=0) +1/2 * Gamma(1;gdt,\infty)
      it1 = it1 + cor_term

      ! 2 order correction
      theta_term = inv_theta2 ! = theta^{-2}
      uig_term = gammauic(-1.0d0, gdt) ! Gamma(-1;gdt,\infty)
      cor_term = coeff2*theta_term*uig_term ! (n=1, k=0) +1/8 * theta^{-2} * Gamma(-1;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term
      
      ! 3 order correction
      theta_term = theta_term*inv_theta2 ! = theta^{-4}
      gdt_term = inv_gdt*inv_gdt2 ! = gdt^{-3}
      uig_term = d2x3*uig_term + (1.0d0 - d2*gdt) * d3*gdt_term * exp_mgdt ! Gamma(-3;gdt,\infty) = 1/6*Gamma(-1;gdt,\infty) + (1-1/2*gdt) 1/3*gdt^{-3} exp(-gdt)
      cor_term = coeff3*theta_term*uig_term ! (n=1, k=0) +1/16 * theta^{-4} * Gamma(-3;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term
      
      ! 4 order correction
      theta_term = theta_term*inv_theta2 ! = theta^{-6}
      gdt_term = gdt_term*inv_gdt2 ! = gdt^{-5}
      uig_term = d4x5*uig_term + (1.0d0 - d4*gdt) * d5*gdt_term * exp_mgdt ! Gamma(-5;gdt,\infty) = 1/20*Gamma(-3;gdt,\infty) + (1-1/4*gdt) 1/5*gdt^{-5} exp(-gdt)
      cor_term = coeff4*theta_term*uig_term ! (n=1, k=0) +5/128 * theta^{-6} * Gamma(-5;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term

      ! 5 order correction
      theta_term = theta_term*inv_theta2 ! = theta^{-8}
      gdt_term = gdt_term*inv_gdt2 ! = gdt^{-7}
      uig_term = d6x7*uig_term + (1.0d0 - d6*gdt) * d7*gdt_term * exp_mgdt ! Gamma(-7;gdt,\infty) = 1/42*Gamma(-5;gdt,\infty) + (1-1/6*gdt) 1/7*gdt^{-7} exp(-gdt)
      cor_term = coeff5*theta_term*uig_term ! (n=1, k=0) +7/256 * theta^{-8} * Gamma(-7;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term

      ! 6 order correction
      theta_term = theta_term*inv_theta2 ! = theta^{-10}
      gdt_term = gdt_term*inv_gdt2 ! = gdt^{-9}
      uig_term = d8x9*uig_term + (1.0d0 - d8*gdt) * d9*gdt_term * exp_mgdt ! Gamma(-9;gdt,\infty) = 1/72*Gamma(-7;gdt,\infty) + (1-1/8*gdt) 1/9*gdt^{-9} exp(-gdt)
      cor_term = coeff6*theta_term*uig_term ! (n=1, k=0) +21/1024 * theta^{-10} * Gamma(-9;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term

      ! 7 order correction
      theta_term = theta_term*inv_theta2 ! = theta^{-12}
      gdt_term = gdt_term*inv_gdt2 ! = gdt^{-11}
      uig_term = d10x11*uig_term + (1.0d0 - d10*gdt) * d11*gdt_term * exp_mgdt ! Gamma(-11;gdt,\infty) = 1/110*Gamma(-9;gdt,\infty) + (1-1/10*gdt) 1/11*gdt^{-11} exp(-gdt)
      cor_term = coeff7*theta_term*uig_term ! (n=1, k=0) +33/2048 * theta^{-12} * Gamma(-11;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term
      
      810 continue
      it1 = it1 * theta
   end subroutine integral_thermal1

   subroutine integral_thermal12(gamma_min, theta, &
      it1, it2)
      implicit none
      real(8), intent(in) :: gamma_min
      real(8), intent(in) :: theta
      real(8), intent(out) :: it1
      real(8), intent(out) :: it2

      real(8) :: theta2, inv_theta, inv_theta2
      real(8) :: gdt2, gdt, inv_gdt, exp_mgdt
      real(8) :: besk0, besk1
      real(8) :: theta_term, gdt_term, uig_term, cxt_term, cor1_term, cor2_term

      if (gamma_min < 1.0d0) then
         error stop "Error in integral_thermal12: gamma_min must be >= 1"
      end if
      if (theta <= 0.0d0) then
         error stop "Error in integral_thermal12: theta must be > 0"
      end if

      theta2 = theta*theta
      inv_theta = 1.0d0/theta
      inv_theta2 = inv_theta*inv_theta
      gdt = gamma_min*inv_theta
      gdt2 = gdt*gdt
      exp_mgdt = exp(-gdt)
      inv_gdt = 1.0d0/gdt

      besk0 = besselk0(inv_theta)
      besk1 = besselk1(inv_theta)
      it1 = 2.0d0*theta*besk1 + besk0
      it2 = (6.0d0*theta + inv_theta)*besk1 + 3.0d0*besk0

      ! 0 order correction
      cxt_term = coeff0*theta2
      uig_term = (6.0d0 + 6.0d0*gdt + 3.0d0*gdt2 + gdt2*gdt)*exp_mgdt ! Gamma(4;gdt,\infty)
      cor2_term = cxt_term*uig_term ! (n=2, k=0) -1 * theta^{2} * Gamma(4;gdt,\infty)
      uig_term = (2.0d0 + 2.0d0*gdt + gdt2)*exp_mgdt ! Gamma(3;gdt,\infty)
      cor1_term = cxt_term*uig_term ! (n=1, k=0) -1 * theta^{2} * Gamma(3;gdt,\infty)
      it1 = it1 + cor1_term
      it2 = it2 + cor2_term
      
      ! 1 order correction
      cxt_term = coeff1
      uig_term = (1.0d0 + gdt)*exp_mgdt ! Gamma(2;gdt,\infty)
      cor2_term = cxt_term*uig_term ! (n=2, k=0) +1/2 * Gamma(2;gdt,\infty)
      uig_term = exp_mgdt ! Gamma(1;gdt,\infty)
      cor1_term = cxt_term*uig_term ! (n=1, k=0) +1/2 * Gamma(1;gdt,\infty)
      it1 = it1 + cor1_term
      it2 = it2 + cor2_term

      ! 2 order correction
      theta_term = inv_theta2 ! = theta^{-2}
      gdt_term = inv_gdt ! = gdt^{-1}
      cxt_term = coeff2*theta_term
      uig_term = gammauic(0.0d0, gdt) ! Gamma(0;gdt,\infty)
      cor2_term = cxt_term*uig_term ! (n=2, k=0) +1/8 * theta^{-2} * Gamma(0;gdt,\infty)
      uig_term = -uig_term + gdt_term*exp_mgdt ! Gamma(-1;gdt,\infty) = -Gamma(0;gdt,\infty) + gdt^{-1} exp(-gdt)
      cor1_term = cxt_term*uig_term ! (n=1, k=0) +1/8 * theta^{-2} * Gamma(-1;gdt,\infty)
      if ( &
         (cor2_term < tol * max(it2, floor_eps)) &
         .and. &
         (cor1_term < tol * max(it1, floor_eps)) &
         ) &
         goto 810
      it1 = it1 + cor1_term
      it2 = it2 + cor2_term
      
      ! 3 order correction
      theta_term = theta_term*inv_theta2 ! = theta^{-4}
      gdt_term = gdt_term*inv_gdt ! = gdt^{-2}
      cxt_term = coeff3*theta_term ! = +1/16 * theta^{-4}
      uig_term = -d2*uig_term + d2*gdt_term*exp_mgdt ! Gamma(-2;gdt,\infty) = -1/2 * Gamma(-1;gdt,\infty) + 1/2*gdt^{-2} exp(-gdt)
      cor2_term = cxt_term*uig_term ! (n=2, k=0) +1/16 * theta^{-4} * Gamma(-2;gdt,\infty)
      gdt_term = gdt_term*inv_gdt ! = gdt^{-3}
      uig_term = -d3*uig_term + d3*gdt_term*exp_mgdt ! Gamma(-3;gdt,\infty) = -1/3 * Gamma(-2;gdt,\infty) + 1/3*gdt^{-3} exp(-gdt)
      cor1_term = cxt_term*uig_term ! (n=1, k=0) +1/16 * theta^{-4} * Gamma(-3;gdt,\infty)
      if ( &
         (cor2_term < tol * max(it2, floor_eps)) &
         .and. &
         (cor1_term < tol * max(it1, floor_eps)) &
         ) &
         goto 810
      it1 = it1 + cor1_term
      it2 = it2 + cor2_term
      
      ! 4 order correction
      theta_term = theta_term*inv_theta2 ! = theta^{-6}
      gdt_term = gdt_term*inv_gdt ! = gdt^{-4}
      cxt_term = coeff4*theta_term ! = +5/128 * theta^{-6}
      uig_term = -d4*uig_term + d4*gdt_term*exp_mgdt ! Gamma(-4;gdt,\infty) = -1/4 * Gamma(-3;gdt,\infty) + 1/4*gdt^{-4} exp(-gdt)
      cor2_term = cxt_term*uig_term ! (n=2, k=0) +5/128 * theta^{-6} * Gamma(-4;gdt,\infty)
      gdt_term = gdt_term*inv_gdt ! = gdt^{-5}
      uig_term = -d5*uig_term + d5*gdt_term*exp_mgdt ! Gamma(-5;gdt,\infty) = -1/5 * Gamma(-4;gdt,\infty) + 1/5*gdt^{-5} exp(-gdt)
      cor1_term = cxt_term*uig_term ! (n=1, k=0) +5/128 * theta^{-6} * Gamma(-5;gdt,\infty)
      if ( &
         (cor2_term < tol * max(it2, floor_eps)) &
         .and. &
         (cor1_term < tol * max(it1, floor_eps)) &
         ) &
         goto 810
      it1 = it1 + cor1_term
      it2 = it2 + cor2_term
         
      ! 5 order correction
      theta_term = theta_term*inv_theta2 ! = theta^{-8}
      gdt_term = gdt_term*inv_gdt ! = gdt^{-6}
      cxt_term = coeff5*theta_term ! = +7/256 * theta^{-8}
      uig_term = -d6*uig_term + d6*gdt_term*exp_mgdt ! Gamma(-6;gdt,\infty) = -1/6 * Gamma(-5;gdt,\infty) + 1/6*gdt^{-6} exp(-gdt)
      cor2_term = cxt_term*uig_term ! (n=2, k=0) +7/256 * theta^{-8} * Gamma(-6;gdt,\infty)
      gdt_term = gdt_term*inv_gdt ! = gdt^{-7}
      uig_term = -d7*uig_term + d7*gdt_term*exp_mgdt ! Gamma(-7;gdt,\infty) = -1/7 * Gamma(-6;gdt,\infty) + 1/7*gdt^{-7} exp(-gdt)
      cor1_term = cxt_term*uig_term ! (n=1, k=0) +7/256 * theta^{-8} * Gamma(-7;gdt,\infty)
      if ( &
         (cor2_term < tol * max(it2, floor_eps)) &
         .and. &
         (cor1_term < tol * max(it1, floor_eps)) &
         ) &
         goto 810
      it1 = it1 + cor1_term
      it2 = it2 + cor2_term
         
      ! 6 order correction
      theta_term = theta_term*inv_theta2 ! = theta^{-10}
      gdt_term = gdt_term*inv_gdt ! = gdt^{-8}
      cxt_term = coeff6*theta_term ! = +21/1024 * theta^{-10}
      uig_term = -d8*uig_term + d8*gdt_term*exp_mgdt ! Gamma(-8;gdt,\infty) = -1/8 * Gamma(-7;gdt,\infty) + 1/8*gdt^{-8} exp(-gdt)
      cor2_term = cxt_term*uig_term ! (n=2, k=0) +21/1024 * theta^{-10} * Gamma(-8;gdt,\infty)
      gdt_term = gdt_term*inv_gdt ! = gdt^{-9}
      uig_term = -d9*uig_term + d9*gdt_term*exp_mgdt ! Gamma(-9;gdt,\infty) = -1/9 * Gamma(-8;gdt,\infty) + 1/9*gdt^{-9} exp(-gdt)
      cor1_term = cxt_term*uig_term ! (n=1, k=0) +21/1024 * theta^{-10} * Gamma(-9;gdt,\infty)
      if ( &
         (cor2_term < tol * max(it2, floor_eps)) &
         .and. &
         (cor1_term < tol * max(it1, floor_eps)) &
         ) &
         goto 810
      it1 = it1 + cor1_term
      it2 = it2 + cor2_term
         
      ! 7 order correction
      theta_term = theta_term*inv_theta2 ! = theta^{-12}
      gdt_term = gdt_term*inv_gdt ! = gdt^{-10}
      cxt_term = coeff7*theta_term ! = +33/2048 * theta^{-12}
      uig_term = -d10*uig_term + d10*gdt_term*exp_mgdt ! Gamma(-10;gdt,\infty) = -1/10 * Gamma(-9;gdt,\infty) + 1/10*gdt^{-10} exp(-gdt)
      cor2_term = cxt_term*uig_term ! (n=2, k=0) +33/2048 * theta^{-12} * Gamma(-10;gdt,\infty)
      gdt_term = gdt_term*inv_gdt ! = gdt^{-11}
      uig_term = -d11*uig_term + d11*gdt_term*exp_mgdt ! Gamma(-11;gdt,\infty) = -1/11 * Gamma(-10;gdt,\infty) + 1/11*gdt^{-11} exp(-gdt)
      cor1_term = cxt_term*uig_term ! (n=1, k=0) +33/2048 * theta^{-12} * Gamma(-11;gdt,\infty)
      it1 = it1 + cor1_term
      it2 = it2 + cor2_term
      
      810 continue
      it1 = it1 * theta
      it2 = it2 * theta2
      
   end subroutine integral_thermal12

   subroutine integral_cpl(p, gamma_min, gamma_max, &
      val)
      implicit none
      real(8), intent(in) :: p, gamma_min, gamma_max
      real(8), intent(out) :: val
      real(8) :: s, a, ln_p1, ln_p2
      
      s = 1.0d0 - p
      a = gamma_min / gamma_max

      ln_p1 = s*log(gamma_max)
      ln_p2 = log(gammauic(s, a))

      val = exp(ln_p1 + ln_p2)

   end subroutine

   subroutine solve_theta(p, gamma_min, gamma_max, xi_e, &
      theta, flag)
      implicit none
      real(8), intent(in) :: p, gamma_min, gamma_max, xi_e
      real(8), intent(out) :: theta
      integer, intent(out) :: flag
      real(8) :: ln_c, int, init_theta, gamma_peak

      flag = 0
      theta = -1.0d0
      if (gamma_max <= gamma_min) then
         flag = flag + 8
         print *, "hybrid_spectrum_kernel.solve_theta: gamma_max <= gamma_min, please check input parameters"
         return
      end if

      call integral_cpl(p, gamma_min, gamma_max, &
         int)
      
      ln_c = log((1.0d0/xi_e - 1.0d0)*sqrt(gamma_min*gamma_min-1.0d0)) &
         + (1.0d0+p)*log(gamma_min) &
         + log(int) &
         + gamma_min/gamma_max
      global_lnc = ln_c
      global_gmin = gamma_min

      init_theta = get_initial_theta()
      if (init_theta < 0.0d0) then
         flag = flag + 2
         return
      end if

      theta = newton_method(init_theta, 1.0d-6, 50)
      if (theta < 0.0d0) then
         flag = flag + 4
         return
      end if

      gamma_peak = get_thermal_peak(theta)
      if (gamma_min < gamma_peak) then
         flag = flag + 1
      end if

   end subroutine solve_theta

   function sfunc(theta) result(val)
      ! the optimizing function
      ! we want to find a theta, which meets
      ! sfunc(theta) == 0.0
      implicit none
      real(8), intent(in) :: theta
      real(8) :: it1, val

      call integral_thermal1(global_gmin, theta, it1)
      val = log(it1) - global_lnc + global_gmin / theta

   end function sfunc

   function get_initial_theta() result(init_theta)
      ! only can be called from solve_theta
      ! return the initial theta for the Newton method
      ! if failed to find an initial value,
      ! the returned the value is -1.0
      implicit none
      real(8) :: init_theta, guess_theta
      integer, parameter :: max_iter = 6
      real(8), parameter :: ctol = 1.0d-6
      integer :: iter, i
      logical :: valid
      real(8) :: k(6), sval

      init_theta = 0.144d0 * global_gmin
      return 
      k = [ &
         4.0d-2, 5.166199d-2, &
         ! 6.672402d-2, 8.617739d-2, &
         1.1130238d-1, 1.4375255d-1, &
         ! 1.8566355d-1, 2.397937d-1, &
         3.0970547d-1, 4.0d-1 &
      ]

      valid = .false.
      
      ! since the sfunc decrease as theta grows
      ! it is convinient to exclude the general interval
      ! by calculating the value of sfunc
      ! at the lower and the upper boundaries of the general interval
      iter_loop: do iter = 1, max_iter
         ! theta is less than the general lower boundary
         guess_theta = k(1) * global_gmin
         sval = sfunc(guess_theta)
         ! extremely lucky
         if (abs(sval) < ctol) then
            init_theta = guess_theta
            valid = .true.
            exit iter_loop
         end if
         if (sval < 0.0d0) then
            k = k * 1.0d-1
            cycle iter_loop
         end if
         
         ! theta is greater than the general upper boundary
         guess_theta = k(6) * global_gmin
         sval = sfunc(guess_theta)
         ! extremely lucky
         if (abs(sval) < ctol) then
            init_theta = guess_theta
            valid = .true.
            exit iter_loop
         end if
         if (sval > 0.0d0) then
            k = k * 1.0d1
            cycle
         end if
         
         ! found the general interval where the theta locates at
         ! scan the grid to find the delicate interval
         grid_loop: do i = 2, 5
            guess_theta = k(i) * global_gmin
            sval = sfunc(guess_theta)
            ! extremely lucky
            if (abs(sval) < ctol) then
               init_theta = guess_theta
               valid = .true.
               exit iter_loop
            end if
            if (sval < 0.0d0) then
               init_theta = 0.5d0 * (guess_theta + k(i - 1) * global_gmin)
               valid = .true.
               exit iter_loop
            end if
            
            guess_theta = k(7 - i) * global_gmin
            sval = sfunc(guess_theta)
            ! extremely lucky
            if (abs(sval) < ctol) then
               init_theta = guess_theta
               valid = .true.
               exit iter_loop
            end if
            if (sval > 0.0d0) then
               init_theta = 0.5d0 * (guess_theta + k(8 - i) * global_gmin)
               valid = .true.
               exit iter_loop
            end if

            ! finally, the theta falls between k(5) and k(6)
            init_theta = 0.5d0 * (k(3) + k(4)) * global_gmin
            valid = .true.
            exit iter_loop
            
         end do grid_loop

      end do iter_loop
      
      if (.not. valid) then
         init_theta = -1.0d0
      end if

   end function get_initial_theta

   function get_initial_theta_10() result(init_theta)
      ! only can be called from solve_theta
      ! return the initial theta for the Newton method
      ! if failed to find an initial value,
      ! the returned the value is -1.0
      implicit none
      real(8) :: init_theta, guess_theta
      integer, parameter :: max_iter = 6
      real(8), parameter :: ctol = 1.0d-6
      integer :: iter, i
      logical :: valid
      real(8) :: k(10), sval

      k = [ &
         4.0d-2, 5.166199d-2, &
         6.672402d-2, 8.617739d-2, &
         1.1130238d-1, 1.4375255d-1, &
         1.8566355d-1, 2.397937d-1, &
         3.0970547d-1, 4.0d-1 &
      ]

      valid = .false.
      
      ! since the sfunc decrease as theta grows
      ! it is convinient to exclude the general interval
      ! by calculating the value of sfunc
      ! at the lower and the upper boundaries of the general interval
      iter_loop: do iter = 1, max_iter
         ! theta is less than the general lower boundary
         guess_theta = k(1) * global_gmin
         sval = sfunc(guess_theta)
         ! extremely lucky
         if (abs(sval) < ctol) then
            init_theta = guess_theta
            valid = .true.
            exit iter_loop
         end if
         if (sval < 0.0d0) then
            k = k * 1.0d-1
            cycle iter_loop
         end if
         
         ! theta is greater than the general upper boundary
         guess_theta = k(10) * global_gmin
         sval = sfunc(guess_theta)
         ! extremely lucky
         if (abs(sval) < ctol) then
            init_theta = guess_theta
            valid = .true.
            exit iter_loop
         end if
         if (sval > 0.0d0) then
            k = k * 1.0d1
            cycle
         end if
         
         ! found the general interval where the theta locates at
         ! scan the grid to find the delicate interval
         grid_loop: do i = 2, 5
            guess_theta = k(i) * global_gmin
            sval = sfunc(guess_theta)
            ! extremely lucky
            if (abs(sval) < ctol) then
               init_theta = guess_theta
               valid = .true.
               exit iter_loop
            end if
            if (sval < 0.0d0) then
               init_theta = 0.5d0 * (guess_theta + k(i - 1) * global_gmin)
               valid = .true.
               exit iter_loop
            end if
            
            guess_theta = k(11 - i) * global_gmin
            sval = sfunc(guess_theta)
            ! extremely lucky
            if (abs(sval) < ctol) then
               init_theta = guess_theta
               valid = .true.
               exit iter_loop
            end if
            if (sval > 0.0d0) then
               init_theta = 0.5d0 * (guess_theta + k(12 - i) * global_gmin)
               valid = .true.
               exit iter_loop
            end if

            ! finally, the theta falls between k(5) and k(6)
            init_theta = 0.5d0 * (k(5) + k(6)) * global_gmin
            valid = .true.
            exit iter_loop
            
         end do grid_loop

      end do iter_loop
      
      if (.not. valid) then
         init_theta = -1.0d0
      end if

   end function get_initial_theta_10

   function newton_method(init_theta, rtol, max_iter) result(best_theta)
      implicit none
      real(8), intent(in) :: init_theta, rtol
      integer, intent(in) :: max_iter
      real(8) :: best_theta, theta
      real(8) :: inv_theta, gdt, it1, it2, rel_shift
      integer :: iter
  
      theta = init_theta
      best_theta = -1.0d0

      ! Newton iteration method
      ! \Theta_{n+1} = \Theta_{n} - f(\Theta_{n}) / f'(\Theta_{n})
      iter_loop: do iter = 1, max_iter
         call integral_thermal12(global_gmin, theta, &
            it1, it2)
         inv_theta = 1.0d0/theta
         gdt = global_gmin*inv_theta

         ! relative shift for each iter
         ! \Theta_{n+1}/\Theta_{n} - 1.0d0 = - f(\Theta_{n}) / f'(\Theta_{n}) / \Theta_{n}
         rel_shift = gdt + log(it1) - global_lnc
         rel_shift = rel_shift / (gdt - inv_theta*it2/it1) 
         if (abs(rel_shift) < rtol) then
            best_theta = theta * (1.0d0 + rel_shift)
            exit iter_loop
         end if

         ! not converged yet
         theta = theta * (1.0d0 + rel_shift)
      end do iter_loop
  
   end function newton_method

   subroutine normalized_hybrid_spec(n_gamma, gamma, p, gamma_min, gamma_max, xi_e, &
      spec)
      integer, intent(in) :: n_gamma
      real(8), intent(in) :: gamma(n_gamma)
      real(8), intent(in) :: p, gamma_min, gamma_max, xi_e
      real(8), intent(out) :: spec(n_gamma)

      real(8) :: theta, it1, int, thermal_constant, cpl_constant, g, ln_val, inv_gmax, inv_theta
      integer :: flag, i
      
      call solve_theta(p, gamma_min, gamma_max, xi_e, &
         theta, flag)

      if (flag > 1) then
         print *, p, gamma_min, gamma_max, xi_e, theta, flag
         spec = -1.0d0
         return 
      end if

      call integral_thermal1(gamma_min, theta, &
         it1)
      call integral_cpl(p, gamma_min, gamma_max, &
         int)
      thermal_constant = log((1.0d0 - xi_e) / it1)
      cpl_constant = log(xi_e/int)
      inv_theta = 1.0d0/theta
      inv_gmax = 1.0d0/gamma_max
      do i = 1, n_gamma
         g = gamma(i)

         if (g > 1.0d0 .and. g < gamma_min) then
            ln_val = log(g * sqrt(g*g - 1.0d0)) - g*inv_theta + thermal_constant
         else if (g >= gamma_min) then
            ln_val = -p * log(g) - g*inv_gmax + cpl_constant
         else
            spec(i) = 0.0d0
            cycle
         end if

         if (ln_val > m700) then
            spec(i) = exp(ln_val)
         else
            spec(i) = 0.0d0
         end if

      end do

   end subroutine normalized_hybrid_spec

   subroutine benchmark_normalized_hybrid_spec(n_gamma, gamma, p, gamma_min, gamma_max, xi_e, &
      spec)
      implicit none
      integer, intent(in) :: n_gamma
      real(8), intent(in) :: gamma(n_gamma)
      real(8), intent(in) :: p, gamma_min, gamma_max, xi_e
      real(8), intent(out) :: spec(n_gamma)
   
      integer i
      real(8) :: t0, t1
   
      call cpu_time(t0)
      do i=1, 10000
         call normalized_hybrid_spec(n_gamma, gamma, p, gamma_min, gamma_max, xi_e, &
            spec)
      end do
      call cpu_time(t1)
      print *, "fast loop=10000, ", t1-t0
   
   end subroutine benchmark_normalized_hybrid_spec

   subroutine normalized_hybrid_spec_lg(n_gamma, gamma, p, gamma_min, gamma_max, xi_e, &
      spec)
      integer, intent(in) :: n_gamma
      real(8), intent(in) :: gamma(n_gamma)
      real(8), intent(in) :: p, gamma_min, gamma_max, xi_e
      real(8), intent(out) :: spec(n_gamma)

      call normalized_hybrid_spec(n_gamma, gamma, p, gamma_min, gamma_max, xi_e, &
         spec)
      
      spec = spec * gamma*ln_10
   end subroutine normalized_hybrid_spec_lg

   function thermal(n_gamma, gamma, gamma_min, theta, xi_e) result(spec)
      implicit none
      integer, intent(in) :: n_gamma
      real(8), intent(in) :: gamma(n_gamma) 
      real(8), intent(in) :: gamma_min, theta, xi_e
      
      real(8) :: spec(n_gamma)
      real(8) :: ln_expr(n_gamma)
      real(8) :: it_1_val
      logical :: mask(n_gamma) 
      real(8), parameter :: tiny_val = 1.0d-30
      real(8) :: scalar_constant
      integer :: i
      real(8) :: g, ln_val

      call integral_thermal1(gamma_min, theta, &
         it_1_val)
      
      spec = 0.0d0
      ln_expr = m700
      mask = (gamma > 1.0d0) .and. (gamma < gamma_min)
      scalar_constant = log((1.0d0 - xi_e) / max(it_1_val, 1d-60))

      !$OMP SIMD
      do i = 1, n_gamma
         g = gamma(i)
         if (g > 1.0d0 .and. g < gamma_min) then
            ln_val = log(g * sqrt(g*g - 1.0d0)) - g/theta + scalar_constant
            if (ln_val > -700.0d0) then
               spec(i) = exp(ln_val)
            else
               spec(i) = 0.0d0
            end if
         else
            spec(i) = 0.0d0
         end if
      end do
      
   end function thermal

   function cpl(n_gamma, gamma, p, gamma_min, gamma_max, xi_e) result(spec)
      implicit none
      integer, intent(in) :: n_gamma
      real(8), intent(in) :: gamma(n_gamma) 
      real(8), intent(in) :: p, gamma_min, gamma_max, xi_e
      
      real(8) :: spec(n_gamma)
      real(8) :: ln_expr(n_gamma)
      real(8) :: int_p_val
      logical :: mask(n_gamma) 
      real(8), parameter :: tiny_val = 1.0d-30
      real(8) :: scalar_constant

      integer :: i
      real(8) :: g, ln_val

      call integral_cpl(p, gamma_min, gamma_max, &
         int_p_val)
      
      spec = 0.0d0
      ln_expr = m700
      mask = gamma >= gamma_min
      scalar_constant = log(xi_e / max(int_p_val, 1d-60))

      !$OMP SIMD
      do i = 1, n_gamma
         g = gamma(i)
         if (g >= gamma_min) then
            ln_val = -p * log(g) - g/gamma_max + scalar_constant
            if (ln_val > -700.0d0) then
               spec(i) = exp(ln_val)
            else
               spec(i) = 0.0d0
            end if
         else
            spec(i) = 0.0d0
         end if
      end do
      
   end function cpl

   function get_thermal_peak(theta) result(gamma_peak)
      implicit none
      real(8), intent(in) :: theta
      real(8) :: gamma_peak, theta2, r, phi

      theta2 = theta * theta
      r = sqrt(3.0 + 4.0 * theta2)
      phi = acos(0.5 * theta * (16.0 * theta2 - 9.0) / (r * r * r))
      gamma_peak = 2.0d0 / 3.0d0 * (r * cos(phi / 3.0) + theta)

   end function get_thermal_peak

end module hybrid_spectrum_kernel_fast