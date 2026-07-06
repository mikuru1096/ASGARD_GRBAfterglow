! f2py: skip
! public: hybrid_coord, hybrid_thermal_coord

module hybrid_spectrum
   use hybrid_special, only: gamma_uic, bessel_k0, bessel_k1
   use electron_coord_common, only: coord_fourvel, coord_from_xg, gamma_from_coord, dxg_dcoord
   implicit none
   private
   public :: hybrid_coord, hybrid_thermal_coord

   real(8), private, parameter :: m700 = -7.0d2

   ! 热核渐近修正系数。
   ! Coefficients for the thermal-kernel asymptotic corrections.
   real(8), private, parameter :: coeff0 = -1.0d0, coeff1 = 1.0d0/2.0d0, coeff2 = 1.0d0/8.0d0
   real(8), private, parameter :: coeff3 = 1.0d0/1.6d1, coeff4 = 5.0d0/1.28d2, coeff5 = 7.0d0/2.56d2
   real(8), private, parameter :: coeff6 = 2.1d1/1.024d3, coeff7 = 3.3d1/2.048d3

   ! Gamma 递推中反复使用的倒数系数。
   ! Reused reciprocal factors in the Gamma-function recurrences.
   real(8), private, parameter :: d2 = 1.0d0/2.0d0, d3 = 1.0d0/3.0d0, d4 = 1.0d0/4.0d0, d5 = 1.0d0/5.0d0
   real(8), private, parameter :: d6 = 1.0d0/6.0d0, d7 = 1.0d0/7.0d0, d8 = 1.0d0/8.0d0, d9 = 1.0d0/9.0d0
   real(8), private, parameter :: d10 = 1.0d0/1.0d1, d11 = 1.0d0/1.1d1, d2x3 = 1.0d0/6.0d0
   real(8), private, parameter :: d4x5 = 1.0d0/2.0d1, d6x7 = 1.0d0/4.2d1, d8x9 = 1.0d0/7.2d1
   real(8), private, parameter :: d10x11 = 1.0d0/1.1d2

   ! 热积分级数截断精度。
   ! Accuracy controls for truncating the thermal integral series.
   real(8), private, parameter :: tol = 1.0d-15, floor_eps = 1.0d-30

   real(8), private, parameter :: ln_10 = log(1.0d1)

   contains

   ! 计算 thermal 分支的 1 阶归一化积分。
   ! Compute the first thermal-branch normalization integral.
   subroutine thermal_int1(gamma_min, theta, &
      it1)
      implicit none
      real(8), intent(in) :: gamma_min, theta
      real(8), intent(out) :: it1

      real(8) :: inv_theta, inv_theta2, theta2, gdt, inv_gdt, inv_gdt2, exp_mgdt, theta_term, gdt_term
      real(8) :: uig_term, cor_term, besk0, besk1

      if (gamma_min < 1.0d0) then
         error stop "Error in thermal_int1: gamma_min must be >= 1"
      end if
      if (theta <= 0.0d0) then
         error stop "Error in thermal_int1: theta must be > 0"
      end if

      theta2 = theta*theta
      inv_theta = 1.0d0/theta
      inv_theta2 = inv_theta*inv_theta
      gdt = gamma_min*inv_theta
      exp_mgdt = exp(-gdt)
      inv_gdt = 1.0d0/gdt
      inv_gdt2 = inv_gdt*inv_gdt

      besk0 = bessel_k0(inv_theta)
      besk1 = bessel_k1(inv_theta)
      it1 = 2.0d0*theta*besk1 + besk0

      ! 0 阶修正。
      ! Zeroth-order correction.
      uig_term = (2.0d0 + 2.0d0*gdt + gdt*gdt)*exp_mgdt ! Gamma(3;gdt,\infty)
      cor_term = coeff0*theta2*uig_term ! (n=1, k=0) -1 * theta^{2} * Gamma(3;gdt,\infty)
      it1 = it1 + cor_term

      ! 1 阶修正。
      ! First-order correction.
      uig_term = exp_mgdt ! Gamma(1;gdt,\infty)
      cor_term = coeff1*uig_term ! (n=1, k=0) +1/2 * Gamma(1;gdt,\infty)
      it1 = it1 + cor_term

      ! 2 阶修正。
      ! Second-order correction.
      theta_term = inv_theta2 ! = theta^{-2}
      uig_term = gamma_uic(-1.0d0, gdt) ! Gamma(-1;gdt,\infty)
      cor_term = coeff2*theta_term*uig_term ! (n=1, k=0) +1/8 * theta^{-2} * Gamma(-1;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term

      ! 3 阶修正。
      ! Third-order correction.
      theta_term = theta_term*inv_theta2 ! = theta^{-4}
      gdt_term = inv_gdt*inv_gdt2 ! = gdt^{-3}
      ! Gamma(-3;gdt,\infty) = 1/6*Gamma(-1;gdt,\infty) + (1-1/2*gdt) 1/3*gdt^{-3} exp(-gdt)
      uig_term = d2x3*uig_term + (1.0d0 - d2*gdt) * d3*gdt_term * exp_mgdt
      cor_term = coeff3*theta_term*uig_term ! (n=1, k=0) +1/16 * theta^{-4} * Gamma(-3;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term

      ! 4 阶修正。
      ! Fourth-order correction.
      theta_term = theta_term*inv_theta2 ! = theta^{-6}
      gdt_term = gdt_term*inv_gdt2 ! = gdt^{-5}
      ! Gamma(-5;gdt,\infty) = 1/20*Gamma(-3;gdt,\infty) + (1-1/4*gdt) 1/5*gdt^{-5} exp(-gdt)
      uig_term = d4x5*uig_term + (1.0d0 - d4*gdt) * d5*gdt_term * exp_mgdt
      cor_term = coeff4*theta_term*uig_term ! (n=1, k=0) +5/128 * theta^{-6} * Gamma(-5;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term

      ! 5 阶修正。
      ! Fifth-order correction.
      theta_term = theta_term*inv_theta2 ! = theta^{-8}
      gdt_term = gdt_term*inv_gdt2 ! = gdt^{-7}
      ! Gamma(-7;gdt,\infty) = 1/42*Gamma(-5;gdt,\infty) + (1-1/6*gdt) 1/7*gdt^{-7} exp(-gdt)
      uig_term = d6x7*uig_term + (1.0d0 - d6*gdt) * d7*gdt_term * exp_mgdt
      cor_term = coeff5*theta_term*uig_term ! (n=1, k=0) +7/256 * theta^{-8} * Gamma(-7;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term

      ! 6 阶修正。
      ! Sixth-order correction.
      theta_term = theta_term*inv_theta2 ! = theta^{-10}
      gdt_term = gdt_term*inv_gdt2 ! = gdt^{-9}
      ! Gamma(-9;gdt,\infty) = 1/72*Gamma(-7;gdt,\infty) + (1-1/8*gdt) 1/9*gdt^{-9} exp(-gdt)
      uig_term = d8x9*uig_term + (1.0d0 - d8*gdt) * d9*gdt_term * exp_mgdt
      cor_term = coeff6*theta_term*uig_term ! (n=1, k=0) +21/1024 * theta^{-10} * Gamma(-9;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term

      ! 7 阶修正。
      ! Seventh-order correction.
      theta_term = theta_term*inv_theta2 ! = theta^{-12}
      gdt_term = gdt_term*inv_gdt2 ! = gdt^{-11}
      ! Gamma(-11;gdt,\infty) = 1/110*Gamma(-9;gdt,\infty) + (1-1/10*gdt) 1/11*gdt^{-11} exp(-gdt)
      uig_term = d10x11*uig_term + (1.0d0 - d10*gdt) * d11*gdt_term * exp_mgdt
      cor_term = coeff7*theta_term*uig_term ! (n=1, k=0) +33/2048 * theta^{-12} * Gamma(-11;gdt,\infty)
      if (cor_term < tol * max(it1, floor_eps)) goto 810
      it1 = it1 + cor_term

      810 continue
      it1 = it1 * theta
   end subroutine thermal_int1

   ! 同时计算 thermal 分支的 1 阶和 2 阶积分，供 Newton 导数使用。
   ! Compute the first and second thermal integrals together for the Newton derivative.
   subroutine thermal_int12(gamma_min, theta, &
      it1, it2)
      implicit none
      real(8), intent(in) :: gamma_min, theta
      real(8), intent(out) :: it1, it2

      real(8) :: theta2, inv_theta, inv_theta2, gdt2, gdt, inv_gdt, exp_mgdt, besk0, besk1, theta_term
      real(8) :: gdt_term, uig_term, cxt_term, cor1_term, cor2_term

      if (gamma_min < 1.0d0) then
         error stop "Error in thermal_int12: gamma_min must be >= 1"
      end if
      if (theta <= 0.0d0) then
         error stop "Error in thermal_int12: theta must be > 0"
      end if

      theta2 = theta*theta
      inv_theta = 1.0d0/theta
      inv_theta2 = inv_theta*inv_theta
      gdt = gamma_min*inv_theta
      gdt2 = gdt*gdt
      exp_mgdt = exp(-gdt)
      inv_gdt = 1.0d0/gdt

      besk0 = bessel_k0(inv_theta)
      besk1 = bessel_k1(inv_theta)
      it1 = 2.0d0*theta*besk1 + besk0
      it2 = (6.0d0*theta + inv_theta)*besk1 + 3.0d0*besk0

      ! 0 阶修正。
      ! Zeroth-order correction.
      cxt_term = coeff0*theta2
      uig_term = (6.0d0 + 6.0d0*gdt + 3.0d0*gdt2 + gdt2*gdt)*exp_mgdt ! Gamma(4;gdt,\infty)
      cor2_term = cxt_term*uig_term ! (n=2, k=0) -1 * theta^{2} * Gamma(4;gdt,\infty)
      uig_term = (2.0d0 + 2.0d0*gdt + gdt2)*exp_mgdt ! Gamma(3;gdt,\infty)
      cor1_term = cxt_term*uig_term ! (n=1, k=0) -1 * theta^{2} * Gamma(3;gdt,\infty)
      it1 = it1 + cor1_term
      it2 = it2 + cor2_term

      ! 1 阶修正。
      ! First-order correction.
      cxt_term = coeff1
      uig_term = (1.0d0 + gdt)*exp_mgdt ! Gamma(2;gdt,\infty)
      cor2_term = cxt_term*uig_term ! (n=2, k=0) +1/2 * Gamma(2;gdt,\infty)
      uig_term = exp_mgdt ! Gamma(1;gdt,\infty)
      cor1_term = cxt_term*uig_term ! (n=1, k=0) +1/2 * Gamma(1;gdt,\infty)
      it1 = it1 + cor1_term
      it2 = it2 + cor2_term

      ! 2 阶修正。
      ! Second-order correction.
      theta_term = inv_theta2 ! = theta^{-2}
      gdt_term = inv_gdt ! = gdt^{-1}
      cxt_term = coeff2*theta_term
      uig_term = gamma_uic(0.0d0, gdt) ! Gamma(0;gdt,\infty)
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

      ! 3 阶修正。
      ! Third-order correction.
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

      ! 4 阶修正。
      ! Fourth-order correction.
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

      ! 5 阶修正。
      ! Fifth-order correction.
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

      ! 6 阶修正。
      ! Sixth-order correction.
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

      ! 7 阶修正。
      ! Seventh-order correction.
      theta_term = theta_term*inv_theta2 ! = theta^{-12}
      gdt_term = gdt_term*inv_gdt ! = gdt^{-10}
      cxt_term = coeff7*theta_term ! = +33/2048 * theta^{-12}
      ! Gamma(-10;gdt,\infty) = -1/10 * Gamma(-9;gdt,\infty) + 1/10*gdt^{-10} exp(-gdt)
      uig_term = -d10*uig_term + d10*gdt_term*exp_mgdt
      cor2_term = cxt_term*uig_term ! (n=2, k=0) +33/2048 * theta^{-12} * Gamma(-10;gdt,\infty)
      gdt_term = gdt_term*inv_gdt ! = gdt^{-11}
      ! Gamma(-11;gdt,\infty) = -1/11 * Gamma(-10;gdt,\infty) + 1/11*gdt^{-11} exp(-gdt)
      uig_term = -d11*uig_term + d11*gdt_term*exp_mgdt
      cor1_term = cxt_term*uig_term ! (n=1, k=0) +33/2048 * theta^{-12} * Gamma(-11;gdt,\infty)
      it1 = it1 + cor1_term
      it2 = it2 + cor2_term

      810 continue
      it1 = it1 * theta
      it2 = it2 * theta2

   end subroutine thermal_int12

   ! 计算 cutoff power-law 分支归一化积分。
   ! Compute the cutoff power-law branch normalization integral.
   subroutine cpl_integral(p, gamma_min, gamma_max, &
      val)
      implicit none
      real(8), intent(in) :: p, gamma_min, gamma_max
      real(8), intent(out) :: val
      real(8) :: s, a, ln_p1, ln_p2

      s = 1.0d0 - p
      a = gamma_min / gamma_max

      ln_p1 = s*log(gamma_max)
      ln_p2 = log(gamma_uic(s, a))

      val = exp(ln_p1 + ln_p2)

   end subroutine

   ! 由非热电子比例求 hybrid thermal 温度。
   ! Solve the hybrid thermal temperature from the nonthermal electron fraction.
   subroutine solve_theta(p, gamma_min, gamma_max, xi_e, &
      theta)
      implicit none
      real(8), intent(in) :: p, gamma_min, gamma_max, xi_e
      real(8), intent(out) :: theta
      real(8) :: ln_c, int, init_theta

      if (gamma_max <= gamma_min) then
         error stop "hybrid_spectrum.solve_theta: gamma_max must exceed gamma_min."
      end if

      call cpl_integral(p, gamma_min, gamma_max, &
         int)

      ln_c = log((1.0d0/xi_e - 1.0d0)*sqrt(gamma_min*gamma_min-1.0d0)) &
         + (1.0d0+p)*log(gamma_min) &
         + log(int) &
         + gamma_min/gamma_max
      init_theta = get_initial_theta(gamma_min, ln_c)
      theta = newton_method(init_theta, 1.0d-6, 50, gamma_min, ln_c)

   end subroutine solve_theta

   function theta_objective(theta, gamma_min, ln_c) result(val)
      implicit none
      real(8), intent(in) :: theta, gamma_min, ln_c
      real(8) :: it1, val

      call thermal_int1(gamma_min, theta, &
         it1)
      val = log(it1) - ln_c + gamma_min/theta

   end function theta_objective

   function get_initial_theta(gamma_min, ln_c) result(init_theta)
      implicit none
      real(8), intent(in) :: gamma_min, ln_c
      real(8) :: init_theta, test_theta, test_val, theta_low, theta_hig, val_low, val_hig
      real(8), dimension(11) :: scale
      integer :: shift, idx_low, idx_hig, idx_mid
      integer, parameter :: max_shift = 8
      real(8), parameter :: close_tol = 1.0d-2

      scale = [ &
         6.0000000d-2, 7.5535520d-2, 9.5093590d-2, &
         1.1971574d-1, 1.5071319d-1, 1.8973666d-1, &
         2.3886430d-1, 3.0071234d-1, 3.7857441d-1, &
         4.7659694d-1, 6.0000000d-1 &
      ]

      theta_low = scale(1)*gamma_min
      val_low = theta_objective(theta_low, gamma_min, ln_c)
      if (abs(val_low) < close_tol) then
         init_theta = theta_low
         return
      end if

      do shift = 1, max_shift
         if (val_low >= 0.0d0) exit
         scale = scale*1.0d-1
         theta_low = scale(1)*gamma_min
         val_low = theta_objective(theta_low, gamma_min, ln_c)
         if (abs(val_low) < close_tol) then
            init_theta = theta_low
            return
         end if
      end do
      if (val_low < 0.0d0) then
         error stop "hybrid_spectrum.get_initial_theta: failed to bracket lower positive-temperature side."
      end if

      theta_hig = scale(11)*gamma_min
      val_hig = theta_objective(theta_hig, gamma_min, ln_c)
      if (abs(val_hig) < close_tol) then
         init_theta = theta_hig
         return
      end if
      do shift = 1, max_shift
         if (val_hig <= 0.0d0) exit
         theta_low = theta_hig
         val_low = val_hig
         scale = scale*1.0d1
         theta_hig = scale(11)*gamma_min
         val_hig = theta_objective(theta_hig, gamma_min, ln_c)
         if (abs(val_hig) < close_tol) then
            init_theta = theta_hig
            return
         end if
      end do
      if (val_hig > 0.0d0) then
         error stop "hybrid_spectrum.get_initial_theta: failed to bracket upper positive-temperature side."
      end if

      idx_low = 1
      idx_hig = 11
      do while (idx_hig - idx_low > 1)
         idx_mid = (idx_hig + idx_low)/2
         test_theta = scale(idx_mid)*gamma_min
         test_val = theta_objective(test_theta, gamma_min, ln_c)
         if (abs(test_val) < close_tol) then
            init_theta = test_theta
            return
         end if
         if (test_val > 0.0d0) then
            idx_low = idx_mid
            theta_low = test_theta
            val_low = test_val
         else
            idx_hig = idx_mid
            theta_hig = test_theta
            val_hig = test_val
         end if
      end do

      init_theta = theta_low - val_low*(theta_hig - theta_low)/(val_hig - val_low)

   end function get_initial_theta

   ! Newton 解 thermal 温度；状态全部经参数传入，保持 reentrant。
   ! Newton solve for the thermal temperature; all state is passed through arguments.
   function newton_method(init_theta, rtol, max_iter, gamma_min, ln_c) result(best_theta)
      implicit none
      real(8), intent(in) :: init_theta, rtol, gamma_min, ln_c
      integer, intent(in) :: max_iter
      real(8) :: best_theta, theta, inv_theta, gdt, it1, it2, rel_shift
      integer :: iter

      theta = init_theta
      best_theta = init_theta

      ! Newton 迭代。
      ! Newton iteration.
      ! \Theta_{n+1} = \Theta_{n} - f(\Theta_{n}) / f'(\Theta_{n})
      iter_loop: do iter = 1, max_iter
         call thermal_int12(gamma_min, theta, &
            it1, it2)
         inv_theta = 1.0d0/theta
         gdt = gamma_min*inv_theta

         ! 每步相对位移。
         ! Relative shift for each iteration.
         ! \Theta_{n+1}/\Theta_{n} - 1.0d0 = - f(\Theta_{n}) / f'(\Theta_{n}) / \Theta_{n}
         rel_shift = gdt + log(it1) - ln_c
         rel_shift = rel_shift / (gdt - inv_theta*it2/it1)
         if (abs(rel_shift) < rtol) then
            best_theta = theta * (1.0d0 + rel_shift)
            return
         end if

         ! 尚未收敛，推进到下一次迭代。
         ! Not converged yet; advance to the next iteration.
         theta = theta * (1.0d0 + rel_shift)
         if (theta <= 0.0d0) then
            error stop "hybrid_spectrum.newton_method: theta left the positive domain."
         end if
      end do iter_loop
      error stop "hybrid_spectrum.newton_method failed to converge."

   end function newton_method

   ! 在四速度坐标单元上构造归一化 hybrid electron spectrum。
   ! Build the normalized hybrid electron spectrum on four-velocity coordinate cells.
   subroutine hybrid_coord(n_gamma, coord_edge, coord_scale, p, gamma_min, gamma_max, xi_e, &
      spec)
      implicit none
      integer, intent(in) :: n_gamma
      real(8), intent(in), dimension(n_gamma+1) :: coord_edge
      real(8), intent(in) :: coord_scale, p, gamma_min, gamma_max, xi_e
      real(8), intent(out), dimension(n_gamma) :: spec

      real(8) :: theta, it1, int, thermal_constant, cpl_constant, inv_gmax, inv_theta
      real(8) :: cell_lo, cell_hi, dy_cell, y_min, seg_sum
      integer :: i

      call solve_theta(p, gamma_min, gamma_max, xi_e, &
         theta)

      call thermal_int1(gamma_min, theta, &
         it1)
      call cpl_integral(p, gamma_min, gamma_max, &
         int)
      thermal_constant = log((1.0d0 - xi_e) / it1)
      cpl_constant = log(xi_e/int)
      inv_theta = 1.0d0/theta
      inv_gmax = 1.0d0/gamma_max
      y_min = coord_from_xg(coord_fourvel, coord_scale, dlog(gamma_min))

      do i = 1, n_gamma
         cell_lo = coord_edge(i)
         cell_hi = coord_edge(i+1)
         dy_cell = cell_hi - cell_lo

         seg_sum = 0.0d0
         if (cell_lo < y_min) call add_segment(cell_lo, min(cell_hi, y_min), .true., seg_sum)
         if (cell_hi > y_min) call add_segment(max(cell_lo, y_min), cell_hi, .false., seg_sum)
         spec(i) = seg_sum/dy_cell
      end do

   contains

   ! 对一个四速度坐标子区间积分 dN/dy。
   ! Integrate dN/dy over one four-velocity coordinate sub-interval.
   subroutine add_segment(y_lo, y_hi, is_thermal, acc)
      implicit none
      logical, intent(in) :: is_thermal
      real(8), intent(in) :: y_lo, y_hi
      real(8), intent(inout) :: acc
      integer :: iq
      real(8) :: half_dy, y_mid, y_eval, gamma, density, jac
      real(8), parameter, dimension(3) :: xi=[-dsqrt(3d0/5d0),0d0,dsqrt(3d0/5d0)], wi=[5d0/9d0,8d0/9d0,5d0/9d0]

      if (y_hi <= y_lo) return

      half_dy = 0.5d0*(y_hi-y_lo)
      y_mid = 0.5d0*(y_hi+y_lo)
      do iq = 1, 3
         y_eval = y_mid + half_dy*xi(iq)
         gamma = gamma_from_coord(coord_fourvel, coord_scale, y_eval)
         density = spec_gamma(gamma, is_thermal)
         jac = gamma*dxg_dcoord(coord_fourvel, coord_scale, y_eval)
         acc = acc + half_dy*wi(iq)*density*jac
      end do
   end subroutine add_segment

   ! 单点物理谱密度 dN/dgamma：低能 thermal，高能 cutoff power law。
   ! Point physical density dN/dgamma: thermal branch at low energy, cutoff power law at high energy.
   real(8) function spec_gamma(g, is_thermal)
      implicit none
      logical, intent(in) :: is_thermal
      real(8), intent(in) :: g
      real(8) :: ln_val

      if (is_thermal) then
         ln_val = log(g * sqrt(g*g - 1.0d0)) - g*inv_theta + thermal_constant
      else
         ln_val = -p * log(g) - g*inv_gmax + cpl_constant
      end if

      if (ln_val <= m700) then
         spec_gamma = 0.0d0
         return
      end if

      spec_gamma = exp(ln_val)
   end function spec_gamma
   end subroutine hybrid_coord

   ! 在 four-velocity 单元上积分 Maxwell-Juttner thermal 分支，归一化为单位热电子数。
   ! Integrate the Maxwell-Juttner thermal branch on four-velocity cells, normalized to unit thermal count.
   subroutine hybrid_thermal_coord(n_gamma, coord_edge, coord_scale, p, gamma_min, gamma_max, xi_e, &
      spec)
      implicit none
      integer, intent(in) :: n_gamma
      real(8), intent(in), dimension(n_gamma+1) :: coord_edge
      real(8), intent(in) :: coord_scale, p, gamma_min, gamma_max, xi_e
      real(8), intent(out), dimension(n_gamma) :: spec

      real(8) :: theta, it1, thermal_constant, inv_theta, y_min, cell_lo, cell_hi, dy_cell, seg_sum
      integer :: i

      call solve_theta(p, gamma_min, gamma_max, xi_e, theta)
      call thermal_int1(gamma_min, theta, it1)
      thermal_constant = -log(it1)
      inv_theta = 1.0d0/theta
      y_min = coord_from_xg(coord_fourvel, coord_scale, dlog(gamma_min))

      do i = 1, n_gamma
         cell_lo = coord_edge(i)
         cell_hi = min(coord_edge(i+1), y_min)
         dy_cell = coord_edge(i+1) - coord_edge(i)
         seg_sum = 0.0d0
         if (cell_hi > cell_lo) call add_segment(cell_lo, cell_hi, seg_sum)
         spec(i) = seg_sum/dy_cell
      end do

   contains

      subroutine add_segment(y_lo, y_hi, acc)
         implicit none
         real(8), intent(in) :: y_lo, y_hi
         real(8), intent(inout) :: acc
         integer :: iq
         real(8) :: half_dy, y_mid, y_eval, gamma, lnval, density, jac
         real(8), parameter, dimension(3) :: xi=[-dsqrt(3d0/5d0),0d0,dsqrt(3d0/5d0)], wi=[5d0/9d0,8d0/9d0,5d0/9d0]

         half_dy = 0.5d0*(y_hi-y_lo)
         y_mid = 0.5d0*(y_hi+y_lo)
         do iq = 1, 3
            y_eval = y_mid + half_dy*xi(iq)
            gamma = gamma_from_coord(coord_fourvel, coord_scale, y_eval)
            lnval = log(gamma * sqrt(gamma*gamma - 1.0d0)) - gamma*inv_theta + thermal_constant
            if (lnval > m700) then
               density = exp(lnval)
               jac = gamma*dxg_dcoord(coord_fourvel, coord_scale, y_eval)
               acc = acc + half_dy*wi(iq)*density*jac
            end if
         end do
      end subroutine add_segment
   end subroutine hybrid_thermal_coord

end module hybrid_spectrum
