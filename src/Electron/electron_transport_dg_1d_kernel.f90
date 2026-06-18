!f2py: skip
module electron_transport_dg_1d_kernel
    use constants
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use electron_injection_profiles, only: electron_dnx_powerlaw_cutoff_value, electron_initial_powerlaw_params
    implicit none
    private

    integer, parameter :: electron_dg1d_poly_order = 12
    integer, parameter :: electron_dg1d_nodes_per_domain = electron_dg1d_poly_order + 1
    integer, parameter :: electron_dg1d_elements_per_segment = 6
    integer, parameter :: electron_dg1d_max_breaks = 3
    integer, parameter :: electron_dg1d_source_quad_order = 4
    real(8), parameter :: electron_dg1d_time_theta = 1.0d0
    real(8), parameter :: electron_dg1d_source_quad_nodes(electron_dg1d_source_quad_order) = &
        (/2.3260475760753783d-5, 1.1852538299698612d-2, 2.0155705894317942d-1, 7.4986401224114997d-1/)
    real(8), parameter :: electron_dg1d_source_quad_weights(electron_dg1d_source_quad_order) = &
        (/5.8639022562157604d-4, 1.1796387950693022d-1, 4.3814275798321245d-1, 4.4330715288489519d-1/)
    real(8), allocatable :: reference_r(:), reference_w(:), reference_bary(:), reference_dmat(:,:)
    real(8), allocatable :: projection_r(:), projection_w(:)

    type, public :: electron_dg1d_mesh
        integer :: ndom = 0
        integer :: nnode = electron_dg1d_nodes_per_domain
        integer :: ntot = 0
        real(8), allocatable :: x_left(:), x_right(:), r(:), w(:), dmat(:,:), bary(:), x(:)
    end type electron_dg1d_mesh

    public :: electron_dg1d_build_mesh, electron_dg1d_initial_state, electron_dg1d_project_state
    public :: electron_dg1d_gamma_nodes, electron_dg1d_source_nodes, electron_dg1d_advance_step
    public :: electron_dg1d_advance_characteristic_step
    public :: electron_dg1d_project_to_grid, electron_dg1d_integral

contains

subroutine electron_dg1d_build_mesh(x_min, x_max, x_break_a, x_break_b, x_break_c, mesh)
    real(8), intent(in) :: x_min, x_max, x_break_a, x_break_b, x_break_c
    type(electron_dg1d_mesh), intent(out) :: mesh
    real(8) :: breaks(electron_dg1d_max_breaks), active(electron_dg1d_max_breaks), min_width
    integer :: n_active

    if (x_max <= x_min) error stop 'electron_dg1d_build_mesh requires x_max > x_min'
    min_width = 1d-6*(x_max - x_min)
    breaks = (/x_break_a, x_break_b, x_break_c/)
    call electron_dg1d_sort_breaks(breaks)
    n_active = 0
    call add_active_break(breaks(1), x_min, x_max, min_width, active, n_active)
    call add_active_break(breaks(2), x_min, x_max, min_width, active, n_active)
    call add_active_break(breaks(3), x_min, x_max, min_width, active, n_active)

    mesh%ndom = (n_active + 1)*electron_dg1d_elements_per_segment
    mesh%nnode = electron_dg1d_nodes_per_domain
    mesh%ntot = mesh%ndom*mesh%nnode
    call allocate_spectral_mesh(mesh)
    call set_domain_bounds(x_min, x_max, active, n_active, mesh)
    call fill_physical_nodes(mesh)
end subroutine electron_dg1d_build_mesh

subroutine electron_dg1d_initial_state(mesh, total_norm, p, gamma_m, gamma_c, gamma_max, state)
    type(electron_dg1d_mesh), intent(in) :: mesh
    real(8), intent(in) :: total_norm, p, gamma_m, gamma_c, gamma_max
    real(8), intent(out) :: state(mesh%ntot)
    logical :: active
    real(8) :: coeff, slope, gamma
    integer :: i

    do i = 1, mesh%ntot
        gamma = ten**mesh%x(i)
        call electron_initial_powerlaw_params(total_norm, p, gamma_m, gamma_c, gamma, active, coeff, slope)
        if (active) then
            state(i) = electron_dnx_powerlaw_cutoff_value(mesh%x(i), coeff, slope, gamma_max)
        else
            state(i) = zero
        endif
    enddo
end subroutine electron_dg1d_initial_state

subroutine electron_dg1d_project_state(old_mesh, old_state, new_mesh, new_state)
    type(electron_dg1d_mesh), intent(in) :: old_mesh, new_mesh
    real(8), intent(in) :: old_state(old_mesh%ntot)
    real(8), intent(out) :: new_state(new_mesh%ntot)
    real(8) :: old_total, new_total
    integer :: k, offset

    do k = 1, new_mesh%ndom
        offset = domain_offset(new_mesh, k)
        call electron_dg1d_project_element(old_mesh, old_state, new_mesh, k, &
                                           new_state(offset + 1:offset + new_mesh%nnode))
    enddo
    call electron_dg1d_integral(old_mesh, old_state, old_total)
    call electron_dg1d_integral(new_mesh, new_state, new_total)
    if (old_total > zero) then
        if (new_total <= zero) error stop 'electron_dg1d_project_state lost positive content'
        new_state = new_state*(old_total/new_total)
    endif
end subroutine electron_dg1d_project_state

subroutine electron_dg1d_gamma_nodes(mesh, gamma_nodes)
    type(electron_dg1d_mesh), intent(in) :: mesh
    real(8), intent(out) :: gamma_nodes(mesh%ntot)

    gamma_nodes = ten**mesh%x
end subroutine electron_dg1d_gamma_nodes

subroutine electron_dg1d_source_nodes(mesh, source_norm, p, gamma_m, gamma_max, source)
    type(electron_dg1d_mesh), intent(in) :: mesh
    real(8), intent(in) :: source_norm, p, gamma_m, gamma_max
    real(8), intent(out) :: source(mesh%ntot)
    real(8) :: gamma
    integer :: i

    do i = 1, mesh%ntot
        gamma = ten**mesh%x(i)
        if (gamma > gamma_m) then
            source(i) = electron_dnx_powerlaw_cutoff_value(mesh%x(i), source_norm, p, gamma_max)
        else
            source(i) = zero
        endif
    enddo
end subroutine electron_dg1d_source_nodes

subroutine electron_dg1d_advance_step(mesh, radius, dr, dloggamma_loss, source, state_in, state_out)
    type(electron_dg1d_mesh), intent(in) :: mesh
    real(8), intent(in) :: radius, dr, dloggamma_loss(mesh%ntot), source(mesh%ntot), state_in(mesh%ntot)
    real(8), intent(out) :: state_out(mesh%ntot)
    real(8) :: block(mesh%nnode,mesh%nnode), local_rhs(mesh%nnode), speed(mesh%ntot)
    real(8) :: dx, coupling, rhs_scale, theta
    integer :: i, j, k, offset, next_offset

    theta = electron_dg1d_time_theta
    rhs_scale = one - theta
    speed = -(dloggamma_loss/dlog(ten) + one/(radius*dlog(ten)))
    state_out = zero
    do k = mesh%ndom, 1, -1
        offset = domain_offset(mesh, k)
        dx = mesh%x_right(k) - mesh%x_left(k)
        block = zero
        local_rhs = state_in(offset + 1:offset + mesh%nnode) + dr*source(offset + 1:offset + mesh%nnode)
        do i = 1, mesh%nnode
            do j = 1, mesh%nnode
                coupling = (two/dx)*(mesh%w(j)/mesh%w(i))*mesh%dmat(j,i)*speed(offset + j)
                block(i,j) = -theta*dr*coupling
                local_rhs(i) = local_rhs(i) + rhs_scale*dr*coupling*state_in(offset + j)
            enddo
            block(i,i) = block(i,i) + one
        enddo
        coupling = (two/dx)*speed(offset + 1)/mesh%w(1)
        block(1,1) = block(1,1) - theta*dr*coupling
        local_rhs(1) = local_rhs(1) + rhs_scale*dr*coupling*state_in(offset + 1)
        if (k < mesh%ndom) then
            next_offset = domain_offset(mesh, k + 1)
            coupling = (two/dx)*speed(offset + mesh%nnode)/mesh%w(mesh%nnode)
            local_rhs(mesh%nnode) = local_rhs(mesh%nnode) - theta*dr*coupling*state_out(next_offset + 1) &
                                    - rhs_scale*dr*coupling*state_in(next_offset + 1)
        endif
        call electron_dg1d_solve_dense(mesh%nnode, block, local_rhs)
        state_out(offset + 1:offset + mesh%nnode) = local_rhs
    enddo
end subroutine electron_dg1d_advance_step

subroutine electron_dg1d_advance_characteristic_step(mesh, dr, a_rad, b_ad, source, state_in, state_out)
    type(electron_dg1d_mesh), intent(in) :: mesh
    real(8), intent(in) :: dr, a_rad, b_ad, source(mesh%ntot), state_in(mesh%ntot)
    real(8), intent(out) :: state_out(mesh%ntot)
    integer :: q

    state_out = zero
    call electron_dg1d_project_characteristic(mesh, state_in, a_rad, b_ad, dr, one, state_out)
    do q = 1, electron_dg1d_source_quad_order
        call electron_dg1d_project_characteristic(mesh, source, a_rad, b_ad, &
                                                  electron_dg1d_source_quad_nodes(q)*dr, &
                                                  electron_dg1d_source_quad_weights(q)*dr, state_out)
    enddo
    call electron_dg1d_zero_negative_cell_means(mesh, state_out)
end subroutine electron_dg1d_advance_characteristic_step

subroutine electron_dg1d_zero_negative_cell_means(mesh, state)
    type(electron_dg1d_mesh), intent(in) :: mesh
    real(8), intent(inout) :: state(mesh%ntot)
    real(8) :: cell_sum
    integer :: k, offset

    do k = 1, mesh%ndom
        offset = domain_offset(mesh, k)
        cell_sum = sum(mesh%w*state(offset + 1:offset + mesh%nnode))
        if (.not. (cell_sum >= zero .and. ieee_is_finite(cell_sum))) then
            state(offset + 1:offset + mesh%nnode) = zero
        endif
    enddo
end subroutine electron_dg1d_zero_negative_cell_means

subroutine electron_dg1d_project_characteristic(mesh, field, a_rad, b_ad, lag, scale, projected)
    type(electron_dg1d_mesh), intent(in) :: mesh
    real(8), intent(in) :: field(mesh%ntot), a_rad, b_ad, lag, scale
    real(8), intent(inout) :: projected(mesh%ntot)
    real(8) :: modal(mesh%nnode), pvals(mesh%nnode)
    real(8) :: dx, lo, hi, mid, half_width, x_old, x_new, r_new, old_value
    real(8) :: x_back_left, x_back_right, old_min, old_max
    integer :: degree, k, k_old, q, m, i, offset

    degree = mesh%nnode - 1
    call ensure_projection_quadrature(mesh%nnode)
    old_min = mesh%x_left(1)
    old_max = mesh%x_right(mesh%ndom)
    do k = 1, mesh%ndom
        offset = domain_offset(mesh, k)
        dx = mesh%x_right(k) - mesh%x_left(k)
        modal = zero
        x_back_left = electron_dg1d_characteristic_back_x(mesh, a_rad, b_ad, lag, mesh%x_left(k))
        x_back_right = electron_dg1d_characteristic_back_x(mesh, a_rad, b_ad, lag, mesh%x_right(k))
        ! The finite DG grid has no represented state/source above old_max; clipped
        ! preimages intentionally discard exponentially small unresolved tails.
        x_back_left = max(old_min, min(old_max, x_back_left))
        x_back_right = max(old_min, min(old_max, x_back_right))
        if (x_back_right > x_back_left) then
            do k_old = 1, mesh%ndom
                if (mesh%x_right(k_old) <= x_back_left) cycle
                if (mesh%x_left(k_old) >= x_back_right) exit
                lo = max(x_back_left, mesh%x_left(k_old))
                hi = min(x_back_right, mesh%x_right(k_old))
                if (hi <= lo) cycle
                mid = 0.5d0*(lo + hi)
                half_width = 0.5d0*(hi - lo)
                do q = 1, mesh%nnode
                    x_old = mid + half_width*projection_r(q)
                    x_new = electron_dg1d_characteristic_forward_x(a_rad, b_ad, lag, x_old)
                    r_new = two*(x_new - mesh%x_left(k))/dx - one
                    old_value = interpolate_domain(mesh, field, k_old, x_old)
                    call legendre_basis_values(degree, r_new, pvals)
                    modal = modal + scale*half_width*projection_w(q)*old_value*pvals
                enddo
            enddo
        endif
        do m = 0, degree
            modal(m + 1) = dble(2*m + 1)*modal(m + 1)/dx
        enddo
        do i = 1, mesh%nnode
            call legendre_basis_values(degree, mesh%r(i), pvals)
            projected(offset + i) = projected(offset + i) + sum(modal*pvals)
        enddo
    enddo
end subroutine electron_dg1d_project_characteristic

real(8) function electron_dg1d_characteristic_back_x(mesh, a_rad, b_ad, lag, x_new) result(x_back)
    type(electron_dg1d_mesh), intent(in) :: mesh
    real(8), intent(in) :: a_rad, b_ad, lag, x_new
    real(8) :: gamma_now, gamma_back, exp_b, denom

    if (lag <= zero) then
        x_back = x_new
        return
    endif
    gamma_now = ten**x_new
    if (abs(b_ad*lag) <= 1d-10) then
        denom = one - a_rad*gamma_now*lag
    else
        exp_b = dexp(-b_ad*lag)
        denom = exp_b - (a_rad/b_ad)*(one - exp_b)*gamma_now
    endif
    if (denom <= zero) then
        x_back = mesh%x_right(mesh%ndom) + one
        return
    endif
    gamma_back = gamma_now/denom
    if (.not. ieee_is_finite(gamma_back)) then
        x_back = mesh%x_right(mesh%ndom) + one
    else
        x_back = dlog10(gamma_back)
    endif
end function electron_dg1d_characteristic_back_x

real(8) function electron_dg1d_characteristic_forward_x(a_rad, b_ad, lag, x_old) result(x_new)
    real(8), intent(in) :: a_rad, b_ad, lag, x_old
    real(8) :: gamma_old, gamma_new, exp_b

    if (lag <= zero) then
        x_new = x_old
        return
    endif
    gamma_old = ten**x_old
    if (abs(b_ad*lag) <= 1d-10) then
        gamma_new = gamma_old/(one + a_rad*gamma_old*lag)
    else
        exp_b = dexp(-b_ad*lag)
        gamma_new = gamma_old*exp_b/(one + (a_rad/b_ad)*gamma_old*(one - exp_b))
    endif
    x_new = dlog10(gamma_new)
end function electron_dg1d_characteristic_forward_x

subroutine electron_dg1d_project_element(old_mesh, old_state, new_mesh, k_new, values)
    type(electron_dg1d_mesh), intent(in) :: old_mesh, new_mesh
    integer, intent(in) :: k_new
    real(8), intent(in) :: old_state(old_mesh%ntot)
    real(8), intent(out) :: values(new_mesh%nnode)
    real(8) :: modal(new_mesh%nnode), pvals(new_mesh%nnode)
    real(8) :: dx_new, lo, hi, mid, half_width, x_eval, r_new, old_value
    integer :: degree, k_old, q, m, i

    degree = new_mesh%nnode - 1
    dx_new = new_mesh%x_right(k_new) - new_mesh%x_left(k_new)
    modal = zero
    call ensure_projection_quadrature(new_mesh%nnode)
    do k_old = 1, old_mesh%ndom
        if (old_mesh%x_right(k_old) <= new_mesh%x_left(k_new)) cycle
        if (old_mesh%x_left(k_old) >= new_mesh%x_right(k_new)) exit
        lo = max(new_mesh%x_left(k_new), old_mesh%x_left(k_old))
        hi = min(new_mesh%x_right(k_new), old_mesh%x_right(k_old))
        if (hi <= lo) cycle
        mid = 0.5d0*(lo + hi)
        half_width = 0.5d0*(hi - lo)
        do q = 1, new_mesh%nnode
            x_eval = mid + half_width*projection_r(q)
            old_value = interpolate_domain(old_mesh, old_state, k_old, x_eval)
            r_new = two*(x_eval - new_mesh%x_left(k_new))/dx_new - one
            call legendre_basis_values(degree, r_new, pvals)
            modal = modal + half_width*projection_w(q)*old_value*pvals
        enddo
    enddo
    do m = 0, degree
        modal(m + 1) = dble(2*m + 1)*modal(m + 1)/dx_new
    enddo
    do i = 1, new_mesh%nnode
        call legendre_basis_values(degree, new_mesh%r(i), pvals)
        values(i) = sum(modal*pvals)
    enddo
end subroutine electron_dg1d_project_element

subroutine electron_dg1d_solve_dense(n, a, b)
    integer, intent(in) :: n
    real(8), intent(inout) :: a(n,n), b(n)
    real(8) :: row_tmp(n), factor, tmp
    integer :: i, j, k, pivot

    do k = 1, n - 1
        pivot = k
        do i = k + 1, n
            if (abs(a(i,k)) > abs(a(pivot,k))) pivot = i
        enddo
        if (.not. (abs(a(pivot,k)) > zero)) error stop 'electron_dg1d_solve_dense singular matrix'
        if (pivot /= k) then
            row_tmp = a(k,:)
            a(k,:) = a(pivot,:)
            a(pivot,:) = row_tmp
            tmp = b(k)
            b(k) = b(pivot)
            b(pivot) = tmp
        endif
        do i = k + 1, n
            factor = a(i,k)/a(k,k)
            a(i,k) = zero
            do j = k + 1, n
                a(i,j) = a(i,j) - factor*a(k,j)
            enddo
            b(i) = b(i) - factor*b(k)
        enddo
    enddo
    if (.not. (abs(a(n,n)) > zero)) error stop 'electron_dg1d_solve_dense singular final pivot'
    do i = n, 1, -1
        tmp = b(i)
        do j = i + 1, n
            tmp = tmp - a(i,j)*b(j)
        enddo
        b(i) = tmp/a(i,i)
    enddo
end subroutine electron_dg1d_solve_dense

subroutine electron_dg1d_project_to_grid(mesh, state, num_gamma, gamma_grid, dndgamma)
    type(electron_dg1d_mesh), intent(in) :: mesh
    integer, intent(in) :: num_gamma
    real(8), intent(in) :: state(mesh%ntot), gamma_grid(num_gamma)
    real(8), intent(out) :: dndgamma(num_gamma)
    real(8) :: x_eval, dnx_eval
    integer :: i, k

    do i = 1, num_gamma
        x_eval = dlog10(gamma_grid(i))
        k = locate_domain(mesh, x_eval)
        dnx_eval = interpolate_domain(mesh, state, k, x_eval)
        dndgamma(i) = dnx_eval/(gamma_grid(i)*dlog(ten))
    enddo
end subroutine electron_dg1d_project_to_grid

subroutine electron_dg1d_integral(mesh, state, total)
    type(electron_dg1d_mesh), intent(in) :: mesh
    real(8), intent(in) :: state(mesh%ntot)
    real(8), intent(out) :: total
    integer :: k, offset
    real(8) :: jac

    total = zero
    do k = 1, mesh%ndom
        offset = domain_offset(mesh, k)
        jac = 0.5d0*(mesh%x_right(k) - mesh%x_left(k))
        total = total + jac*sum(mesh%w*state(offset + 1:offset + mesh%nnode))
    enddo
end subroutine electron_dg1d_integral

subroutine add_active_break(x_break, x_min, x_max, min_width, active, n_active)
    real(8), intent(in) :: x_break, x_min, x_max, min_width
    real(8), intent(inout) :: active(electron_dg1d_max_breaks)
    integer, intent(inout) :: n_active

    if (x_break <= x_min + min_width .or. x_break >= x_max - min_width) return
    if (n_active >= 1) then
        if (minval(abs(x_break - active(1:n_active))) <= min_width) return
    endif
    n_active = n_active + 1
    active(n_active) = x_break
end subroutine add_active_break

subroutine electron_dg1d_sort_breaks(breaks)
    real(8), intent(inout) :: breaks(electron_dg1d_max_breaks)
    real(8) :: tmp
    integer :: i, j

    do i = 1, electron_dg1d_max_breaks - 1
        do j = i + 1, electron_dg1d_max_breaks
            if (breaks(j) < breaks(i)) then
                tmp = breaks(i)
                breaks(i) = breaks(j)
                breaks(j) = tmp
            endif
        enddo
    enddo
end subroutine electron_dg1d_sort_breaks

subroutine allocate_spectral_mesh(mesh)
    type(electron_dg1d_mesh), intent(inout) :: mesh

    allocate(mesh%x_left(mesh%ndom), mesh%x_right(mesh%ndom), &
             mesh%r(mesh%nnode), mesh%w(mesh%nnode), mesh%dmat(mesh%nnode,mesh%nnode), &
             mesh%bary(mesh%nnode), mesh%x(mesh%ntot))
    call ensure_reference_spectral(mesh%nnode)
    mesh%r = reference_r
    mesh%w = reference_w
    mesh%bary = reference_bary
    mesh%dmat = reference_dmat
end subroutine allocate_spectral_mesh

subroutine ensure_reference_spectral(nnode)
    integer, intent(in) :: nnode

    if (allocated(reference_r)) return
    allocate(reference_r(nnode), reference_w(nnode), reference_bary(nnode), reference_dmat(nnode,nnode))
    call lgl_nodes_weights(nnode, reference_r, reference_w)
    call barycentric_weights(nnode, reference_r, reference_bary)
    call differentiation_matrix(nnode, reference_r, reference_bary, reference_dmat)
end subroutine ensure_reference_spectral

subroutine ensure_projection_quadrature(nnode)
    integer, intent(in) :: nnode

    if (allocated(projection_r)) return
    allocate(projection_r(nnode), projection_w(nnode))
    call gauss_legendre_nodes_weights(nnode, projection_r, projection_w)
end subroutine ensure_projection_quadrature

subroutine set_domain_bounds(x_min, x_max, active, n_active, mesh)
    real(8), intent(in) :: x_min, x_max, active(electron_dg1d_max_breaks)
    integer, intent(in) :: n_active
    type(electron_dg1d_mesh), intent(inout) :: mesh
    real(8) :: segment_left(electron_dg1d_max_breaks + 1), segment_right(electron_dg1d_max_breaks + 1), dx
    integer :: element_index, i, iseg, nseg

    nseg = n_active + 1
    segment_left = zero
    segment_right = zero
    segment_left(1) = x_min
    segment_right(nseg) = x_max
    if (n_active > 0) then
        segment_left(2:nseg) = active(1:n_active)
        segment_right(1:n_active) = active(1:n_active)
    endif

    element_index = 0
    do iseg = 1, nseg
        dx = (segment_right(iseg) - segment_left(iseg))/dble(electron_dg1d_elements_per_segment)
        do i = 1, electron_dg1d_elements_per_segment
            element_index = element_index + 1
            mesh%x_left(element_index) = segment_left(iseg) + dble(i - 1)*dx
            mesh%x_right(element_index) = segment_left(iseg) + dble(i)*dx
        enddo
    enddo
end subroutine set_domain_bounds

subroutine fill_physical_nodes(mesh)
    type(electron_dg1d_mesh), intent(inout) :: mesh
    real(8) :: dx
    integer :: k, i, idx

    do k = 1, mesh%ndom
        dx = mesh%x_right(k) - mesh%x_left(k)
        do i = 1, mesh%nnode
            idx = domain_offset(mesh, k) + i
            mesh%x(idx) = mesh%x_left(k) + 0.5d0*(mesh%r(i) + one)*dx
        enddo
    enddo
end subroutine fill_physical_nodes

subroutine lgl_nodes_weights(nnode, r, w)
    integer, intent(in) :: nnode
    real(8), intent(out) :: r(nnode), w(nnode)
    integer :: degree, i, iter
    real(8) :: x, pn, dpn, ddpn

    degree = nnode - 1
    if (degree < 1) error stop 'lgl_nodes_weights requires nnode >= 2'
    r(1) = -one
    r(nnode) = one
    w(1) = two/(dble(degree)*dble(degree + 1))
    w(nnode) = w(1)
    do i = 2, nnode - 1
        x = -dcos(pi*dble(i - 1)/dble(degree))
        do iter = 1, 50
            call legendre_value_derivative(degree, x, pn, dpn)
            ddpn = (two*x*dpn - dble(degree)*dble(degree + 1)*pn)/(one - x*x)
            x = x - dpn/ddpn
            if (abs(dpn) <= 1d-14) exit
        enddo
        r(i) = x
        call legendre_value_derivative(degree, r(i), pn, dpn)
        w(i) = two/(dble(degree)*dble(degree + 1)*pn*pn)
    enddo
end subroutine lgl_nodes_weights

subroutine legendre_value_derivative(degree, x, pn, dpn)
    integer, intent(in) :: degree
    real(8), intent(in) :: x
    real(8), intent(out) :: pn, dpn
    real(8) :: pnm2, pnm1, pnow
    integer :: n

    if (degree == 0) then
        pn = one
        dpn = zero
        return
    endif
    pnm2 = one
    pnm1 = x
    if (degree == 1) then
        pn = pnm1
    else
        do n = 2, degree
            pnow = (dble(2*n - 1)*x*pnm1 - dble(n - 1)*pnm2)/dble(n)
            pnm2 = pnm1
            pnm1 = pnow
        enddo
        pn = pnm1
    endif
    dpn = dble(degree)*(pnm2 - x*pn)/(one - x*x)
end subroutine legendre_value_derivative

subroutine legendre_basis_values(degree, x, pvals)
    integer, intent(in) :: degree
    real(8), intent(in) :: x
    real(8), intent(out) :: pvals(degree + 1)
    integer :: n

    pvals(1) = one
    if (degree == 0) return
    pvals(2) = x
    do n = 2, degree
        pvals(n + 1) = (dble(2*n - 1)*x*pvals(n) - dble(n - 1)*pvals(n - 1))/dble(n)
    enddo
end subroutine legendre_basis_values

subroutine gauss_legendre_nodes_weights(nnode, r, w)
    integer, intent(in) :: nnode
    real(8), intent(out) :: r(nnode), w(nnode)
    real(8) :: x, pn, dpn
    integer :: i, iter

    do i = 1, nnode
        x = -dcos(pi*(dble(i) - 0.25d0)/(dble(nnode) + 0.5d0))
        do iter = 1, 50
            call legendre_value_derivative(nnode, x, pn, dpn)
            x = x - pn/dpn
            if (abs(pn) <= 1d-14) exit
        enddo
        r(i) = x
        call legendre_value_derivative(nnode, r(i), pn, dpn)
        w(i) = two/((one - r(i)*r(i))*dpn*dpn)
    enddo
end subroutine gauss_legendre_nodes_weights

subroutine barycentric_weights(nnode, nodes, bary)
    integer, intent(in) :: nnode
    real(8), intent(in) :: nodes(nnode)
    real(8), intent(out) :: bary(nnode)
    integer :: i, j

    bary = one
    do i = 1, nnode
        do j = 1, nnode
            if (j /= i) bary(i) = bary(i)*(nodes(i) - nodes(j))
        enddo
        bary(i) = one/bary(i)
    enddo
end subroutine barycentric_weights

subroutine differentiation_matrix(nnode, nodes, bary, dmat)
    integer, intent(in) :: nnode
    real(8), intent(in) :: nodes(nnode), bary(nnode)
    real(8), intent(out) :: dmat(nnode,nnode)
    integer :: i, j

    dmat = zero
    do i = 1, nnode
        do j = 1, nnode
            if (j /= i) dmat(i,j) = bary(j)/(bary(i)*(nodes(i) - nodes(j)))
        enddo
        dmat(i,i) = -sum(dmat(i,:))
    enddo
end subroutine differentiation_matrix

integer function domain_offset(mesh, k) result(offset)
    type(electron_dg1d_mesh), intent(in) :: mesh
    integer, intent(in) :: k

    offset = (k - 1)*mesh%nnode
end function domain_offset

integer function locate_domain(mesh, x_eval) result(k_found)
    type(electron_dg1d_mesh), intent(in) :: mesh
    real(8), intent(in) :: x_eval
    integer :: k

    if (x_eval <= mesh%x_left(1)) then
        k_found = 1
        return
    endif
    do k = 1, mesh%ndom
        if (x_eval <= mesh%x_right(k)) then
            k_found = k
            return
        endif
    enddo
    k_found = mesh%ndom
end function locate_domain

real(8) function interpolate_domain(mesh, state, k, x_eval) result(value)
    type(electron_dg1d_mesh), intent(in) :: mesh
    real(8), intent(in) :: state(mesh%ntot), x_eval
    integer, intent(in) :: k
    real(8) :: r_eval, numerator, denominator, diff
    integer :: i, idx, offset

    offset = domain_offset(mesh, k)
    r_eval = two*(x_eval - mesh%x_left(k))/(mesh%x_right(k) - mesh%x_left(k)) - one
    numerator = zero
    denominator = zero
    do i = 1, mesh%nnode
        idx = offset + i
        diff = r_eval - mesh%r(i)
        if (abs(diff) <= 1d-14) then
            value = state(idx)
            return
        endif
        numerator = numerator + mesh%bary(i)*state(idx)/diff
        denominator = denominator + mesh%bary(i)/diff
    enddo
    value = numerator/denominator
end function interpolate_domain

end module electron_transport_dg_1d_kernel
