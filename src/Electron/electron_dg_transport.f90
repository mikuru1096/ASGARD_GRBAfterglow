!f2py: skip
module electron_dg_transport
    use constants
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use electron_coord_common, only: coord_loggamma, coord_fourvel, &
                                                 coord_from_xg, xg_from_coord, &
                                                 gamma_from_coord, dxg_dcoord
    use electron_injection_profiles, only: dnx_cutoff, exp_cutoff, &
                                           pl_params
    implicit none
    private

    ! DG 谱元离散参数：LGL 节点、分段断点和源项投影积分阶数。
    ! DG spectral-element setup: LGL nodes, mesh breaks, and source quadrature order.
    integer, parameter :: dg_order = 12, dg_nodes = dg_order + 1, dg_segments = 6, dg_breaks = 3
    integer, parameter :: dg_quadn = 4, dg_loggamma = coord_loggamma, dg_fourvel = coord_fourvel
    real(8), parameter, dimension(dg_quadn) :: dg_qnodes= [2.3260475760753783d-5, &
        1.1852538299698612d-2, 2.0155705894317942d-1, 7.4986401224114997d-1]
    real(8), parameter, dimension(dg_quadn) :: dg_qweights= [5.8639022562157604d-4, &
        1.1796387950693022d-1, 4.3814275798321245d-1, 4.4330715288489519d-1]
    real(8), allocatable, dimension(:) :: reference_r,reference_w,reference_bary
    real(8), allocatable, dimension(:,:) :: reference_dmat,reference_transport,reference_basis
    real(8), allocatable, dimension(:) :: projection_r,projection_w
    real(8), allocatable, dimension(:,:) :: projection_basis
    !$omp threadprivate(reference_r,reference_w,reference_bary,reference_dmat, &
    !$omp& reference_transport,reference_basis,projection_r,projection_w,projection_basis)

    ! DG 网格记录：只保存坐标、节点、权重和导数矩阵，不承载推进逻辑。
    ! DG mesh record: stores coordinates, nodes, weights, and matrices without solver behavior.
    type, public :: dg_mesh
        integer :: ndom = 0, nnode = dg_nodes, ntot = 0, coord_kind = dg_loggamma
        real(8) :: coord_scale = 1d0
        real(8), allocatable, dimension(:) :: x_left,x_right,r,w,bary,x
        real(8), allocatable, dimension(:,:) :: dmat
        real(8), allocatable, dimension(:) :: x_gamma,gamma,dxgamma_dcoord,dln_dcoord
    end type dg_mesh

    public :: dg_initial_state, dg_project_state, dg_build_mesh, dg_project_source, dg_kinetic_source
    public :: dg_advance_step, dg_scale_content, dg_limit_positive, dg_char_step
    public :: dg_project_cells, dg_integral, dg_tail_fraction

contains

! 在 four-velocity 坐标上构造 DG 网格，用 gamma 断点标出注入和冷却特征位置。
! Build the DG mesh in four-velocity coordinates, with gamma breaks near injection/cooling features.
subroutine dg_build_mesh(x_min, x_max, break_a, break_b, break_c, gamma_scale, mesh)
    real(8), intent(in) :: x_min, x_max, break_a, break_b, break_c, gamma_scale
    type(dg_mesh), intent(out) :: mesh
    real(8), dimension(dg_breaks) :: breaks
    real(8) :: fourvel2

    if (gamma_scale <= 1d0) error stop 'dg_build_mesh requires gamma_scale > 1'
    fourvel2 = gamma_scale*gamma_scale - 1d0
    breaks = [break_a, break_b, break_c]
    call dg_build_coord(dg_fourvel, fourvel2, &
                                        x_min, x_max, breaks, mesh)
end subroutine dg_build_mesh

! 从给定坐标类型和断点生成分段 LGL 单元。
! Build segmented LGL elements for the selected coordinate and active break set.
subroutine dg_build_coord(coord_kind, coord_scale, x_min, x_max, x_breaks, mesh)
    integer, intent(in) :: coord_kind
    real(8), intent(in), dimension(dg_breaks) :: x_breaks
    real(8), intent(in) :: coord_scale,x_min,x_max
    type(dg_mesh), intent(out) :: mesh
    real(8), dimension(dg_breaks) :: breaks,active
    real(8) :: min_width, coord_min, coord_max
    integer :: n_active

    coord_min = coord_from_xg(coord_kind, coord_scale, x_min)
    coord_max = coord_from_xg(coord_kind, coord_scale, x_max)
    min_width = 1d-6*(coord_max - coord_min)
    breaks = x_breaks
    do n_active = 1, dg_breaks
        breaks(n_active) = coord_from_xg(coord_kind, coord_scale, breaks(n_active))
    enddo
    call sort_breaks(breaks)
    n_active = 0
    call add_active_break(breaks(1), coord_min, coord_max, min_width, active, n_active)
    call add_active_break(breaks(2), coord_min, coord_max, min_width, active, n_active)
    call add_active_break(breaks(3), coord_min, coord_max, min_width, active, n_active)

    mesh%ndom = (n_active + 1)*dg_segments
    mesh%nnode = dg_nodes
    mesh%ntot = mesh%ndom*mesh%nnode
    mesh%coord_kind = coord_kind
    mesh%coord_scale = coord_scale
    call allocate_spectral_mesh(mesh)
    call set_domain_bounds(coord_min, coord_max, active, n_active, mesh)
    call fill_physical_nodes(mesh)
end subroutine dg_build_coord

! 把初始 power-law cutoff 分布投影到 DG 节点。
! Project the initial power-law cutoff distribution onto DG nodes.
subroutine dg_initial_state(mesh, total_norm, p, gamma_m, gamma_c, gamma_max, state)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in) :: total_norm, p, gamma_m, gamma_c, gamma_max
    real(8), intent(out), dimension(mesh%ntot) :: state
    logical :: active
    real(8) :: coeff, slope, gamma
    integer :: i

    do i = 1, mesh%ntot
        gamma = mesh%gamma(i)
        call pl_params(total_norm, p, gamma_m, gamma_c, gamma, active, coeff, slope)
        if (active) then
            state(i) = dnx_cutoff(mesh%x_gamma(i), coeff, slope, gamma_max)* &
                       mesh%dxgamma_dcoord(i)
        else
            state(i) = 0d0
        endif
    enddo
end subroutine dg_initial_state

! 在网格重建后把旧 DG 状态守恒投影到新网格。
! Conservatively project a DG state after the mesh is rebuilt.
subroutine dg_project_state(old_mesh, old_state, new_mesh, new_state)
    type(dg_mesh), intent(in) :: old_mesh, new_mesh
    real(8), intent(in), dimension(old_mesh%ntot) :: old_state
    real(8), intent(out), dimension(new_mesh%ntot) :: new_state
    integer :: k, offset

    do k = 1, new_mesh%ndom
        offset = (k - 1)*new_mesh%nnode
        call dg_project_element(old_mesh, old_state, new_mesh, k, &
                                           new_state(offset + 1:offset + new_mesh%nnode))
    enddo
end subroutine dg_project_state

! 将热后非热电子注入源项投影到 DG 空间。
! Project the nonthermal electron injection source into the DG space.
subroutine dg_project_source(mesh, source_norm, p, gamma_m, gamma_max, source)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in) :: source_norm, p, gamma_m, gamma_max
    real(8), intent(out), dimension(mesh%ntot) :: source
    real(8), dimension(mesh%nnode) :: modal,pvals
    real(8) :: dx, mid, half_width, x_eval, gamma, value
    integer :: degree, k, q, m, i, offset

    source = 0d0
    if (source_norm <= 0d0) return
    degree = mesh%nnode - 1
    call ensure_projection_quadrature(mesh%nnode)
    do k = 1, mesh%ndom
        offset = (k - 1)*mesh%nnode
        dx = mesh%x_right(k) - mesh%x_left(k)
        mid = 0.5d0*(mesh%x_left(k) + mesh%x_right(k))
        half_width = 0.5d0*dx
        modal = 0d0
        do q = 1, mesh%nnode
            x_eval = mid + half_width*projection_r(q)
            gamma = gamma_from_coord(mesh%coord_kind, mesh%coord_scale, x_eval)
            value = 0d0
            if (gamma > gamma_m) then
                value = dnx_cutoff( &
                    xg_from_coord(mesh%coord_kind, mesh%coord_scale, x_eval), &
                    source_norm, p, gamma_max)* &
                    dxg_dcoord(mesh%coord_kind, mesh%coord_scale, x_eval)
            endif
            pvals = projection_basis(:,q)
            modal = modal + half_width*projection_w(q)*value*pvals
        enddo
        do m = 0, degree
            modal(m + 1) = dble(2*m + 1)*modal(m + 1)/dx
        enddo
        do i = 1, mesh%nnode
            pvals = reference_basis(:,i)
            source(offset + i) = sum(modal*pvals)
        enddo
    enddo
    call dg_limit_positive(mesh, source)
end subroutine dg_project_source

! 将 reverse-shock kinetic 源项投影到 DG 空间。
! Project the reverse-shock kinetic source into the DG space.
subroutine dg_kinetic_source(mesh, source_norm, p, gamma_m, gamma_max, source)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in) :: source_norm, p, gamma_m, gamma_max
    real(8), intent(out), dimension(mesh%ntot) :: source
    real(8), dimension(mesh%nnode) :: modal,pvals
    real(8) :: dx, mid, half_width, x_eval, gamma, value
    integer :: degree, k, q, m, i, offset

    source = 0d0
    if (source_norm <= 0d0) return
    degree = mesh%nnode - 1
    call ensure_projection_quadrature(mesh%nnode)
    do k = 1, mesh%ndom
        offset = (k - 1)*mesh%nnode
        dx = mesh%x_right(k) - mesh%x_left(k)
        mid = 0.5d0*(mesh%x_left(k) + mesh%x_right(k))
        half_width = 0.5d0*dx
        modal = 0d0
        do q = 1, mesh%nnode
            x_eval = mid + half_width*projection_r(q)
            gamma = gamma_from_coord(mesh%coord_kind, mesh%coord_scale, x_eval)
            value = 0d0
            if (gamma > gamma_m) value = source_norm*gamma*(gamma - 1d0)**(-p)* &
                                         exp_cutoff(gamma, gamma_max)* &
                                         dxg_dcoord(mesh%coord_kind, mesh%coord_scale, x_eval)
            pvals = projection_basis(:,q)
            modal = modal + half_width*projection_w(q)*value*pvals
        enddo
        do m = 0, degree
            modal(m + 1) = dble(2*m + 1)*modal(m + 1)/dx
        enddo
        do i = 1, mesh%nnode
            pvals = reference_basis(:,i)
            source(offset + i) = sum(modal*pvals)
        enddo
    enddo
    call dg_limit_positive(mesh, source)
end subroutine dg_kinetic_source

! 按目标总数缩放 DG 状态，供网格投影和源项模板复用。
! Rescale a DG state to the requested total content for mesh projection or source templates.
subroutine dg_scale_content(mesh, target_content, state)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in) :: target_content
    real(8), intent(inout), dimension(mesh%ntot) :: state
    real(8) :: dg_content

    if (target_content > 0d0) then
        call dg_integral(mesh, state, dg_content)
        if (dg_content <= 0d0) error stop 'dg_scale_content requires positive DG content'
        state = state*(target_content/dg_content)
    endif
end subroutine dg_scale_content

! 对每个单元执行保持平均值的正性限制。
! Apply a cell-average-preserving positivity limiter inside each element.
subroutine dg_limit_positive(mesh, state)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(inout), dimension(mesh%ntot) :: state
    real(8) :: cell_average, cell_min, theta, x_eval
    integer :: k, offset, q

    call ensure_projection_quadrature(mesh%nnode)
    do k = 1, mesh%ndom
        offset = (k - 1)*mesh%nnode
        cell_average = 0.5d0*sum(mesh%w*state(offset + 1:offset + mesh%nnode))
        cell_min = minval(state(offset + 1:offset + mesh%nnode))
        do q = 1, mesh%nnode
            x_eval = 0.5d0*((1d0 - projection_r(q))*mesh%x_left(k) + (1d0 + projection_r(q))*mesh%x_right(k))
            cell_min = min(cell_min, interpolate_domain(mesh, state, k, x_eval))
        enddo
        if (cell_min < 0d0) then
            if (cell_average > 0d0) then
                theta = min(1d0, cell_average/(cell_average - cell_min))
                state(offset + 1:offset + mesh%nnode) = cell_average + &
                    theta*(state(offset + 1:offset + mesh%nnode) - cell_average)
            else
                state(offset + 1:offset + mesh%nnode) = 0d0
            end if
        end if
    enddo
end subroutine dg_limit_positive

! 隐式 DG 输运步：冷却速度单向时用逐单元回代，符号混合时切到 dense solve。
! Implicit DG transport: use element back substitution for single-direction cooling, dense solve for mixed signs.
subroutine dg_advance_step(mesh, adiabatic_rate, dr, dloggamma_loss, source, state_in, state_out)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in), dimension(mesh%ntot) :: dloggamma_loss,source,state_in
    real(8), intent(in) :: adiabatic_rate,dr
    real(8), intent(out), dimension(mesh%ntot) :: state_out
    real(8), dimension(dg_nodes,dg_nodes) :: block
    real(8), dimension(dg_nodes) :: local_rhs
    real(8), dimension(mesh%ntot) :: speed
    real(8) :: dx, coupling, scale, column_coeff
    integer :: i, j, k, offset, next_offset

    speed = -(dloggamma_loss + adiabatic_rate)/mesh%dln_dcoord
    if (maxval(speed) > 0d0) then
        call dg_dense_step(mesh, dr, speed, source, state_in, state_out)
        return
    endif
    state_out = 0d0
    do k = mesh%ndom, 1, -1
        offset = (k - 1)*mesh%nnode
        dx = mesh%x_right(k) - mesh%x_left(k)
        scale = 2d0/dx
        block = 0d0
        local_rhs = state_in(offset + 1:offset + mesh%nnode) + dr*source(offset + 1:offset + mesh%nnode)
        do j = 1, mesh%nnode
            column_coeff = -dr*scale*speed(offset + j)
            do i = 1, mesh%nnode
                block(i,j) = column_coeff*reference_transport(i,j)
            enddo
        enddo
        do i = 1, mesh%nnode
            block(i,i) = block(i,i) + 1d0
        enddo
        if (k > 1 .or. mesh%coord_kind /= dg_fourvel) then
            coupling = scale*speed(offset + 1)/mesh%w(1)
            block(1,1) = block(1,1) - dr*coupling
        endif
        if (k < mesh%ndom) then
            next_offset = k*mesh%nnode
            coupling = scale*speed(offset + mesh%nnode)/mesh%w(mesh%nnode)
            local_rhs(mesh%nnode) = local_rhs(mesh%nnode) - dr*coupling*state_out(next_offset + 1)
        endif
        call dg_solve_block(block, local_rhs)
        state_out(offset + 1:offset + mesh%nnode) = local_rhs
    enddo
end subroutine dg_advance_step

! 组装全局 dense DG 输运矩阵，用于速度变号的冷却/加热情形。
! Assemble and solve the global dense DG transport system for sign-changing speeds.
subroutine dg_dense_step(mesh, dr, speed, source, state_in, state_out)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in), dimension(mesh%ntot) :: speed,source,state_in
    real(8), intent(in) :: dr
    real(8), intent(out), dimension(mesh%ntot) :: state_out
    real(8), allocatable, dimension(:,:) :: kmat,amat
    real(8), allocatable, dimension(:) :: rhs
    integer :: ntot, row

    ntot = mesh%ntot
    allocate(kmat(ntot,ntot), amat(ntot,ntot), rhs(ntot))
    kmat = 0d0
    call dg_transport_matrix(mesh, speed, kmat)
    amat = -dr*kmat
    do row = 1, ntot
        amat(row,row) = amat(row,row) + 1d0
    enddo
    rhs = state_in + dr*source
    call dg_solve_dense(ntot, amat, rhs)
    state_out = rhs
    call dg_limit_positive(mesh, state_out)
    deallocate(kmat, amat, rhs)
end subroutine dg_dense_step

! 组装 DG 通量矩阵，显式处理单元面上游耦合。
! Assemble the DG flux matrix with explicit upwind coupling across element faces.
subroutine dg_transport_matrix(mesh, speed, kmat)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in), dimension(mesh%ntot) :: speed
    real(8), intent(inout), dimension(mesh%ntot,mesh%ntot) :: kmat
    real(8) :: dx, face_speed, coeff
    integer :: i, j, k, offset, row, col, neighbor_offset

    do k = 1, mesh%ndom
        offset = (k - 1)*mesh%nnode
        dx = mesh%x_right(k) - mesh%x_left(k)
        do i = 1, mesh%nnode
            row = offset + i
            do j = 1, mesh%nnode
                col = offset + j
                kmat(row,col) = kmat(row,col) + (2d0/dx)*reference_transport(i,j)*speed(offset + j)
            enddo
        enddo

        if (k > 1) then
            row = offset + 1
            face_speed = speed(row)
            coeff = (2d0/dx)*face_speed/mesh%w(1)
            if (face_speed > 0d0) then
                neighbor_offset = (k - 2)*mesh%nnode
                col = neighbor_offset + mesh%nnode
            else
                col = row
            endif
            kmat(row,col) = kmat(row,col) + coeff
        endif

        if (k < mesh%ndom) then
            row = offset + mesh%nnode
            face_speed = speed(row)
            coeff = -(2d0/dx)*face_speed/mesh%w(mesh%nnode)
            if (face_speed > 0d0) then
                col = row
            else
                neighbor_offset = k*mesh%nnode
                col = neighbor_offset + 1
            endif
            kmat(row,col) = kmat(row,col) + coeff
        else
            row = offset + mesh%nnode
            face_speed = speed(row)
            if (face_speed > 0d0) then
                coeff = -(2d0/dx)*face_speed/mesh%w(mesh%nnode)
                kmat(row,row) = kmat(row,row) + coeff
            endif
        endif
    enddo
end subroutine dg_transport_matrix

! 特征线 DG 步：分别投影并限制旧谱响应与源响应，再线性叠加。
! DG characteristic step: project and limit homogeneous and source responses separately, then add them.
subroutine dg_char_step(mesh, dr, a_rad, b_ad, source, state_in, state_out)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in), dimension(mesh%ntot) :: source,state_in
    real(8), intent(in) :: dr,a_rad,b_ad
    real(8), intent(out), dimension(mesh%ntot) :: state_out
    real(8), dimension(mesh%ntot) :: source_response
    integer :: q

    state_out = 0d0
    call dg_project_char(mesh, state_in, a_rad, b_ad, dr, 1d0, state_out)
    call dg_zero_bad(mesh, state_out)
    call dg_limit_positive(mesh, state_out)
    source_response = 0d0
    do q = 1, dg_quadn
        call dg_project_char(mesh, source, a_rad, b_ad, &
                            dg_qnodes(q)*dr,dg_qweights(q)*dr, source_response)
    enddo
    call dg_zero_bad(mesh, source_response)
    call dg_limit_positive(mesh, source_response)
    state_out = state_out + source_response
end subroutine dg_char_step

! 清空非有限或负平均的单元，避免坏单元参与后续投影。
! Clear cells with non-finite or negative averages before later DG projection.
subroutine dg_zero_bad(mesh, state)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(inout), dimension(mesh%ntot) :: state
    real(8) :: cell_sum
    integer :: k, offset

    do k = 1, mesh%ndom
        offset = (k - 1)*mesh%nnode
        cell_sum = sum(mesh%w*state(offset + 1:offset + mesh%nnode))
        if (.not. (cell_sum >= 0d0 .and. ieee_is_finite(cell_sum))) then
            state(offset + 1:offset + mesh%nnode) = 0d0
        endif
    enddo
end subroutine dg_zero_bad

! 将场沿冷却特征线投影到当前 DG 网格。
! Project a field onto the current DG grid along cooling characteristics.
subroutine dg_project_char(mesh, field, a_rad, b_ad, lag, scale, projected)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in), dimension(mesh%ntot) :: field
    real(8), intent(in) :: a_rad,b_ad,lag,scale
    real(8), intent(inout), dimension(mesh%ntot) :: projected
    real(8), dimension(mesh%nnode) :: modal,pvals
    real(8) :: dx, lo, hi, mid, half_width, x_old, x_new, r_new, old_value
    real(8) :: back_l, back_r, old_min, old_max, low_content
    integer :: degree, k, k_old, q, m, i, offset

    degree = mesh%nnode - 1
    call ensure_projection_quadrature(mesh%nnode)
    old_min = mesh%x_left(1)
    old_max = mesh%x_right(mesh%ndom)
    low_content = 0d0
    call dg_low_content(mesh, field, a_rad, b_ad, lag, scale, &
                                                   old_min, old_max, low_content)
    do k = 1, mesh%ndom
        offset = (k - 1)*mesh%nnode
        dx = mesh%x_right(k) - mesh%x_left(k)
        modal = 0d0
        back_l = dg_back_x(mesh, a_rad, b_ad, lag, mesh%x_left(k))
        back_r = dg_back_x(mesh, a_rad, b_ad, lag, mesh%x_right(k))
        ! The finite DG grid has no represented state/source above old_max; clipped
        ! preimages intentionally discard exponentially small unresolved tails.
        back_l = max(old_min, min(old_max, back_l))
        back_r = max(old_min, min(old_max, back_r))
        if (back_r > back_l) then
            do k_old = 1, mesh%ndom
                if (mesh%x_right(k_old) <= back_l) cycle
                if (mesh%x_left(k_old) >= back_r) exit
                lo = max(back_l, mesh%x_left(k_old))
                hi = min(back_r, mesh%x_right(k_old))
                if (hi <= lo) cycle
                mid = 0.5d0*(lo + hi)
                half_width = 0.5d0*(hi - lo)
                do q = 1, mesh%nnode
                    x_old = mid + half_width*projection_r(q)
                    x_new = dg_forward_x(mesh, a_rad, b_ad, lag, x_old)
                    r_new = 2d0*(x_new - mesh%x_left(k))/dx - 1d0
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
    if (low_content /= 0d0) then
        dx = mesh%x_right(1) - mesh%x_left(1)
        projected(1:mesh%nnode) = projected(1:mesh%nnode) + low_content/dx
    endif
end subroutine dg_project_char

! four-velocity 坐标低能闭边界会收集越界回溯质量。
! Closed low boundary in four-velocity coordinates gathers traced-back mass outside the grid.
subroutine dg_low_content(mesh, field, a_rad, b_ad, lag, scale, &
                                                     old_min, old_max, content)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in), dimension(mesh%ntot) :: field
    real(8), intent(in) :: a_rad,b_ad,lag,scale,old_min,old_max
    real(8), intent(out) :: content
    real(8) :: x_cut, lo, hi, interval_content
    integer :: k_old

    content = 0d0
    if (mesh%coord_kind /= dg_fourvel .or. lag <= 0d0) return
    x_cut = dg_back_x(mesh, a_rad, b_ad, lag, old_min)
    if (x_cut <= old_min) return
    x_cut = min(x_cut, old_max)
    do k_old = 1, mesh%ndom
        if (mesh%x_left(k_old) >= x_cut) exit
        lo = max(old_min, mesh%x_left(k_old))
        hi = min(x_cut, mesh%x_right(k_old))
        if (hi <= lo) cycle
        call dg_interval_int(mesh, field, k_old, lo, hi, interval_content)
        content = content + scale*interval_content
    enddo
end subroutine dg_low_content

! 积分单个旧单元内的 DG 多项式质量。
! Integrate DG polynomial content over an interval inside an old element.
subroutine dg_interval_int(mesh, field, k, lo, hi, content)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in), dimension(mesh%ntot) :: field
    real(8), intent(in) :: lo,hi
    integer, intent(in) :: k
    real(8), intent(out) :: content
    real(8) :: mid, half_width, x_eval
    integer :: q

    mid = 0.5d0*(lo + hi)
    half_width = 0.5d0*(hi - lo)
    content = 0d0
    do q = 1, mesh%nnode
        x_eval = mid + half_width*projection_r(q)
        content = content + half_width*projection_w(q)*interpolate_domain(mesh, field, k, x_eval)
    enddo
end subroutine dg_interval_int

! 解析冷却轨道的反向映射 x_new -> x_back。
! Backward map along the analytic cooling characteristic, x_new to x_back.
real(8) function dg_back_x(mesh, a_rad, b_ad, lag, x_new) result(x_back)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in) :: a_rad, b_ad, lag, x_new
    real(8) :: gamma_now, gamma_back, exp_b, denom

    if (lag <= 0d0) then
        x_back = x_new
        return
    endif
    gamma_now = gamma_from_coord(mesh%coord_kind, mesh%coord_scale, x_new)
    if (abs(b_ad*lag) <= 1d-10) then
        denom = 1d0 - a_rad*gamma_now*lag
    else
        exp_b = dexp(-b_ad*lag)
        denom = exp_b - (a_rad/b_ad)*(1d0 - exp_b)*gamma_now
    endif
    if (denom <= 0d0) then
        x_back = mesh%x_right(mesh%ndom) + 1d0
        return
    endif
    gamma_back = gamma_now/denom
    if (.not. ieee_is_finite(gamma_back)) then
        x_back = mesh%x_right(mesh%ndom) + 1d0
    else
        x_back = coord_from_xg(mesh%coord_kind, mesh%coord_scale, dlog(gamma_back))
    endif
end function dg_back_x

! 解析冷却轨道的正向映射 x_old -> x_new。
! Forward map along the analytic cooling characteristic, x_old to x_new.
real(8) function dg_forward_x(mesh, a_rad, b_ad, lag, x_old) result(x_new)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in) :: a_rad, b_ad, lag, x_old
    real(8) :: gamma_old, gamma_new, exp_b

    if (lag <= 0d0) then
        x_new = x_old
        return
    endif
    gamma_old = gamma_from_coord(mesh%coord_kind, mesh%coord_scale, x_old)
    if (abs(b_ad*lag) <= 1d-10) then
        gamma_new = gamma_old/(1d0 + a_rad*gamma_old*lag)
    else
        exp_b = dexp(-b_ad*lag)
        gamma_new = gamma_old*exp_b/(1d0 + (a_rad/b_ad)*gamma_old*(1d0 - exp_b))
    endif
    x_new = coord_from_xg(mesh%coord_kind, mesh%coord_scale, dlog(gamma_new))
end function dg_forward_x

! 将旧网格单元贡献投影到新网格的目标单元。
! Project old-mesh element content into a target element on the new mesh.
subroutine dg_project_element(old_mesh, old_state, new_mesh, k_new, values)
    type(dg_mesh), intent(in) :: old_mesh, new_mesh
    integer, intent(in) :: k_new
    real(8), intent(in), dimension(old_mesh%ntot) :: old_state
    real(8), intent(out), dimension(new_mesh%nnode) :: values
    real(8), dimension(new_mesh%nnode) :: modal,pvals,emom
    real(8) :: dx_new, xgl, xgr, old_l, old_r, lo, hi, mid, half_width
    real(8) :: old_coord, new_coord, x_gamma, r_new, old_value, contribution
    real(8) :: target_n,target_e
    integer :: degree, k_old, q, m, i

    degree = new_mesh%nnode - 1
    dx_new = new_mesh%x_right(k_new) - new_mesh%x_left(k_new)
    xgl = xg_from_coord(new_mesh%coord_kind, new_mesh%coord_scale, new_mesh%x_left(k_new))
    xgr = xg_from_coord(new_mesh%coord_kind, new_mesh%coord_scale, new_mesh%x_right(k_new))
    old_l = coord_from_xg(old_mesh%coord_kind, old_mesh%coord_scale, xgl)
    old_r = coord_from_xg(old_mesh%coord_kind, old_mesh%coord_scale, xgr)
    modal = 0d0
    target_n = 0d0
    target_e = 0d0
    call ensure_projection_quadrature(new_mesh%nnode)
    ! Split every target cell at old DG faces and integrate in the old coordinate.
    ! The P0/P1 coefficients are constrained by number and kinetic-energy moments.
    ! 在旧 DG 单元面处分段，并以零阶粒子数和一阶动能矩约束 P0/P1 系数。
    do k_old = 1, old_mesh%ndom
        lo = max(old_l, old_mesh%x_left(k_old))
        hi = min(old_r, old_mesh%x_right(k_old))
        if (hi <= lo) cycle
        mid = 0.5d0*(lo + hi)
        half_width = 0.5d0*(hi - lo)
        do q = 1, old_mesh%nnode
            old_coord = mid + half_width*projection_r(q)
            x_gamma = xg_from_coord(old_mesh%coord_kind, old_mesh%coord_scale, old_coord)
            new_coord = coord_from_xg(new_mesh%coord_kind, new_mesh%coord_scale, x_gamma)
            r_new = 2d0*(new_coord - new_mesh%x_left(k_new))/dx_new - 1d0
            old_value = interpolate_domain(old_mesh, old_state, k_old, old_coord)
            call legendre_basis_values(degree, r_new, pvals)
            contribution = half_width*projection_w(q)*old_value
            modal = modal + contribution*pvals
            target_n = target_n + contribution
            target_e = target_e + contribution*(dexp(x_gamma) - 1d0)
        enddo
    enddo
    do m = 0, degree
        modal(m + 1) = dble(2*m + 1)*modal(m + 1)/dx_new
    enddo
    emom = 0d0
    mid = 0.5d0*(new_mesh%x_left(k_new) + new_mesh%x_right(k_new))
    half_width = 0.5d0*dx_new
    do q = 1, new_mesh%nnode
        new_coord = mid + half_width*projection_r(q)
        x_gamma = xg_from_coord(new_mesh%coord_kind, new_mesh%coord_scale, new_coord)
        call legendre_basis_values(degree, projection_r(q), pvals)
        emom = emom + half_width*projection_w(q)*(dexp(x_gamma) - 1d0)*pvals
    enddo
    modal(1) = target_n/dx_new
    modal(2) = (target_e - modal(1)*emom(1) - sum(modal(3:new_mesh%nnode)* &
                 emom(3:new_mesh%nnode)))/emom(2)
    do i = 1, new_mesh%nnode
        call legendre_basis_values(degree, new_mesh%r(i), pvals)
        values(i) = sum(modal*pvals)
    enddo
end subroutine dg_project_element

! 求解固定 LGL 块矩阵，供单向隐式输运回代使用。
! Solve a fixed-size LGL block system for the single-direction implicit transport path.
subroutine dg_solve_block(a, b)
    real(8), intent(inout), dimension(dg_nodes,dg_nodes) :: a
    real(8), intent(inout), dimension(dg_nodes) :: b
    real(8) :: factor, inv_pivot, tmp
    integer :: i, j, k, pivot

    do k = 1, dg_nodes - 1
        pivot = k
        do i = k + 1, dg_nodes
            if (abs(a(i,k)) > abs(a(pivot,k))) pivot = i
        enddo
        if (.not. (abs(a(pivot,k)) > 0d0)) error stop 'dg_solve_block singular matrix'
        if (pivot /= k) then
            do j = k, dg_nodes
                tmp = a(k,j)
                a(k,j) = a(pivot,j)
                a(pivot,j) = tmp
            enddo
            tmp = b(k)
            b(k) = b(pivot)
            b(pivot) = tmp
        endif
        inv_pivot = 1d0/a(k,k)
        do i = k + 1, dg_nodes
            a(i,k) = a(i,k)*inv_pivot
            b(i) = b(i) - a(i,k)*b(k)
        enddo
        do j = k + 1, dg_nodes
            factor = a(k,j)
            do i = k + 1, dg_nodes
                a(i,j) = a(i,j) - a(i,k)*factor
            enddo
        enddo
        do i = k + 1, dg_nodes
            a(i,k) = 0d0
        enddo
    enddo
    if (.not. (abs(a(dg_nodes,dg_nodes)) > 0d0)) &
        error stop 'dg_solve_block singular final pivot'
    do i = dg_nodes, 1, -1
        tmp = b(i)
        do j = i + 1, dg_nodes
            tmp = tmp - a(i,j)*b(j)
        enddo
        b(i) = tmp/a(i,i)
    enddo
end subroutine dg_solve_block

! 求解全局 dense 线性系统，供速度变号情形使用。
! Solve the global dense linear system for sign-changing transport.
subroutine dg_solve_dense(n, a, b)
    integer, intent(in) :: n
    real(8), intent(inout), dimension(n,n) :: a
    real(8), intent(inout), dimension(n) :: b
    real(8) :: factor, inv_pivot, tmp
    integer :: i, j, k, pivot

    do k = 1, n - 1
        pivot = k
        do i = k + 1, n
            if (abs(a(i,k)) > abs(a(pivot,k))) pivot = i
        enddo
        if (.not. (abs(a(pivot,k)) > 0d0)) error stop 'dg_solve_dense singular matrix'
        if (pivot /= k) then
            do j = k, n
                tmp = a(k,j)
                a(k,j) = a(pivot,j)
                a(pivot,j) = tmp
            enddo
            tmp = b(k)
            b(k) = b(pivot)
            b(pivot) = tmp
        endif
        inv_pivot = 1d0/a(k,k)
        do i = k + 1, n
            factor = a(i,k)*inv_pivot
            a(i,k) = 0d0
            do j = k + 1, n
                a(i,j) = a(i,j) - factor*a(k,j)
            enddo
            b(i) = b(i) - factor*b(k)
        enddo
    enddo
    if (.not. (abs(a(n,n)) > 0d0)) error stop 'dg_solve_dense singular final pivot'
    do i = n, 1, -1
        tmp = b(i)
        do j = i + 1, n
            tmp = tmp - a(i,j)*b(j)
        enddo
        b(i) = tmp/a(i,i)
    enddo
end subroutine dg_solve_dense

! 将 DG 状态投影回外部能量坐标单元平均量。
! Project the DG state back to external coordinate-cell averages.
subroutine dg_project_cells(mesh, state, num_gamma, coord_edge, dN_coord)
    type(dg_mesh), intent(in) :: mesh
    integer, intent(in) :: num_gamma
    real(8), intent(in), dimension(mesh%ntot) :: state
    real(8), intent(in), dimension(num_gamma+1) :: coord_edge
    real(8), intent(out), dimension(num_gamma) :: dN_coord
    real(8) :: lo, hi, mid, half_width, x_eval, cell_content
    integer :: i, k, k_start, q

    call ensure_projection_quadrature(mesh%nnode)
    k_start = 1
    do i = 1, num_gamma
        cell_content = 0d0
        do while (k_start < mesh%ndom .and. mesh%x_right(k_start) <= coord_edge(i))
            k_start = k_start + 1
        enddo
        k = k_start
        do while (k <= mesh%ndom)
            if (mesh%x_left(k) >= coord_edge(i+1)) exit
            lo = max(coord_edge(i), mesh%x_left(k))
            hi = min(coord_edge(i+1), mesh%x_right(k))
            if (hi <= lo) then
                k = k + 1
                cycle
            endif
            mid = 0.5d0*(lo + hi)
            half_width = 0.5d0*(hi - lo)
            do q = 1, mesh%nnode
                x_eval = mid + half_width*projection_r(q)
                cell_content = cell_content + &
                               half_width*projection_w(q)*interpolate_domain(mesh, state, k, x_eval)
            enddo
            k = k + 1
        enddo
        dN_coord(i) = cell_content/(coord_edge(i+1) - coord_edge(i))
    enddo
end subroutine dg_project_cells

! 计算 DG 状态总粒子数。
! Integrate total particle content in the DG state.
subroutine dg_integral(mesh, state, total)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in), dimension(mesh%ntot) :: state
    real(8), intent(out) :: total
    integer :: k, offset
    real(8) :: jac

    total = 0d0
    do k = 1, mesh%ndom
        offset = (k - 1)*mesh%nnode
        jac = 0.5d0*(mesh%x_right(k) - mesh%x_left(k))
        total = total + jac*sum(mesh%w*state(offset + 1:offset + mesh%nnode))
    enddo
end subroutine dg_integral

! 计算高能尾部矩，用于判定是否需要扩展活动网格。
! Compute a high-energy tail moment to decide whether the active mesh must expand.
subroutine dg_tail_fraction(mesh, state, gamma_cut, moment_power, fraction)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in), dimension(mesh%ntot) :: state
    real(8), intent(in) :: gamma_cut,moment_power
    real(8), intent(out) :: fraction
    real(8) :: total, tail, jac, contribution
    integer :: k, i, idx, offset

    total = 0d0
    tail = 0d0
    do k = 1, mesh%ndom
        offset = (k - 1)*mesh%nnode
        jac = 0.5d0*(mesh%x_right(k) - mesh%x_left(k))
        do i = 1, mesh%nnode
            idx = offset + i
            contribution = jac*mesh%w(i)*state(idx)*mesh%gamma(idx)**moment_power
            total = total + contribution
            if (mesh%gamma(idx) >= gamma_cut) tail = tail + contribution
        enddo
    enddo
    if (total > 0d0) then
        fraction = tail/total
    else
        fraction = 0d0
    endif
end subroutine dg_tail_fraction

! 将用户给出的断点加入活动断点集合。
! Add a user-provided break into the active break set.
subroutine add_active_break(x_break, x_min, x_max, min_width, active, n_active)
    real(8), intent(in) :: x_break, x_min, x_max, min_width
    real(8), intent(inout), dimension(dg_breaks) :: active
    integer, intent(inout) :: n_active

    if (x_break <= x_min + min_width .or. x_break >= x_max - min_width) return
    if (n_active >= 1) then
        if (minval(abs(x_break - active(1:n_active))) <= min_width) return
    endif
    n_active = n_active + 1
    active(n_active) = x_break
end subroutine add_active_break

! 对活动断点排序，保证分段网格从低到高生成。
! Sort active breaks so segmented elements are built from low to high coordinate.
subroutine sort_breaks(breaks)
    real(8), intent(inout), dimension(dg_breaks) :: breaks
    real(8) :: tmp
    integer :: i, j

    do i = 1, dg_breaks - 1
        do j = i + 1, dg_breaks
            if (breaks(j) < breaks(i)) then
                tmp = breaks(i)
                breaks(i) = breaks(j)
                breaks(j) = tmp
            endif
        enddo
    enddo
end subroutine sort_breaks

! 分配网格数组并装载参考 LGL 矩阵。
! Allocate mesh arrays and load the reference LGL matrices.
subroutine allocate_spectral_mesh(mesh)
    type(dg_mesh), intent(inout) :: mesh

    allocate(mesh%x_left(mesh%ndom), mesh%x_right(mesh%ndom), &
             mesh%r(mesh%nnode), mesh%w(mesh%nnode), mesh%dmat(mesh%nnode,mesh%nnode), &
             mesh%bary(mesh%nnode), mesh%x(mesh%ntot), mesh%x_gamma(mesh%ntot), mesh%gamma(mesh%ntot), &
             mesh%dxgamma_dcoord(mesh%ntot), mesh%dln_dcoord(mesh%ntot))
    call ensure_reference_spectral(mesh%nnode)
    mesh%r = reference_r
    mesh%w = reference_w
    mesh%bary = reference_bary
    mesh%dmat = reference_dmat
end subroutine allocate_spectral_mesh

! 初始化参考 LGL 节点、权重、重心权重和导数矩阵。
! Initialize reference LGL nodes, weights, barycentric weights, and differentiation matrices.
subroutine ensure_reference_spectral(nnode)
    integer, intent(in) :: nnode
    real(8), dimension(nnode) :: pvals
    integer :: degree, i, j

    if (allocated(reference_r)) return
    allocate(reference_r(nnode), reference_w(nnode), reference_bary(nnode), &
             reference_dmat(nnode,nnode), reference_transport(nnode,nnode), &
             reference_basis(nnode,nnode))
    call lgl_nodes_weights(nnode, reference_r, reference_w)
    degree = nnode - 1
    do i = 1, nnode
        call legendre_basis_values(degree, reference_r(i), pvals)
        reference_basis(:,i) = pvals
    enddo
    call barycentric_weights(nnode, reference_r, reference_bary)
    call differentiation_matrix(nnode, reference_r, reference_bary, reference_dmat)
    do j = 1, nnode
        do i = 1, nnode
            reference_transport(i,j) = (reference_w(j)/reference_w(i))*reference_dmat(j,i)
        enddo
    enddo
end subroutine ensure_reference_spectral

! 初始化投影用 Gauss-Legendre 积分节点。
! Initialize Gauss-Legendre quadrature nodes for projection.
subroutine ensure_projection_quadrature(nnode)
    integer, intent(in) :: nnode
    real(8), dimension(nnode) :: pvals
    integer :: degree, q

    if (allocated(projection_r)) return
    allocate(projection_r(nnode), projection_w(nnode), projection_basis(nnode,nnode))
    call gauss_nodes_weights(nnode, projection_r, projection_w)
    degree = nnode - 1
    do q = 1, nnode
        call legendre_basis_values(degree, projection_r(q), pvals)
        projection_basis(:,q) = pvals
    enddo
end subroutine ensure_projection_quadrature

! 根据活动断点设置每个 DG 单元的坐标边界。
! Set coordinate bounds for every DG element from the active breaks.
subroutine set_domain_bounds(x_min, x_max, active, n_active, mesh)
    real(8), intent(in), dimension(dg_breaks) :: active
    real(8), intent(in) :: x_min,x_max
    integer, intent(in) :: n_active
    type(dg_mesh), intent(inout) :: mesh
    real(8), dimension(dg_breaks+1) :: segment_left,segment_right
    real(8) :: dx
    integer :: element_index, i, iseg, nseg

    nseg = n_active + 1
    segment_left = 0d0
    segment_right = 0d0
    segment_left(1) = x_min
    segment_right(nseg) = x_max
    if (n_active > 0) then
        segment_left(2:nseg) = active(1:n_active)
        segment_right(1:n_active) = active(1:n_active)
    endif

    element_index = 0
    do iseg = 1, nseg
        dx = (segment_right(iseg) - segment_left(iseg))/dble(dg_segments)
        do i = 1, dg_segments
            element_index = element_index + 1
            mesh%x_left(element_index) = segment_left(iseg) + dble(i - 1)*dx
            mesh%x_right(element_index) = segment_left(iseg) + dble(i)*dx
        enddo
    enddo
end subroutine set_domain_bounds

! 填充 DG 节点的物理 gamma 坐标和坐标 Jacobian。
! Fill physical gamma coordinates and coordinate Jacobians at DG nodes.
subroutine fill_physical_nodes(mesh)
    type(dg_mesh), intent(inout) :: mesh
    real(8) :: dx
    integer :: k, i, idx

    do k = 1, mesh%ndom
        dx = mesh%x_right(k) - mesh%x_left(k)
        do i = 1, mesh%nnode
            idx = (k - 1)*mesh%nnode + i
            mesh%x(idx) = mesh%x_left(k) + 0.5d0*(mesh%r(i) + 1d0)*dx
            mesh%x_gamma(idx) = xg_from_coord(mesh%coord_kind, mesh%coord_scale, mesh%x(idx))
            mesh%gamma(idx) = dexp(mesh%x_gamma(idx))
            mesh%dxgamma_dcoord(idx) = dxg_dcoord(mesh%coord_kind, mesh%coord_scale, mesh%x(idx))
            mesh%dln_dcoord(idx) = mesh%dxgamma_dcoord(idx)
        enddo
    enddo
end subroutine fill_physical_nodes

! 生成 Legendre-Gauss-Lobatto 节点和权重。
! Generate Legendre-Gauss-Lobatto nodes and weights.
subroutine lgl_nodes_weights(nnode, r, w)
    integer, intent(in) :: nnode
    real(8), intent(out), dimension(nnode) :: r,w
    integer :: degree, i, iter
    real(8) :: x, pn, dpn, ddpn

    degree = nnode - 1
    if (degree < 1) error stop 'lgl_nodes_weights requires nnode >= 2'
    r(1) = -1d0
    r(nnode) = 1d0
    w(1) = 2d0/(dble(degree)*dble(degree + 1))
    w(nnode) = w(1)
    do i = 2, nnode - 1
        x = -dcos(pi*dble(i - 1)/dble(degree))
        do iter = 1, 50
            call legendre_value_derivative(degree, x, pn, dpn)
            ddpn = (2d0*x*dpn - dble(degree)*dble(degree + 1)*pn)/(1d0 - x*x)
            x = x - dpn/ddpn
            if (abs(dpn) <= 1d-14) exit
        enddo
        r(i) = x
        call legendre_value_derivative(degree, r(i), pn, dpn)
        w(i) = 2d0/(dble(degree)*dble(degree + 1)*pn*pn)
    enddo
end subroutine lgl_nodes_weights

! 计算 Legendre 多项式及其导数。
! Evaluate a Legendre polynomial and its derivative.
subroutine legendre_value_derivative(degree, x, pn, dpn)
    integer, intent(in) :: degree
    real(8), intent(in) :: x
    real(8), intent(out) :: pn, dpn
    real(8) :: pnm2, pnm1, pnow
    integer :: n

    if (degree == 0) then
        pn = 1d0
        dpn = 0d0
        return
    endif
    pnm2 = 1d0
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
    dpn = dble(degree)*(pnm2 - x*pn)/(1d0 - x*x)
end subroutine legendre_value_derivative

! 计算从 0 到 degree 的 Legendre 基函数值。
! Evaluate Legendre basis values from degree 0 through the requested degree.
subroutine legendre_basis_values(degree, x, pvals)
    integer, intent(in) :: degree
    real(8), intent(in) :: x
    real(8), intent(out), dimension(degree+1) :: pvals
    integer :: n

    pvals(1) = 1d0
    if (degree == 0) return
    pvals(2) = x
    do n = 2, degree
        pvals(n + 1) = (dble(2*n - 1)*x*pvals(n) - dble(n - 1)*pvals(n - 1))/dble(n)
    enddo
end subroutine legendre_basis_values

! 生成 Gauss-Legendre 投影积分节点和权重。
! Generate Gauss-Legendre projection nodes and weights.
subroutine gauss_nodes_weights(nnode, r, w)
    integer, intent(in) :: nnode
    real(8), intent(out), dimension(nnode) :: r,w
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
        w(i) = 2d0/((1d0 - r(i)*r(i))*dpn*dpn)
    enddo
end subroutine gauss_nodes_weights

! 计算 Lagrange 插值的重心权重。
! Compute barycentric weights for Lagrange interpolation.
subroutine barycentric_weights(nnode, nodes, bary)
    integer, intent(in) :: nnode
    real(8), intent(in), dimension(nnode) :: nodes
    real(8), intent(out), dimension(nnode) :: bary
    integer :: i, j

    bary = 1d0
    do i = 1, nnode
        do j = 1, nnode
            if (j /= i) bary(i) = bary(i)*(nodes(i) - nodes(j))
        enddo
        bary(i) = 1d0/bary(i)
    enddo
end subroutine barycentric_weights

! 由重心权重生成谱导数矩阵。
! Build the spectral differentiation matrix from barycentric weights.
subroutine differentiation_matrix(nnode, nodes, bary, dmat)
    integer, intent(in) :: nnode
    real(8), intent(in), dimension(nnode) :: nodes,bary
    real(8), intent(out), dimension(nnode,nnode) :: dmat
    integer :: i, j

    dmat = 0d0
    do i = 1, nnode
        do j = 1, nnode
            if (j /= i) dmat(i,j) = bary(j)/(bary(i)*(nodes(i) - nodes(j)))
        enddo
        dmat(i,i) = -sum(dmat(i,:))
    enddo
end subroutine differentiation_matrix

! 在指定 DG 单元内做重心插值。
! Interpolate within a DG element using barycentric interpolation.
real(8) function interpolate_domain(mesh, state, k, x_eval) result(value)
    type(dg_mesh), intent(in) :: mesh
    real(8), intent(in), dimension(mesh%ntot) :: state
    real(8), intent(in) :: x_eval
    integer, intent(in) :: k
    real(8) :: r_eval, numerator, denominator, diff
    integer :: i, idx, offset

    offset = (k - 1)*mesh%nnode
    r_eval = 2d0*(x_eval - mesh%x_left(k))/(mesh%x_right(k) - mesh%x_left(k)) - 1d0
    numerator = 0d0
    denominator = 0d0
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

end module electron_dg_transport
