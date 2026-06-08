program cpu_wavefront_routes_probe
    use omp_lib
    implicit none
    integer, parameter :: dp=kind(1.0d0)
    integer :: n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,n_repeat
    real(dp), allocatable :: face(:,:,:),source(:,:,:),input(:,:),reference(:,:,:),candidate(:,:,:)

    call read_args(n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,n_repeat)
    allocate(face(n_gamma-1,n_step,n_items),source(n_gamma,n_step,n_items),input(n_gamma,n_items), &
             reference(n_gamma,n_step,n_items),candidate(n_gamma,n_step,n_items))
    call fill_problem(n_items,n_gamma,n_step,face,source,input)
    call print_header()
    call benchmark_routes(n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,n_repeat, &
                          face,source,input,reference,candidate)
contains

subroutine read_args(n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,n_repeat)
    integer, intent(out) :: n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,n_repeat
    character(len=64) :: arg

    n_threads=8; n_items=300; n_gamma=96; n_step=300; tile_gamma=32; tile_step=32; n_repeat=3
    call get_command_argument(1,arg); if (len_trim(arg) > 0) read(arg,*) n_threads
    call get_command_argument(2,arg); if (len_trim(arg) > 0) read(arg,*) n_items
    call get_command_argument(3,arg); if (len_trim(arg) > 0) read(arg,*) n_gamma
    call get_command_argument(4,arg); if (len_trim(arg) > 0) read(arg,*) n_step
    call get_command_argument(5,arg); if (len_trim(arg) > 0) read(arg,*) tile_gamma
    call get_command_argument(6,arg); if (len_trim(arg) > 0) read(arg,*) tile_step
    call get_command_argument(7,arg); if (len_trim(arg) > 0) read(arg,*) n_repeat
end subroutine read_args

subroutine fill_problem(n_items,n_gamma,n_step,face,source,input)
    integer, intent(in) :: n_items,n_gamma,n_step
    real(dp), intent(out) :: face(n_gamma-1,n_step,n_items),source(n_gamma,n_step,n_items),input(n_gamma,n_items)
    integer :: item,i_gamma,i_step
    real(dp) :: item_scale,gamma_x,step_x

    do item=1,n_items
        item_scale=1.0_dp+0.0007_dp*dble(item-1)
        do i_gamma=1,n_gamma
            gamma_x=dble(i_gamma)/dble(n_gamma)
            input(i_gamma,item)=item_scale*exp(-3.0_dp*gamma_x)
            do i_step=1,n_step
                step_x=dble(i_step)/dble(n_step)
                source(i_gamma,i_step,item)=item_scale*1.0e-5_dp*(1.0_dp+gamma_x)*(1.0_dp+0.25_dp*step_x)
            end do
        end do
        do i_step=1,n_step
            step_x=dble(i_step)/dble(n_step)
            do i_gamma=1,n_gamma-1
                gamma_x=dble(i_gamma)/dble(n_gamma)
                face(i_gamma,i_step,item)=0.005_dp*item_scale*(1.0_dp+0.5_dp*gamma_x+0.25_dp*step_x)
            end do
        end do
    end do
end subroutine fill_problem

subroutine print_header()
    print '(A)', 'route,threads,items,gamma,steps,tile_gamma,tile_step,seconds,baseline_seconds,speedup,max_abs_err,max_rel_err'
end subroutine print_header

subroutine benchmark_routes(n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,n_repeat, &
                            face,source,input,reference,candidate)
    integer, intent(in) :: n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,n_repeat
    real(dp), intent(in) :: face(n_gamma-1,n_step,n_items),source(n_gamma,n_step,n_items),input(n_gamma,n_items)
    real(dp), intent(out) :: reference(n_gamma,n_step,n_items),candidate(n_gamma,n_step,n_items)
    real(dp) :: baseline_s,route_s,max_abs,max_rel

    call omp_set_num_threads(n_threads)
    call solve_many_sequence(n_items,n_gamma,n_step,face,source,input,reference)
    baseline_s=time_sequence(n_repeat,n_items,n_gamma,n_step,face,source,input,reference)
    route_s=time_route3(n_repeat,n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,face,source,input,candidate)
    call max_error(n_items,n_gamma,n_step,reference,candidate,max_abs,max_rel)
    call print_row('3_task_depend',n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,route_s,baseline_s,max_abs,max_rel)
    route_s=time_route2(n_repeat,n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,face,source,input,candidate)
    call max_error(n_items,n_gamma,n_step,reference,candidate,max_abs,max_rel)
    call print_row('2_tile_wave',n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,route_s,baseline_s,max_abs,max_rel)
    route_s=time_route1(n_repeat,n_threads,n_items,n_gamma,n_step,face,source,input,candidate)
    call max_error(n_items,n_gamma,n_step,reference,candidate,max_abs,max_rel)
    call print_row('1_item_batch',n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,route_s,baseline_s,max_abs,max_rel)
end subroutine benchmark_routes

function time_sequence(n_repeat,n_items,n_gamma,n_step,face,source,input,out) result(seconds)
    integer, intent(in) :: n_repeat,n_items,n_gamma,n_step
    real(dp), intent(in) :: face(n_gamma-1,n_step,n_items),source(n_gamma,n_step,n_items),input(n_gamma,n_items)
    real(dp), intent(out) :: out(n_gamma,n_step,n_items)
    integer :: i_repeat
    real(dp) :: seconds,t0

    t0=omp_get_wtime()
    do i_repeat=1,n_repeat
        call solve_many_sequence(n_items,n_gamma,n_step,face,source,input,out)
    end do
    seconds=(omp_get_wtime()-t0)/dble(n_repeat)
end function time_sequence

function time_route3(n_repeat,n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,face,source,input,out) result(seconds)
    integer, intent(in) :: n_repeat,n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step
    real(dp), intent(in) :: face(n_gamma-1,n_step,n_items),source(n_gamma,n_step,n_items),input(n_gamma,n_items)
    real(dp), intent(out) :: out(n_gamma,n_step,n_items)
    integer :: i_repeat
    real(dp) :: seconds,t0

    t0=omp_get_wtime()
    do i_repeat=1,n_repeat
        call solve_many_task_depend(n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,face,source,input,out)
    end do
    seconds=(omp_get_wtime()-t0)/dble(n_repeat)
end function time_route3

function time_route2(n_repeat,n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,face,source,input,out) result(seconds)
    integer, intent(in) :: n_repeat,n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step
    real(dp), intent(in) :: face(n_gamma-1,n_step,n_items),source(n_gamma,n_step,n_items),input(n_gamma,n_items)
    real(dp), intent(out) :: out(n_gamma,n_step,n_items)
    integer :: i_repeat
    real(dp) :: seconds,t0

    t0=omp_get_wtime()
    do i_repeat=1,n_repeat
        call solve_many_tile_wave(n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,face,source,input,out)
    end do
    seconds=(omp_get_wtime()-t0)/dble(n_repeat)
end function time_route2

function time_route1(n_repeat,n_threads,n_items,n_gamma,n_step,face,source,input,out) result(seconds)
    integer, intent(in) :: n_repeat,n_threads,n_items,n_gamma,n_step
    real(dp), intent(in) :: face(n_gamma-1,n_step,n_items),source(n_gamma,n_step,n_items),input(n_gamma,n_items)
    real(dp), intent(out) :: out(n_gamma,n_step,n_items)
    integer :: i_repeat
    real(dp) :: seconds,t0

    t0=omp_get_wtime()
    do i_repeat=1,n_repeat
        call solve_many_item_batch(n_threads,n_items,n_gamma,n_step,face,source,input,out)
    end do
    seconds=(omp_get_wtime()-t0)/dble(n_repeat)
end function time_route1

subroutine solve_many_sequence(n_items,n_gamma,n_step,face,source,input,out)
    integer, intent(in) :: n_items,n_gamma,n_step
    real(dp), intent(in) :: face(n_gamma-1,n_step,n_items),source(n_gamma,n_step,n_items),input(n_gamma,n_items)
    real(dp), intent(out) :: out(n_gamma,n_step,n_items)
    integer :: item

    do item=1,n_items
        call solve_item_sequence(n_gamma,n_step,face(:,:,item),source(:,:,item),input(:,item),out(:,:,item))
    end do
end subroutine solve_many_sequence

subroutine solve_many_item_batch(n_threads,n_items,n_gamma,n_step,face,source,input,out)
    integer, intent(in) :: n_threads,n_items,n_gamma,n_step
    real(dp), intent(in) :: face(n_gamma-1,n_step,n_items),source(n_gamma,n_step,n_items),input(n_gamma,n_items)
    real(dp), intent(out) :: out(n_gamma,n_step,n_items)
    integer :: item

    !$omp parallel do private(item) schedule(dynamic) num_threads(n_threads)
    do item=1,n_items
        call solve_item_sequence(n_gamma,n_step,face(:,:,item),source(:,:,item),input(:,item),out(:,:,item))
    end do
    !$omp end parallel do
end subroutine solve_many_item_batch

subroutine solve_item_sequence(n_gamma,n_step,face,source,input,out)
    integer, intent(in) :: n_gamma,n_step
    real(dp), intent(in) :: face(n_gamma-1,n_step),source(n_gamma,n_step),input(n_gamma)
    real(dp), intent(out) :: out(n_gamma,n_step)
    integer :: i_gamma,i_step
    real(dp) :: left_face,right_face,previous_step,upper_neighbor,rhs,diag

    do i_step=1,n_step
        do i_gamma=n_gamma,1,-1
            left_face=0.0_dp; right_face=0.0_dp
            if (i_gamma > 1) left_face=face(i_gamma-1,i_step)
            if (i_gamma < n_gamma) right_face=face(i_gamma,i_step)
            diag=1.0_dp+left_face
            if (i_gamma == 1) diag=1.0_dp+right_face
            previous_step=0.0_dp; upper_neighbor=0.0_dp
            if (i_step > 1) previous_step=out(i_gamma,i_step-1)
            if (i_gamma < n_gamma) upper_neighbor=out(i_gamma+1,i_step)
            rhs=source(i_gamma,i_step)+previous_step+right_face*upper_neighbor
            if (i_step == 1) rhs=rhs+input(i_gamma)
            out(i_gamma,i_step)=rhs/diag
        end do
    end do
end subroutine solve_item_sequence

subroutine solve_many_tile_wave(n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,face,source,input,out)
    integer, intent(in) :: n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step
    real(dp), intent(in) :: face(n_gamma-1,n_step,n_items),source(n_gamma,n_step,n_items),input(n_gamma,n_items)
    real(dp), intent(out) :: out(n_gamma,n_step,n_items)
    integer :: n_gamma_block,n_step_block,wave,item,block_step,j_gamma,block_gamma

    out=0.0_dp
    n_gamma_block=(n_gamma+tile_gamma-1)/tile_gamma
    n_step_block=(n_step+tile_step-1)/tile_step
    !$omp parallel num_threads(n_threads) private(wave,item,block_step,j_gamma,block_gamma)
    do wave=2,n_gamma_block+n_step_block
        !$omp do collapse(2) schedule(dynamic)
        do item=1,n_items
            do block_step=1,n_step_block
                j_gamma=wave-block_step
                if (j_gamma >= 1 .and. j_gamma <= n_gamma_block) then
                    block_gamma=n_gamma_block-j_gamma+1
                    call solve_tile(n_gamma,n_step,tile_gamma,tile_step,block_gamma,block_step, &
                                    face(:,:,item),source(:,:,item),input(:,item),out(:,:,item))
                end if
            end do
        end do
        !$omp end do
    end do
    !$omp end parallel
end subroutine solve_many_tile_wave

subroutine solve_many_task_depend(n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,face,source,input,out)
    integer, intent(in) :: n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step
    real(dp), intent(in) :: face(n_gamma-1,n_step,n_items),source(n_gamma,n_step,n_items),input(n_gamma,n_items)
    real(dp), intent(out) :: out(n_gamma,n_step,n_items)
    integer :: n_gamma_block,n_step_block,item,block_step,j_gamma,block_gamma
    integer, allocatable :: token(:,:,:)

    out=0.0_dp
    n_gamma_block=(n_gamma+tile_gamma-1)/tile_gamma
    n_step_block=(n_step+tile_step-1)/tile_step
    allocate(token(0:n_gamma_block+1,0:n_step_block,n_items))
    token=0
    !$omp parallel num_threads(n_threads)
    !$omp single
    do item=1,n_items
        do block_step=1,n_step_block
            do j_gamma=1,n_gamma_block
                block_gamma=n_gamma_block-j_gamma+1
                !$omp task firstprivate(item,block_gamma,block_step) &
                !$omp depend(in: token(block_gamma+1,block_step,item), token(block_gamma,block_step-1,item)) &
                !$omp depend(out: token(block_gamma,block_step,item))
                call solve_tile(n_gamma,n_step,tile_gamma,tile_step,block_gamma,block_step, &
                                face(:,:,item),source(:,:,item),input(:,item),out(:,:,item))
                token(block_gamma,block_step,item)=1
                !$omp end task
            end do
        end do
    end do
    !$omp end single
    !$omp end parallel
    deallocate(token)
end subroutine solve_many_task_depend

subroutine solve_tile(n_gamma,n_step,tile_gamma,tile_step,block_gamma,block_step,face,source,input,out)
    integer, intent(in) :: n_gamma,n_step,tile_gamma,tile_step,block_gamma,block_step
    real(dp), intent(in) :: face(n_gamma-1,n_step),source(n_gamma,n_step),input(n_gamma)
    real(dp), intent(inout) :: out(n_gamma,n_step)
    integer :: i_gamma,i_step,gamma_lo,gamma_hi,step_lo,step_hi
    real(dp) :: left_face,right_face,previous_step,upper_neighbor,rhs,diag

    gamma_lo=(block_gamma-1)*tile_gamma+1
    gamma_hi=min(block_gamma*tile_gamma,n_gamma)
    step_lo=(block_step-1)*tile_step+1
    step_hi=min(block_step*tile_step,n_step)
    do i_step=step_lo,step_hi
        do i_gamma=gamma_hi,gamma_lo,-1
            left_face=0.0_dp; right_face=0.0_dp
            if (i_gamma > 1) left_face=face(i_gamma-1,i_step)
            if (i_gamma < n_gamma) right_face=face(i_gamma,i_step)
            diag=1.0_dp+left_face
            if (i_gamma == 1) diag=1.0_dp+right_face
            previous_step=0.0_dp; upper_neighbor=0.0_dp
            if (i_step > 1) previous_step=out(i_gamma,i_step-1)
            if (i_gamma < n_gamma) upper_neighbor=out(i_gamma+1,i_step)
            rhs=source(i_gamma,i_step)+previous_step+right_face*upper_neighbor
            if (i_step == 1) rhs=rhs+input(i_gamma)
            out(i_gamma,i_step)=rhs/diag
        end do
    end do
end subroutine solve_tile

subroutine max_error(n_items,n_gamma,n_step,reference,candidate,max_abs,max_rel)
    integer, intent(in) :: n_items,n_gamma,n_step
    real(dp), intent(in) :: reference(n_gamma,n_step,n_items),candidate(n_gamma,n_step,n_items)
    real(dp), intent(out) :: max_abs,max_rel
    integer :: item,i_gamma,i_step
    real(dp) :: diff,scale

    max_abs=0.0_dp; max_rel=0.0_dp
    do item=1,n_items
        do i_step=1,n_step
            do i_gamma=1,n_gamma
                diff=abs(reference(i_gamma,i_step,item)-candidate(i_gamma,i_step,item))
                scale=max(abs(reference(i_gamma,i_step,item)),1.0e-300_dp)
                max_abs=max(max_abs,diff)
                max_rel=max(max_rel,diff/scale)
            end do
        end do
    end do
end subroutine max_error

subroutine print_row(route,n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,route_s,baseline_s,max_abs,max_rel)
    character(len=*), intent(in) :: route
    integer, intent(in) :: n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step
    real(dp), intent(in) :: route_s,baseline_s,max_abs,max_rel

    print '(A,",",I0,",",I0,",",I0,",",I0,",",I0,",",I0,",",ES14.6,",",ES14.6,",",F10.4,",",ES14.6,",",ES14.6)', &
          trim(route),n_threads,n_items,n_gamma,n_step,tile_gamma,tile_step,route_s,baseline_s,baseline_s/route_s,max_abs,max_rel
end subroutine print_row

end program cpu_wavefront_routes_probe
