program DNS
    implicit none

    ! Declare all the parameters which will be used 
    integer, parameter :: nx = 400, ny = 200, step_x = 20, step_y = 20, max_iter = 5000, Re = 100
    real(8), parameter :: Lx = 5.0d0, Ly = 1.0d0, rho = 1.25d0, u_inlet = 1.0d0, CFL = 0.9d0
    real(8) :: dx, dy, dt, nu, t = 0, t_final = 5.0d0
    logical, dimension(nx,ny) :: mask = .true.
    
    ! Declare the looping variables
    integer :: i, j
    
    ! Declare the variables which will be used
    real(8), dimension(nx, ny) :: u, v, p, vn, un, RHS
    real(8), dimension(nx, ny) :: u_star, v_star

    ! Declare the step geometry
    do i = 1, step_x
        do j = 1, step_y
            mask(i,j) = .false.
        end do
    end do 

    ! Initialise the variables
    dx = Lx/real(nx - 1,8)
    dy = Ly/real(ny - 1,8)

    u = 0.0d0
    v = 0.0d0
    p = 0.0d0
    u_star = 0.0d0  
    v_star = 0.0d0 
    nu = (u_inlet * Ly)/Re
    RHS = 0.0d0 

    print *, "========================================"
    print *, "Starting backward-facing step simulation"
    print *, "USING CORRECTED FRACTIONAL STEP METHOD"
    print *, "========================================"
    print *, "Grid: ", nx, " x ", ny
    print *, "Reynolds number: ", Re
    print *, "Domain: ", Lx, " x ", Ly
    print *, "dx = ", dx, ", dy = ", dy
    print *, "========================================", new_line('a')

    ! Main time loop
    i = 0
    do while (t < t_final)
        un = u
        vn = v
        
        call Set_Dt(u, v, dt, dx, dy, CFL, nx, ny, nu)
        
        call Compute_Intermediate_Velocity(u_star, v_star, u, v, mask, dx, dy, dt, nu, nx, ny)
        
        call Apply_Vel_BC(u_star, v_star, u_inlet, step_x, step_y, nx, ny)
        
        call Build_Pressure_RHS(u_star, v_star, mask, dx, dy, RHS, rho, dt, nx, ny)
        
        call Solve_Pressure_Multigrid(p, RHS, mask, dx, dy, nx, ny)
        
        call Correct_Velocity(u, v, u_star, v_star, p, mask, dx, dy, dt, rho, nx, ny)
        
        call Apply_Vel_BC(u, v, u_inlet, step_x, step_y, nx, ny)

        ! Diagnostics
        if (mod(i,50) == 0) then
            call Check_Divergence(u, v, mask, dx, dy, nx, ny, t)
        end if
        
        t = t + dt
        i = i + 1
    end do

    print *, new_line('a'), "========================================", new_line('a')
    print *, "Simulation completed!"
    print *, "Final time: ", t
    print *, "Total time steps: ", i
    print *, "========================================", new_line('a')

    call Save_Solution(u, v, p, mask, nx, ny, dx, dy)
    print *, "Solution saved to solution.dat"

contains
    subroutine Apply_Vel_BC(u, v, u_inlet,step_x, step_y, nx, ny)
        implicit none    
        integer, intent(in) :: nx, ny, step_x, step_y
        real(8), intent(inout) :: u(nx,ny), v(nx,ny) 
        real(8), intent(in) :: u_inlet
        integer :: j
        
        ! Inlet (above the step)
        do j = step_y + 1, ny
            u(1,j) = u_inlet
            v(1,j) = 0.0d0
        end do

        ! Outlet (zero gradient)
        do j = 1, ny
            u(nx,j) = u(nx-1,j)
            v(nx,j) = v(nx-1,j)
        end do

        ! Top Wall
        do j = 1, nx
            u(j,ny) = 0.0d0
            v(j,ny) = 0.0d0
        end do

        ! Bottom Wall (after the step)
        do j = step_x + 1, nx
            u(j,1) = 0.0d0
            v(j,1) = 0.0d0
        end do

        ! Step wall (Vertical face)
        do j = 1, step_y
            u(step_x + 1,j) = 0.0d0
            v(step_x + 1,j) = 0.0d0
        end do

        ! Step wall (Horizontal top face)
        do j = 1, step_x
            u(j,step_y + 1) = 0.0d0
            v(j,step_y + 1) = 0.0d0
        end do

    end subroutine Apply_Vel_BC

    subroutine Set_Dt(u, v, dt, dx, dy, CFL, nx, ny, nu)
        implicit none
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: u(nx, ny), v(nx,ny), dx, dy, CFL, nu
        real(8), intent(inout) :: dt

        ! Local Variables
        real(8)  :: dt_cfl, dt_diff

        dt_cfl = min(CFL * dx / max(maxval(abs(u)), 1e-8), &
                     CFL * dy / max(maxval(abs(v)), 1e-8))
        dt_diff = min(0.25d0 * dx**2/nu, 0.25d0 * dy**2/nu)

        dt = min(dt_cfl, dt_diff)
    end subroutine Set_Dt

    subroutine Check_Divergence(u, v, mask, dx, dy, nx, ny, t)
        implicit none
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: u(nx,ny), v(nx,ny), dx, dy, t
        logical, intent(in) :: mask(nx,ny)
        integer :: i, j
        real(8) :: div, max_div, avg_div, sum_div
        integer :: count_pts
        
        max_div = 0.0d0
        sum_div = 0.0d0
        count_pts = 0
        
        do i = 2, nx-1
            do j = 2, ny-1
                if (.not. mask(i,j)) cycle
                div = abs((u(i+1,j) - u(i-1,j))/(2.0d0*dx) + &
                         (v(i,j+1) - v(i,j-1))/(2.0d0*dy))
                max_div = max(max_div, div)
                sum_div = sum_div + div
                count_pts = count_pts + 1
            end do
        end do
        
        avg_div = sum_div / real(count_pts, 8)
        print '(A,F8.4,A,ES10.3,A,ES10.3)', &
            "  t=", t, "  Max Div=", max_div, "  Avg Div=", avg_div
    end subroutine Check_Divergence

    subroutine Save_Solution(u, v, p, mask, nx, ny, dx, dy)
        implicit none
        
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: u(nx,ny), v(nx,ny), p(nx,ny), dx, dy
        logical, intent(in) :: mask(nx,ny)
        
        integer :: i, j, mask_val
        real(8) :: x, y

        open(unit = 100, file = 'solution.dat', status = 'replace')
        write(100,*) '# x y u v p mask' 

        do j = 1, ny
            do i = 1, nx
                x = dble(i-1)*dx
                y = dble(j-1)*dy

                if (mask(i,j)) then
                    mask_val = 1
                else
                    mask_val = 0
                end if
                
                write(100,'(6(ES14.6, 1x))') x, y, u(i,j), v(i,j), p(i,j), dble(mask_val)
            end do
            write(100,*) ! Newline for gnuplot
        end do
        
        close(100)
        
        print *, "  Solution data format: x, y, u, v, p, mask"
        print *, "  Use gnuplot or Python to visualize"

    end subroutine Save_Solution

    subroutine Compute_Intermediate_Velocity(u_star, v_star, u, v, mask, dx, dy, dt, nu, nx, ny)
        implicit none
        integer, intent(in) :: nx, ny
        real(8), intent(out) :: u_star(nx,ny), v_star(nx,ny)
        real(8), intent(in) :: u(nx,ny), v(nx,ny)
        logical, intent(in) :: mask(nx,ny)
        real(8), intent(in) :: dx, dy, dt, nu
        
        integer :: i, j
        real(8) :: du_dx, du_dy, dv_dx, dv_dy
        real(8) :: d2u_dx2, d2u_dy2, d2v_dx2, d2v_dy2
        
        ! Initialize with current velocities
        u_star = u
        v_star = v
        
        do i = 2, nx - 1
            do j = 2, ny - 1
                if (.not. mask(i,j)) cycle
                
                ! Convective terms (central difference)
                du_dx = (u(i+1,j) - u(i-1,j)) / (2.0d0*dx)
                du_dy = (u(i,j+1) - u(i,j-1)) / (2.0d0*dy)
                
                dv_dx = (v(i+1,j) - v(i-1,j)) / (2.0d0*dx)
                dv_dy = (v(i,j+1) - v(i,j-1)) / (2.0d0*dy)
                
                ! Diffusive terms
                d2u_dx2 = (u(i+1,j) - 2.0d0*u(i,j) + u(i-1,j)) / (dx*dx)
                d2u_dy2 = (u(i,j+1) - 2.0d0*u(i,j) + u(i,j-1)) / (dy*dy)
                
                d2v_dx2 = (v(i+1,j) - 2.0d0*v(i,j) + v(i-1,j)) / (dx*dx)
                d2v_dy2 = (v(i,j+1) - 2.0d0*v(i,j) + v(i,j-1)) / (dy*dy)
                
                ! Update intermediate velocities (NO PRESSURE)
                u_star(i,j) = u(i,j) + dt * ( &
                    - (u(i,j)*du_dx + v(i,j)*du_dy) &  ! Convection
                    + nu * (d2u_dx2 + d2u_dy2) )        ! Diffusion
                
                v_star(i,j) = v(i,j) + dt * ( &
                    - (u(i,j)*dv_dx + v(i,j)*dv_dy) &  ! Convection
                    + nu * (d2v_dx2 + d2v_dy2) )        ! Diffusion
            end do
        end do
    end subroutine Compute_Intermediate_Velocity

    subroutine Build_Pressure_RHS(u_star, v_star, mask, dx, dy, RHS, rho, dt, nx, ny)
        implicit none
        integer, intent(in) :: nx, ny
        real(8), dimension(nx, ny), intent(in) :: u_star, v_star
        logical, dimension(nx,ny), intent(in) :: mask
        real(8), intent(in) :: dy, dx, rho, dt
        real(8), dimension(nx, ny), intent(out) :: RHS
        
        real(8) :: du_dx, dv_dy
        integer :: i, j
        
        RHS = 0.0d0
        
        ! RHS = (rho/dt) * div(u_star)
        do i = 2, nx - 1
            do j = 2, ny - 1
                if (.not. mask(i,j)) cycle
                
                ! Divergence of intermediate velocity
                du_dx = (u_star(i+1,j) - u_star(i-1,j)) / (2.0d0*dx)
                dv_dy = (v_star(i,j+1) - v_star(i,j-1)) / (2.0d0*dy)
                
                RHS(i,j) = (rho / dt) * (du_dx + dv_dy)
            end do
        end do
        
        ! Boundary conditions for RHS
        RHS(1,:) = 0.0d0
        RHS(nx,:) = 0.0d0
        RHS(:,1) = 0.0d0
        RHS(:,ny) = 0.0d0
        
    end subroutine Build_Pressure_RHS

    subroutine Correct_Velocity(u, v, u_star, v_star, p, mask, dx, dy, dt, rho, nx, ny)
        implicit none
        integer, intent(in) :: nx, ny
        real(8), intent(out) :: u(nx,ny), v(nx,ny)
        real(8), intent(in) :: u_star(nx,ny), v_star(nx,ny), p(nx,ny)
        logical, intent(in) :: mask(nx,ny)
        real(8), intent(in) :: dx, dy, dt, rho
        
        integer :: i, j
        real(8) :: dp_dx, dp_dy
        
        ! Initialize with intermediate velocities
        u = u_star
        v = v_star
        
        ! Apply pressure correction
        do i = 2, nx - 1
            do j = 2, ny - 1
                if (.not. mask(i,j)) cycle
                
                ! Pressure gradients
                dp_dx = (p(i+1,j) - p(i-1,j)) / (2.0d0*dx)
                dp_dy = (p(i,j+1) - p(i,j-1)) / (2.0d0*dy)
                
                ! Correct velocities to enforce incompressibility
                u(i,j) = u_star(i,j) - (dt/rho) * dp_dx
                v(i,j) = v_star(i,j) - (dt/rho) * dp_dy
            end do
        end do
    end subroutine Correct_Velocity

    subroutine Solve_Pressure_Multigrid(p, RHS, mask, dx, dy, nx, ny)
        implicit none
        
        real(8), dimension(nx, ny), intent(inout) :: p, RHS
        logical, dimension(nx,ny), intent(in) :: mask
        real(8), intent(in) :: dy, dx
        integer, intent(in) :: nx, ny
        
        ! Local variables
        integer :: v1 = 2, v2 = 1  
        integer :: cycle, max_cycles = 10
        real(8) :: tolerance = 1.0e-4
        real(8), dimension(nx,ny) :: p_old
        real(8) :: error
        
        ! Initialize pressure BC
        call Apply_Pressure_BC(p, mask, nx, ny)
        
        ! Do multiple V-cycles until convergence
        do cycle = 1, max_cycles
            p_old = p
            
            ! One V-cycle
            call V_Cycle(p, RHS, mask, dx, dy, nx, ny, v1, v2, 1)
            
            ! Apply BC after V-cycle
            call Apply_Pressure_BC(p, mask, nx, ny)
            
            ! Check convergence
            error = maxval(abs(p - p_old), mask=mask)
            
            ! if (mod(cycle, 2) == 0 .or. cycle == 1) then
            !     print '(A,I3,A,ES10.3)', "  MG cycle ", cycle, "  error=", error
            ! end if
            
            ! if (error < tolerance .and. cycle > 2) then
            !     print *, "  Converged in ", cycle, " cycles"
            !     exit
            ! end if
        end do
        
    end subroutine Solve_Pressure_Multigrid

    recursive subroutine V_Cycle(p, RHS, mask, dx, dy, nx, ny, v1, v2, level)
        implicit none
        integer, intent(in) :: nx, ny, v1, v2, level
        real(8), dimension(nx, ny), intent(inout) :: p
        real(8), dimension(nx, ny), intent(in) :: RHS
        logical, dimension(nx,ny), intent(in) :: mask
        real(8), intent(in) :: dx, dy
        
        ! Local variables
        integer :: nx_coarse, ny_coarse, i, j
        real(8), allocatable :: residual(:,:), error_coarse(:,:), RHS_coarse(:,:)
        real(8), allocatable :: error_fine(:,:)
        logical, allocatable :: mask_coarse(:,:)
        real(8) :: dx_coarse, dy_coarse
        real(8), parameter :: damping = 0.8d0  ! Damping factor for stability
        
        ! Base case: solve directly on coarsest grid
        if (nx <= 25 .or. ny <= 25) then
            call Smooth_GS(p, RHS, mask, dx, dy, nx, ny, 50)
            return
        end if
        
        ! Pre-smoothing
        call Smooth_GS(p, RHS, mask, dx, dy, nx, ny, v1)
        
        ! Compute residual
        allocate(residual(nx, ny))
        call Compute_Residual(p, RHS, residual, mask, dx, dy, nx, ny)
        
        ! Restrict to coarse grid
        nx_coarse = (nx + 1) / 2
        ny_coarse = (ny + 1) / 2
        allocate(RHS_coarse(nx_coarse, ny_coarse))
        allocate(mask_coarse(nx_coarse, ny_coarse))
        allocate(error_coarse(nx_coarse, ny_coarse))
        
        ! Use simple injection for stability with complex geometry
        call Restrict_Injection(residual, RHS_coarse, mask, mask_coarse, &
                                nx, ny, nx_coarse, ny_coarse)
        
        error_coarse = 0.0d0
        dx_coarse = dx * 2.0d0
        dy_coarse = dy * 2.0d0
        
        ! Recursive call on coarse grid
        call V_Cycle(error_coarse, RHS_coarse, mask_coarse, dx_coarse, dy_coarse, &
                    nx_coarse, ny_coarse, v1, v2, level + 1)
        
        ! Prolongate error back to fine grid
        allocate(error_fine(nx, ny))
        call Prolongate_Simple(error_coarse, error_fine, mask_coarse, &
                            nx_coarse, ny_coarse, nx, ny)
        
        ! Correct solution with damping for stability
        do i = 1, nx
            do j = 1, ny
                if (mask(i,j)) then
                    p(i,j) = p(i,j) + damping * error_fine(i,j)
                end if
            end do
        end do
        
        ! Post-smoothing
        call Smooth_GS(p, RHS, mask, dx, dy, nx, ny, v2)
        
        ! Cleanup
        deallocate(residual, RHS_coarse, mask_coarse, error_coarse, error_fine)
        
    end subroutine V_Cycle

    subroutine Smooth_GS(p, RHS, mask, dx, dy, nx, ny, iterations)
        implicit none
        integer, intent(in) :: nx, ny, iterations
        real(8), dimension(nx, ny), intent(inout) :: p
        real(8), dimension(nx, ny), intent(in) :: RHS
        logical, dimension(nx,ny), intent(in) :: mask
        real(8), intent(in) :: dx, dy
        
        integer :: iter, i, j
        real(8) :: dx2, dy2, coeff
        
        dx2 = dx * dx
        dy2 = dy * dy
        coeff = 1.0d0 / (2.0d0 * (dx2 + dy2))
        
        do iter = 1, iterations
            ! Gauss-Seidel sweep
            do i = 2, nx - 1
                do j = 2, ny - 1
                    if (.not. mask(i,j)) cycle
                    
                    p(i,j) = coeff * ( &
                        dy2 * (p(i+1,j) + p(i-1,j)) + &
                        dx2 * (p(i,j+1) + p(i,j-1)) - &
                        dx2 * dy2 * RHS(i,j) )
                end do
            end do
            
            ! Apply BC every few iterations
            if (mod(iter, 5) == 0 .or. iter == iterations) then
                call Apply_Pressure_BC(p, mask, nx, ny)
            end if
        end do
    end subroutine Smooth_GS

    subroutine Apply_Pressure_BC(p, mask, nx, ny)
        implicit none
        integer, intent(in) :: nx, ny
        real(8), dimension(nx, ny), intent(inout) :: p
        logical, dimension(nx,ny), intent(in) :: mask
        integer :: i, j
        
        ! Neumann BC: dp/dn = 0
        do j = 1, ny
            if (mask(1,j)) p(1,j) = p(2,j)      ! Inlet
            if (mask(nx,j)) p(nx,j) = 0         ! Outlet
        end do
        
        do i = 1, nx
            if (mask(i,1)) p(i,1) = p(i,2)      ! Bottom
            if (mask(i,ny)) p(i,ny) = p(i,ny-1) ! Top
        end do
        
        ! Handle masked regions (set to zero or neighbor average)
        do i = 1, nx
            do j = 1, ny
                if (.not. mask(i,j)) p(i,j) = 0.0d0
            end do
        end do
    end subroutine Apply_Pressure_BC

    subroutine Compute_Residual(p, RHS, residual, mask, dx, dy, nx, ny)
        implicit none
        integer, intent(in) :: nx, ny
        real(8), dimension(nx, ny), intent(in) :: p, RHS
        real(8), dimension(nx, ny), intent(out) :: residual
        logical, dimension(nx,ny), intent(in) :: mask
        real(8), intent(in) :: dx, dy
        
        integer :: i, j
        real(8) :: dx2, dy2, laplacian
        
        dx2 = dx * dx
        dy2 = dy * dy
        residual = 0.0d0
        
        do i = 2, nx - 1
            do j = 2, ny - 1
                if (.not. mask(i,j)) cycle
                
                ! Laplacian
                laplacian = (p(i+1,j) - 2.0d0*p(i,j) + p(i-1,j)) / dx2 + &
                        (p(i,j+1) - 2.0d0*p(i,j) + p(i,j-1)) / dy2
                
                ! Residual = RHS - Laplacian(p)
                residual(i,j) = RHS(i,j) - laplacian
            end do
        end do
        
        ! Zero at boundaries
        residual(1,:) = 0.0d0
        residual(nx,:) = 0.0d0
        residual(:,1) = 0.0d0
        residual(:,ny) = 0.0d0
        
    end subroutine Compute_Residual

    subroutine Restrict_Injection(fine, coarse, mask_fine, mask_coarse, &
                                nx_fine, ny_fine, nx_coarse, ny_coarse)
        implicit none
        integer, intent(in) :: nx_fine, ny_fine, nx_coarse, ny_coarse
        real(8), dimension(nx_fine, ny_fine), intent(in) :: fine
        real(8), dimension(nx_coarse, ny_coarse), intent(out) :: coarse
        logical, dimension(nx_fine, ny_fine), intent(in) :: mask_fine
        logical, dimension(nx_coarse, ny_coarse), intent(out) :: mask_coarse
        
        integer :: i, j, if, jf
        
        coarse = 0.0d0
        mask_coarse = .false.
        
        ! Simple injection: take every other point
        do i = 1, nx_coarse
            do j = 1, ny_coarse
                if = 2*i - 1
                jf = 2*j - 1
                
                if (if <= nx_fine .and. jf <= ny_fine) then
                    mask_coarse(i,j) = mask_fine(if, jf)
                    if (mask_coarse(i,j)) then
                        coarse(i,j) = fine(if, jf)
                    end if
                end if
            end do
        end do
    end subroutine Restrict_Injection

    subroutine Prolongate_Simple(coarse, fine, mask_coarse, nx_coarse, ny_coarse, nx_fine, ny_fine)
        implicit none
        integer, intent(in) :: nx_coarse, ny_coarse, nx_fine, ny_fine
        real(8), dimension(nx_coarse, ny_coarse), intent(in) :: coarse
        real(8), dimension(nx_fine, ny_fine), intent(out) :: fine
        logical, dimension(nx_coarse, ny_coarse), intent(in) :: mask_coarse
        
        integer :: i, j, ic, jc
        integer :: i_fine, j_fine
        
        fine = 0.0d0
        
        ! Method: Direct injection + simple averaging
        ! This is more stable than full bilinear near boundaries
        
        ! Step 1: Inject coarse grid points directly
        do ic = 1, nx_coarse
            do jc = 1, ny_coarse
                if (.not. mask_coarse(ic, jc)) cycle
                
                i_fine = 2*ic - 1
                j_fine = 2*jc - 1
                
                if (i_fine <= nx_fine .and. j_fine <= ny_fine) then
                    fine(i_fine, j_fine) = coarse(ic, jc)
                end if
            end do
        end do
        
        ! Step 2: Interpolate horizontal edges
        do ic = 1, nx_coarse - 1
            do jc = 1, ny_coarse
                if (.not. (mask_coarse(ic,jc) .and. mask_coarse(ic+1,jc))) cycle
                
                i_fine = 2*ic
                j_fine = 2*jc - 1
                
                if (i_fine <= nx_fine .and. j_fine <= ny_fine) then
                    fine(i_fine, j_fine) = 0.5d0 * (coarse(ic,jc) + coarse(ic+1,jc))
                end if
            end do
        end do
        
        ! Step 3: Interpolate vertical edges
        do ic = 1, nx_coarse
            do jc = 1, ny_coarse - 1
                if (.not. (mask_coarse(ic,jc) .and. mask_coarse(ic,jc+1))) cycle
                
                i_fine = 2*ic - 1
                j_fine = 2*jc
                
                if (i_fine <= nx_fine .and. j_fine <= ny_fine) then
                    fine(i_fine, j_fine) = 0.5d0 * (coarse(ic,jc) + coarse(ic,jc+1))
                end if
            end do
        end do
        
        ! Step 4: Interpolate center points (bilinear)
        do ic = 1, nx_coarse - 1
            do jc = 1, ny_coarse - 1
                if (.not. (mask_coarse(ic,jc) .and. mask_coarse(ic+1,jc) .and. &
                        mask_coarse(ic,jc+1) .and. mask_coarse(ic+1,jc+1))) cycle
                
                i_fine = 2*ic
                j_fine = 2*jc
                
                if (i_fine <= nx_fine .and. j_fine <= ny_fine) then
                    fine(i_fine, j_fine) = 0.25d0 * ( &
                        coarse(ic,jc) + coarse(ic+1,jc) + &
                        coarse(ic,jc+1) + coarse(ic+1,jc+1) )
                end if
            end do
        end do
        
    end subroutine Prolongate_Simple

end program DNS