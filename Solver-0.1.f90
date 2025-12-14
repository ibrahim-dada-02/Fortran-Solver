program DNS
    implicit none

    ! Declare all the parameters which will be used 
    integer, parameter :: nx = 200, ny = 100, step_x = 50, step_y = 50, max_iter = 5000, Re = 100
    real(8), parameter :: Lx = 5.0d0, Ly = 1.0d0, rho = 1.25d0, u_inlet = 1.0d0, CFL = 0.9d0
    real(8) :: dx, dy, dt, nu, t = 0, t_final = 5.0d0
    logical, dimension(nx,ny) :: mask = .true.
    
    ! Declare the looping variables
    integer :: i, j
    
    ! Declare the variables which will be used
    real(8), dimension(nx, ny) :: u, v, p, vn, un, RHS

    ! Declare the step
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
    nu = (u_inlet * Ly)/Re
    RHS = 0.0d0 

    print *, "Starting backward-facing step simulation"
    print *, "Grid: ", nx, " x ", ny
    print *, "Reynolds number: ", Re

    ! Loop
    i = 0
    do while (t < t_final)
        un = u
        vn = v
        ! Set dt
        call Set_Dt(u, v, dt, dx, dy, CFL, nx, ny, nu)

        ! Build RHS for pressure Poisson equation
        call Build_Pressure_RHS(u, v, mask, dx, dy,  RHS, rho, dt,  nx, ny)
        
        ! Solve pressure Poisson equation
        call Solve_Pressure_SOR(p, RHS, mask, dx, dy, max_iter, nx, ny)
        
        ! Update velocity field
        call Update_and_Correct_Velocity(u, v, un, vn, p, mask, dx, dy, dt, nu, nx, ny, rho)
        
        ! Apply boundary conditions
        call Apply_Vel_BC(u, v, u_inlet,step_x, step_y, nx, ny)

        if (mod(i,5) == 0) call Check_Divergence(u, v, mask, dx, dy, nx, ny, t)
        
        t = t + dt
        i = i + 1
    end do

contains
    ! Solving Pressure has 2 parts, building RHS and using a numerical scheme to calculate using interations
    subroutine Build_Pressure_RHS(u, v, mask, dx, dy, RHS, rho, dt, nx, ny) 
        implicit none

        ! Declare the dummy variables
        integer, intent(in) :: nx, ny
        real(8), dimension(nx, ny), intent(in) :: u, v
        logical, dimension(nx,ny), intent(in) :: mask
        real(8), intent(in) :: dy, dx, rho, dt
        real(8), dimension(nx, ny), intent(out) :: RHS

        ! Declare the the local variables
        real(8) :: du_dx, dv_dy

        ! Declare the looping variables
        integer :: i , j
 
        ! RHS = rho/dt * div
        ! div = du/dx + dv/dy

        ! Calculate the derivative using central difference method
        do i = 2, nx - 1
            do j = 2, ny - 1
                if (.not. mask(i,j)) cycle

                ! Calculate Div
                du_dx = (u(i + 1,j) - u(i - 1,j))/ (2 * dx)
                dv_dy = (v(i,j + 1) - v(i,j - 1))/ (2 * dy)

                ! Calculate RHS
                RHS(i,j) = (rho / dt) * (du_dx + dv_dy) 

            end do
        end do
    end subroutine Build_Pressure_RHS

    subroutine Solve_Pressure_SOR(p, RHS, mask, dx, dy, max_iter, nx, ny)
        implicit none
        
        ! Declare the dummy variables
        real(8), dimension(nx, ny), intent(inout) :: p, RHS
        logical, dimension(nx,ny), intent(in) :: mask
        real(8), intent(in) :: dy, dx
        integer, intent(in) :: max_iter, nx, ny

        ! New Pressure
        real(8), dimension(nx,ny) :: pn
        real(8) :: dx2, dy2, coeff, omega = 1.7, p_gs

        ! Declare the looping variables
        integer :: i , j, iter

        ! Find squares
        dx2 = dx ** 2
        dy2 = dy ** 2

        ! Calculate Coeff
        coeff = 1.0d0/(2.0d0*(dx2 + dy2))

        ! Calculate the double derivative:
        do iter = 1, max_iter
            ! Pressure Iterations
            pn = p
            do i = 2, nx - 1
                do j = 2, ny - 1
                    if (.not. mask(i,j)) cycle
                    p_gs = coeff * ( dy2*(p(i+1,j) + p(i-1,j)) &
                        + dx2*(p(i,j+1) + p(i,j-1)) &
                        - dx2*dy2*RHS(i,j) )

                    p(i,j) = (1.0d0 - omega)*p(i,j) + omega*p_gs
                end do
            end do
            
            ! Fix Boundary conditions
            do j = 1, ny
                if (mask(1,j)) p(1,j) = p(2,j)
                if (mask(nx,j)) p(nx,j) = p(nx - 1,j)
            end do
            do i = 1, nx
                if (mask(i,1)) p(i,1) = p(i,2)
                if (mask(i,ny)) p(i,ny) = p(i,ny-1)
            end do

            if (maxval(abs(p - pn)) < 1e-4 .and. iter > 5) exit          
        end do
    end subroutine Solve_Pressure_SOR

    subroutine Update_and_Correct_Velocity(u, v, un, vn, p, mask, dx, dy, dt, nu, nx, ny, rho)
        implicit none
        integer, intent(in) :: nx, ny
        real(8), intent(inout) :: u(nx,ny), v(nx,ny)
        real(8), intent(in) :: un(nx,ny), vn(nx,ny), p(nx,ny)
        logical, intent(in) :: mask(nx,ny)
        real(8), intent(in) :: dx, dy, dt, nu, rho
        integer :: i, j
        real(8) :: du_dx, du_dy, dv_dx, dv_dy, d2u_dx2, d2u_dy2, d2v_dx2, d2v_dy2
        
        do i = 2, nx - 1
            do j = 2, ny - 1
                ! Mask Check
                if (.not. mask(i,j))cycle

                ! Convective Terms
                du_dx = (un(i+1,j) - un(i-1,j))/(2.0d0*dx)
                du_dy = (un(i,j+1) - un(i,j-1))/(2.0d0*dy)
                    
                dv_dx = (vn(i+1,j) - vn(i-1,j))/(2.0d0*dx)
                dv_dy = (vn(i,j+1) - vn(i,j-1))/(2.0d0*dy)
                
                ! Diffusive terms
                d2u_dx2 = (un(i+1,j) - 2.0d0*un(i,j) + un(i-1,j))/(dx*dx)
                d2u_dy2 = (un(i,j+1) - 2.0d0*un(i,j) + un(i,j-1))/(dy*dy)
                    
                d2v_dx2 = (vn(i+1,j) - 2.0d0*vn(i,j) + vn(i-1,j))/(dx*dx)
                d2v_dy2 = (vn(i,j+1) - 2.0d0*vn(i,j) + vn(i,j-1))/(dy*dy)

                ! Update velocities
                u(i,j) = un(i,j) & 
                    + dt * ( &
                    - (un(i,j) * du_dx + vn(i,j) * du_dy) & ! Convection Term
                    + (nu * (d2u_dx2 + d2u_dy2)) & ! Diffusion Term
                    - (1.0d0/rho) * ((p(i + 1,j) - p(i - 1,j)) / (2.0d0 * dx)) & ! Pressure Term
                    )

                v(i,j) = vn(i,j) & 
                    + dt * ( &
                    - (un(i,j) * dv_dx + vn(i,j) * dv_dy) & ! Convection Term
                    + (nu * (d2v_dx2 + d2v_dy2)) & ! Diffusion Term
                    - (1.0d0/rho) * ((p(i,j + 1) - p(i,j - 1)) / (2.0d0 * dy)) & ! Pressure Term
                    )
            end do
        end do
    end subroutine Update_and_Correct_Velocity
    
    subroutine Apply_Vel_BC(u, v, u_inlet,step_x, step_y, nx, ny)
        implicit none    
        integer, intent(in) :: nx, ny, step_x, step_y
        real(8), intent(inout) :: u(nx,ny), v(nx,ny) 
        real(8), intent(in) :: u_inlet
        integer :: j
        
        ! Inlet
        do j = step_y + 1, ny
            u(1,j) = u_inlet
            v(1,j) = 0
        end do

        ! Outlet (zero gradient)
        do j = 1, ny
            u(nx,j) = u(nx-1,j)
            v(nx,j) = v(nx-1,j)
        end do

        ! Top Wall
        do j = 1, nx
            u(j,ny) = 0
            v(j,ny) = 0
        end do

        ! Bottom Wall
        do j = step_x + 1, nx
            u(j,1) = 0
            v(j,1) = 0
        end do

        ! Step wall (Vertical)
        do j = 1, step_y
            u(step_x + 1,j) = 0.0d0
            v(step_x + 1,j) = 0.0d0
        end do

        ! Step wall (Horizontal)
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

        dt_cfl = min(CFL * dx / max(maxval(u), 1e-8), CFL * dy / max(maxval(v), 1e-8))
        dt_diff = min(0.25 * dx ** 2/nu, 0.25 * dy ** 2/nu)

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
                div = abs((u(i+1,j) - u(i-1,j))/(2.0d0*dx) + (v(i,j+1) - v(i,j-1))/(2.0d0*dy))
                max_div = max(max_div, div)
                sum_div = sum_div + div
                count_pts = count_pts + 1
            end do
        end do
        
        avg_div = sum_div / real(count_pts, 8)
        print *, "t=", t, "  Max Div=", max_div, "  Avg Div=", avg_div
    end subroutine Check_Divergence

    subroutine Save_Solution(u, v, p, mask, nx, ny, dx, dy)
        implicit none
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: u(nx,ny), v(nx,ny), p(nx,ny), dx, dy
        integer, intent(in) :: mask(nx,ny)
        integer :: i, j
        real(8) :: x, y

        open(unit = 10, file = 'solution.dat', status = 'replace')
        write(10,*) '# x y u v p mask'

        do j = 1, ny
            do i = 1, nx
                x = dble(i-1)*dx
                y = dble(j-1)*dy
                write(10,'(6E16.8)') x, y, u(i,j), v(i,j), p(i,j), dble(mask(i,j))
            end do
            write(10,*)
        end do
        
        close(10)

    end subroutine Save_Solution
end program DNS
