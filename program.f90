PROGRAM platePrandtlNavierStokes
    IMPLICIT NONE
    INTEGER(2), PARAMETER :: io = 12
    INTEGER(4) :: ni, nj
    INTEGER(4) :: i, j, s_max
    REAL(8) :: l, h, dx, dy, u_0, nu, eps, dt, a, cfl
    REAL(8), ALLOCATABLE :: x_node(:,:), y_node(:,:)
    REAL(8), ALLOCATABLE :: x_cell(:,:), y_cell(:,:)
    REAL(8), ALLOCATABLE :: u_c(:,:), v_c(:,:), p_c(:,:)
    REAL(8), ALLOCATABLE :: u_n(:,:), v_n(:,:), p_n(:,:)


    CALL DataInput(io, l, h, ni, nj, u_0, nu, s_max, eps, cfl)

    ALLOCATE(x_node(ni,nj))
    ALLOCATE(y_node(ni,nj))
    ALLOCATE(x_cell(0:ni, 0:nj))
    ALLOCATE(y_cell(0:ni, 0:nj))

    ALLOCATE(u_c(0:ni, 0:nj))
    ALLOCATE(v_c(0:ni, 0:nj))
    ALLOCATE(p_c(0:ni, 0:nj))

    ALLOCATE(u_n(ni,nj))   
    ALLOCATE(v_n(ni,nj))   
    ALLOCATE(p_n(ni,nj))

    CALL MeshMaking(ni, nj, l, h, dx, dy, x_node, y_node, x_cell, y_cell, u_0, a, cfl, dt, nu)

    CALL InitialConditionsPrandtl(ni, nj, u_0, u_n, p_n)

    CALL InitialConditionsNavierStokes(ni, nj, u_0, u_c, v_c, p_c)

    CALL BoundaryConditionsPrandtl(ni, nj, u_0, u_n, v_n)
    
    CALL BoundaryConditionsPrandtl(ni, nj, u_0, u_n, v_n)    

    CALL SolverPrandtl(ni, nj, s_max, dx, dy, nu, eps, u_0, u_n, v_n)

    CALL SolverNavierStokes(ni, nj, s_max, dx, dy, nu, eps, u_0, u_c, v_c, p_c, a, dt)

    CALL OutputFieldsNode(io, ni, nj, x_node, y_node, u_n, v_n, p_n)

    CALL OutputFieldsCell(io, ni, nj, x_cell, y_cell, u_c, v_c, p_c)

    DEALLOCATE(x_node)
    DEALLOCATE(y_node)
    DEALLOCATE(x_cell)
    DEALLOCATE(y_cell)

    DEALLOCATE(u_c)
    DEALLOCATE(v_c)
    DEALLOCATE(p_c)
    
    DEALLOCATE(u_n)   
    DEALLOCATE(v_n)   
    DEALLOCATE(p_n)

END PROGRAM


SUBROUTINE DataInput(io, l, h, ni, nj, u_0, nu, s_max, eps, cfl)
    ! Takes input data from file input.txt
    IMPLICIT NONE
    INTEGER(2) :: io
    INTEGER(4) :: ni, nj, s_max
    REAL(8) :: l, h, u_0, nu, eps, cfl
    INTENT(IN) io
    INTENT(OUT) l, h, ni, nj, u_0, nu, s_max, eps, cfl

    WRITE(*,*) 'READING INPUT FILE'
    OPEN(io,FILE='INPUT.TXT')
    READ(io,*) l
    READ(io,*) h
    READ(io,*) ni
    READ(io,*) nj
    READ(io,*) u_0
    READ(io,*) nu
    READ(io,*) s_max
    READ(io,*) eps
    READ(io,*) cfl

    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE


SUBROUTINE MeshMaking(ni, nj, l, h, dx, dy, x_node, y_node, x_cell, y_cell, u_0, a, cfl, dt, nu)
    ! Makes mesh for numerical solution Prandtl (node) and 
    ! Navier-Stokes (cell) systems of equations
    IMPLICIT NONE
    INTEGER(4) :: ni, nj, i, j
    REAL(8) :: l, h, dx, dy, u_0, a, cfl, dt, nu
    REAL(8), DIMENSION(ni,nj) :: x_node, y_node
    REAL(8), DIMENSION(0:ni, 0:nj) :: x_cell, y_cell
    INTENT(IN) l, h, ni, nj, u_0, cfl, nu
    INTENT(OUT) dx, dy, x_node, y_node, x_cell, y_cell, a, dt

    WRITE(*,*) 'MESH MAKING'

    dx = l / (ni - 1)
    dy = h / (nj - 1)
    dt = cfl * MIN(0.5D0 * dx * dx / nu, 0.5D0 * dy * dy / nu, dx / u_0)
    a = 1D0 / (u_0 * u_0)

    DO i = 1, ni
        DO j = 1, nj
            x_node(i,j) = (i - 1) * dx
            y_node(i,j) = (j - 1) * dy
        END DO
    END DO

    DO i = 0, ni
        DO j = 0, nj
            x_cell(i,j) = (i - 5D-1) * dx
            y_cell(i,j) = (j - 5D-1) * dy
        END DO
    END DO

    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE


SUBROUTINE InitialConditionsPrandtl(ni, nj, u_0, u, p)
    ! Initial uniform velocity condition in the inlet
    IMPLICIT NONE
    INTEGER(4) :: ni, nj
    REAL(8) :: u_0
    REAL(8), DIMENSION(ni,nj) :: u, p
    INTENT(IN) ni, nj, u_0
    INTENT(OUT) u, p

    WRITE(*,*) 'INITIAL CONDITIONS APPLYING (PRANDTL)'
    
    p(1:ni, 1:nj) = 0D0

    u(1, 1:nj) = u_0
    

    WRITE(*,*) 'SUCCESS'
    
    END SUBROUTINE


SUBROUTINE InitialConditionsNavierStokes(ni, nj, u_0, u, v, p)
    ! Initial velocities and pressure conditions in area
    IMPLICIT NONE
    INTEGER(4) :: ni, nj
    REAL(8) :: u_0
    REAL(8), DIMENSION(0:ni,0:nj) :: u, v, p
    INTENT(IN) ni, nj, u_0
    INTENT(OUT) u, v, p

    WRITE(*,*) 'INITIAL CONDITIONS APPLYING (NAVIER-STOKES)'
    
    u(0:ni, 0:nj) = u_0

    v(0:ni, 0:nj) = 0D0       

    p(0:ni, 0:nj) = 0D0   

    WRITE(*,*) 'SUCCESS'
    
    END SUBROUTINE


SUBROUTINE BoundaryConditionsPrandtl(ni, nj, u_0, u, v)
    ! Boundary no-slip and uniform conditions for velocity
    IMPLICIT NONE
    INTEGER(4) :: ni, nj
    REAL(8) :: u_0
    REAL(8), DIMENSION(ni,nj) :: u, v
    INTENT(IN) ni, nj, u_0
    INTENT(OUT) u, v
   
    u(1:ni, 1) = 0D0
    v(1:ni, 1) = 0D0
    u(1:ni, nj) = u_0
  
    END SUBROUTINE


SUBROUTINE BoundaryConditionsNavierStokes(ni, nj, u_0, u, v, p)
    ! Boundary conditions for velocities and pressure
    IMPLICIT NONE
    INTEGER(4) :: ni, nj
    REAL(8) :: u_0
    REAL(8), DIMENSION(0:ni,0:nj) :: u, v, p
    INTENT(IN) ni, nj, u_0
    INTENT(OUT) u, v, p
    
    u(0, 1:nj) = u_0
    v(0, 1:nj) = 0D0
    p(0, 1:nj) = p(1, 1:nj)

    u(ni, 1:nj) = u(ni - 1, 1:nj)
    v(ni, 1:nj) = v(ni - 1, 1:nj)
    p(ni, 1:nj) = 0D0

    u(1:ni, 0) = - u(1:ni, 1)
    v(1:ni, 0) = - v(1:ni, 1)
    p(1:ni, 0) = p(1:ni, 1)

    u(1:ni, nj) = u(1:ni, nj - 1)
    v(1:ni, nj) = v(1:ni, nj - 1)
    p(1:ni, nj) = 0D0
   
    END SUBROUTINE


SUBROUTINE ThomasAlgorithm(k_max, a, b, c, d, res)
    ! Solution of tridiagonal system 
    ! a_{k} * x_{k - 1} + b_{k} * x_{k} + c_{k} * x_{k + 1} = d_{k}
    ! with k_max unknowns (a_{1} = 0, c_{k_max} = 0)
    IMPLICIT NONE
    INTEGER(4) :: k, k_max
    REAL(8), DIMENSION(k_max) :: a, b, c, d, alpha, beta, res
    INTENT(IN) k_max, a, b, c, d
    INTENT(OUT) res
   
    alpha(2) = - c(1) / b(1)
    beta(2) = d(1) / b(1)

    DO k = 2, k_max - 1
        alpha(k + 1) = - c(k) / (b(k) + a(k) * alpha(k))
        beta(k + 1) = (d(k) - a(k) * beta(k)) / (b(k) + a(k) * alpha(k))
    END DO

    res(k_max) = (d(k) - a(k) * beta(k)) / (b(k) + a(k) * alpha(k))

    DO k = k_max - 1, 1, -1
        res(k) = alpha(k + 1) * res(k + 1) + beta(k + 1)
    END DO

    END SUBROUTINE


SUBROUTINE SolverPrandtl(ni, nj, s_max, dx, dy, nu, eps, u_0, u, v)
    ! Solver for Prandtl (Simplified Navier-Stokes) equations system
    IMPLICIT NONE
    LOGICAL(1), EXTERNAL :: ConvergenceCheckPrandtl
    INTEGER(4) :: i, j, s, ni, nj, s_max
    REAL(8) :: dx, dy, nu, eps, u_0
    REAL(8), DIMENSION(ni,nj) :: u, v
    REAL(8), DIMENSION(nj) :: u_temp, v_temp, a, b, c, d

    WRITE(*,*) 'SOLVING EQUATIONS (PRANDTL)'

    DO i = 2, ni
            
        u_temp = u(i - 1, :)
        v_temp = v(i - 1, :)

        DO s = 1, s_max

            a(1) = 0D0
            b(1) = 1D0
            c(1) = 0D0
            d(1) = 0D0

            DO j = 2, nj - 1
                a(j) = - v_temp(j - 1) / (2D0 * dy) - nu / dy**2
                b(j) = u_temp(j) / dx + 2D0 * nu / dy**2
                c(j) = v_temp(j + 1) / (2D0 * dy) - nu / dy**2
                d(j) = u(i - 1, j)**2D0 / dx
            END DO

            a(nj) = 0D0
            b(nj) = 1D0
            c(nj) = 0D0
            d(nj) = u_0

            CALL ThomasAlgorithm(nj, a, b, c, d, u(i, :))

            DO j = 2, nj
                v(i,j) = v(i, j - 1) - dy / (2 * dx) * (u(i,j) - u(i - 1, j) + u(i, j - 1) - u(i - 1, j - 1))
            END DO    
            
            IF ((ConvergenceCheckPrandtl(u(i, :), u_temp, nj, eps)) .AND. (ConvergenceCheckPrandtl(v(i, :), v_temp, nj, eps))) THEN
                WRITE(*,*) 'SOLUTION CONVERGED BY RESIDUALS, NODE №', I, ', s = ', s
                EXIT 
            END IF

            u_temp = u(i, :)
            v_temp = v(i, :)

            IF (s == s_max) THEN
                WRITE(*,*) 'SOLUTION CONVERGED BY ITERATIONS BOUNDARY, NODE №', I
            END IF

        END DO

    END DO

    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE


SUBROUTINE SolverNavierStokes(ni, nj, s_max, dx, dy, nu, eps, u_0, u, v, p, a, dt)
    !Solver for Navier-Stokes system of equations
    IMPLICIT NONE
    REAL(8), EXTERNAL :: HalfIndexValue, ResidualNavierStokes
    LOGICAL(1), EXTERNAL :: ConvergenceCheckNavierStokes
    INTEGER(4) :: i, j, ni, nj, s, s_max
    REAL(8) :: u_hat_left, u_hat_right
    REAL(8) :: v_hat_top, v_hat_bot
    REAL(8) :: u_left, u_right, u_top, u_bot
    REAL(8) :: v_left, v_right, v_top, v_bot
    REAL(8) :: p_left, p_right, p_top, p_bot
    REAL(8) :: dx, dy, nu, eps, u_0, a, dt, res_u, res_v, res_p 
    REAL(8), DIMENSION(0:ni,0:nj) :: u_old, v_old, p_old, u, v, p
    INTENT(IN) ni, nj, s_max, dx, dy, nu, eps, u_0, a, dt
    INTENT(OUT) u, v, p

    WRITE(*,*) 'SOLVING EQUATIONS (NAVIER-STOKES)'

    DO s = 1, s_max

        u_old = u
        v_old = v
        p_old = p

        DO i = 1, ni - 1
            DO j = 1, nj - 1

                u_hat_left = 5D-1 * (u(i - 1, j) + u(i,j))
                u_hat_right = 5D-1 * (u(i, j) + u(i + 1,j))
                
                v_hat_top = 5D-1 * (v(i, j) + v(i,j + 1))
                v_hat_bot = 5D-1 * (v(i, j - 1) + v(i,j))

                u_left = HalfIndexValue(u_hat_left, u(i, j), u(i - 1,j))
                v_left = HalfIndexValue(u_hat_left, v(i, j), v(i - 1,j))
                p_left = HalfIndexValue(u_hat_left, p(i - 1, j), p(i,j))
                
                u_right = HalfIndexValue(u_hat_right, u(i + 1, j), u(i,j))
                v_right = HalfIndexValue(u_hat_right, v(i + 1, j), v(i,j))
                p_right = HalfIndexValue(u_hat_right, p(i, j), p(i + 1,j))
                
                u_top = HalfIndexValue(v_hat_top, u(i, j + 1), u(i,j))
                v_top = HalfIndexValue(v_hat_top, v(i, j + 1), v(i,j))
                p_top = HalfIndexValue(v_hat_top, p(i, j), p(i,j + 1))
                
                u_bot = HalfIndexValue(v_hat_bot, u(i, j), u(i,j - 1))
                v_bot = HalfIndexValue(v_hat_bot, v(i, j), v(i,j - 1))
                p_bot = HalfIndexValue(v_hat_bot, p(i, j - 1), p(i,j))

                p(i,j) = p(i,j) - dt * ((u_right - u_left) / dx + (v_top - v_bot) / dy) / a
                u(i,j) = u(i,j) - dt * ((u_hat_right * u_right - u_hat_left * u_left) / dx &
                    + (v_hat_top * u_top - v_hat_bot * u_bot) / dy &
                    + (p_right - p_left) / dx &
                    - nu * (u(i + 1, j) - 2 * u(i,j) + u(i - 1, j)) / dx ** 2 &
                    - nu * (u(i, j + 1) - 2 * u(i,j) + u(i, j - 1)) / dy ** 2 )
                v(i,j) = v(i,j) - dt * ((u_hat_right * v_right - u_hat_left * v_left) / dx &
                    + (v_hat_top * v_top - v_hat_bot * v_bot) / dy &
                    + (p_top - p_bot) / dy &
                    - nu * (v(i + 1, j) - 2 * v(i,j) + v(i - 1, j)) / dx ** 2 &
                    - nu * (v(i, j + 1) - 2 * v(i,j) + v(i, j - 1)) / dy ** 2 )

            END DO
        END DO

        CALL BoundaryConditionsNavierStokes(ni, nj, u_0, u, v, p)

        res_u = ResidualNavierStokes(u, u_old, ni, nj, dt)
        res_v = ResidualNavierStokes(v, v_old, ni, nj, dt)
        res_p = ResidualNavierStokes(p, p_old, ni, nj, dt)

        ! CALL ResidualLogs(1, s, res_u, res_v, res_p)

        IF (ConvergenceCheckNavierStokes(u, u_old, ni, nj, eps, dt) &
            .OR. ConvergenceCheckNavierStokes(v, v_old, ni, nj, eps, dt) &
            .OR. ConvergenceCheckNavierStokes(p, p_old, ni, nj, eps, dt)) THEN
                WRITE(*,*) 'SOLUTION CONVERGED BY RESIDUALS'
                EXIT
        END IF

        IF (MOD(s,100) == 0) THEN
            WRITE(*,*) s, 'ITERATIONS MADE'
        END IF

        IF (s == s_max) THEN
            WRITE(*,*) 'SOLUTION CONVERGED BY ITERATIONS BOUNDARY'
        END IF

    END DO

    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE


SUBROUTINE OutputFieldsCell(io, ni, nj, x, y, u, v, p)
    ! Cells-based results output
    IMPLICIT NONE
    INTEGER(2) :: io
    INTEGER(4) :: ni, nj
    REAL(8), DIMENSION(0:ni, 0:nj) :: x, y
    REAL(8), DIMENSION(0:ni, 0:nj) :: u, v, p
    INTENT(IN) io, ni, nj, x, y, u, v, p
        
    WRITE(*,*) 'RESULTS OUTPUT (NAVIER-STOKES)' 
    OPEN(io, FILE='RES_NS.PLT')
    WRITE(io,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
    WRITE(io,*) 'ZONE I=', ni, ', J=', nj, ', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
    WRITE(io,'(100E25.16)') x(1:ni, 1:nj) 
    WRITE(io,*) 
    WRITE(io,'(100E25.16)') y(1:ni, 1:nj)
    WRITE(io,*) 
    WRITE(io,'(100E25.16)') u(1:ni - 1, 1:nj - 1)
    WRITE(io,*) 
    WRITE(io,'(100E25.16)') v(1:ni - 1, 1:nj - 1)
    WRITE(io,*) 
    WRITE(io,'(100E25.16)') p(1:ni - 1, 1:nj - 1)
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE 


SUBROUTINE OutputFieldsNode(io, ni, nj, x, y, u, v, p)
    ! Nodes-based results output
    IMPLICIT NONE
    INTEGER(2) :: io
    INTEGER(4) :: ni, nj
    REAL(8), DIMENSION(ni,nj) :: x, y
    REAL(8), DIMENSION(ni,nj) :: u, v, p
    INTENT(IN) io, ni, nj, x, y, u, v, p
    
    WRITE(*,*) 'RESULTS OUTPUT (PRANDTL) ' 
    OPEN(io,FILE='RES_PR.PLT')
    WRITE(io,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
    WRITE(io,*) 'ZONE I=', ni, ', J=', nj, ', DATAPACKING=BLOCK'
    WRITE(io,'(100E25.16)') x(1:ni, 1:nj)
    WRITE(io,*) 
    WRITE(io,'(100E25.16)') y(1:ni, 1:nj)
    WRITE(io,*)
    WRITE(io,'(100E25.16)') u(1:ni, 1:nj)
    WRITE(io,*)
    WRITE(io,'(100E25.16)') v(1:ni, 1:nj)
    WRITE(io,*)
    WRITE(io,'(100E25.16)') p(1:ni, 1:nj)
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE 


SUBROUTINE ResidualLogs(io, s, res_u, res_v, res_p)
    ! Nodes-based results output
    IMPLICIT NONE
    INTEGER(2) :: io
    REAL(8) :: s, res_u, res_v, res_p
    INTENT(IN) io, res_u, res_v, res_p, s
     
    OPEN(io,FILE='RESIDUALS.PLT')
    WRITE(io,*) 'VARIABLES = "ITER", "RES_U", "RES_V", "RES_P"' 
    WRITE(io,'(100E25.16)') s, res_u, res_v, res_p
    CLOSE(io)

    END SUBROUTINE


LOGICAL(1) FUNCTION ConvergenceCheckPrandtl(a, b, n, eps)
    IMPLICIT NONE
    REAL(8), EXTERNAL :: ResidualPrandtl
    REAL(8), DIMENSION(n) :: a, b
    REAL(8) :: eps
    INTEGER(4) :: n

    ConvergenceCheckPrandtl = (ResidualPrandtl(a, b, n) < eps)

    END FUNCTION


LOGICAL(1) FUNCTION ConvergenceCheckNavierStokes(a, b, ni, nj, eps, dt)
    IMPLICIT NONE
    REAL(8), EXTERNAL :: ResidualNavierStokes
    REAL(8), DIMENSION(ni, nj) :: a, b, dif
    REAL(8) :: eps, dt
    INTEGER(4) :: ni, nj

    ConvergenceCheckNavierStokes = (ResidualNavierStokes(a, b, ni, nj, dt) < eps)

    END FUNCTION


REAL(8) FUNCTION HalfIndexValue(arg, minus_res, plus_res)
    IMPLICIT NONE
    REAL(8) :: arg, minus_res, plus_res

    IF (arg < 0D0) THEN
        HalfIndexValue = minus_res
    ELSE
        HalfIndexValue = plus_res
    END IF

    END FUNCTION


REAL(8) FUNCTION ResidualPrandtl(a, b, n)
    IMPLICIT NONE
    REAL(8), DIMENSION(n) :: a, b, dif
    INTEGER(4) :: i, n

    DO i = 1, n
        dif(i) = abs(a(i) - b(i)) / abs(a(i))
    END DO

    ResidualPrandtl = MAXVAL(dif)

    END FUNCTION


REAL(8) FUNCTION ResidualNavierStokes(a, b, ni, nj, dt)
    IMPLICIT NONE
    REAL(8), DIMENSION(ni, nj) :: a, b, dif
    REAL(8) :: dt
    INTEGER(4) :: i, j, ni, nj

    DO i = 1, ni
        DO j = 1, nj
            dif(i,j) = abs(a(i,j) - b(i,j)) / abs(a(i,j))
        END DO
    END DO

    ResidualNavierStokes = MAXVAL(dif) / dt 

    END FUNCTION
