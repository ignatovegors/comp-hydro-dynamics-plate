PROGRAM  platePrandtlNavierStokes
    IMPLICIT NONE
    LOGICAL(1), EXTERNAL :: ConvergenceCheck
    INTEGER(2), PARAMETER :: io = 12
    INTEGER(4) :: ni, nj, niter
    INTEGER(4) :: i, j, s, s_max
    REAL(8) :: l, h, dx, dy, u_0, nu, eps
    REAL(8), ALLOCATABLE :: x_node(:,:), y_node(:,:)
    REAL(8), ALLOCATABLE :: x_cell(:,:), y_cell(:,:)
    REAL(8), ALLOCATABLE :: u_c(:,:), v_c(:,:), p_c(:,:)
    REAL(8), ALLOCATABLE :: u_n(:,:), v_n(:,:), p_n(:,:)
    REAL(8), ALLOCATABLE :: u_temp(:), v_temp(:), a(:), b(:), c(:), d(:)

    CALL DataInput(io, l, h, ni, nj, u_0, nu, s_max, eps)

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

    CALL MeshMaking(ni, nj, l, h, dx, dy, x_node, y_node, x_cell, y_cell)

    Call InitialConditionsPrandtl(ni, nj, u_0, u_n, p_n)

    CALL BoundaryConditionsPrandtl(ni, nj, u_0, u_n, v_n)


    !****************** Solve equations ********************       
    ALLOCATE(u_temp(nj))
    ALLOCATE(v_temp(nj))
    ALLOCATE(a(nj))
    ALLOCATE(b(nj))
    ALLOCATE(c(nj))
    ALLOCATE(d(nj))
    
    DO i = 2, ni
            
        u_temp = u_n(i - 1, :)
        v_temp = v_n(i - 1, :)

        DO s = 1, s_max

            a(1) = 0D0
            b(1) = 1D0
            c(1) = 0D0
            d(1) = 0D0

            DO j = 2, nj - 1
                a(j) = - v_temp(j - 1) / (2D0 * dy) - nu / dy**2
                b(j) = u_temp(j) / dx + 2D0 * nu / dy**2
                c(j) = v_temp(j + 1) / (2D0 * dy) - nu / dy**2
                d(j) = u_n(i - 1, j)**2D0 / dx
            END DO

            a(nj) = 0D0
            b(nj) = 1D0
            c(nj) = 0D0
            d(nj) = u_0

            CALL ThomasAlgorithm(nj, a, b, c, d, u_n(i, :))

            DO j = 2, nj
                v_n(i,j) = v_n(i, j - 1) - dy / (2 * dx) * (u_n(i,j) - u_n(i - 1, j) + u_n(i, j - 1) - u_n(i - 1, j - 1))
            END DO    
            
            IF ((ConvergenceCheck(u_n(i, :), u_temp, nj, eps)) .AND. (ConvergenceCheck(v_n(i, :), v_temp, nj, eps))) THEN
                WRITE(*,*) 'SOLUTION CONVERGED, NODE №', I 
                EXIT 
            END IF

            u_temp = u_n(i, :)
            v_temp = v_n(i, :)

        END DO

    END DO

    !*******************************************************

    ! Navier - Stokes

    !*******************************************************

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


SUBROUTINE DataInput(io, l, h, ni, nj, u_0, nu, s_max, eps)
    ! Takes input data from file input.txt
    IMPLICIT NONE
    INTEGER(2) :: io
    INTEGER(4) :: ni, nj, s_max
    REAL(8) :: l, h, u_0, nu, eps
    INTENT(IN) io
    INTENT(OUT) l, h, ni, nj, u_0, nu, s_max, eps

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
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE


SUBROUTINE MeshMaking(ni, nj, l, h, dx, dy, x_node, y_node, x_cell, y_cell)
    ! Makes mesh for numerical solution Prandtl (node) and 
    ! Navier-Stokes (cell) systems of equations
    IMPLICIT NONE
    INTEGER(4) :: ni, nj, i, j
    REAL(8) :: l, h, dx, dy
    REAL(8), DIMENSION(ni,nj) :: x_node, y_node
    REAL(8), DIMENSION(0:ni, 0:nj) :: x_cell, y_cell
    INTENT(IN) l, h, ni, nj
    INTENT(OUT) dx, dy, x_node, y_node, x_cell, y_cell

    WRITE(*,*) 'MESH MAKING'

    dx = l / (ni - 1)
    dy = h / (nj - 1)

    DO i = 1, ni
        DO j = 1, nj
            x_node(i,j) = (i - 1) * dx
            y_node(i,j) = (j - 1) * dy
        END DO
    END DO

    x_cell(0, 1:nj) = - dx / 2D0
    y_cell(0, 1:nj) = y_node(1, 1:nj) + dy / 2D0
    x_cell(1:ni, 0) = x_node(1:ni, 1) + dx / 2D0
    y_cell(1:ni, 0) = - dy / 2D0

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


SUBROUTINE BoundaryConditionsPrandtl(ni, nj, u_0, u, v)
    ! Boundary no-slip and uniform conditions for velocity
    IMPLICIT NONE
    INTEGER(4) :: ni, nj
    REAL(8) :: u_0
    REAL(8), DIMENSION(ni,nj) :: u, v
    INTENT(IN) ni, nj, u_0
    INTENT(OUT) u, v

    WRITE(*,*) 'BOUNDARY CONDITIONS APPLYING (PRANDTL)'
    
    u(1:ni, 1) = 0D0
    v(1:ni, 1) = 0D0
    u(1:ni, nj) = u_0

    WRITE(*,*) 'SUCCESS'
    
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


SUBROUTINE OutputFieldsCell(io, ni, nj, x, y, u, v, p)
    ! Cells-based results output
    IMPLICIT NONE
    INTEGER(2) :: io
    INTEGER(4) :: ni, nj
    REAL(8), DIMENSION(ni,nj) :: x, y
    REAL(8), DIMENSION(0:ni, 0:nj) :: u, v, p
    INTENT(IN) io, ni, nj, x, y, u, v, p
        
    WRITE(*,*) 'RESULTS OUTPUT (NAVIER-STOKES)' 
    OPEN(io, FILE='RES_NS.PLT')
    WRITE(io,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
    WRITE(io,*) 'ZONE I=', ni, ', J=', nj, ', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
    WRITE(io,'(100E25.16)') x(1:ni, 1:nj) 
    WRITE(io,'(100E25.16)') y(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') u(1:ni - 1, 1:nj - 1)
    WRITE(io,'(100E25.16)') v(1:ni - 1, 1:nj - 1)
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
    WRITE(io,'(100E25.16)') y(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') u(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') v(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') p(1:ni, 1:nj)
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    END  SUBROUTINE 


LOGICAL(1) FUNCTION ConvergenceCheck(a, b, n, eps)
    IMPLICIT NONE
    REAL(8), DIMENSION(n) :: a, b, dif
    REAL(8) :: eps
    INTEGER(4) :: i, n

    DO i = 1, n
        dif(i) = abs(a(i) - b(i)) / abs(a(i))
    END DO

    ConvergenceCheck = (MAXVAL(dif) < eps)

    END FUNCTION
