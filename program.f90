PROGRAM  platePrandtl
IMPLICIT NONE
INTEGER(2), PARAMETER :: io = 12
INTEGER(4) :: ni, nj, niter
INTEGER(4) :: i, j
REAL(8) :: l, h, dx, dy, u_0
REAL(8), ALLOCATABLE :: x_node(:,:), y_node(:,:)
REAL(8), ALLOCATABLE :: x_cell(:,:), y_cell(:,:)
REAL(8), ALLOCATABLE :: u_c(:,:), v_c(:,:), p_c(:,:)
REAL(8), ALLOCATABLE :: u_n(:,:), v_n(:,:), p_n(:,:)

    CALL DataInput(io, l, h, ni, nj, u_0)

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

    u_n = 5
    v_n = 7
    p_n = 3

    CALL BoundaryConditionsPrandtl(ni, nj, u_0, u_n, v_n, p_n)


    !****************** Solve equations ********************       

    ! Prandtl

    ! Navier - Stokes

    !*******************************************************

    WRITE(*,*) 'RESULTS OUTPUT (PRANDTL)' 
    OPEN(io, FILE='RES_PR.PLT')
    CALL OutputFieldsNode(io, ni, nj, x_node, y_node, u_n, v_n, p_n)
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    WRITE(*,*) 'RESULTS OUTPUT (NAVIER-STOKES) ' 
    OPEN(io,FILE='RES_NS.PLT')
    CALL OutputFieldsCell(io, ni, nj, x_node, y_node, u_c, v_c, p_c)
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

END PROGRAM


SUBROUTINE DataInput(io, l, h, ni, nj, u_0)
IMPLICIT NONE
INTEGER(2) :: io
INTEGER(4) :: ni, nj
REAL(8) :: l, h, u_0
INTENT(IN) io
INTENT(OUT) l, h, ni, nj, u_0

    WRITE(*,*) 'READING INPUT FILE'
    OPEN(io,FILE='INPUT.TXT')
    READ(io,*) l
    READ(io,*) h
    READ(io,*) ni
    READ(io,*) nj
    READ(io,*) u_0
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

END SUBROUTINE


SUBROUTINE MeshMaking(ni, nj, l, h, dx, dy, x_node, y_node, x_cell, y_cell)
IMPLICIT NONE
INTEGER(4) :: ni, nj, i, j
REAL(8) :: l, h, dx, dy
REAL(8), DIMENSION(ni,nj) :: x_node, y_node
REAL(8), DIMENSION(0:ni, 0:nj) :: x_cell, y_cell
INTENT(IN) l, h, ni, nj
INTENT(OUT) dx, dy, x_node, y_node, x_cell, y_cell

    dx = l / (ni - 1)
    dy = h / (nj - 1)

    DO i = 1, ni
        DO j = 1, nj
            x_node(i,j) = (i - 1) * dx
            y_node(i,j) = (j - 1) * dy
        END DO
    END DO

    x_cell(0, 1:nj) = - dx / 2
    y_cell(0, 1:nj) = y_node(1, 1:nj) + dy / 2
    x_cell(1:ni, 0) = x_node(1:ni, 1) + dx / 2
    y_cell(1:ni, 0) = - dy / 2

END SUBROUTINE


SUBROUTINE BoundaryConditionsPrandtl(ni, nj, u_0, u, v, p)
IMPLICIT NONE
INTEGER(4) :: ni, nj
REAL(8) :: u_0
REAL(8), DIMENSION(ni,nj) :: u, v, p
INTENT(IN) ni, nj, u_0
INTENT(OUT) u, v, p
    
    u(1:ni, 1) = 0
    v(1:ni, 1) = 0
    u(1:ni, nj) = u_0
    p(1:ni, 1:nj) = 0
    
END SUBROUTINE


SUBROUTINE OutputFieldsCell(io, ni, nj, x, y, u, v, p)
IMPLICIT NONE
INTEGER(2) :: io
INTEGER(4) :: ni, nj
REAL(8), DIMENSION(ni,nj) :: x, y
REAL(8), DIMENSION(0:ni, 0:nj) :: u, v, p
INTENT(IN) io, ni, nj, x, y, u, v, p

    WRITE(io,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
    WRITE(io,*) 'ZONE I=',ni,', J=',nj,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
    WRITE(io,'(100E25.16)') x(1:ni, 1:nj) 
    WRITE(io,'(100E25.16)') y(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') u(1:ni - 1, 1:nj - 1)
    WRITE(io,'(100E25.16)') v(1:ni - 1, 1:nj - 1)
    WRITE(io,'(100E25.16)') p(1:ni - 1, 1:nj - 1)

END SUBROUTINE 


SUBROUTINE OutputFieldsNode(io, ni, nj, x, y, u, v, p)
IMPLICIT NONE
INTEGER(2) :: io
INTEGER(4) :: ni, nj
REAL(8), DIMENSION(ni,nj) :: x, y
REAL(8), DIMENSION(ni,nj) :: u, v, p
INTENT(IN) io, ni, nj, x, y, u, v, p

    WRITE(io,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
    WRITE(io,*) 'ZONE I=', ni, ', J=', nj, ', DATAPACKING=BLOCK'
    WRITE(io,'(100E25.16)') x(1:ni, 1:nj) 
    WRITE(io,'(100E25.16)') y(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') u(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') v(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') p(1:ni, 1:nj)

END  SUBROUTINE 
