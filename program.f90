PROGRAM  platePrandtl
IMPLICIT NONE
INTEGER(2), PARAMETER :: io = 12
INTEGER(4) :: ni, nj, niter
INTEGER(4) :: i, j
REAL(8) :: l, h, dx, dy
REAL(8), ALLOCATABLE :: x_node(:,:), y_node(:,:)
REAL(8), ALLOCATABLE :: x_cell(:,:), y_cell(:,:)
REAL(8), ALLOCATABLE :: u_c(:,:), v_c(:,:), p_c(:,:)
REAL(8), ALLOCATABLE :: u_n(:,:), v_n(:,:), p_n(:,:)

    WRITE(*,*) 'READING INPUT FILE'
    OPEN(io,FILE='INPUT.TXT')
    READ(io,*) l
    READ(io,*) h
    READ(io,*) ni
    READ(io,*) nj
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

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

    DO i = 1, ni
        DO j = 1, nj
            x_cell(i,j) = x_node(i,j) + dx / 2
            y_cell(i,j) = y_node(i,j) + dy / 2
        END DO
    END DO

    DO i = 0, ni
        DO j = 0, nj
            u_c(i,j) = x_cell(i,j)
            v_c(i,j) = y_cell(i,j)
            p_c(i,j) = x_cell(i,j) + y_cell(i,j)
        END DO
    END DO

    DO i = 1, ni
        DO j = 1, nj
            u_n(i,j) = x_node(i,j)
            v_n(i,j) = y_node(i,j)
            p_n(i,j) = x_node(i,j) + y_node(i,j)       
        END DO
    END DO


    !****************** Solve equations ********************       

    ! Prandtl

    ! Navier - Stokes

    !*******************************************************

    WRITE(*,*) 'RESULTS OUTPUT (PRANDTL)' 
    OPEN(io, FILE='RES_PR.PLT')
    Call OutputFieldsNode(io, ni, nj, x_node, y_node, u_n, v_n, p_n)
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    WRITE(*,*) 'RESULTS OUTPUT (NAVIER-STOKES) ' 
    OPEN(io,FILE='RES_NS.PLT')
    Call OutputFieldsCell(io, ni, nj, x_node, y_node, u_c, v_c, p_c)
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

END PROGRAM


SUBROUTINE OutputFieldsCell(io, ni, nj, x, y, u, v, p)
IMPLICIT NONE
INTEGER(2) :: io
INTEGER(4) :: ni, nj
REAL(8), DIMENSION(ni,nj) :: x, y
REAL(8), DIMENSioN(0:ni, 0:nj) :: u, v, p

    WRITE(io,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
    WRITE(io,*) 'ZONE I=',ni,', J=',nj,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
    WRITE(io,'(100E25.16)') x(1:ni, 1:nj) 
    WRITE(io,'(100E25.16)') y(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') u(1:ni-1, 1:nj-1)
    WRITE(io,'(100E25.16)') v(1:ni-1, 1:nj-1)
    WRITE(io,'(100E25.16)') p(1:ni-1, 1:nj-1)

END SUBROUTINE 


SUBROUTINE OutputFieldsNode(io, ni, nj, x, y, u, v, p)
IMPLICIT NONE
INTEGER(2) :: io
INTEGER(4) :: ni, nj
REAL(8), DIMENSION(ni,nj) :: x,y
REAL(8), DIMENSION(ni,nj) :: u,v,p

    WRITE(io,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
    WRITE(io,*) 'ZONE I=', ni, ', J=', nj, ', DATAPACKING=BLOCK'
    WRITE(io,'(100E25.16)') x(1:ni, 1:nj) 
    WRITE(io,'(100E25.16)') y(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') u(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') v(1:ni, 1:nj)
    WRITE(io,'(100E25.16)') p(1:ni, 1:nj)

END  SUBROUTINE 
