        PROGRAM  Pr
        IMPLICIT NONE
 
        INTEGER, PARAMETER:: IO = 12 ! input-output unit
        INTEGER NI, NJ, NITER
        INTEGER I,J
        REAL L,H,dx,dy
        REAL,ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
        REAL,ALLOCATABLE :: X_Cell(:,:),Y_Cell(:,:)
        REAL,ALLOCATABLE :: U_c(:,:),V_c(:,:),P_c(:,:)
        REAL,ALLOCATABLE :: U_n(:,:),V_n(:,:),P_n(:,:)

        write(*,*) 'Read input file'
        open(IO,FILE='Input.txt')
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        close(IO)
   
        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates
        allocate(X_Cell(0:NI,0:NJ)) ! Cell Centers
        allocate(Y_Cell(0:NI,0:NJ)) ! Cell Centers

!*******************  Cell-centered variables **********
        allocate(U_c(0:NI,0:NJ))   ! Velocity U
        allocate(V_c(0:NI,0:NJ))   ! Velocity V
        allocate(P_c(0:NI,0:NJ))   ! Pressure

!*******************  Node variables ******************
        allocate(U_n(NI,NJ))   ! Velocity U
        allocate(V_n(NI,NJ))   ! Velocity V
        allocate(P_n(NI,NJ))   ! Pressure

        dx=L/(NI-1)
        dy=H/(NJ-1)

        do I=1,NI
          do J=1,NJ
            X_Node(I,J)=(I-1)*dx
            Y_Node(I,J)=(J-1)*dy
          end do
        end do

        X_Cell(0,1:NJ)=-dx/2
        Y_Cell(0,1:NJ)=Y_Node(1,1:NJ)+dy/2
        X_Cell(1:NI,0)=X_Node(1:NI,1)+dx/2
        Y_Cell(1:NI,0)=-dy/2
        do I=1,NI
          do J=1,NJ
            X_Cell(I,J)=X_Node(I,J)+dx/2
            Y_Cell(I,J)=Y_Node(I,J)+dy/2
          end do
        end do

!************************* INITIAL FIELD ********************* 

       do I=0,NI
         do J=0,NJ
           U_c(I,J)=X_Cell(I,J)
           V_c(I,J)=Y_Cell(I,J)
           P_c(I,J)=X_Cell(I,J)+Y_Cell(I,J)
          end do
       end do

       do I=1,NI
         do J=1,NJ
           U_n(I,J)=X_Node(I,J)
           V_n(I,J)=Y_Node(I,J)
           P_n(I,J)=X_Node(I,J)+Y_Node(I,J)       
         end do
       end do


!****************** Solve equations ********************       

! Prandtl

! Navier - Stokes


!****************** Output Results ********************           
  
       write(*,*) 'Output data node (Prandtl)' 
       Open(IO,FILE='data_pr.plt')
       Call OutputFields_Node(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P_n)
       Close(IO)

       write(*,*) 'Output data cell (Navier - Stokes) ' 
       Open(IO,FILE='data_ns.plt')
       Call OutputFields_Cell(IO,NI,NJ,X_Node,Y_Node,U_c,V_c,P_c)
       Close(IO)

       End program


       SUBROUTINE OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
         IMPLICIT NONE

         INTEGER NI,NJ,IO
         REAL, DIMENSION(NI,NJ):: X,Y
         REAL, DIMENSION(0:NI,0:NJ)::U,V,P
       
         Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
         Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
         Write(IO,'(100E25.16)') X(1:NI,1:NJ) 
         Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
         Write(IO,'(100E25.16)') U(1:NI-1,1:NJ-1)
         Write(IO,'(100E25.16)') V(1:NI-1,1:NJ-1)
         Write(IO,'(100E25.16)') P(1:NI-1,1:NJ-1)

       END SUBROUTINE 


       SUBROUTINE OutputFields_Node(IO,NI,NJ,X,Y,U,V,P)
         IMPLICIT NONE
 
         INTEGER NI,NJ,IO
         REAL, DIMENSION(NI,NJ):: X,Y
         REAL, DIMENSION(NI,NJ):: U,V,P

         Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
         Write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
         Write(IO,'(100E25.16)') X(1:NI,1:NJ) 
         Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
         Write(IO,'(100E25.16)') U(1:NI,1:NJ)
         Write(IO,'(100E25.16)') V(1:NI,1:NJ)
         Write(IO,'(100E25.16)') P(1:NI,1:NJ)

       END  SUBROUTINE 
      