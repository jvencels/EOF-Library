! *****************************************************************************/
! * License
! *    This file is part of EOF-Library.
! *
! *    EOF-Library is free software: you can redistribute it and/or modify it
! *    under the terms of the GNU General Public License as published by
! *    the Free Software Foundation, either version 3 of the License, or
! *    (at your option) any later version.
! *
! *    EOF-Library is distributed in the hope that it will be useful, but WITHOUT
! *    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! *    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! *    for more details.
! *
! *    You should have received a copy of the GNU General Public License
! *    along with EOF-Library.  If not, see <http://www.gnu.org/licenses/>.
! *
! *****************************************************************************/
! *
! *  Authors: Juris Vencels (EOF Consulting, Latvia)
! *           Peter Råback (CSC - IT Center for Science, Finland)
! *  
! *  Web:     https://eof-library.com
! *
! *  Original Date: 21.12.2016
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Module for receiving interpolated fields from OpenFOAM.
!------------------------------------------------------------------------------
MODULE OpenFOAM2ElmerSolverUtils

  USE DefUtils
  USE MeshUtils
  USE ElementUtils
  USE ParticleUtils

  IMPLICIT NONE

! derived types
!------------------------------------------------------------------------------  
  TYPE OFproc_t
    INTEGER :: reqSend, reqRecv, globalRank
    TYPE(Variable_t), POINTER :: OFVar
    TYPE(Mesh_t), POINTER :: OFMesh
    LOGICAL,POINTER :: foundCells(:)
    INTEGER,POINTER :: foundCellsIndx(:)
    INTEGER :: nElements
    INTEGER :: nFoundElements
    INTEGER,POINTER :: foundElementIndx(:)
    REAL(KIND=dp),POINTER :: recvValues(:)
    LOGICAL :: boxOverlap
  END TYPE OFproc_t

  TYPE(VariableTable_t), ALLOCATABLE :: OFVarTable(:)
!  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE ::   OFVarNames(:)

  INTEGER :: totOFRanks, OFRanksStart, ElmerRanksStart, myGlobalRank, &
             totGlobalRanks, myLocalRank, totLocalRanks, nVars, nElements, &
             totElementsFound, nBodiesToComm

!  INTEGER, ALLOCATABLE :: commElementId(:), elemPerNode(:)
!  LOGICAL, ALLOCATABLE :: commMaterial(:), commBody(:), commBodyForce(:)
  REAL(KIND=dp), POINTER :: commElementX(:), commElementY(:), commElementZ(:)
  TYPE(OFproc_t), ALLOCATABLE, TARGET :: OFp(:,:) ! [OFrank][body]
  REAL(KIND=dp) :: myBoundBox(3,2) ! [x,y,z][min,max]
  REAL(KIND=dp), POINTER :: ELboundBoxes(:,:,:,:) ! [x,y,z][min,max][rank][body]
  INTEGER, POINTER :: OF_EL_overlap(:,:,:) ! [ELrank][OFrank][body]

END MODULE OpenFOAM2ElmerSolverUtils

!------------------------------------------------------------------------------
SUBROUTINE MPI_TEST_SLEEP( req, ierr )

  USE ISO_C_BINDING
  USE OpenFOAM2ElmerSolverUtils

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE usleep(n) bind(C)
      USE ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(c_int32_t), VALUE :: n
    END SUBROUTINE usleep
  END INTERFACE

  !------------------------------------------------------------------------------
  INTEGER :: req, ierr
  LOGICAL :: Flag

  DO WHILE ( .TRUE. )
    CALL MPI_TEST( req, Flag, MPI_STATUS_IGNORE, ierr )
    IF (Flag) EXIT
    CALL usleep(1000_c_int32_t)
  END DO

END SUBROUTINE MPI_TEST_SLEEP

!------------------------------------------------------------------------------
SUBROUTINE findOverlappingBoxes(s)

  USE OpenFOAM2ElmerSolverUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  INTEGER :: ierr,  i, s
  INTEGER :: status(MPI_STATUS_SIZE)

  CALL MPI_ALLGATHER(myBoundBox, 6, MPI_DOUBLE, ELboundBoxes(:,:,:,s), 6, MPI_DOUBLE, ELMER_COMM_WORLD, ierr)

  IF ( myLocalRank==0 ) THEN
    CALL MPI_SEND(ELboundBoxes(:,:,:,s), totLocalRanks*2*3, MPI_DOUBLE, &
                  OFp(0,s) % globalRank, 1001, MPI_COMM_WORLD, ierr)
    CALL MPI_RECV(OF_EL_overlap(:,:,s), totOFRanks*totLocalRanks, MPI_INTEGER, &
                  OFp(0,s) % globalRank, 1002, MPI_COMM_WORLD, status, ierr)
  END IF

  CALL MPI_Bcast(OF_EL_overlap(:,:,s), totOFRanks*totLocalRanks, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

  DO i=0,totOFRanks-1
    OFp(i,s) % boxOverlap = (OF_EL_overlap(myLocalRank,i,s)==1)
  END DO

END SUBROUTINE findOverlappingBoxes

!------------------------------------------------------------------------------
SUBROUTINE OpenFOAM2ElmerSolver( Model,Solver,dt,TransientSimulation )
  
  USE DefUtils
  USE MeshUtils
  USE ElementUtils
  USE ParticleUtils
  USE OpenFOAM2ElmerSolverUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt, exitcond
  LOGICAL :: TransientSimulation

  ! local variables
  !------------------------------------------------------------------------------  
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var, VarCoord 
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  INTEGER :: i, j, k, l, n, ierr, OFstatus, s
  INTEGER :: status(MPI_STATUS_SIZE)
  LOGICAL :: Found
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp) :: DgScale, InvDgScale, ElemCenter(3), VarCenter
  INTEGER, POINTER :: Blist(:)

  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: Visited = .FALSE., UserDefinedCoordinates

  SAVE Visited, Mesh
  
  
  !------------------------------------------------------------------------------  

  CALL Info('OpenFOAM2ElmerSolver','-----------------------------------------', Level=4 )

  ! The variable containing the field contributions
  !--------------------------------------------------------------------------
  Params => GetSolverParams()
  Mesh => GetMesh()

  IF(.NOT. Visited ) THEN
    nVars = 0
    DO i=1,100
      VarName = ListGetString( Params, 'Target Variable '//TRIM(I2S(i)), Found )
      IF(.NOT. Found ) THEN
        IF (i==1) CALL Fatal('OpenFOAM2ElmerSolver','> Target Variable 1 < must exist for the solver!')
        EXIT
      ELSE
        nVars = nVars + 1
      END IF
    END DO
    
    ALLOCATE(OFVarTable(nVars))
    DO i=1,nVars
      VarName = ListGetString( Params, 'Target Variable '//TRIM(I2S(i)), Found )
      Var => VariableGet( Mesh % Variables, varName, ThisOnly = .TRUE. )
      IF(.NOT. ASSOCIATED( Var ) ) THEN
        CALL Fatal('OpenFOAM2ElmerSolver','Variable does not exist in Elmer mesh: '//TRIM(VarName))
      END IF
      OFVarTable(i) % Variable => Var
    END DO
    
    CALL Info('OpenFOAM2ElmerSolver','Number of target variables: '//TRIM(I2S(nVars)),Level=5)

    ! We can map variables that exist on nodal points, center points, or gauss points.
    ! However, the size and type of variables must be the same
    Var => OFVarTable(1) % Variable
    DO j=1,nVars
      IF( OFVarTable(j) % Variable % TYPE /= Var % TYPE ) THEN
        CALL Fatal('OpenFOAM2ElmerSolver','Currently it is assumed that all variables are of same type!')   
      END IF
      IF( SIZE( OFVarTable(j) % Variable % Values ) /= SIZE( Var % Values ) ) THEN
        CALL Fatal('OpenFOAM2ElmerSolver','Currently it is assumed that all variables are of same size!')   
      END IF
    END DO

    IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN
      DgScale = ListGetCReal( Params,'DG Interpolation Shrink Factor', Found )
      IF(.NOT. Found ) DgScale = 1.0_dp / SQRT( 3.0_dp )
      InvDgScale = 1.0_dp / DgScale
      WRITE( Message,'(A,ES12.3)') 'Using DG interpolation shrink factor:',DgScale
      CALL Info('OpenFOAM2ElmerSolver', Message )
    END IF

    ! MPI coupling
    !------------------------------------------------------------------------
    myLocalRank   = ParEnv % MyPE
    totLocalRanks = ParEnv % PEs

    IF( ListCheckPresent( Params,'Bodies') ) THEN
      BList => ListGetIntegerArray( Params, 'Bodies', Found ) 
      nBodiesToComm = SIZE(BList)
    ELSE
      nBodiesToComm = 1
    END IF

    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myGlobalRank, ierr )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, totGlobalRanks, ierr )

    totOFRanks = totGlobalRanks - totLocalRanks

    IF(totOFRanks==0) THEN
      CALL Fatal('OpenFOAM2ElmerSolver','MPI communicator does not have OpenFOAM procs!')
    END IF

    IF (myLocalRank == 0) ElmerRanksStart = myGlobalRank
    CALL MPI_BCAST(ElmerRanksStart, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

    IF (ElmerRanksStart==0) THEN
      OFRanksStart = totLocalRanks
    ELSE
      OFRanksStart = 0
    END IF

    CALL Info('OpenFOAM2ElmerSolver', 'Allocating OpenFOAM data structures',Level=3)
    ALLOCATE( OF_EL_overlap(0:totLocalRanks-1,0:totOFRanks-1,nBodiesToComm) )
    ALLOCATE( ELboundBoxes(3,2,0:totLocalRanks-1,nBodiesToComm) )
    ALLOCATE( OFp(0:totOFRanks-1,nBodiesToComm) )

    ! Get getboundBox
    myBoundBox(1,1) = MINVAL(Mesh % Nodes % x)
    myBoundBox(1,2) = MAXVAL(Mesh % Nodes % x)
    myBoundBox(2,1) = MINVAL(Mesh % Nodes % y)
    myBoundBox(2,2) = MAXVAL(Mesh % Nodes % y)
    myBoundBox(3,1) = MINVAL(Mesh % Nodes % z)
    myBoundBox(3,2) = MAXVAL(Mesh % Nodes % z)

    ! Here we assume that the variable that we are looking for is always of proper size
    ! specified by the 1st variable to be interpolated.
    nElements = size( Var % Values )

    ALLOCATE( commElementX(nElements),  &
              commElementY(nElements),  &
              commElementZ(nElements) )

    VarName = ListGetString( Params,'Interpolate Coordinate', UserDefinedCoordinates )
    
    IF( UserDefinedCoordinates ) THEN
      CALL Info('OpenFOAM2ElmerSolver','Using precomputed coordinates for interpolation')

      VarCoord => VariableGet( Mesh % Variables, VarName )
      IF( .NOT. ASSOCIATED( VarCoord ) ) THEN
        VarCoord => VariableGet( Mesh % Variables,"ip coordinate" )
      END IF
      IF( .NOT. ASSOCIATED( VarCoord ) ) THEN
        CALL Fatal('OpenFOAM2ElmerSolver','Missing >Interpolate Coordinate< needed for IPs')
      END IF
      IF( SIZE( VarCoord % Values ) / 3  /= nElements ) THEN
        CALL Fatal('OpenFOAM2ElmerSolver','Size of coordinates should match size of variable x 3!')
      END IF
      commElementX = VarCoord % Values(1::3)
      commElementY = VarCoord % Values(2::3)
      commElementZ = VarCoord % Values(3::3)

    ELSE IF (nElements > 0) THEN
      CALL Info('OpenFOAM2ElmerSolver','Allocating coordinates for interpolation')

      IF( Var % TYPE == Variable_on_nodes ) THEN

        ! Set coordinates for nodal interpolation
        ! Note that these are ordered such that f_i related to (x,y,z)_i 
        !------------------------------------------------------------------------
        DO i = 1, Mesh % NumberOfNodes 
          j = Var % Perm( i )
          IF( j > 0 ) THEN
            commElementX(j) = Mesh % Nodes % x(i)
            commElementY(j) = Mesh % Nodes % y(i)
            commElementZ(j) = Mesh % Nodes % z(i)            
          END IF
        END DO

      ELSE IF( Var % TYPE == Variable_on_elements ) THEN

        ! Set coordinates for elemental interpolation
        ! Note that these are ordered such that f_i related to (x,y,z)_i 
        !------------------------------------------------------------------------
        DO i = 1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
          Element => Mesh % Elements(i)
          n = Element % TYPE % NumberOfNodes
          j = Var % Perm( Element % ElementIndex )          
          IF( j > 0 ) THEN
            commElementX(j) = SUM( Mesh % Nodes % x( Element % NodeIndexes ) ) / n
            commElementY(j) = SUM( Mesh % Nodes % y( Element % NodeIndexes ) ) / n
            commElementZ(j) = SUM( Mesh % Nodes % z( Element % NodeIndexes ) ) / n
          END IF
        END DO

      ELSE IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN

        ! Set coordinates for elemental interpolation using DG
        !------------------------------------------------------------------------
        DO i = 1, Mesh % NumberOfBulkElements
          Element => Mesh % Elements(i)
          n = Element % TYPE % NumberOfNodes
          IF( ANY( Var % Perm( Element % DGIndexes ) == 0 ) ) CYCLE

          ! Compute element center
          ElemCenter(1) = SUM( Mesh % Nodes % x( Element % NodeIndexes ) ) / n
          ElemCenter(2) = SUM( Mesh % Nodes % y( Element % NodeIndexes ) ) / n
          ElemCenter(3) = SUM( Mesh % Nodes % z( Element % NodeIndexes ) ) / n

          ! Scale the element to be smaller
          commElementX(Var % Perm(Element % DGIndexes)) = ElemCenter(1) + &
              DgScale * ( Mesh % Nodes % x( Element % NodeIndexes ) - ElemCenter(1) )

          commElementY(Var % Perm(Element % DGIndexes)) = ElemCenter(2) + &
              DgScale * ( Mesh % Nodes % y( Element % NodeIndexes ) - ElemCenter(2) )

          commElementZ(Var % Perm(Element % DGIndexes)) = ElemCenter(3) + &
              DgScale * ( Mesh % Nodes % z( Element % NodeIndexes ) - ElemCenter(3) )
        END DO

      ELSE IF( Var % TYPE == Variable_on_gauss_points ) THEN
        CALL Fatal('OpenFOAM2ElmerSolver','Coordinates for gauss point interpolation not given')

      ELSE 
        CALL Fatal('OpenFOAM2ElmerSolver','Unknown variable type to create interpolation coordinates')      
      END IF
    END IF

    CALL Info('OpenFOAM2ElmerSolver','Total number of bodiesto communicate: '//TRIM(I2S(nBodiesToComm)), Level=3 )

    totElementsFound = 0

    DO s=1,nBodiesToComm
      DO i=0,totOFRanks-1
        ALLOCATE( OFp(i,s) % OFMesh )
        OFp(i,s) % globalRank = i + OFRanksStart
      END DO
      CALL findOverlappingBoxes(s)

      ! Starting communication
      !------------------------------------------------------------------------
      CALL Info('OpenFOAM2ElmerSolver','Starting communication for Body #'//TRIM(I2S(s)), Level=3 )

      DO i = 0, totOFRanks - 1
        IF(.NOT.OFp(i,s) % boxOverlap) CYCLE
        ! Number of Elmer elements
        CALL MPI_ISEND(nElements, 1, MPI_INTEGER, &
            OFp(i,s) % globalRank, 899, MPI_COMM_WORLD, OFp(i,s) % reqSend, ierr)
      END DO

      CALL Info('OpenFOAM2ElmerSolver','Sending coordinates', Level=3 )
      DO i = 0, totOFRanks - 1
        IF(.NOT.OFp(i,s) % boxOverlap) CYCLE
        ! Number of Elmer elements
        CALL MPI_TEST_SLEEP(OFp(i,s) % reqSend, ierr)
        ! Element X coords
        CALL MPI_ISEND(commElementX, nElements, MPI_DOUBLE, &
            OFp(i,s) % globalRank, 897, MPI_COMM_WORLD, OFp(i,s) % reqSend, ierr)
        CALL MPI_REQUEST_FREE(OFp(i,s) % reqSend, ierr)
        CALL MPI_ISEND(commElementY, nElements, MPI_DOUBLE, &
            OFp(i,s) % globalRank, 897, MPI_COMM_WORLD, OFp(i,s) % reqSend, ierr)
        CALL MPI_REQUEST_FREE(OFp(i,s) % reqSend, ierr)
        CALL MPI_ISEND(commElementZ, nElements, MPI_DOUBLE, &
            OFp(i,s) % globalRank, 897, MPI_COMM_WORLD, OFp(i,s) % reqSend, ierr)
      END DO

      DO i = 0, totOFRanks - 1
        OFp(i,s) % nFoundElements = 0 ! keep this
        IF(.NOT.OFp(i,s) % boxOverlap) CYCLE
        ! Element Z coords
        CALL MPI_TEST_SLEEP(OFp(i,s) % reqSend, ierr)
        ! Number of found elements
        CALL MPI_IRECV(OFp(i,s) % nFoundElements, 1, MPI_INTEGER, &
            OFp(i,s) % globalRank, 895, MPI_COMM_WORLD, OFp(i,s) % reqRecv, ierr)
      END DO

      DO i = 0, totOFRanks - 1
        IF(.NOT.OFp(i,s) % boxOverlap) CYCLE
        ! Number of found elements
        CALL MPI_TEST_SLEEP(OFp(i,s) % reqRecv, ierr)

        IF (OFp(i,s) % nFoundElements == 0) CYCLE
        totElementsFound = totElementsFound + OFp(i,s) % nFoundElements

        ALLOCATE ( OFp(i,s) % foundElementIndx(OFp(i,s) % nFoundElements), &
                   OFp(i,s) % recvValues(OFp(i,s) % nFoundElements))

        ! Indexes of found elements
        CALL MPI_IRECV(OFp(i,s) % foundElementIndx, OFp(i,s) % nFoundElements, MPI_INTEGER, &
            OFp(i,s) % globalRank, 894, MPI_COMM_WORLD, OFp(i,s) % reqRecv, ierr)
      END DO

      DO i = 0, totOFRanks - 1
        IF (OFp(i,s) % nFoundElements == 0) CYCLE
        ! Indexes of found elements
        CALL MPI_TEST_SLEEP(OFp(i,s) % reqRecv, ierr)
        OFp(i,s) % foundElementIndx = OFp(i,s) % foundElementIndx + 1
      END DO

    END DO ! nBodiesToComm

    IF (totElementsFound < nElements) THEN
      CALL Fatal('OpenFOAM2ElmerSolver','Elmer #'//TRIM(I2S(myLocalRank))//' has ' &
                 //TRIM(I2S(nElements))//' elements, OpenFOAM found '//TRIM(I2S(totElementsFound)))
    END IF

  END IF ! Visited

  DO s=1,nBodiesToComm
    ! Receive simulation status
    CALL MPI_IRECV( OFstatus, 1, MPI_INTEGER, OFp(0,s) % globalRank, 799, MPI_COMM_WORLD, OFp(0,s) % reqRecv, ierr)
    CALL MPI_TEST_SLEEP(OFp(0,s) % reqRecv, ierr)

    IF (OFstatus.NE.1) THEN
      CALL Info('OpenFOAM2ElmerSolver','Elmer has last iteration!', Level=3 )
      exitcond = ListGetCReal( CurrentModel % Simulation,'Exit Condition',Found)
      IF(.NOT.Found) CALL ListAddConstReal(CurrentModel % Simulation,'Exit Condition',1.0_dp)
    END IF

    ! Receive fields
    DO j=1,nVars
      Var => OFVarTable(j) % Variable 
      CALL Info('OpenFOAM2ElmerSolver','Receiving: '//TRIM(Var % name),Level=5)

      DO i = 0, totOFRanks - 1
        IF (OFp(i,s) % nFoundElements == 0) CYCLE
        CALL MPI_IRECV( OFp(i,s) % recvValues, OFp(i,s) % nFoundElements, MPI_DOUBLE, &
                        OFp(i,s) % globalRank, 900, MPI_COMM_WORLD, OFp(i,s) % reqRecv, ierr)
      END DO

      DO i = 0, totOFRanks - 1
        IF (OFp(i,s) % nFoundElements == 0) CYCLE
        CALL MPI_TEST_SLEEP(OFp(i,s) % reqRecv, ierr)

#ifndef SKIP_FP_FIX
        ! Fix floating point exception
        WHERE (OFp(i,s) % recvValues<EPSILON(1.e0) &
                         .AND. OFp(i,s) % recvValues > -EPSILON(1.e0))
          OFp(i,s) % recvValues = 0
        END WHERE
#else
#warning "Dirty fix for floating point exception is disabled!"
#endif

        ! Directly interpolate to the points needed, don't reinterpolate
        DO k = 1, OFp(i,s) % nFoundElements
          l = OFp(i,s) % foundElementIndx(k)
          Var % Values(l) = OFp(i,s) % recvValues(k)
        END DO
      END DO

      ! We used shrinked version of DG element. Now extrapolate the interpolated values
      ! to the nodes.
      IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN

        DO i = 1, Mesh % NumberOfBulkElements
          Element => Mesh % Elements(i)
          n = Element % TYPE % NumberOfNodes
          IF( ANY( Var % Perm( Element % DGIndexes ) == 0 ) ) CYCLE

          ! Compute element center
          VarCenter = SUM( Var % Values(Var % Perm(Element % DGIndexes ) ) ) / n

          Var % Values(Var % Perm(Element % DGIndexes ) ) = &
              VarCenter + InvDgScale * ( Var % Values(Var % Perm(Element % DGIndexes ) ) - VarCenter )
        END DO

      END IF ! Variable_on_nodes_on_elements

    END DO ! nVars
  END DO ! nBodiesToComm

  VISITED = .TRUE.

  CALL Info('OpenFOAM2ElmerSolver','All done', Level=4 )
  CALL Info('OpenFOAM2ElmerSolver','-----------------------------------------', Level=4 )
  
!------------------------------------------------------------------------------

END SUBROUTINE OpenFOAM2ElmerSolver
