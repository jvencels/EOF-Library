! *****************************************************************************/
! * Copyright (C) 2016-2018 University of Latvia
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
! ******************************************************************************
! *
! *  Authors: Juris Vencels (University of Latvia)
! *           Peter Råback (CSC - IT Center for Science, Finland)
! *  
! *  Email:   juris.vencels@lu.lv
! *  
! *  Web:     https://eof-library.com
! *
! *  Original Date: 29.09.2016
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Module for sending interpolated fields to OpenFOAM.
!------------------------------------------------------------------------------
MODULE Elmer2OpenFOAMSolverUtils

  USE DefUtils
  USE MeshUtils
  USE ElementUtils
  USE ParticleUtils

  IMPLICIT NONE

  TYPE OFproc_t
    INTEGER :: reqSend, reqRecv, globalRank
    TYPE(Variable_t), POINTER :: OFVar
    TYPE(Mesh_t), POINTER :: OFMesh
    LOGICAL,POINTER :: foundCells(:)
    INTEGER,POINTER :: foundCellsIndx(:)
    INTEGER :: nFoundCells
    LOGICAL :: boxOverlap
  END TYPE OFproc_t

  TYPE(OFproc_t), ALLOCATABLE, TARGET :: OFp(:)
  INTEGER :: totOFRanks, OFRanksStart, ElmerRanksStart, myGlobalRank, &
                   totGlobalRanks, myLocalRank, totLocalRanks, nVars
  REAL(KIND=dp) :: myBoundBox(3,2) ! [x,y,z][min,max]
  REAL(KIND=dp), POINTER :: ELboundBoxes(:,:,:) ! [x,y,z][min,max][rank]
  INTEGER, POINTER :: OF_EL_overlap(:,:) ! [ELrank][OFrank]

END MODULE Elmer2OpenFOAMSolverUtils

!------------------------------------------------------------------------------
SUBROUTINE MPI_TEST_SLEEP( req, ierr )

  USE ISO_C_BINDING, ONLY : C_LONG
  USE Elmer2OpenFOAMSolverUtils

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
SUBROUTINE findOverlappingBoxes()

  USE Elmer2OpenFOAMSolverUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  INTEGER :: ierr,  i
  INTEGER :: status(MPI_STATUS_SIZE)

  CALL MPI_ALLGATHER(myBoundBox, 6, MPI_DOUBLE, ELboundBoxes, 6, MPI_DOUBLE, ELMER_COMM_WORLD, ierr)

  IF ( myLocalRank==0 ) THEN
    CALL MPI_SEND(ELboundBoxes, totLocalRanks*2*3, MPI_DOUBLE, OFp(0) % globalRank, 1001, MPI_COMM_WORLD, ierr)
    CALL MPI_RECV(OF_EL_overlap, totOFRanks*totLocalRanks, MPI_INTEGER, OFp(0) % globalRank, 1002, MPI_COMM_WORLD, status, ierr)
  END IF

  CALL MPI_Bcast(OF_EL_overlap, totOFRanks*totLocalRanks, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

  DO i=0,totOFRanks-1
    OFp(i) % boxOverlap = (OF_EL_overlap(myLocalRank,i)==1)
  END DO

END SUBROUTINE findOverlappingBoxes

!------------------------------------------------------------------------------
SUBROUTINE Elmer2OpenFOAMSolver( Model,Solver,dt,TransientSimulation )
  
  USE DefUtils
  USE Interpolation
  USE MeshUtils
  USE ElementUtils
  USE ParticleUtils
  USE Elmer2OpenFOAMSolverUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt, exitcond
  LOGICAL :: TransientSimulation

  ! local variables
  !------------------------------------------------------------------------------  
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  INTEGER :: i, j, ierr, OFstatus
  INTEGER :: status(MPI_STATUS_SIZE)
  LOGICAL :: Found
  REAL(KIND=dp) :: commTime
  CHARACTER(LEN=15) :: timeStr
  
  INTERFACE
    SUBROUTINE InterpolateMeshToMeshQ( OldMesh, NewMesh, OldVariables, NewVariables, &
        UseQuadrantTree, Projector, MaskName, FoundNodes, NewMaskPerm, KeepUnfoundNodes )
      USE Types
      TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
      TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
      LOGICAL, OPTIONAL :: UseQuadrantTree,FoundNodes(:)
      CHARACTER(LEN=*),OPTIONAL :: MaskName
      TYPE(Projector_t), POINTER, OPTIONAL :: Projector
      INTEGER, OPTIONAL, POINTER :: NewMaskPerm(:)  !< Mask the new variable set by the given MaskName when trying to define the interpolation.
       LOGICAL, OPTIONAL :: KeepUnfoundNodes  !< Do not disregard unfound nodes from projector
    END SUBROUTINE InterpolateMeshToMeshQ
  END INTERFACE

  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL, SAVE :: VISITED = .FALSE.

  !------------------------------------------------------------------------------  

  CALL Info('Elmer2OpenFOAMSolver','-----------------------------------------', Level=4 )

  commTime = MPI_WTIME()

  ! The variable containing the field contributions
  !--------------------------------------------------------------------------
  Params => GetSolverParams()
  Mesh => GetMesh()

  IF (.NOT. VISITED) THEN
    nVars = 0

    DO i=1,100
      VarName = ListGetString( Params, 'Target Variable '//TRIM(I2S(i)), Found )
      IF(.NOT. Found ) THEN
        IF (i==1) CALL Fatal('Elmer2OpenFOAMSolver','> Target Variable 1 < must exist for the solver!')
        EXIT
      ELSE
       ! Test that the variable exists in the primary mesh
        Var => VariableGet(Mesh % Variables, VarName )
        IF(.NOT. ASSOCIATED( Var ) ) THEN
          CALL Fatal('Elmer2OpenFOAMSolver','Variable '//TRIM(VarName)//' does not exist in Elmer mesh!')
        ELSE
          nVars = nVars + 1
        END IF
      END IF
    END DO

    CALL Info('Elmer2OpenFOAMSolver','Number of target variables: '//TRIM(I2S(nVars)),Level=3)

    ! MPI coupling
    !------------------------------------------------------------------------
    myLocalRank   = ParEnv % MyPE
    totLocalRanks = ParEnv % PEs
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myGlobalRank, ierr )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, totGlobalRanks, ierr )

    totOFRanks = totGlobalRanks - totLocalRanks

    IF(totOFRanks==0) THEN
      CALL Fatal('Elmer2OpenFOAMSolver','MPI communicator does not have OpenFOAM procs!')
    END IF

    IF (myLocalRank == 0) ElmerRanksStart = myGlobalRank
    CALL MPI_BCAST(ElmerRanksStart, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

    IF (ElmerRanksStart==0) THEN
       OFRanksStart = totLocalRanks
    ELSE
       OFRanksStart = 0
    END IF

    CALL Info('Elmer2OpenFOAMSolver', 'Allocating OpenFOAM data structures',Level=3)
    ALLOCATE( OFp(0:totOFRanks-1) )

    DO i=0,totOFRanks-1
      ALLOCATE( OFp(i) % OFMesh )
      OFp(i) % globalRank = i + OFRanksStart
    END DO

    ! Get getboundBox
    myBoundBox(1,1) = MINVAL(Mesh % Nodes % x)
    myBoundBox(1,2) = MAXVAL(Mesh % Nodes % x)
    myBoundBox(2,1) = MINVAL(Mesh % Nodes % y)
    myBoundBox(2,2) = MAXVAL(Mesh % Nodes % y)
    myBoundBox(3,1) = MINVAL(Mesh % Nodes % z)
    myBoundBox(3,2) = MAXVAL(Mesh % Nodes % z)

    ALLOCATE( OF_EL_overlap(0:totLocalRanks-1,0:totOFRanks-1) )
    ALLOCATE( ELboundBoxes(3,2,0:totLocalRanks-1) )

    CALL findOverlappingBoxes()

    ! Starting communication
    !------------------------------------------------------------------------

    DO i = 0, totOFRanks - 1
      IF(.NOT.OFp(i) % boxOverlap) CYCLE
      ! Number of OpenFOAM cells
      CALL MPI_IRECV(OFp(i) % OFMesh % NumberOfNodes, 1, MPI_INTEGER, &
                      OFp(i) % globalRank, 999, MPI_COMM_WORLD, OFp(i) % reqRecv, ierr)
    END DO

    DO i = 0, totOFRanks - 1
      IF(.NOT.OFp(i) % boxOverlap) CYCLE
      CALL MPI_TEST_SLEEP(OFp(i) % reqRecv, ierr)

      ALLOCATE( OFp(i) % OFMesh % Nodes, &
                OFp(i) % OFMesh % Variables, &
                OFp(i) % OFVar ) 
      ALLOCATE( OFp(i) % OFMesh % Nodes % x( OFp(i) % OFMesh % NumberOfNodes ), &
                OFp(i) % OFMesh % Nodes % y( OFp(i) % OFMesh % NumberOfNodes ), &
                OFp(i) % OFMesh % Nodes % z( OFp(i) % OFMesh % NumberOfNodes ), &
                OFp(i) % foundCells( OFp(i) % OFMesh % NumberOfNodes ) )

      OFp(i) % OFMesh % NumberOfBulkElements = 0
      OFp(i) % OFMesh % NumberOfBoundaryElements = 0
      OFp(i) % OFMesh % Projector => NULL()
      OFp(i) % foundCells = .FALSE.

      ! Cell x coordinates
      CALL MPI_IRECV(OFp(i) % OFMesh % Nodes % x, OFp(i) % OFMesh % NumberOfNodes, MPI_DOUBLE, &
                     OFp(i) % globalRank, 997, MPI_COMM_WORLD, OFp(i) % reqRecv, ierr)
      CALL MPI_REQUEST_FREE(OFp(i) % reqRecv, ierr)
      CALL MPI_IRECV(OFp(i) % OFMesh % Nodes % y, OFp(i) % OFMesh % NumberOfNodes, MPI_DOUBLE, &
                     OFp(i) % globalRank, 997, MPI_COMM_WORLD, OFp(i) % reqRecv, ierr)
      CALL MPI_REQUEST_FREE(OFp(i) % reqRecv, ierr)
      CALL MPI_IRECV(OFp(i) % OFMesh % Nodes % z, OFp(i) % OFMesh % NumberOfNodes, MPI_DOUBLE, &
                     OFp(i) % globalRank, 997, MPI_COMM_WORLD, OFp(i) % reqRecv, ierr)
    END DO

    CALL Info('Elmer2OpenFOAMSolver','Projecting field to OpenFOAM cell centers',Level=10) 
    DO i = 0, totOFRanks - 1
      OFp(i) % nFoundCells = 0 ! keep this
      IF(.NOT.OFp(i) % boxOverlap) CYCLE
      ! wait for z coordinates
      CALL MPI_TEST_SLEEP(OFp(i) % reqRecv, ierr)

      IF ( CoordinateSystemDimension() == 2 ) THEN
        IF ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
             CurrentCoordinateSystem() == CylindricSymmetric ) THEN
          OFp(i) % OFMesh % Nodes % x = SQRT(OFp(i) % OFMesh % Nodes % x**2 + OFp(i) % OFMesh % Nodes % z**2)
        END IF
        OFp(i) % OFMesh % Nodes % z = 0
      END IF

      CALL InterpolateMeshToMeshQ( OldMesh          = Mesh, &
                                   NewMesh          = OFp(i) % OFMesh, &
                                   UseQuadrantTree  = .TRUE., &
                                   Projector        = OFp(i) % OFMesh % Projector, &
                                   FoundNodes       = OFp(i) % foundCells, &
                                   KeepUnfoundNodes = .FALSE.)

      OFp(i) % nFoundCells = COUNT(OFp(i) % foundCells)

      ! Number of cells found in each Elmer process
      CALL MPI_ISEND( OFp(i) % nFoundCells, 1, MPI_INTEGER, &
                      OFp(i) % globalRank, 995, MPI_COMM_WORLD, OFp(i) % reqSend, ierr)
    END DO

    DO i = 0, totOFRanks - 1
      IF(.NOT.OFp(i) % boxOverlap) CYCLE
      ! wait for nFoundCells
      CALL MPI_TEST_SLEEP(OFp(i) % reqSend, ierr)

      IF ( OFp(i) % nFoundCells == 0 ) CYCLE
      ALLOCATE( OFp(i) % OFVar % Values( OFp(i) % nFoundCells ), &
                OFp(i) % foundCellsIndx( OFp(i) % nFoundCells ), &
                OFp(i) % OFVar % Perm( OFp(i) % nFoundCells ) )

      OFp(i) % OFVar % Perm = (/ (j, j = 1, OFp(i) % nFoundCells) /)
      OFp(i) % foundCellsIndx = PACK((/ (j, j = 0, OFp(i) % OFMesh % NumberOfNodes-1) /),OFp(i) % foundCells)

      ! Indexes for cells that were found on this piece of Elmer mesh
      CALL MPI_ISEND( OFp(i) % foundCellsIndx, OFp(i) % nFoundCells, MPI_INTEGER, &
                      OFp(i) % globalRank, 994, MPI_COMM_WORLD, OFp(i) % reqSend, ierr)
    END DO

    DO i = 0, totOFRanks - 1
      IF ( OFp(i) % nFoundCells == 0 ) CYCLE
      CALL MPI_TEST_SLEEP(OFp(i) % reqSend, ierr)
    END DO

  END IF ! .NOT. VISITED

  ! Receive simulation status
  CALL MPI_IRECV( OFstatus, 1, MPI_INTEGER, OFp(0) % globalRank, 799, MPI_COMM_WORLD, OFp(0) % reqRecv, ierr)
  CALL MPI_TEST_SLEEP(OFp(0) % reqRecv, ierr)

  IF (OFstatus.NE.1) THEN
    CALL Info('Elmer2OpenFOAM','Elmer has last iteration!', Level=3 )
    exitcond = ListGetCReal( CurrentModel % Simulation,'Exit Condition',Found)
	  IF(.NOT.Found) CALL ListAddConstReal(CurrentModel % Simulation,'Exit Condition',1.0_dp)
  END IF

  ! Send fields
  DO j=1,nVars
    VarName = ListGetString( Params, 'Target Variable '//TRIM(I2S(j)), Found )
    Var => VariableGet( Mesh % Variables, VarName )
    IF(.NOT. ASSOCIATED( Var ) ) THEN
      CALL Fatal('Elmer2OpenFOAMSolver','Variable '//TRIM(VarName)//' does not exist in Elmer mesh!')
    END IF

    DO i = 0, totOFRanks - 1
      IF ( OFp(i) % nFoundCells == 0 ) CYCLE
      OFp(i) % OFVar % Values = 0
      CALL CRS_ApplyProjector( OFp(i) % OFMesh % Projector % Matrix, Var % Values, &
                   Var % Perm, OFp(i) % OFVar % Values, OFp(i) % OFVar % Perm )

      CALL MPI_ISEND( OFp(i) % OFVar % Values, OFp(i) % nFoundCells, MPI_DOUBLE, &
                      OFp(i) % globalRank, 1000, MPI_COMM_WORLD, OFp(i) % reqSend, ierr)
    END DO

    DO i = 0, totOFRanks - 1
      IF ( OFp(i) % nFoundCells == 0 ) CYCLE
      CALL MPI_TEST_SLEEP(OFp(i) % reqSend, ierr)
    END DO
  END DO

  VISITED = .TRUE.
  write(timeStr , '(F9.5)') MPI_WTIME() - commTime

  CALL Info('Elmer2OpenFOAM',' = '//TRIM(timeStr)//' s', Level=3 )

  CALL Info('Elmer2OpenFOAMSolver','All done', Level=4 )
  CALL Info('Elmer2OpenFOAMSolver','-----------------------------------------', Level=4 )
  
END SUBROUTINE Elmer2OpenFOAMSolver
