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
! *  Web:     http://eof-library.com
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
  END TYPE OFproc_t

  TYPE(OFproc_t), ALLOCATABLE, TARGET :: OFp(:)
  INTEGER :: totOFRanks, OFRanksStart, ElmerRanksStart, myGlobalRank, &
                   totGlobalRanks, myLocalRank, totLocalRanks, nVars

END MODULE Elmer2OpenFOAMSolverUtils

!------------------------------------------------------------------------------
SUBROUTINE MPI_TEST_SLEEP( req, ierr )

  USE ISO_C_BINDING, ONLY : C_LONG
  USE Elmer2OpenFOAMSolverUtils

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE nanosleep(n) BIND(C,name="nanosleep")
      USE ISO_C_BINDING, ONLY : C_LONG
      INTEGER(C_LONG), VALUE :: n
    END SUBROUTINE nanosleep
  END INTERFACE

  !------------------------------------------------------------------------------
  INTEGER :: req, ierr
  LOGICAL :: Flag

  DO WHILE ( .TRUE. )
    CALL MPI_TEST( req, Flag, MPI_STATUS_IGNORE, ierr )
    IF (Flag) EXIT
    CALL nanosleep(10000000_C_LONG)
  END DO

END SUBROUTINE MPI_TEST_SLEEP

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

  LOGICAL, SAVE :: VISITED = .FALSE.

  !------------------------------------------------------------------------------  

  CALL Info('Elmer2OpenFOAMSolver','-----------------------------------------', Level=4 )

  IF (.NOT. VISITED) THEN
    ! The variable containing the field contributions
    !--------------------------------------------------------------------------
    Params => GetSolverParams()
    nVars = 0

    DO i=1,100
      VarName = ListGetString( Params, 'Target Variable '//TRIM(I2S(i)), Found )
      IF(.NOT. Found ) THEN
        IF (i==1) CALL Fatal('Elmer2OpenFOAMSolver','> Target Variable 1 < must exist for the solver!')
        EXIT
      ELSE
       ! Test that the variable exists in the primary mesh
        Var => VariableGet(CurrentModel % Mesh % Variables, VarName )
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

    CALL MPI_ALLREDUCE(myGlobalRank, ElmerRanksStart, 1, MPI_INTEGER, MPI_MIN, ELMER_COMM_WORLD, ierr)

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

    ! Starting communication
    !------------------------------------------------------------------------

    DO i = 0, totOFRanks - 1
      ! Number of OpenFOAM cells
      CALL MPI_IRECV(OFp(i) % OFMesh % NumberOfNodes, 1, MPI_INTEGER, &
                      OFp(i) % globalRank, 999, MPI_COMM_WORLD, OFp(i) % reqRecv, ierr)
    END DO

    DO i = 0, totOFRanks - 1
      CALL MPI_TEST_SLEEP(OFp(i) % reqRecv, ierr)

      CALL Info('Elmer2OpenFOAMSolver','Receiving '//TRIM(I2S(OFp(i) % OFMesh % NumberOfNodes))// &
                ' cells from OpenFOAM proc #'//TRIM(I2S(i)),Level=3)

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
      CALL MPI_RECV(OFp(i) % OFMesh % Nodes % x, OFp(i) % OFMesh % NumberOfNodes, MPI_DOUBLE, &
                     OFp(i) % globalRank, 998, MPI_COMM_WORLD, status, ierr)
    END DO
    DO i = 0, totOFRanks - 1
      ! wait for x coordinates
      !CALL MPI_TEST_SLEEP(OFp(i) % reqRecv, ierr)

      ! Cell y coordinates
      CALL MPI_RECV(OFp(i) % OFMesh % Nodes % y, OFp(i) % OFMesh % NumberOfNodes, MPI_DOUBLE, &
                     OFp(i) % globalRank, 997, MPI_COMM_WORLD, status, ierr)
    END DO

    DO i = 0, totOFRanks - 1
      ! wait for y coordinates
      !CALL MPI_TEST_SLEEP(OFp(i) % reqRecv, ierr)

      ! Cell z coordinates
      CALL MPI_RECV(OFp(i) % OFMesh % Nodes % z, OFp(i) % OFMesh % NumberOfNodes, MPI_DOUBLE, &
                     OFp(i) % globalRank, 996, MPI_COMM_WORLD, status, ierr)
    END DO

    CALL Info('Elmer2OpenFOAMSolver','Projecting field to OpenFOAM cell centers',Level=10) 
    DO i = 0, totOFRanks - 1
      ! wait for z coordinates
      !CALL MPI_TEST_SLEEP(OFp(i) % reqRecv, ierr)
      IF ( CoordinateSystemDimension() == 2 ) OFp(i) % OFMesh % Nodes % z = 0

      CALL InterpolateMeshToMeshQ( OldMesh          = CurrentModel % Mesh, &
                                   NewMesh          = OFp(i) % OFMesh, &
                                   UseQuadrantTree  = .TRUE., &
                                   Projector        = OFp(i) % OFMesh % Projector, &
                                   FoundNodes       = OFp(i) % foundCells, &
                                   KeepUnfoundNodes = .FALSE.)

      OFp(i) % nFoundCells = COUNT(OFp(i) % foundCells)

      IF ( OFp(i) % nFoundCells > 0 ) THEN
        ALLOCATE( OFp(i) % OFVar % Values( OFp(i) % nFoundCells ), &
                  OFp(i) % foundCellsIndx( OFp(i) % nFoundCells ), &
                  OFp(i) % OFVar % Perm( OFp(i) % nFoundCells ) )

        OFp(i) % OFVar % Perm = (/ (j, j = 1, OFp(i) % nFoundCells) /)
      END IF

      ! Number of cells found in each Elmer process
      CALL MPI_ISEND( OFp(i) % nFoundCells, 1, MPI_INTEGER, &
                      OFp(i) % globalRank, 995, MPI_COMM_WORLD, OFp(i) % reqSend, ierr)
    END DO

    DO i = 0, totOFRanks - 1
      ! wait for nFoundCells
      CALL MPI_TEST_SLEEP(OFp(i) % reqSend, ierr)
      IF ( OFp(i) % nFoundCells > 0 ) THEN
        OFp(i) % foundCellsIndx = PACK((/ (j, j = 0, OFp(i) % OFMesh % NumberOfNodes-1) /),OFp(i) % foundCells)

        ! Indexes for cells that were found on this piece of Elmer mesh
        CALL MPI_ISEND( OFp(i) % foundCellsIndx, OFp(i) % nFoundCells, MPI_INTEGER, &
                        OFp(i) % globalRank, 994, MPI_COMM_WORLD, OFp(i) % reqSend, ierr)
      END IF
    END DO

    DO i = 0, totOFRanks - 1
      IF ( OFp(i) % nFoundCells > 0 ) THEN
        ! wait for foundCells
        CALL MPI_TEST_SLEEP(OFp(i) % reqSend, ierr)
      END IF
    END DO

  END IF ! .NOT. VISITED

  Params => GetSolverParams()

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
    Var => VariableGet( CurrentModel % Mesh % Variables, VarName )
    IF(.NOT. ASSOCIATED( Var ) ) THEN
      CALL Fatal('Elmer2OpenFOAMSolver','Variable '//TRIM(VarName)//' does not exist in Elmer mesh!')
    END IF

    DO i = 0, totOFRanks - 1
      IF ( OFp(i) % nFoundCells > 0 ) THEN
        OFp(i) % OFVar % Values = 0
        CALL CRS_ApplyProjector( OFp(i) % OFMesh % Projector % Matrix, Var % Values, &
                     Var % Perm, OFp(i) % OFVar % Values, OFp(i) % OFVar % Perm )

        CALL MPI_ISEND( OFp(i) % OFVar % Values, OFp(i) % nFoundCells, MPI_DOUBLE, &
                        OFp(i) % globalRank, 1000, MPI_COMM_WORLD, OFp(i) % reqSend, ierr)
      END IF
    END DO

    DO i = 0, totOFRanks - 1
      IF ( OFp(i) % nFoundCells > 0 ) THEN
        CALL MPI_TEST_SLEEP(OFp(i) % reqSend, ierr)
      END IF
    END DO
  END DO

  VISITED = .TRUE.

  CALL Info('Elmer2OpenFOAMSolver','All done', Level=4 )
  CALL Info('Elmer2OpenFOAMSolver','-----------------------------------------', Level=4 )
  
END SUBROUTINE Elmer2OpenFOAMSolver
