!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juris Vencels, Peter Råback 
! *  Email:   juris.vencels@gmail.com
! *  Web:     http://vencels.com
! *  Address: University of Latvia
! *           Laboratory for mathematical modelling of 
! *               environmental and technological processes
! *           Zellu Str. 23, Riga, LV-1002, Latvia
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
  END TYPE OFproc_t

  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE ::   OFVarNames(:)

  INTEGER :: totOFRanks, OFRanksStart, ElmerRanksStart, myGlobalRank, &
             totGlobalRanks, myLocalRank, totLocalRanks, nVars, nElements, &
             totElementsFound

  INTEGER, ALLOCATABLE :: commElementId(:), elemPerNode(:)
  LOGICAL, ALLOCATABLE :: commMaterial(:), commBody(:)
  REAL(KIND=dp), ALLOCATABLE :: commElementX(:), commElementY(:), commElementZ(:)
  TYPE(OFproc_t), ALLOCATABLE, TARGET :: OFp(:)

END MODULE OpenFOAM2ElmerSolverUtils

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
  TYPE(Variable_t), POINTER :: Var
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  INTEGER :: i, j, k, n, ierr, OFstatus
  LOGICAL :: Found
  TYPE(ValueList_t), POINTER :: Material
  TYPE(ValueListEntry_t), POINTER :: ptrVar
  TYPE(Element_t), POINTER :: Element
  LOGICAL, SAVE :: VISITED = .FALSE.

  !------------------------------------------------------------------------------  

  CALL Info('OpenFOAM2ElmerSolver','-----------------------------------------', Level=4 )

  IF (.NOT. VISITED) THEN
    ! The variable containing the field contributions
    !--------------------------------------------------------------------------
    Params => GetSolverParams()

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

    ALLOCATE(OFVarNames(nVars))
    ALLOCATE(commMaterial(CurrentModel % NumberOfMaterials))
    ALLOCATE(commBody(CurrentModel % NumberOfBodies))

    commMaterial = .FALSE.
    commBody = .FALSE.

    DO i=1,nVars
      VarName = ListGetString( Params, 'Target Variable '//TRIM(I2S(i)), Found )

      ! Find variables, materials and bodies that need to be communicated with OF
      DO j=1,CurrentModel % NumberOfMaterials
        Material => CurrentModel % Materials(j) % Values
        ptrVar => ListFind( Material, VarName, Found )
        IF(Found) THEN
          IF (ptrVar%DepNameLen==0) CYCLE
          Var => VariableGet(CurrentModel % Mesh % Variables, ptrVar%DependName  )
          IF(ASSOCIATED( Var ) ) THEN
            ! Name of communicated variable
            OFVarNames(i) = ptrVar%DependName
            CALL Info('OpenFOAM2ElmerSolver','Interpolated variable: '//OFVarNames(i),Level=3)

            ! Material that has to be communicated
            commMaterial(j) = .TRUE.
            CALL Info('OpenFOAM2ElmerSolver','Material: '//TRIM(I2S(j)),Level=3)

            ! Bodies that have to be communicated
            DO k=1,CurrentModel % NumberOfBodies
              IF (j==ListGetInteger( CurrentModel % Bodies(k) % Values, 'Material')) THEN
                commBody(k) = .TRUE.
                CALL Info('OpenFOAM2ElmerSolver','Body: '//TRIM(I2S(k)),Level=3)
              END IF
            END DO ! NumberOfBodies
            EXIT
          END IF
        END IF
      END DO ! NumberOfMaterials

      ! Test that the variable exists in the primary mesh
      IF(.NOT. ASSOCIATED( Var ) ) THEN
        CALL Fatal('OpenFOAM2ElmerSolver','Variable does not exist in Elmer mesh!')
      END IF
    END DO ! nVars

    CALL Info('OpenFOAM2ElmerSolver','Number of target variables: '//TRIM(I2S(nVars)),Level=3)

    ! MPI coupling
    !------------------------------------------------------------------------
    myLocalRank   = ParEnv % MyPE
    totLocalRanks = ParEnv % PEs
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myGlobalRank, ierr )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, totGlobalRanks, ierr )

    totOFRanks = totGlobalRanks - totLocalRanks

    IF(totOFRanks==0) THEN
      CALL Fatal('OpenFOAM2ElmerSolver','MPI communicator does not have OpenFOAM procs!')
    END IF

    CALL MPI_ALLREDUCE(myGlobalRank, ElmerRanksStart, 1, MPI_INTEGER, MPI_MIN, ELMER_COMM_WORLD, ierr)

    IF (ElmerRanksStart==0) THEN
       OFRanksStart = totLocalRanks
    ELSE
       OFRanksStart = 0
    END IF

    CALL Info('OpenFOAM2ElmerSolver', 'Allocating OpenFOAM data structures',Level=3)
    ALLOCATE( OFp(0:totOFRanks-1) )

    DO i=0,totOFRanks-1
      ALLOCATE( OFp(i) % OFMesh )
      OFp(i) % globalRank = i + OFRanksStart
    END DO


    ! Count elements that need to be communicated
    !------------------------------------------------------------------------
    nElements = 0
    VarName = ListGetString( Params, 'Target Variable 1', Found )
    ptrVar => ListFind( Material, VarName, Found )
    DO i = 1, GetNOFActive()
      Element => GetActiveElement(i)
      IF (commBody(Element % BodyId)) nElements = nElements + 1
    END DO

    ! Get element Id's and cell centres
    !------------------------------------------------------------------------
    IF (nElements>0) THEN
      ALLOCATE( commElementId(nElements), &
                commElementX(nElements),  &
                commElementY(nElements),  &
                commElementZ(nElements), &
                elemPerNode(CurrentModel % Mesh % NumberOfNodes) )
      elemPerNode = 0

      j = 0
      DO i = 1, GetNOFActive()
        Element => GetActiveElement(i)
        IF (commBody(Element % BodyId)) THEN
          j = j+1
          commElementId(j) = i
          n = Element % TYPE % NumberOfNodes
          commElementX(j) = SUM( CurrentModel % Mesh % Nodes % x( Element % NodeIndexes ) ) / n
          commElementY(j) = SUM( CurrentModel % Mesh % Nodes % y( Element % NodeIndexes ) ) / n
          commElementZ(j) = SUM( CurrentModel % Mesh % Nodes % z( Element % NodeIndexes ) ) / n
          elemPerNode( Element % NodeIndexes ) = elemPerNode( Element % NodeIndexes ) +1
        END IF
      END DO
    END IF

    ! Starting communication
    !------------------------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

    DO i = 0, totOFRanks - 1
      ! Number of Elmer elements
       CALL MPI_ISEND(nElements, 1, MPI_INTEGER, &
                      OFp(i) % globalRank, 899, MPI_COMM_WORLD, OFp(i) % reqSend, ierr)
    END DO

    DO i = 0, totOFRanks - 1
      ! Number of Elmer elements
      CALL MPI_WAIT( OFp(i) % reqSend, MPI_STATUS_IGNORE, ierr )
      ! Element X coords
      CALL MPI_ISEND(commElementX, nElements, MPI_DOUBLE, &
                     OFp(i) % globalRank, 899, MPI_COMM_WORLD, OFp(i) % reqSend, ierr)
    END DO

    DO i = 0, totOFRanks - 1
      ! Element X coords
      CALL MPI_WAIT( OFp(i) % reqSend, MPI_STATUS_IGNORE, ierr )
      ! Element Y coords
      CALL MPI_ISEND(commElementY, nElements, MPI_DOUBLE, &
                     OFp(i) % globalRank, 899, MPI_COMM_WORLD, OFp(i) % reqSend, ierr)
    END DO

    DO i = 0, totOFRanks - 1
      ! Element Y coords
      CALL MPI_WAIT( OFp(i) % reqSend, MPI_STATUS_IGNORE, ierr )
      ! Element Z coords
      CALL MPI_ISEND(commElementZ, nElements, MPI_DOUBLE, &
                     OFp(i) % globalRank, 899, MPI_COMM_WORLD, OFp(i) % reqSend, ierr)
    END DO

    DO i = 0, totOFRanks - 1
      ! Element Z coords
      CALL MPI_WAIT( OFp(i) % reqSend, MPI_STATUS_IGNORE, ierr )
      ! Number of found elements
      CALL MPI_IRECV(OFp(i) % nFoundElements, 1, MPI_INTEGER, &
                     OFp(i) % globalRank, 899, MPI_COMM_WORLD, OFp(i) % reqRecv, ierr)
    END DO

    totElementsFound = 0
    DO i = 0, totOFRanks - 1
      ! Number of found elements
      CALL MPI_WAIT( OFp(i) % reqRecv, MPI_STATUS_IGNORE, ierr )
      IF (OFp(i) % nFoundElements>0) THEN
        totElementsFound = totElementsFound + OFp(i) % nFoundElements
        ALLOCATE ( OFp(i) % foundElementIndx(OFp(i) % nFoundElements), &
                   OFp(i) % recvValues(OFp(i) % nFoundElements))
        ! Indexes of found elements
        CALL MPI_IRECV(OFp(i) % foundElementIndx, OFp(i) % nFoundElements, MPI_INTEGER, &
                       OFp(i) % globalRank, 899, MPI_COMM_WORLD, OFp(i) % reqRecv, ierr)
      END IF
    END DO

    IF (totElementsFound .NE. nElements) THEN
      CALL Fatal('OpenFOAM2ElmerSolver','Elmer #'//TRIM(I2S(myLocalRank))//' has ' &
                 //TRIM(I2S(nElements))//' elements, OpenFOAM found '//TRIM(I2S(totElementsFound)))
    END IF

    DO i = 0, totOFRanks - 1
      IF (OFp(i) % nFoundElements>0) THEN
        ! Indexes of found elements
        CALL MPI_WAIT( OFp(i) % reqRecv, MPI_STATUS_IGNORE, ierr )
        OFp(i) % foundElementIndx = OFp(i) % foundElementIndx + 1
      END IF
    END DO

  END IF ! .NOT. VISITED

  ! Receive simulation status
  CALL MPI_IRECV( OFstatus, 1, MPI_INTEGER, OFp(0) % globalRank, 799, MPI_COMM_WORLD, OFp(0) % reqRecv, ierr)
  CALL MPI_WAIT( OFp(0) % reqRecv, MPI_STATUS_IGNORE, ierr )
  IF (OFstatus.NE.1) THEN
    CALL Info('OpenFOAM2ElmerSolver','Elmer has last iteration!', Level=3 )
    exitcond = ListGetCReal( CurrentModel % Simulation,'Exit Condition',Found)
	IF(.NOT.Found) CALL ListAddConstReal(CurrentModel % Simulation,'Exit Condition',1.0_dp)
  END IF

  ! Receive fields
  DO j=1,nVars
    VarName = OFVarNames(j)
    Var => VariableGet(CurrentModel % Mesh % Variables, VarName )
    Var % Values = 0

    DO i = 0, totOFRanks - 1
      IF ( OFp(i) % nFoundElements > 0 ) THEN
        CALL MPI_IRECV( OFp(i) % recvValues, OFp(i) % nFoundElements, MPI_DOUBLE, &
                        OFp(i) % globalRank, 900, MPI_COMM_WORLD, OFp(i) % reqRecv, ierr)
      END IF
    END DO

    DO i = 0, totOFRanks - 1
      IF ( OFp(i) % nFoundElements > 0 ) THEN
        CALL MPI_WAIT( OFp(i) % reqRecv, MPI_STATUS_IGNORE, ierr)
        DO k = 1, OFp(i) % nFoundElements
          Element => GetActiveElement(commElementId(OFp(i) % foundElementIndx(k)))
          Var % Values(Var % Perm(Element % NodeIndexes)) = Var % Values(Var % Perm(Element % NodeIndexes))&
                 + OFp(i) % recvValues(k) / elemPerNode(Element % NodeIndexes)
        END DO
      END IF
    END DO
  END DO

  VISITED = .TRUE.

  CALL Info('OpenFOAM2ElmerSolver','All done', Level=4 )
  CALL Info('OpenFOAM2ElmerSolver','-----------------------------------------', Level=4 )
  
!------------------------------------------------------------------------------  

END SUBROUTINE OpenFOAM2ElmerSolver
