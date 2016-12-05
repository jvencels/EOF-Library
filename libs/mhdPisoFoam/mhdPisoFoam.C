/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pisoFoam

Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "PstreamGlobals.H"
#include <new>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    int nCells;
    nCells = mesh.cells().size();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // OpenFOAM communicator
    int totGlobalRanks, myGlobalRank;
    int totLocalRanks, myLocalRank;
    int totElmerRanks, ElmerRanksStart;

    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &totGlobalRanks);
    MPI_Comm_size(Foam::PstreamGlobals::MPI_OF_WORLD, &totLocalRanks);
    MPI_Comm_rank(Foam::PstreamGlobals::MPI_OF_WORLD, &myLocalRank);

    totElmerRanks = totGlobalRanks-totLocalRanks;
    if (myGlobalRank<totLocalRanks) {
        ElmerRanksStart = totLocalRanks;
    } else {
        ElmerRanksStart = 0;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    MPI_Request reqSend, reqRecv;
    int indx;
    int * nCellsPerElmerProc;
    double *cellCentres_x, *cellCentres_y, *cellCentres_z, *recvBuff;

    nCellsPerElmerProc = new (std::nothrow) int[totElmerRanks];
    if (nCellsPerElmerProc == nullptr) {
        FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError); 
    }

    cellCentres_x = new (std::nothrow) double[nCells];
    cellCentres_y = new (std::nothrow) double[nCells];
    cellCentres_z = new (std::nothrow) double[nCells];
    recvBuff      = new (std::nothrow) double[nCells];
    if (cellCentres_x == nullptr || cellCentres_y == nullptr || cellCentres_z == nullptr || recvBuff == nullptr) {
        FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError); 
    }
            
    MPI_Barrier(MPI_COMM_WORLD);

    Info<< "Sending data to Elmer.." << endl;

    for ( indx=0; indx<totElmerRanks; indx++ ) {
        MPI_Isend(&nCells, 1, MPI_INTEGER, indx+ElmerRanksStart, 999, MPI_COMM_WORLD, &reqSend);
    }

    Info<< "Preparing cell centres.." << endl;
    forAll(mesh.cells(),cellI) 
    { 
      cellCentres_x[cellI] = mesh.C()[cellI].component(0);
      cellCentres_y[cellI] = mesh.C()[cellI].component(1);
      cellCentres_z[cellI] = mesh.C()[cellI].component(2);
    }

    Info<< "Sending cell centres to Elmer.." << endl;
    for ( indx=0; indx<totElmerRanks; indx++ ) {
        MPI_Isend(cellCentres_x, nCells, MPI_DOUBLE, indx+ElmerRanksStart, 1000, MPI_COMM_WORLD, &reqSend);
        MPI_Wait(&reqSend, MPI_STATUS_IGNORE);
    }
    for ( indx=0; indx<totElmerRanks; indx++ ) {
        MPI_Isend(cellCentres_y, nCells, MPI_DOUBLE, indx+ElmerRanksStart, 1001, MPI_COMM_WORLD, &reqSend);
        MPI_Wait(&reqSend, MPI_STATUS_IGNORE);
    }
    for ( indx=0; indx<totElmerRanks; indx++ ) {
        MPI_Isend(cellCentres_z, nCells, MPI_DOUBLE, indx+ElmerRanksStart, 1002, MPI_COMM_WORLD, &reqSend);
        MPI_Wait(&reqSend, MPI_STATUS_IGNORE);
    }

    Info<< "Sending to Elmer is done.." << endl;
/*
    for ( indx=0; indx<totElmerRanks; indx++ ) {
        MPI_Isend(&nCells, 3*nCells, MPI_DOUBLE, indx+ElmerRanksStart, 1, MPI_COMM_WORLD, &reqSend);
    }

    for ( indx=0; indx<totElmerRanks; indx++ ) {
        nCellsPerElmerProc[indx] = 0;
        MPI_Irecv(&nCellsPerElmerProc[indx], 1, MPI_INTEGER, indx+ElmerRanksStart, 0, MPI_COMM_WORLD, &reqRecv);
        MPI_Wait(&reqRecv, MPI_STATUS_IGNORE);
    }
*/

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "Receiving JxB from Elmer.." << nl << endl;
        MPI_Barrier(MPI_COMM_WORLD);

        for ( indx=0; indx<totElmerRanks; indx++ ) {
            MPI_Irecv(recvBuff, nCells, MPI_DOUBLE, indx+ElmerRanksStart, 1003, MPI_COMM_WORLD, &reqRecv);
            MPI_Wait(&reqRecv, MPI_STATUS_IGNORE);
            forAll(mesh.cells(),cellI) 
            { 
                JxB[cellI].component(0) = recvBuff[cellI];
            }
        }
        for ( indx=0; indx<totElmerRanks; indx++ ) {
            MPI_Irecv(recvBuff, nCells, MPI_DOUBLE, indx+ElmerRanksStart, 1004, MPI_COMM_WORLD, &reqRecv);
            MPI_Wait(&reqRecv, MPI_STATUS_IGNORE);
            forAll(mesh.cells(),cellI) 
            { 
                JxB[cellI].component(1) = recvBuff[cellI];
            }
        }
        for ( indx=0; indx<totElmerRanks; indx++ ) {
            MPI_Irecv(recvBuff, nCells, MPI_DOUBLE, indx+ElmerRanksStart, 1005, MPI_COMM_WORLD, &reqRecv);
            MPI_Wait(&reqRecv, MPI_STATUS_IGNORE);
            forAll(mesh.cells(),cellI) 
            { 
                JxB[cellI].component(2) = recvBuff[cellI];
            }
        }
        for ( indx=0; indx<totElmerRanks; indx++ ) {
            MPI_Irecv(recvBuff, nCells, MPI_DOUBLE, indx+ElmerRanksStart, 1006, MPI_COMM_WORLD, &reqRecv);
            MPI_Wait(&reqRecv, MPI_STATUS_IGNORE);
            forAll(mesh.cells(),cellI) 
            { 
                JH[cellI] = recvBuff[cellI];
            }
        }
        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            #include "UEqn.H"
            #include "TEqn.H"

            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    delete[] nCellsPerElmerProc;
    delete[] cellCentres_x, cellCentres_y, cellCentres_z;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
