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

\*---------------------------------------------------------------------------*/

#include "Elmer.H"
#include <new>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Elmer::Elmer(const fvMesh& mesh)
:
mesh_(mesh)
{
    int i;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &totGlobalRanks);

    MPI_Comm_size(Foam::PstreamGlobals::MPI_OF_WORLD, &totLocalRanks);
    MPI_Comm_rank(Foam::PstreamGlobals::MPI_OF_WORLD, &myLocalRank);

    totElmerRanks = totGlobalRanks-totLocalRanks;

    MPI_Allreduce(&myGlobalRank, &OFRanksStart, 1, MPI_INTEGER, MPI_MIN, Foam::PstreamGlobals::MPI_OF_WORLD);

    if (OFRanksStart==0) {
        ElmerRanksStart = totLocalRanks;
    } else {
        ElmerRanksStart = 0;
    }

    nCells = mesh_.cells().size();

    ELp = new (std::nothrow) ElmerProc_t[totElmerRanks];

    if (ELp == nullptr) {
        FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError); 
    }

    cellCentres_x = new (std::nothrow) double[nCells];
    cellCentres_y = new (std::nothrow) double[nCells];
    cellCentres_z = new (std::nothrow) double[nCells];

    if (cellCentres_x == nullptr || cellCentres_y == nullptr || cellCentres_z == nullptr) {
        FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError); 
    } else {
        forAll(mesh_.cells(),cellI) 
        { 
            cellCentres_x[cellI] = mesh_.C()[cellI].component(0);
            cellCentres_y[cellI] = mesh_.C()[cellI].component(1);
            cellCentres_z[cellI] = mesh_.C()[cellI].component(2);
        }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    MPI_Barrier(MPI_COMM_WORLD);

    Info<< "Sending data to Elmer.." << endl;

    for ( i=0; i<totElmerRanks; i++ ) {
        ELp[i].globalRank = i+ElmerRanksStart;
        MPI_Send(&nCells, 1, MPI_INTEGER, ELp[i].globalRank, 999, MPI_COMM_WORLD);
    }

    for ( i=0; i<totElmerRanks; i++ ) {
        MPI_Send(cellCentres_x, nCells, MPI_DOUBLE, ELp[i].globalRank, 999, MPI_COMM_WORLD);
    }
    for ( i=0; i<totElmerRanks; i++ ) {
        MPI_Send(cellCentres_y, nCells, MPI_DOUBLE, ELp[i].globalRank, 999, MPI_COMM_WORLD);
    }
    for ( i=0; i<totElmerRanks; i++ ) {
        MPI_Send(cellCentres_z, nCells, MPI_DOUBLE, ELp[i].globalRank, 999, MPI_COMM_WORLD);
    }

    int totCellsFound = 0;
    for ( i=0; i<totElmerRanks; i++ ) {
        MPI_Recv(&ELp[i].nFoundCells, 1, MPI_INTEGER, ELp[i].globalRank, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        totCellsFound += ELp[i].nFoundCells;
    }

    if (totCellsFound != nCells) {
        FatalErrorInFunction << "Wrong number of received cells" << Foam::abort(FatalError); 
    }

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundCells > 0 ) {
            ELp[i].foundCellsIndx = new (std::nothrow) int[ELp[i].nFoundCells];
            ELp[i].recvValues = new (std::nothrow) double[ELp[i].nFoundCells];

            if (ELp[i].foundCellsIndx == nullptr || ELp[i].recvValues == nullptr) {
                FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError); 
            }

            MPI_Recv(ELp[i].foundCellsIndx, ELp[i].nFoundCells, MPI_INTEGER,
                      ELp[i].globalRank, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}


void Foam::Elmer::getScalar(volScalarField& field)
{
    int i, j;

    Info<< "Receiving scalar field from Elmer.." << endl;

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundCells > 0 ) {
            MPI_Recv(ELp[i].recvValues, ELp[i].nFoundCells, MPI_DOUBLE, ELp[i].globalRank, 
                      1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (j=0; j<ELp[i].nFoundCells; j++ ) {
                field[ELp[i].foundCellsIndx[j]] = ELp[i].recvValues[j];
            }
        }
    }
}


void Foam::Elmer::getVector(volVectorField& field)
{
    int i, j, dim;

    Info<< "Receiving vector field from Elmer.." << endl;

    for (dim=0; dim<3; dim++) { 
        for ( i=0; i<totElmerRanks; i++ ) {
            if ( ELp[i].nFoundCells > 0 ) {
                MPI_Recv(ELp[i].recvValues, ELp[i].nFoundCells, MPI_DOUBLE, ELp[i].globalRank, 
                          1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (j=0; j<ELp[i].nFoundCells; j++ ) {
                    field[ELp[i].foundCellsIndx[j]].component(dim) = ELp[i].recvValues[j];
                }
            }
        }
    }
}

// ************************************************************************* //
