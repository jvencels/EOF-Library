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
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Elmer::Elmer(const fvMesh& mesh, int mode)
:
mesh_(mesh),
mode_(mode)
{
    int i, j, k, flag;

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

    ELp = new (std::nothrow) ElmerProc_t[totElmerRanks];

    if (ELp == nullptr) {
        FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError); 
    }

    for ( i=0; i<totElmerRanks; i++ ) {
        ELp[i].globalRank = i+ElmerRanksStart;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Receiving from Elmer
    if (mode_==-1) {
        nCells = mesh_.cells().size();

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

        Info<< "Sending data to Elmer.." << endl;

        for ( i=0; i<totElmerRanks; i++ ) {
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
            while ( true ) {
                MPI_Iprobe(ELp[i].globalRank, 999, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
                if (flag) break;
                nanosleep((const struct timespec[]){{0, 500000000L}}, NULL);
            }
            MPI_Recv(&ELp[i].nFoundCells, 1, MPI_INTEGER, ELp[i].globalRank, 999, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            totCellsFound += ELp[i].nFoundCells;
        }

        if (totCellsFound < nCells) {
            FatalErrorInFunction << "OpenFOAM #" << myLocalRank << " has " << nCells
                                 << " cells, Elmer found " << totCellsFound << Foam::abort(FatalError); 
        }

        for ( i=0; i<totElmerRanks; i++ ) {
            if ( ELp[i].nFoundCells > 0 ) {
                ELp[i].foundCellsIndx = new (std::nothrow) int[ELp[i].nFoundCells];
                ELp[i].recvBuffer0 = new (std::nothrow) double[ELp[i].nFoundCells];

                if (ELp[i].foundCellsIndx == nullptr || ELp[i].recvBuffer0 == nullptr) {
                    FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError); 
                }

                MPI_Recv(ELp[i].foundCellsIndx, ELp[i].nFoundCells, MPI_INTEGER,
                          ELp[i].globalRank, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Sending fields to Elmer
    if (mode_==1) {

        //interpolationDict = mesh.solutionDict().subDict("interpolationSchemes");
        //fvSchemes = mesh.lookupObject<IOdictionary>("fvSchemes");

        // Extract the dictionary from the database
        const dictionary& fvSchemes = mesh.lookupObject<IOdictionary>
        (
           "fvSchemes"
        );

        // Exctract subdictionary from the main dictionary
        interpolationDict = fvSchemes.subDict("interpolationSchemes");

        for ( i=0; i<totElmerRanks; i++ ) {
            while ( true ) {
                MPI_Iprobe(ELp[i].globalRank, 899, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
                if (flag) break;
                nanosleep((const struct timespec[]){{0, 500000000L}}, NULL);
            }
            MPI_Recv(&ELp[i].nElem, 1, MPI_INTEGER, ELp[i].globalRank, 899, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            Pout<< "Received " << ELp[i].nElem << " elements from Elmer #" << i << endl;

            ELp[i].sendBuffer0 = new (std::nothrow) double[ELp[i].nElem];
            ELp[i].sendBuffer1 = new (std::nothrow) double[ELp[i].nElem];
            ELp[i].sendBuffer2 = new (std::nothrow) double[ELp[i].nElem];
            ELp[i].foundElement = new (std::nothrow) label[ELp[i].nElem];
            ELp[i].positions = new (std::nothrow) point[ELp[i].nElem];
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            MPI_Recv(ELp[i].sendBuffer0, ELp[i].nElem, MPI_DOUBLE, ELp[i].globalRank, 899, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            MPI_Recv(ELp[i].sendBuffer1, ELp[i].nElem, MPI_DOUBLE, ELp[i].globalRank, 899, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            MPI_Recv(ELp[i].sendBuffer2, ELp[i].nElem, MPI_DOUBLE, ELp[i].globalRank, 899, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        Info<< "Searching for cells.." << endl;
        for ( i=0; i<totElmerRanks; i++ ) {
            ELp[i].nFoundElements = 0;
            for ( j=0; j<ELp[i].nElem; j++ ) {
                ELp[i].positions[j].x() = ELp[i].sendBuffer0[j];
                ELp[i].positions[j].y() = ELp[i].sendBuffer1[j];
                ELp[i].positions[j].z() = ELp[i].sendBuffer2[j];
                ELp[i].foundElement[j] = mesh.findCell(ELp[i].positions[j]);
                if (ELp[i].foundElement[j] != -1) ELp[i].nFoundElements++;
            }
            Pout<< "Found " << ELp[i].nFoundElements << " elements from Elmer #" << i << endl;
            ELp[i].foundElementIndx = new (std::nothrow) int[ELp[i].nFoundElements];
            MPI_Send(&ELp[i].nFoundElements, 1, MPI_INTEGER, ELp[i].globalRank, 899, MPI_COMM_WORLD);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            if (ELp[i].nFoundElements > 0) {
                k = 0;
                for ( j=0; j<ELp[i].nElem; j++ ) {
                    if (ELp[i].foundElement[j] != -1) {
                        ELp[i].foundElementIndx[k] = j;
                        k++;
                    }
                }
                MPI_Send(ELp[i].foundElementIndx, ELp[i].nFoundElements, MPI_INTEGER, 
                         ELp[i].globalRank, 899, MPI_COMM_WORLD);
            }
        }
    }
}


void Foam::Elmer::recvScalar(volScalarField& field)
{
    int i, j, flag;

    Info<< "Receiving scalar field from Elmer.." << endl;

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundCells > 0 ) {
            while ( true ) {
                MPI_Iprobe(ELp[i].globalRank, 1000, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
                if (flag) break;
                nanosleep((const struct timespec[]){{0, 500000000L}}, NULL);
            }
            MPI_Recv(ELp[i].recvBuffer0, ELp[i].nFoundCells, MPI_DOUBLE, ELp[i].globalRank, 
                      1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (j=0; j<ELp[i].nFoundCells; j++ ) {
                field[ELp[i].foundCellsIndx[j]] = ELp[i].recvBuffer0[j];
            }
        }
    }
}


void Foam::Elmer::recvVector(volVectorField& field)
{
    int i, j, dim, flag;

    Info<< "Receiving vector field from Elmer.." << endl;

    for (dim=0; dim<3; dim++) { 
        for ( i=0; i<totElmerRanks; i++ ) {
            if ( ELp[i].nFoundCells > 0 ) {
                while ( true ) {
                    MPI_Iprobe(ELp[i].globalRank, 1000, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
                    if (flag) break;
                    nanosleep((const struct timespec[]){{0, 500000000L}}, NULL);
                }
                MPI_Recv(ELp[i].recvBuffer0, ELp[i].nFoundCells, MPI_DOUBLE, ELp[i].globalRank, 
                          1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (j=0; j<ELp[i].nFoundCells; j++ ) {
                    field[ELp[i].foundCellsIndx[j]].component(dim) = ELp[i].recvBuffer0[j];
                }
            }
        }
    }
}


void Foam::Elmer::sendScalar(volScalarField& field)
{
    int i, j;

    autoPtr<interpolation<scalar> > interpField = interpolation<scalar>::New(interpolationDict, field);

    Info<< "Sending scalar field to Elmer.." << endl;

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundElements > 0 ) {
            for (j=0; j<ELp[i].nFoundElements; j++) {
                ELp[i].sendBuffer0[j] = interpField->
                       interpolate(ELp[i].positions[j], ELp[i].foundElement[ELp[i].foundElementIndx[j]]);
            }
            MPI_Send(ELp[i].sendBuffer0, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank, 
                      900, MPI_COMM_WORLD);
        }
    }
}


void Foam::Elmer::sendVector(volVectorField& field)
{
    int i, j;

    autoPtr<interpolation<vector> > interpField = interpolation<vector>::New(interpolationDict, field);

    Info<< "Sending vector field to Elmer.." << endl;

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundElements > 0 ) {
            for (j=0; j<ELp[i].nFoundElements; j++) {
                vector tmpVector = interpField->
                       interpolate(ELp[i].positions[j], ELp[i].foundElement[ELp[i].foundElementIndx[j]]);
                ELp[i].sendBuffer0[j] = tmpVector.component(0);
                ELp[i].sendBuffer1[j] = tmpVector.component(1);
                ELp[i].sendBuffer2[j] = tmpVector.component(2);
            }
        }
    }

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundElements > 0 ) {
            MPI_Send(ELp[i].sendBuffer0, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank,
                       900, MPI_COMM_WORLD);
            MPI_Send(ELp[i].sendBuffer1, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank,
                       900, MPI_COMM_WORLD);
            MPI_Send(ELp[i].sendBuffer2, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank,
                       900, MPI_COMM_WORLD);
        }
    }
}


void Foam::Elmer::sendStatus(int status)
{
    int i;

    if (myLocalRank==0) {
        for ( i=0; i<totElmerRanks; i++ ) {
            MPI_Send(&status, 1, MPI_INTEGER, ELp[i].globalRank,799, MPI_COMM_WORLD);
        }
    }
}
// ************************************************************************* //
