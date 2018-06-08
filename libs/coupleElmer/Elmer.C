/*---------------------------------------------------------------------------*\
Copyright (C) 2016-2018 University of Latvia
-------------------------------------------------------------------------------
License
    This file is part of EOF-Library.

    EOF-Library is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EOF-Library is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with EOF-Library.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Elmer

Description
    Coupling with FEM open source multiphysical simulation software Elmer

SourceFiles
    Elmer.C

Authors:   Juris Vencels
Email:     juris.vencels@lu.lv
Web:       http://eof-library.com
Address:   University of Latvia
           Laboratory for mathematical modelling of 
               environmental and technological processes
           Zellu Str. 23, Riga, LV-1002, Latvia

\*---------------------------------------------------------------------------*/

#include "Elmer.H"
#include <new>
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Elmer::Elmer(const fvMesh& mesh, int mode)
:
mesh_(mesh),
mode_(mode),
myBoundBox(mesh.points(),false)
{
    initialize();
}


Foam::Elmer::Elmer(const fvMesh& mesh, int mode, bool init)
:
mesh_(mesh),
mode_(mode),
myBoundBox(mesh.points(),false)
{
    if(init) initialize();
}


void Foam::Elmer::initialize()
{
    int i, j, k, flag;
    double commTime = MPI_Wtime();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &totGlobalRanks);

    MPI_Comm_size(PstreamGlobals::MPI_COMM_FOAM, &totLocalRanks);
    MPI_Comm_rank(PstreamGlobals::MPI_COMM_FOAM, &myLocalRank);

    totElmerRanks = totGlobalRanks-totLocalRanks;

    if (myLocalRank==0) OFRanksStart = myGlobalRank;

    MPI_Bcast(&OFRanksStart, 1, MPI_INT, 0, PstreamGlobals::MPI_COMM_FOAM);

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

    OFboundBoxes = new (std::nothrow) double[totLocalRanks*2*3];
    ELboundBoxes = new (std::nothrow) double[totElmerRanks*2*3];
    OF_EL_overlap = new (std::nothrow) int[totLocalRanks*totElmerRanks];
    findOverlappingBoxes();

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
            if (!ELp[i].boxOverlap) continue;
            MPI_Isend(&nCells, 1, MPI_INT, ELp[i].globalRank, 999, MPI_COMM_WORLD, &ELp[i].reqSend);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqSend);
            MPI_Isend(cellCentres_x, nCells, MPI_DOUBLE, ELp[i].globalRank, 997, MPI_COMM_WORLD, &ELp[i].reqSend);
            MPI_Request_free(&ELp[i].reqSend);
            MPI_Isend(cellCentres_y, nCells, MPI_DOUBLE, ELp[i].globalRank, 997, MPI_COMM_WORLD, &ELp[i].reqSend);
            MPI_Request_free(&ELp[i].reqSend);
            MPI_Isend(cellCentres_z, nCells, MPI_DOUBLE, ELp[i].globalRank, 997, MPI_COMM_WORLD, &ELp[i].reqSend);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            ELp[i].nFoundCells = 0; // keep this
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqSend);
            MPI_Irecv(&ELp[i].nFoundCells, 1, MPI_INT, ELp[i].globalRank, 995,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
        }

        int totCellsFound = 0;
        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqRecv);
            totCellsFound += ELp[i].nFoundCells;
        }

        if (totCellsFound < nCells) {
            FatalErrorInFunction << "OpenFOAM #" << myLocalRank << " has " << nCells
                                 << " cells, Elmer found " << totCellsFound << Foam::abort(FatalError); 
        }

        for ( i=0; i<totElmerRanks; i++ ) {
            if ( ELp[i].nFoundCells == 0 ) continue;
            ELp[i].foundCellsIndx = new (std::nothrow) int[ELp[i].nFoundCells];
            ELp[i].recvBuffer0 = new (std::nothrow) double[ELp[i].nFoundCells];

            if (ELp[i].foundCellsIndx == nullptr || ELp[i].recvBuffer0 == nullptr) {
                FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
            }

            MPI_Irecv(ELp[i].foundCellsIndx, ELp[i].nFoundCells, MPI_INT, ELp[i].globalRank, 994,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
        }

        for ( i=0; i<totElmerRanks; i++ ) {
            if ( ELp[i].nFoundCells == 0 ) continue;
            MPI_Test_Sleep(ELp[i].reqRecv);
        }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Sending fields to Elmer
    if (mode_==1) {

        // Extract the dictionary from the database
        const dictionary& fvSchemes = mesh_.lookupObject<IOdictionary>
        (
           "fvSchemes"
        );

        // Exctract subdictionary from the main dictionary
        interpolationDict = fvSchemes.subDict("interpolationSchemes");

        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Irecv(&ELp[i].nElem, 1, MPI_INT, ELp[i].globalRank, 899, MPI_COMM_WORLD, &ELp[i].reqRecv);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqRecv);
            ELp[i].sendBuffer0 = new (std::nothrow) double[ELp[i].nElem];
            ELp[i].sendBuffer1 = new (std::nothrow) double[ELp[i].nElem];
            ELp[i].sendBuffer2 = new (std::nothrow) double[ELp[i].nElem];
            ELp[i].foundElement = new (std::nothrow) label[ELp[i].nElem];

            if (ELp[i].sendBuffer0 == nullptr || ELp[i].sendBuffer1 == nullptr || 
                ELp[i].sendBuffer2 == nullptr || ELp[i].foundElement == nullptr) {
                FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
            }

            MPI_Irecv(ELp[i].sendBuffer0, ELp[i].nElem, MPI_DOUBLE, ELp[i].globalRank, 897,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
            MPI_Request_free(&ELp[i].reqRecv);
            MPI_Irecv(ELp[i].sendBuffer1, ELp[i].nElem, MPI_DOUBLE, ELp[i].globalRank, 897,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
            MPI_Request_free(&ELp[i].reqRecv);
            MPI_Irecv(ELp[i].sendBuffer2, ELp[i].nElem, MPI_DOUBLE, ELp[i].globalRank, 897,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
        }

        Info<< "Searching for cells.." << endl;
        for ( i=0; i<totElmerRanks; i++ ) {
            ELp[i].nFoundElements = 0; // keep this
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqRecv);
            for ( j=0; j<ELp[i].nElem; j++ ) {
                point tmpPoint(ELp[i].sendBuffer0[j],ELp[i].sendBuffer1[j],ELp[i].sendBuffer2[j]);

                ELp[i].foundElement[j] = mesh_.findCell(tmpPoint);
                if (ELp[i].foundElement[j] > -1) ELp[i].nFoundElements++;
            }
            //Pout<< "Found " << ELp[i].nFoundElements << " elements from Elmer #" << i << endl;
            MPI_Isend(&ELp[i].nFoundElements, 1, MPI_INT, ELp[i].globalRank, 895,
                      MPI_COMM_WORLD, &ELp[i].reqSend);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqSend);

            if (ELp[i].nFoundElements == 0) continue;
            ELp[i].foundElementIndx = new (std::nothrow) int[ELp[i].nFoundElements];
            ELp[i].foundElementCellIndx = new (std::nothrow) int[ELp[i].nFoundElements];
            ELp[i].positions = new (std::nothrow) point[ELp[i].nFoundElements];

            if (ELp[i].foundElementIndx == nullptr ||
                ELp[i].foundElementCellIndx == nullptr ||
                ELp[i].positions == nullptr) {
                FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
            }

            k = 0;
            for ( j=0; j<ELp[i].nElem; j++ ) {
                if (ELp[i].foundElement[j] > -1) {
                    ELp[i].foundElementIndx[k] = j;
                    ELp[i].foundElementCellIndx[k] = ELp[i].foundElement[j];
                    ELp[i].positions[k].x() = ELp[i].sendBuffer0[j];
                    ELp[i].positions[k].y() = ELp[i].sendBuffer1[j];
                    ELp[i].positions[k].z() = ELp[i].sendBuffer2[j];
                    k++;
                }
            }
            MPI_Isend(ELp[i].foundElementIndx, ELp[i].nFoundElements, MPI_INT, ELp[i].globalRank, 894,
                      MPI_COMM_WORLD, &ELp[i].reqSend);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            if (ELp[i].nFoundElements == 0) continue;
            MPI_Test_Sleep(ELp[i].reqSend);
        }

        Info<< "OpenFOAM2Elmer Init = " << MPI_Wtime()-commTime << " s" << nl << endl;
    }
}


void Foam::Elmer::recvScalar(volScalarField& field)
{
    int i, j;

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundCells == 0 ) continue;
        MPI_Irecv(ELp[i].recvBuffer0, ELp[i].nFoundCells, MPI_DOUBLE, ELp[i].globalRank, 1000,
                  MPI_COMM_WORLD, &ELp[i].reqRecv);
    }
    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundCells == 0 ) continue;
        MPI_Test_Sleep(ELp[i].reqRecv);
        for (j=0; j<ELp[i].nFoundCells; j++ ) {
            field[ELp[i].foundCellsIndx[j]] = ELp[i].recvBuffer0[j];
        }
    }
}


void Foam::Elmer::recvVector(volVectorField& field)
{
    int i, j, dim;

    for (dim=0; dim<3; dim++) { 
        for ( i=0; i<totElmerRanks; i++ ) {
            if ( ELp[i].nFoundCells == 0 ) continue;
            MPI_Irecv(ELp[i].recvBuffer0, ELp[i].nFoundCells, MPI_DOUBLE, ELp[i].globalRank, 1000,
                      MPI_COMM_WORLD, &ELp[i].reqRecv);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            if ( ELp[i].nFoundCells == 0 ) continue;
            MPI_Test_Sleep(ELp[i].reqRecv);
            for (j=0; j<ELp[i].nFoundCells; j++ ) {
                field[ELp[i].foundCellsIndx[j]].component(dim) = ELp[i].recvBuffer0[j];
            }
        }
    }
}


void Foam::Elmer::sendScalar(volScalarField& field)
{
    int i, j;
    double commTime = MPI_Wtime();

    autoPtr<interpolation<scalar> > interpField = interpolation<scalar>::New(interpolationDict, field);

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundElements == 0 ) continue;
        for (j=0; j<ELp[i].nFoundElements; j++) {
            ELp[i].sendBuffer0[j] = interpField->
                   interpolate(ELp[i].positions[j], ELp[i].foundElementCellIndx[j]);
        }
        MPI_Isend(ELp[i].sendBuffer0, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank, 900,
                  MPI_COMM_WORLD, &ELp[i].reqSend);
    }
    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundElements == 0 ) continue;
        MPI_Test_Sleep(ELp[i].reqSend);
    }

    Info<< "OpenFOAM2Elmer: scalar = " << MPI_Wtime()-commTime << " s" << nl << endl;
}


void Foam::Elmer::sendVector(volVectorField& field)
{
    int i, j;
    double commTime = MPI_Wtime();

    autoPtr<interpolation<vector> > interpField = interpolation<vector>::New(interpolationDict, field);

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundElements == 0 ) continue;
        for (j=0; j<ELp[i].nFoundElements; j++) {
            vector tmpVector = interpField->
                   interpolate(ELp[i].positions[j], ELp[i].foundElementCellIndx[j]);
            ELp[i].sendBuffer0[j] = tmpVector.component(0);
            ELp[i].sendBuffer1[j] = tmpVector.component(1);
            ELp[i].sendBuffer2[j] = tmpVector.component(2);
        }
        MPI_Isend(ELp[i].sendBuffer0, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank, 900,
                 MPI_COMM_WORLD, &ELp[i].reqSend);
        MPI_Request_free(&ELp[i].reqSend);
        MPI_Isend(ELp[i].sendBuffer1, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank, 900,
                 MPI_COMM_WORLD, &ELp[i].reqSend);
        MPI_Request_free(&ELp[i].reqSend);
        MPI_Isend(ELp[i].sendBuffer2, ELp[i].nFoundElements, MPI_DOUBLE, ELp[i].globalRank, 900,
                 MPI_COMM_WORLD, &ELp[i].reqSend);
    }

    for ( i=0; i<totElmerRanks; i++ ) {
        if ( ELp[i].nFoundElements == 0 ) continue;
        MPI_Test_Sleep(ELp[i].reqSend);
    }

    Info<< "OpenFOAM2Elmer: vector = " << MPI_Wtime()-commTime << " s" << nl << endl;
}


void Foam::Elmer::sendStatus(int status)
{
    int i;

    if (myLocalRank==0) {
        for ( i=0; i<totElmerRanks; i++ ) {
            MPI_Isend(&status, 1, MPI_INT, ELp[i].globalRank,799, MPI_COMM_WORLD, &ELp[i].reqSend);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            MPI_Test_Sleep(ELp[i].reqSend);
        }
    }
}


void Foam::Elmer::findOverlappingBoxes()
{
    int OFoffset = myLocalRank*totElmerRanks;

    if ( myLocalRank==0 ) {
        MPI_Recv(ELboundBoxes, totElmerRanks*2*3, MPI_DOUBLE, ELp[0].globalRank, 1001,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Bcast(ELboundBoxes, totElmerRanks*2*3, MPI_DOUBLE, 0, PstreamGlobals::MPI_COMM_FOAM);

    for (int i=0; i<totElmerRanks; i++ ) {
        point tmpPointMin;
        point tmpPointMax;

        int ELoffset = i*2*3;
        tmpPointMin.x() = ELboundBoxes[ELoffset];
        tmpPointMin.y() = ELboundBoxes[ELoffset+1];
        tmpPointMin.z() = ELboundBoxes[ELoffset+2];
        tmpPointMax.x() = ELboundBoxes[ELoffset+3];
        tmpPointMax.y() = ELboundBoxes[ELoffset+4];
        tmpPointMax.z() = ELboundBoxes[ELoffset+5];
        boundBox tmpBoundBox(tmpPointMin, tmpPointMax);

        ELp[i].boxOverlap = myBoundBox.overlaps(tmpBoundBox);        
        OF_EL_overlap[OFoffset+i] = int(ELp[i].boxOverlap);
    }

    MPI_Allgather(&OF_EL_overlap[OFoffset], totElmerRanks, MPI_INT,
        OF_EL_overlap, totElmerRanks, MPI_INT, PstreamGlobals::MPI_COMM_FOAM);

    if ( myLocalRank==0 ) {
        MPI_Send(OF_EL_overlap, totLocalRanks*totElmerRanks, MPI_INT, ELp[0].globalRank, 1002,
                 MPI_COMM_WORLD);
    }
}


void Foam::Elmer::MPI_Test_Sleep(MPI_Request& req)
{
    int flag;

    while ( true ) {
        MPI_Test( &req, &flag, MPI_STATUS_IGNORE);
        if (flag) break;
        nanosleep((const struct timespec[]){{0, 1000000L}}, NULL);
    }
}

// ************************************************************************* //
