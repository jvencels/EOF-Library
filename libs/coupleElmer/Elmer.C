/*---------------------------------------------------------------------------*\
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

Authors:   Juris Vencels (EOF Consulting, Latvia)

Web:       http://eof-library.com

\*---------------------------------------------------------------------------*/

#include "Elmer.H"
#include <new>
#include <iostream>
#include <fstream>
#include "interpolation.H"
#if (FOAM_MAJOR_VERSION == v1812 || FOAM_MAJOR_VERSION < 10)
#include "dynamicFvMesh.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class meshT>
Foam::Elmer<meshT>::Elmer(const meshT& mesh, int mode, bool init, bool multiregion)
:
mesh_(mesh),
mode_(mode),
myBoundBox(mesh.points(),false)
{
    multiregion_ = multiregion;
    if(init) initialize();
}


template <class meshT>
void Foam::Elmer<meshT>::initialize()
{
    int i, j, k;
    double commTime = MPI_Wtime();
    double localTime = 0;

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

    ELp = new ElmerProc_t[totElmerRanks];

    if (ELp == nullptr) {
        FatalErrorInFunction << "Failed to allocate memory" << Foam::abort(FatalError);
    }

    for ( i=0; i<totElmerRanks; i++ ) {
        ELp[i].globalRank = i+ElmerRanksStart;
    }

    OFboundBoxes = new double[totLocalRanks*2*3];
    ELboundBoxes = new double[totElmerRanks*2*3];
    OF_EL_overlap = new int[totLocalRanks*totElmerRanks];
    findOverlappingBoxes();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Receiving from Elmer
    if (mode_==-1) {
        nCells = mesh_.cells().size();

        cellCentres_x = new double[nCells];
        cellCentres_y = new double[nCells];
        cellCentres_z = new double[nCells];

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
            ELp[i].foundCellsIndx = new int[ELp[i].nFoundCells];
            ELp[i].recvBuffer0 = new double[ELp[i].nFoundCells];

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
        const dictionary& fvSchemes = mesh_.template lookupObject<IOdictionary>
        (
           "fvSchemes"
        );

        // Exctract subdictionary from the main dictionary
        interpolationDict = fvSchemes.subDict("interpolationSchemes");

        int nCommElem = 0;

        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Irecv(&ELp[i].nElem, 1, MPI_INT, ELp[i].globalRank, 899, MPI_COMM_WORLD, &ELp[i].reqRecv);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqRecv);

            nCommElem += ELp[i].nElem;

            ELp[i].sendBuffer0 = new double[ELp[i].nElem];
            ELp[i].sendBuffer1 = new double[ELp[i].nElem];
            ELp[i].sendBuffer2 = new double[ELp[i].nElem];
            ELp[i].foundElement = new label[ELp[i].nElem];

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

        int nElemDone = 0;
        int nElemDonePrev = 0;

        Info<< "Searching for cells.." << endl;
        for ( i=0; i<totElmerRanks; i++ ) {
            ELp[i].nFoundElements = 0; // keep this
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqRecv);

            string fname = "O2E-" + std::to_string(myLocalRank) + "-" + std::to_string(i) + ".out";
            std::ifstream f(fname);

            if(f.good() && !multiregion_) {
                Info<< "O2E pair file exist.." << endl;
                for ( j=0; j<ELp[i].nElem; j++ ) {
                    f >> ELp[i].foundElement[j];
                    if (ELp[i].foundElement[j] > -1) ELp[i].nFoundElements++;
                }
            }
            else {
                Info<< "O2E pair file does not exist, creating.." << endl;

                f.close();
                std::ofstream f(fname);

                for ( j=0; j<ELp[i].nElem; j++ ) {
                    point tmpPoint(ELp[i].sendBuffer0[j],ELp[i].sendBuffer1[j],ELp[i].sendBuffer2[j]);

                    nElemDone++;
                    if(MPI_Wtime()-localTime>30 || nElemDone==nCommElem) {
                        Pout << 100.0*nElemDone/nCommElem << "% done, search speed "
                             << 2*(nElemDone-nElemDonePrev) << " points/min, remaining "
                             << nCommElem-nElemDone << " points" << endl;
                        localTime = MPI_Wtime();
                        nElemDonePrev = nElemDone;
                    }

#if (FOAM_MAJOR_VERSION == 2)
#warning "You are using OF v2.x.x. findCell uses slow search algorithm, be patient!"
                    ELp[i].foundElement[j] = mesh_.findCell(tmpPoint,polyMesh::FACEPLANES);
#else
                    ELp[i].foundElement[j] = mesh_.findCell(tmpPoint);
#endif

                    if (ELp[i].foundElement[j] > -1) ELp[i].nFoundElements++;

                    if (!multiregion_) f << ELp[i].foundElement[j] << std::endl;
                }
            }

            f.close();

            //Pout<< "Found " << ELp[i].nFoundElements << " elements from Elmer #" << i << endl;
            MPI_Isend(&ELp[i].nFoundElements, 1, MPI_INT, ELp[i].globalRank, 895,
                      MPI_COMM_WORLD, &ELp[i].reqSend);
        }
        for ( i=0; i<totElmerRanks; i++ ) {
            if (!ELp[i].boxOverlap) continue;
            MPI_Test_Sleep(ELp[i].reqSend);

            if (ELp[i].nFoundElements == 0) continue;
            ELp[i].foundElementIndx = new int[ELp[i].nFoundElements];
            ELp[i].foundElementCellIndx = new int[ELp[i].nFoundElements];
            ELp[i].positions = new point[ELp[i].nFoundElements];

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

        MPI_Barrier(PstreamGlobals::MPI_COMM_FOAM);
        Info<< "OpenFOAM2Elmer Init = " << MPI_Wtime()-commTime << " s" << nl << endl;
    }
}


template <class meshT>
void Foam::Elmer<meshT>::recvScalar(volScalarField& field)
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


template <class meshT>
void Foam::Elmer<meshT>::recvVector(volVectorField& field)
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


template <class meshT>
void Foam::Elmer<meshT>::recvSymmTensor(volSymmTensorField& field)
{
    int i, j, dim;

    for (dim=0; dim<6; dim++) { 
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


template <class meshT>
void Foam::Elmer<meshT>::recvDiagTensor(volTensorField& field)
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
                field[ELp[i].foundCellsIndx[j]].component(4*dim) = ELp[i].recvBuffer0[j];
            }
        }
    }
}


template <class meshT>
void Foam::Elmer<meshT>::sendScalar(volScalarField& field)
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


template <class meshT>
void Foam::Elmer<meshT>::sendVector(volVectorField& field)
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


template <class meshT>
void Foam::Elmer<meshT>::sendStatus(int status)
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


template <class meshT>
void Foam::Elmer<meshT>::findOverlappingBoxes()
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


template <class meshT>
void Foam::Elmer<meshT>::MPI_Test_Sleep(MPI_Request& req)
{
    int flag;

    while ( true ) {
        MPI_Test( &req, &flag, MPI_STATUS_IGNORE);
        if (flag) break;
        //nanosleep((const struct timespec[]){{0, 1000000L}}, NULL);
        struct timespec tim;
        tim.tv_sec = 0;
        tim.tv_nsec = 1000000L;
        nanosleep(&tim, NULL);
    }
}

// explicit instantiations
template class Elmer<fvMesh>;
#if (FOAM_MAJOR_VERSION == v1812 || FOAM_MAJOR_VERSION < 10)
template class Elmer<dynamicFvMesh>;
#endif

// ************************************************************************* //
