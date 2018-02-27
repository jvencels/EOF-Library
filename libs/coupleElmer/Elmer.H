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

Original Date: 12.12.2016.

\*---------------------------------------------------------------------------*/

#ifndef Elmer_H
#define Elmer_H

// Stuff for compatibility with OF2.4
#ifndef FatalErrorInFunction
#ifdef __GNUC__
   #define FUNCTION_NAME __PRETTY_FUNCTION__
#else
   #define FUNCTION_NAME __func__
#endif
 //- Report an error message using Foam::FatalError
 //  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define FatalErrorInFunction FatalErrorIn(FUNCTION_NAME)
#endif

#include "fvCFD.H"
#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class Elmer Declaration
\*---------------------------------------------------------------------------*/

class Elmer
{

public:

    // Public data

        //- Construct and initialize
        Elmer(const fvMesh& mesh, int mode);

        //- Send interpolated scalar field at element centres to Elmer
        void sendScalar(volScalarField& field);

        //- Send interpolated vector field at element centres to Elmer
        void sendVector(volVectorField& field);

        //- Get interpolated scalar field at cell centres from Elmer
        void recvScalar(volScalarField& field);

        //- Get interpolated vector field at cell centres from Elmer
        void recvVector(volVectorField& field);

        //- Simulation status (1=ok, 0=lastIter, -1=error)
        void sendStatus(int status);

        //- Get bound box for local mesh partition
        void getboundBox(const fvMesh& mesh);

        //- Compare bounding box of this OF partition with all EL ranks
        void findOverlappingBoxes();

        //- Test for incoming message that consumes less CPU
        void MPI_Test_Sleep(MPI_Request& req);

        //- Data structure for each Elmer MPI rank
        class ElmerProc_t
        {
            public:
                // MPI stuff
		        MPI_Request reqSend, reqRecv;
		        int globalRank;

                // Receiving fields from Elmer
		        int *foundCellsIndx;
		        int nFoundCells;
				double *recvBuffer0;

                // Sending fields to Elmer
                int nElem;
                double *sendBuffer0, *sendBuffer1, *sendBuffer2;
                label *foundElement;
                int *foundElementIndx;
                int *foundElementCellIndx;
		        int nFoundElements;
                point *positions;

                // Bounding box
                bool boxOverlap;
        };

private:

    // Private data

        //- Number of cells for this OpenFOAM MPI rank
        const fvMesh& mesh_;

        //- Coupling mode: 1=send, -1=receive, 0=both
        const int mode_;

        //- Number of cells for this OpenFOAM MPI rank
        int nCells;

        //- My rank's cell centres
        double *cellCentres_x, *cellCentres_y, *cellCentres_z;

        //- Total global MPI ranks
        int totGlobalRanks;

        //- My rank at MPI_COMM_WORLD
        int myGlobalRank;

        //- Total OpenFOAM ranks
        int totLocalRanks;

        //- My rank at MPI_OF_WORLD
        int myLocalRank;

        //- Total ELmer ranks
        int totElmerRanks;
    
        //- Elmer ranks start in MPI_COMM_WORLD
        int ElmerRanksStart;

        //- OpenFOAM ranks start in MPI_COMM_WORLD
        int OFRanksStart;

        //- Data structure for each Elmer MPI rank
        ElmerProc_t *ELp;

        // Interpolation types used for interpolating Elmer's elemental field
        dictionary interpolationDict;

        // Bounding box for this OF rank
        boundBox myBoundBox;

        //- Bounding boxes for all OF ranks [rank][min,max][x,y,z]
        double *OFboundBoxes;

        //- Bounding boxes for all EL ranks [rank][min,max][x,y,z]
        double *ELboundBoxes;

        //- Boolean matrix for storing overlappings [OFranks][ELranks]
        int *OF_EL_overlap;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
