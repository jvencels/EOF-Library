/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2025 OpenCFD Ltd.
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

#include "Pstream.H"
#include "Map.H"
#include "UPstreamWrapping.H"

#include <cinttypes>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  Pstream_CommonRoutines
#define Pstream_CommonRoutines(Native, TaggedType)                            \
void Foam::UPstream::allToAll                                                 \
(                                                                             \
    const UList<Native>& sendData,                                            \
    UList<Native>& recvData,                                                  \
    const int communicator                                                    \
)                                                                             \
{                                                                             \
    PstreamDetail::allToAll                                                   \
    (                                                                         \
        sendData, recvData, TaggedType, communicator                          \
    );                                                                        \
}


Pstream_CommonRoutines(int32_t, MPI_INT32_T);
Pstream_CommonRoutines(int64_t, MPI_INT64_T);
// Future?
// Pstream_CommonRoutines(uint32_t, MPI_UINT32_T);
// Pstream_CommonRoutines(uint64_t, MPI_UINT64_T);

#undef Pstream_CommonRoutines


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  Pstream_CommonRoutines
#define Pstream_CommonRoutines(Native, TaggedType)                            \
void Foam::UPstream::allToAllConsensus                                        \
(                                                                             \
    const UList<Native>& sendData,                                            \
    UList<Native>& recvData,                                                  \
    const int tag,                                                            \
    const int communicator                                                    \
)                                                                             \
{                                                                             \
    PstreamDetail::allToAllConsensus                                          \
    (                                                                         \
        sendData, recvData, TaggedType, tag, communicator                     \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::UPstream::allToAllConsensus                                        \
(                                                                             \
    const Map<Native>& sendData,                                              \
    Map<Native>& recvData,                                                    \
    const int tag,                                                            \
    const int communicator                                                    \
)                                                                             \
{                                                                             \
    PstreamDetail::allToAllConsensus                                          \
    (                                                                         \
        sendData, recvData, TaggedType, tag, communicator                     \
    );                                                                        \
}


Pstream_CommonRoutines(int32_t, MPI_INT32_T);
Pstream_CommonRoutines(int64_t, MPI_INT64_T);

#undef Pstream_CommonRoutines


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  Pstream_CommonRoutines
#define Pstream_CommonRoutines(Native, TaggedType)                            \
void Foam::UPstream::allToAllv                                                \
(                                                                             \
    const Native* sendData,                                                   \
    const UList<int>& sendCounts,                                             \
    const UList<int>& sendOffsets,                                            \
    Native* recvData,                                                         \
    const UList<int>& recvCounts,                                             \
    const UList<int>& recvOffsets,                                            \
    const int communicator                                                    \
)                                                                             \
{                                                                             \
    PstreamDetail::allToAllv                                                  \
    (                                                                         \
        sendData, sendCounts, sendOffsets,                                    \
        recvData, recvCounts, recvOffsets,                                    \
        TaggedType, communicator                                              \
    );                                                                        \
}

// Unused: Pstream_CommonRoutines(char, MPI_BYTE);

#undef Pstream_CommonRoutines

// ************************************************************************* //
