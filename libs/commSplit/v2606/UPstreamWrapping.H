/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

InNamespace
    Foam::PstreamDetail

Description
    Functions to wrap MPI_Bcast, MPI_Allreduce, MPI_Iallreduce etc.

SourceFiles
    UPstreamWrapping.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_UPstreamWrapping_H
#define Foam_UPstreamWrapping_H

#include "openfoam_mpi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

// MPI_Bcast (default from root=0)
// No fail/abort handling
template<class Type>
bool broadcast
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    const int communicator,
    int root = 0
);

// MPI_Reduce, using root=0
template<class Type>
void reduce
(
    const Type* sendData,   // Use nullptr for in-place operation
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const int communicator,             // Communicator
    UPstream::Request* req = nullptr    // Non-null for non-blocking
);

// MPI_Allreduce or MPI_Iallreduce : in-place operation
template<class Type>
void allReduce
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const int communicator,             // Communicator
    UPstream::Request* req = nullptr    // Non-null for non-blocking
);

// MPI_Scan or MPI_Exscan
template<class Type>
void scanReduce
(
    const Type* sendData,
    Type* recvData,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const int communicator,             // Communicator
    const int exclusive                 // Use Exscan instead of Scan
);

// MPI_Alltoall or MPI_Ialltoall with one element per rank
template<class Type>
void allToAll
(
    const UList<Type>& sendData,
    UList<Type>& recvData,
    MPI_Datatype datatype,
    const int communicator,             // Communicator
    UPstream::Request* req = nullptr    // Non-null for non-blocking
);


// MPI_Alltoallv or MPI_Ialltoallv
template<class Type>
void allToAllv
(
    const Type* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,

    Type* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,

    MPI_Datatype datatype,
    const int communicator,             // Communicator
    UPstream::Request* req = nullptr    // Non-null for non-blocking
);


// Non-blocking consensual integer (size) exchange
template<class Type>
void allToAllConsensus
(
    const UList<Type>& sendData,
    UList<Type>& recvData,
    MPI_Datatype datatype,
    const int tag,                      // Message tag
    const int communicator              // Communicator
);


// Non-blocking consensual integer (size) exchange
template<class Type>
void allToAllConsensus
(
    const Map<Type>& sendData,
    Map<Type>& recvData,
    MPI_Datatype datatype,
    const int tag,                      // Message tag
    const int communicator              // Communicator
);


// MPI_Gather or MPI_Igather
// Uses recvData as send/recv when sendData is nullptr
template<class Type>
void gather
(
    const Type* sendData,   // Local send value
    Type* recvData,         // On master: recv buffer. Ignored elsewhere
    int count,              // Per rank send/recv count. Globally consistent!
    MPI_Datatype datatype,              // The send/recv data type
    const int communicator,             // Communicator
    UPstream::Request* req = nullptr    // Non-null for non-blocking
);


// MPI_Scatter or MPI_Iscatter
// Uses recvData as send/recv when sendData is nullptr
template<class Type>
void scatter
(
    const Type* sendData,   // On master: send buffer. Ignored elsewhere
    Type* recvData,         // Local recv value
    int count,              // Per rank send/recv count. Globally consistent!
    MPI_Datatype datatype,              // The send/recv data type
    const int communicator,             // Communicator
    UPstream::Request* req = nullptr    // Non-null for non-blocking
);


// MPI_Gatherv or MPI_Igatherv
template<class Type>
void gatherv
(
    const Type* sendData,
    int sendCount,              // Ignored on master if recvCounts[0] == 0

    Type* recvData,                     // Ignored on non-root rank
    const UList<int>& recvCounts,       // Ignored on non-root rank
    const UList<int>& recvOffsets,      // Ignored on non-root rank

    MPI_Datatype datatype,              // The send/recv data type
    const int communicator,             // Communicator
    UPstream::Request* req = nullptr    // Non-null for non-blocking
);


// MPI_Scatterv or MPI_Iscatterv
template<class Type>
void scatterv
(
    const Type* sendData,               // Ignored on non-root rank
    const UList<int>& sendCounts,       // Ignored on non-root rank
    const UList<int>& sendOffsets,      // Ignored on non-root rank

    Type* recvData,
    int recvCount,

    MPI_Datatype datatype,              // The send/recv data type
    const int communicator,             // Communicator
    UPstream::Request* req = nullptr    // Non-null for non-blocking
);


// MPI_Allgather or MPI_Iallgather
template<class Type>
void allGather
(
    Type* allData,                      // The send/recv data
    int count,                          // The send/recv count per element

    MPI_Datatype datatype,              // The send/recv data type
    const int communicator,             // Communicator
    UPstream::Request* req = nullptr    // Non-null for non-blocking
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PstreamDetail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UPstreamWrapping.txx"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
