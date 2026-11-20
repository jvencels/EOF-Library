/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
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

Namespace
    Foam::PstreamGlobals

Description
    Global functions and variables for working with parallel streams,
    but principally for MPI.

SourceFiles
    PstreamGlobals.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PstreamGlobals_H
#define Foam_PstreamGlobals_H

#include "DynamicList.H"
#include "FixedList.H"
// UPstream, UPstream::Request, UPstream::Window, ...
#include "openfoam_mpi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamGlobals
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//- Communicator of the OpenFOAM ranks only (commSplit override).
//  Equal to the world communicator after UPstream::init() has split
//  MPI_COMM_WORLD by colour; MPI_COMM_WORLD keeps all programs of an
//  MPMD job (eg, a coupled FEM solver started in the same mpirun).
extern MPI_Comm MPI_COMM_FOAM;

// Track if MPI_Comm_free is needed for communicator index in MPICommunicators_
extern DynamicList<bool> pendingMPIFree_;

// Current communicators, which may be allocated or predefined
// (eg, MPI_COMM_SELF, MPI_COMM_WORLD)
extern DynamicList<MPI_Comm> MPICommunicators_;

//- Outstanding non-blocking operations.
extern DynamicList<MPI_Request> outstandingRequests_;

// The fundamental count for each UPstream::dataTypes entry
typedef Foam::FixedList<int, 17> DataTypeCountLookupTable;

//- Fundamental count for each valid UPstream::dataTypes entry
//- Indexed by UPstream::dataTypes enum
extern DataTypeCountLookupTable dataTypesCount_;

// For UPstream::dataTypes lookup, includes space for last 'invalid' entry
typedef Foam::FixedList<MPI_Datatype, 17> DataTypeLookupTable;

//- MPI data types corresponding to fundamental and OpenFOAM types.
//- Indexed by UPstream::dataTypes enum
extern DataTypeLookupTable MPIdataTypes_;

// For UPstream::opCodes lookup, includes space for last 'invalid' entry
typedef Foam::FixedList<MPI_Op, 13> OpCodesLookupTable;

//- MPI operation types, indexed by UPstream::opCodes enum
extern OpCodesLookupTable MPIopCodes_;



// * * * * * * * * * * * * * * * Communicators * * * * * * * * * * * * * * * //

//- Initialize bookkeeping for MPI communicator index
void initCommunicator(const label index);

//- Fatal if communicator is outside the allocated range
inline void checkCommunicator(int comm, int rank)
{
    if (FOAM_UNLIKELY(comm < 0 || comm >= MPICommunicators_.size()))
    {
        FatalErrorInFunction
            << "rank:" << rank << " : illegal communicator "
            << comm << nl
            << "Communicator should be within range [0,"
            << MPICommunicators_.size()
            << ')' << Foam::abort(FatalError);
    }
}

//- True if warn communicator is active and not equal to given communicator
inline bool warnCommunicator(int comm) noexcept
{
    return (UPstream::warnComm >= 0 && comm != UPstream::warnComm);
}


// * * * * * * * * * * * * * * * * Data Types  * * * * * * * * * * * * * * * //

//- Create mapping into MPIdataTypes_ and define user data types
void initDataTypes();

//- Free any user data types
void deinitDataTypes();

//- Debugging only: check if data type mappings are non-null
bool checkDataTypes();

//- Debugging only: print data type names (all or just user-defined)
void printDataTypes(bool all = false);

//- Lookup of dataTypes enumeration as an MPI_Datatype
inline MPI_Datatype getDataType(UPstream::dataTypes id)
{
    return MPIdataTypes_[static_cast<int>(id)];
}

//- Fatal if data type is not valid
inline void checkDataType(UPstream::dataTypes id)
{
    if (id == UPstream::dataTypes::invalid)
    {
        FatalErrorInFunction
            << "Invalid data type"
            << Foam::abort(FatalError);
    }
}

//- Return MPI internal name for specified MPI_Datatype
std::string dataType_name(MPI_Datatype datatype);

//- Return MPI internal name for dataTypes enumeration
inline std::string dataType_name(UPstream::dataTypes id)
{
    return dataType_name(MPIdataTypes_[static_cast<int>(id)]);
}


// * * * * * * * * * * * * * * * * Op Codes  * * * * * * * * * * * * * * * * //

//- Create mapping into MPIopCodes_
void initOpCodes();

//- Free any user-defined op codes
void deinitOpCodes();

//- Debugging only: check if op code mappings are non-null
bool checkOpCodes();

//- Lookup of opCodes enumeration as an MPI_Op
inline MPI_Op getOpCode(UPstream::opCodes id)
{
    return MPIopCodes_[static_cast<int>(id)];
}

//- Fatal if opcode is not valid
inline void checkOpCode(UPstream::opCodes id)
{
    if (id == UPstream::opCodes::invalid)
    {
        FatalErrorInFunction
            << "Invalid operation code"
            << Foam::abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Requests  * * * * * * * * * * * * * * * * //

//- Reset UPstream::Request to MPI_REQUEST_NULL
//  Does not affect the stack of outstanding requests
inline void reset_request(UPstream::Request* req) noexcept
{
    if (req) *req = UPstream::Request(MPI_REQUEST_NULL);
}


//- Transcribe MPI_Request to UPstream::Request
//- (does not affect the stack of outstanding requests)
//- or else push onto list of outstanding requests
inline void push_request
(
    MPI_Request request,
    UPstream::Request* req = nullptr
)
{
    if (req)
    {
        // Transcribe as UPstream::Request
        *req = UPstream::Request(request);
    }
    else if (MPI_REQUEST_NULL != request)
    {
        // Push onto list of requests
        PstreamGlobals::outstandingRequests_.push_back(request);
    }
}


// * * * * * * * * * * * * * * Convenience Methods * * * * * * * * * * * * * //

//- Broadcast a single int64 value.
//
//  Ensures consistent data types. Used within the following:
//  - UIPBstream::bufferIPCrecv()
//  - UOPBstream::bufferIPCsend()
//  - UOPBstream::send(Foam::zero, ...)

inline bool broadcast_int64(int64_t& value, int comm, int root = 0)
{
    return
    (
        MPI_SUCCESS
     == MPI_Bcast(&value, 1, MPI_INT64_T, root, MPICommunicators_[comm])
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PstreamGlobals
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
