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
#include "PstreamGlobals.H"
#include "UPstreamWrapping.H"
#include "vector.H"  // for debugging

#undef STRINGIFY
#undef STRING_QUOTE

#define STRINGIFY(content) #content
#define STRING_QUOTE(input) STRINGIFY(input)

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

inline bool is_nonAggregate(Foam::UPstream::dataTypes id) noexcept
{
    return
    (
        int(id) >= int(Foam::UPstream::dataTypes::Basic_begin)
     && int(id)  < int(Foam::UPstream::dataTypes::Basic_end)
    )
    ||
    (
        int(id) >= int(Foam::UPstream::dataTypes::User_begin)
     && int(id)  < int(Foam::UPstream::dataTypes::User_end)
    );
}

// Local function to print some error information
inline void printErrorNonIntrinsic
(
    const char* context,
    Foam::UPstream::dataTypes dataTypeId
)
{
    using namespace Foam;

    FatalError
        << "Bad input for " << context << ": likely a programming problem\n"
        << "    Non-intrinsic/non-user data (type:" << int(dataTypeId) << ")\n"
        << Foam::endl;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_gather
(
    const void* sendData,       // Type checking done by caller
    void* recvData,             // Type checking done by caller
    int count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller

    const int communicator,     // Index into MPICommunicators_
    UPstream::Request* req
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_gather] :";

        // Appears to be an in-place request
        if
        (
            UPstream::master(communicator)
         && (!sendData || (sendData == recvData))
        )
        {
            Perr<< " (inplace)";
        }

        Perr<< " type:" << int(dataTypeId) << " count:" << count
            << " comm:" << communicator
            << Foam::endl;
    }

    {
        // Regular gather

        PstreamDetail::gather
        (
            sendData,
            recvData,
            count,
            datatype,
            communicator,
            req
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_scatter
(
    const void* sendData,       // Type checking done by caller
    void* recvData,             // Type checking done by caller
    int count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller

    const int communicator,     // Index into MPICommunicators_
    UPstream::Request* req
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_scatter] :";

        // Appears to be an in-place request
        if
        (
            UPstream::master(communicator)
         && (!recvData || (sendData == recvData))
        )
        {
            Perr<< " (inplace)";
        }

        Perr<< " type:" << int(dataTypeId) << " count:" << count
            << " comm:" << communicator
            << Foam::endl;
    }

    {
        // Regular scatter

        PstreamDetail::scatter
        (
            sendData,
            recvData,
            count,
            datatype,
            communicator,
            req
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_allgather
(
    void* allData,        // Type checking done by caller
    int count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller

    const int communicator,     // Index into MPICommunicators_
    UPstream::Request* req
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_allgather] :"
            << " type:" << int(dataTypeId) << " count:" << count
            << " comm:" << communicator
            << Foam::endl;
    }

    {
        // Regular all gather

        PstreamDetail::allGather
        (
            allData,
            count,
            datatype,
            communicator,
            req
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_gatherv
(
    const void* sendData,
    int sendCount,
    void* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,

    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const int communicator
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    // Runtime assert that we are not using aggregated data types
    if (FOAM_UNLIKELY(!is_nonAggregate(dataTypeId)))
    {
        FatalErrorInFunction;
        printErrorNonIntrinsic("MPI_Gatherv()", dataTypeId);
        FatalError << Foam::abort(FatalError);
    }

    const label np = UPstream::nProcs(communicator);

    // For total-size calculation,
    // don't rely on recvOffsets being (np+1)
    const int totalSize =
    (
        (UPstream::master(communicator) && np > 1)
      ? (recvOffsets[np-1] + recvCounts[np-1])
      : 0
    );

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_gatherv] :"
            << " type:" << int(dataTypeId)
            << " count:" << sendCount
            << " total:" << totalSize
            << " comm:" << communicator
            << " recvCounts:" << flatOutput(recvCounts)
            << " recvOffsets:" << flatOutput(recvOffsets)
            << Foam::endl;
    }

    {
        PstreamDetail::gatherv
        (
            sendData, sendCount,
            recvData, recvCounts, recvOffsets,
            datatype, communicator
        );
    }

    // Extended debugging. Limit to master:

    #if 0
    if (FOAM_UNLIKELY(UPstream::debug))
    {
        if (UPstream::master(communicator))
        {
            switch (dataTypeId)
            {
                #undef  dataPrinter
                #define dataPrinter(enumType, nativeType)           \
                case UPstream::dataTypes::enumType :                \
                {                                                   \
                    UList<nativeType> combined                      \
                    (                                               \
                        static_cast<nativeType*>(recvData),         \
                        totalSize                                   \
                    );                                              \
                                                                    \
                    Info<< "[mpi_gatherv] => "                      \
                    "List<" STRING_QUOTE(nativeType) "> ";          \
                    combined.writeList(Info) << Foam::endl;         \
                                                                    \
                    break;                                          \
                }

                // Some common types
                dataPrinter(type_int32, int32_t);
                dataPrinter(type_int64, int64_t);
                dataPrinter(type_float, float);
                dataPrinter(type_double, double);
                dataPrinter(type_3float, floatVector);
                dataPrinter(type_3double, doubleVector);

                // Some other type
                default: break;
                #undef dataPrinter
            }
        }
    }
    #endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_scatterv
(
    const void* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,

    void* recvData,
    int recvCount,

    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const int communicator
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    // Runtime assert that we are not using aggregated data types
    if (FOAM_UNLIKELY(!is_nonAggregate(dataTypeId)))
    {
        FatalErrorInFunction;
        printErrorNonIntrinsic("MPI_Scatterv()", dataTypeId);
        FatalError << Foam::abort(FatalError);
    }

    {
        PstreamDetail::scatterv
        (
            sendData, sendCounts, sendOffsets,
            recvData, recvCount,
            datatype, communicator
        );
    }
}


// ************************************************************************* //
