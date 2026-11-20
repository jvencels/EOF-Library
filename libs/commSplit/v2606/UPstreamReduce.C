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
#include "PstreamReduceOps.H"
#include "UPstreamWrapping.H"

#include <cinttypes>

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Special reductions for bool

void Foam::UPstream::reduceAnd(bool& value, const int communicator)
{
    PstreamDetail::allReduce(&value, 1, MPI_C_BOOL, MPI_LAND, communicator);
}


void Foam::UPstream::reduceOr(bool& value, const int communicator)
{
    PstreamDetail::allReduce(&value, 1, MPI_C_BOOL, MPI_LOR, communicator);
}




// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static inline bool is_basic_dataType(Foam::UPstream::dataTypes id) noexcept
{
    return
    (
        int(id) >= int(Foam::UPstream::dataTypes::Basic_begin)
     && int(id)  < int(Foam::UPstream::dataTypes::Basic_end)
    );
}

static inline bool is_reduce_opCode(Foam::UPstream::opCodes id) noexcept
{
    return
    (
        int(id) >= int(Foam::UPstream::opCodes::Basic_begin)
     && int(id)  < int(Foam::UPstream::opCodes::Basic_end)
    );
}


namespace
{

using namespace Foam;

// Local function to print some error information
void printErrorMessage
(
    const void* values,
    const UPstream::dataTypes datatype_id,
    const UPstream::opCodes opcode_id
)
{
    FatalError
        << "Bad input for reduce(): likely a programming problem\n";

    if (!is_basic_dataType(datatype_id))
    {
        FatalError<< "    Non-basic data tyoe (" << int(datatype_id) << ")\n";
    }

    if (!is_reduce_opCode(opcode_id))
    {
        FatalError<< "    Invalid reduce op (" << int(opcode_id) << ")\n";
    }

    if (values == nullptr)
    {
        FatalError<< "    nullptr for values\n";
    }
    FatalError<< Foam::endl;
}

} // End anonymous namespace


// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

// The intel-mpi version of MPI_Reduce() does not accept IN_PLACE
// operations (issue #3331)
//
// The open-mpi version (tested up to 4.1) accepts IN_PLACE but fails
// with an MPI_ARG_ERR message.
//
// Do not assume that anyone actually supports this!

#undef Foam_vendor_supports_INPLACE_REDUCE

void Foam::UPstream::mpi_reduce
(
    void* values,                          // Type checking done by caller
    int count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const UPstream::opCodes opCodeId,      // Proper code passed by caller
    const int communicator,                // Index into MPICommunicators_
    UPstream::Request* req
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);
    MPI_Op optype = PstreamGlobals::getOpCode(opCodeId);

    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do - ignore
        return;
    }
    if
    (
        FOAM_UNLIKELY
        (
            !is_basic_dataType(dataTypeId)
         || !is_reduce_opCode(opCodeId)
         || (values == nullptr)
        )
    )
    {
        FatalErrorInFunction;
        printErrorMessage(values, dataTypeId, opCodeId);
        FatalError << Foam::abort(FatalError);
    }

    const bool withTopo =
    (
        (req == nullptr)
     && UPstream::usingTopoControl(UPstream::topoControls::reduce)
     && UPstream::usingNodeComms(communicator)
    );

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_reduce] : (inplace)"
            << " op:" << int(opCodeId)
            << " type:" << int(dataTypeId) << " count:" << count
            << " comm:" << communicator
            << " topo:" << withTopo << Foam::endl;
        // error::printStack(Perr);
    }

    // Workaround for missing/broken in-place handling.
    // Use a local buffer to send the data from.
    // - probably not thread-safe

    #ifndef Foam_vendor_supports_INPLACE_REDUCE
    static std::unique_ptr<char[]> work;
    static int work_len(0);

    const int num_bytes = [datatype,count]()
    {
        int size = 1;
        MPI_Type_size(datatype, &size);
        return (size * count);
    }();

    if (work_len < num_bytes)
    {
        // Min length to avoid many initial re-allocations
        work_len = std::max(256, num_bytes);
        work.reset();
        work = std::make_unique<char[]>(work_len);
    }
    void* send_buffer = work.get();

    std::memcpy(send_buffer, values, num_bytes);
    #else
    void* send_buffer = values;  // ie, in-place
    #endif

    if (withTopo)
    {
        // Topological reduce

        // Stage 1: local reduction within a node -> onto the node leader
        if (UPstream::is_parallel(UPstream::commLocalNode_))
        {
            if (FOAM_UNLIKELY(UPstream::debug))
            {
                Perr<< "[mpi_reduce] : (inplace)"
                    << " op:" << int(opCodeId)
                    << " type:" << int(dataTypeId) << " count:" << count
                    << " comm:" << UPstream::commLocalNode_
                    << " stage-1" << Foam::endl;
            }

            PstreamDetail::reduce
            (
                send_buffer,
                values,
                count,
                datatype,
                optype,
                UPstream::commLocalNode_
            );
        }

        // Stage 2: reduce between node leaders -> world leader
        if (UPstream::is_parallel(UPstream::commInterNode_))
        {
            // Transcribe the previous results as input for this stage
            #ifndef Foam_vendor_supports_INPLACE_REDUCE
            std::memcpy(send_buffer, values, num_bytes);
            #endif

            if (FOAM_UNLIKELY(UPstream::debug))
            {
                Perr<< "[mpi_reduce] : (inplace)"
                    << " op:" << int(opCodeId)
                    << " type:" << int(dataTypeId) << " count:" << count
                    << " comm:" << UPstream::commInterNode_
                    << " stage-2" << Foam::endl;
            }

            PstreamDetail::reduce
            (
                send_buffer,
                values,
                count,
                datatype,
                optype,
                UPstream::commInterNode_
            );
        }
    }
    else
    {
        // Regular reduce

        PstreamDetail::reduce
        (
            send_buffer,
            values,
            count,
            datatype,
            optype,
            communicator,
            req
        );
    }
}


void Foam::UPstream::mpi_allreduce
(
    void* values,                          // Type checking done by caller
    int count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const UPstream::opCodes opCodeId,      // Proper code passed by caller
    const int communicator,                // Index into MPICommunicators_
    UPstream::Request* req
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);
    MPI_Op optype = PstreamGlobals::getOpCode(opCodeId);

    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do - ignore
        return;
    }
    if
    (
        FOAM_UNLIKELY
        (
            !is_basic_dataType(dataTypeId)
         || !is_reduce_opCode(opCodeId)
         || (values == nullptr)
        )
    )
    {
        FatalErrorInFunction;
        printErrorMessage(values, dataTypeId, opCodeId);
        FatalError << Foam::abort(FatalError);
    }

    const bool withTopo =
    (
        (req == nullptr)
     && UPstream::usingTopoControl(UPstream::topoControls::reduce)
     && UPstream::usingNodeComms(communicator)
    );

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_allreduce] :"
            << " op:" << int(opCodeId)
            << " type:" << int(dataTypeId) << " count:" << count
            << " comm:" << communicator
            << " topo:" << withTopo << Foam::endl;
    }

    if (withTopo)
    {
        // Topological allReduce

        // Stage 1: local reduction within a node -> onto the node leader
        if (UPstream::is_parallel(UPstream::commLocalNode_))
        {
            if (FOAM_UNLIKELY(UPstream::debug))
            {
                Perr<< "[mpi_allreduce] :"
                    << " op:" << int(opCodeId)
                    << " type:" << int(dataTypeId) << " count:" << count
                    << " comm:" << UPstream::commLocalNode_
                    << " stage-1:reduce" << Foam::endl;
            }

            UPstream::mpi_reduce
            (
                values,
                count,
                dataTypeId,
                opCodeId,
                UPstream::commLocalNode_
            );
        }

        // Stage 2: all-reduce between node leaders
        if (UPstream::is_parallel(UPstream::commInterNode_))
        {
            if (FOAM_UNLIKELY(UPstream::debug))
            {
                Perr<< "[mpi_allreduce] :"
                    << " op:" << int(opCodeId)
                    << " type:" << int(dataTypeId) << " count:" << count
                    << " comm:" << UPstream::commInterNode_
                    << " stage-2:allreduce" << Foam::endl;
            }

            PstreamDetail::allReduce
            (
                values,
                count,
                datatype,
                optype,
                UPstream::commInterNode_
            );
        }

        // Finally, broadcast the information from each local node leader
        if (UPstream::is_parallel(UPstream::commLocalNode_))
        {
            if (FOAM_UNLIKELY(UPstream::debug))
            {
                Perr<< "[mpi_allreduce] :"
                    << " op:" << int(opCodeId)
                    << " type:" << int(dataTypeId) << " count:" << count
                    << " comm:" << UPstream::commLocalNode_
                    << " stage-3:broadcast" << Foam::endl;
            }

            PstreamDetail::broadcast
            (
                values,
                count,
                datatype,
                UPstream::commLocalNode_
            );
        }
    }
    else
    {
        // Regular allReduce

        PstreamDetail::allReduce
        (
            values,
            count,
            datatype,
            optype,
            communicator,
            req
        );
    }
}


// Based on the experience with MPI_Reduce (see above),
// do not assume that anyone actually supports the inplace option

void Foam::UPstream::mpi_scan_reduce
(
    void* values,                          // Type checking done by caller
    int count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const UPstream::opCodes opCodeId,      // Proper code passed by caller
    const int communicator,                // Index into MPICommunicators_
    const bool exclusive                   // Exclusive scan
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);
    MPI_Op optype = PstreamGlobals::getOpCode(opCodeId);

    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do - ignore
        return;
    }
    if
    (
        FOAM_UNLIKELY
        (
            !is_basic_dataType(dataTypeId)
         || !is_reduce_opCode(opCodeId)
         || (values == nullptr)
        )
    )
    {
        FatalErrorInFunction;
        printErrorMessage(values, dataTypeId, opCodeId);
        FatalError << Foam::abort(FatalError);
    }

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_scan_reduce] : (inplace)"
            << " op:" << int(opCodeId)
            << " type:" << int(dataTypeId) << " count:" << count
            << " comm:" << communicator
            << " excl:" << exclusive << Foam::endl;
        // error::printStack(Perr);
    }

    // Avoid the possibility of missing/broken in-place handling.
    // Use a local buffer to send the data from.
    // - probably not thread-safe

    static std::unique_ptr<char[]> work;
    static int work_len(0);

    const int num_bytes = [datatype,count]()
    {
        int size = 1;
        MPI_Type_size(datatype, &size);
        return (size * count);
    }();

    if (work_len < num_bytes)
    {
        // Min length to avoid many initial re-allocations
        work_len = std::max(256, num_bytes);
        work.reset();
        work = std::make_unique<char[]>(work_len);
    }
    void* send_buffer = work.get();

    std::memcpy(send_buffer, values, num_bytes);

    PstreamDetail::scanReduce
    (
        send_buffer,
        values,
        count,
        datatype,
        optype,
        communicator,
        exclusive
    );

    // In exclusive mode, the value for rank=0 is degenerate and the
    // standard doesn't say which value (if any) it will have.
    //
    // Overwrite with the original value so we are certain to have
    // a valid but perhaps meaningless value.

    if (exclusive && UPstream::master(communicator))
    {
        std::memcpy(values, send_buffer, num_bytes);
    }
}


// ************************************************************************* //
