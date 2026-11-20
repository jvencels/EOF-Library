/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
    Copyright (C) 2026 Keysight Technologies
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

#include "UIPstream.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"
#include "IOstreams.H"

// FUTURE? probe and receive message
// - as of 2023-06 appears to be broken with INTELMPI + PMI-2 (slurm)
//   and perhaps other places so currently avoid

// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

// General blocking/non-blocking MPI receive
std::streamsize Foam::UPstream::mpi_receive
(
    const UPstream::commsTypes commsType,
    void* buf,                      // Type checking done by caller
    std::streamsize count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const int fromProcNo,
    const int tag,
    const int communicator,
    UPstream::Request* req
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    PstreamGlobals::reset_request(req);

    // Could check if nonBlocking and request are consistently specified...


    // TODO: some corrective action, at least when not nonBlocking
    #if 0
    // No warnings here, just on the sender side.
    if (count > std::streamsize(INT_MAX))
    {
        Perr<< "[mpi_recv] from rank " << fromProcNo
            << " exceeds INT_MAX values of "
            << PstreamGlobals::dataType_name(datatype)
            << Foam::endl;

        error::printStack(Perr);
    }
    #endif

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        Perr<< "[mpi_recv] : starting recv from:" << fromProcNo
            << " type:" << int(dataTypeId)
            << " count:" << label(count)
            << " tag:" << tag << " comm:" << communicator
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << " warnComm:" << UPstream::warnComm
            << Foam::endl;
        error::printStack(Perr);
    }
    else if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_recv] : starting recv from:" << fromProcNo
            << " type:" << int(dataTypeId)
            << " count:" << label(count)
            << " tag:" << tag << " comm:" << communicator
            << " commsType:" << UPstream::commsTypeNames[commsType]
            << Foam::endl;
    }

    int returnCode = MPI_ERR_UNKNOWN;

    profilingPstream::beginTiming();

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::scheduled
    )
    {
        // Not UPstream::commsTypes::nonBlocking

        MPI_Status status;

        {
            returnCode = MPI_Recv
            (
                buf,
                count,
                datatype,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
               &status
            );
        }

        profilingPstream::addGatherTime();

        if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
        {
            FatalErrorInFunction
                << "[mpi_recv] : cannot receive message from:"
                << fromProcNo
                << " type:" << int(dataTypeId)
                << " count:" << label(count) << " tag:" << tag
                << Foam::abort(FatalError);
            return 0;
        }
        else if (FOAM_UNLIKELY(UPstream::debug))
        {
            Perr<< "[mpi_recv] : finished recv from:"
                << fromProcNo
                << " type:" << int(dataTypeId)
                << " count:" << label(count) << " tag:" << tag
                << Foam::endl;
        }

        // Check size of message read (number of basic elements)
        MPI_Count num_recv(0);
        MPI_Get_elements_x(&status, datatype, &num_recv);

        // Errors
        if (FOAM_UNLIKELY(num_recv == MPI_UNDEFINED || int64_t(num_recv) < 0))
        {
            FatalErrorInFunction
                << "[mpi_recv] : receive from:" << fromProcNo
                << " type:" << int(dataTypeId)
                << " received count is undefined or negative value"
                << Foam::abort(FatalError);
        }
        else
        {
            // From number of basic elements to number of 'datatype'
            num_recv /= PstreamGlobals::dataTypesCount_[int(dataTypeId)];
        }

        if (FOAM_UNLIKELY(int64_t(num_recv) > int64_t(UList<char>::max_size())))
        {
            FatalErrorInFunction
                << "[mpi_recv] : receive from:" << fromProcNo
                << " type:" << int(dataTypeId)
                << " received count is larger than UList<T>::max_size()"
                << Foam::abort(FatalError);
        }
        else if (FOAM_UNLIKELY(count < std::streamsize(num_recv)))
        {
            FatalErrorInFunction
                << "[mpi_recv] : receive from:" << fromProcNo
                << " type:" << int(dataTypeId)
                << " count:" << label(count)
                << " buffer is too small for incoming message ("
                << label(num_recv) << ')'
                << Foam::abort(FatalError);
        }

        return std::streamsize(num_recv);
    }
    else if (commsType == UPstream::commsTypes::nonBlocking)
    {
        MPI_Request request;

        {
            returnCode = MPI_Irecv
            (
                buf,
                count,
                datatype,
                fromProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );
        }

        if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
        {
            FatalErrorInFunction
                << "[mpi_recv] : cannot start non-blocking receive from:"
                << fromProcNo
                << " type:" << int(dataTypeId)
                << " count:" << label(count)
                << Foam::abort(FatalError);

            return 0;
        }

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();


        if (FOAM_UNLIKELY(UPstream::debug))
        {
            Perr<< "[mpi_recv] : started non-blocking recv from:"
                << fromProcNo
                << " type:" << int(dataTypeId)
                << " count:" << label(count) << " tag:" << tag
                << " request:" <<
                (req ? label(-1) : PstreamGlobals::outstandingRequests_.size())
                << Foam::endl;
        }

        // Assume the message will be completely received.
        return count;
    }

    FatalErrorInFunction
        << "Unsupported communications type " << int(commsType)
        << Foam::abort(FatalError);

    return 0;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UIPstream::bufferIPCrecv()
{
    // Called by constructor
    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "UIPstream IPC read buffer :"
            << " from:" << fromProcNo_
            << " tag:" << tag_ << " comm:" << comm_
            << " wanted size:" << recvBuf_.capacity()
            << Foam::endl;
    }

    // Fallback value
    messageSize_ = recvBuf_.capacity();

    if (commsType() == UPstream::commsTypes::nonBlocking)
    {
        // Non-blocking
        // ~~~~~~~~~~~~
        // No chance of probing for size nor relying on the returned message
        // size (since it returns immediately without any further checks)
        //
        // Fortunately there are not many (any?) places that are using
        // a non-blocking IPstream with streaming anyhow.

        messageSize_ = recvBuf_.size();
    }
    else if (!recvBuf_.capacity())
    {
        // No buffer size allocated/specified - probe size of incoming message
        profilingPstream::beginTiming();

        const int source = (fromProcNo_ < 0) ? MPI_ANY_SOURCE : fromProcNo_;
        // Supporting MPI_ANY_TAG is not particularly useful...

        MPI_Status status;

        MPI_Probe
        (
            source,
            tag_,
            PstreamGlobals::MPICommunicators_[comm_],
           &status
        );

        if (fromProcNo_ < 0)
        {
            // Use the correct (probed) source
            fromProcNo_ = status.MPI_SOURCE;
        }

        profilingPstream::addProbeTime();

        // Buffer of characters (bytes)
        MPI_Count num_recv(0);
        MPI_Get_elements_x(&status, MPI_BYTE, &num_recv);

        // Errors
        if (FOAM_UNLIKELY(num_recv == MPI_UNDEFINED || int64_t(num_recv) < 0))
        {
            FatalErrorInFunction
                << "UIPstream IPC read buffer from:" << fromProcNo_
                << " received count is undefined or negative value"
                << Foam::abort(FatalError);
        }

        // Count is already in basic elements, no need to scale the result

        if (FOAM_UNLIKELY(int64_t(num_recv) > int64_t(UList<char>::max_size())))
        {
            FatalErrorInFunction
                << "UIPstream IPC read buffer from:" << fromProcNo_
                << " received count is larger than UList<T>::max_size()"
                << Foam::abort(FatalError);
        }

        if (FOAM_UNLIKELY(UPstream::debug))
        {
            Perr<< "UIPstream::bufferIPCrecv : probed size:"
                << label(num_recv) << Foam::endl;
        }

        recvBuf_.resize(label(num_recv));
        messageSize_ = label(num_recv);
    }

    // Errors
    if (FOAM_UNLIKELY(fromProcNo_ < 0))
    {
        FatalErrorInFunction
            << "UIPstream IPC failed with ANY_SOURCE probing?"
            << " cannot have proc=" << fromProcNo_
            << Foam::abort(FatalError);
    }

    std::streamsize count = UPstream::mpi_receive
    (
        commsType(),
        recvBuf_.data(),        // buffer
        messageSize_,           // expected size
        UPstream::dataTypes::type_byte,  // MPI_BYTE
        fromProcNo_,
        tag_,
        comm_,
        nullptr   // UPstream::Request
    );

    // Errors
    if (FOAM_UNLIKELY(count < 0))
    {
        FatalErrorInFunction
            << "UIPstream IPC read buffer from:" << fromProcNo_
            << " with negative size?"
            << Foam::abort(FatalError);
    }
    else if (FOAM_UNLIKELY(int64_t(count) > int64_t(UList<char>::max_size())))
    {
        FatalErrorInFunction
            << "UIPstream IPC read buffer from:" << fromProcNo_
            << " received size is larger than UList<T>::max_size()"
            << Foam::abort(FatalError);
    }

    // Set addressed size. Leave actual allocated memory intact.
    recvBuf_.resize(label(count));
    messageSize_ = label(count);

    if (recvBuf_.empty())
    {
        setEof();
    }
}


// ************************************************************************* //
