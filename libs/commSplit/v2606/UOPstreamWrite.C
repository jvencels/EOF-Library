/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

#include "UOPstream.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::UOPstream::bufferIPCsend()
{
    return UOPstream::write
    (
        commsType(),
        toProcNo_,
        sendBuf_.cdata(),
        sendBuf_.size(),
        tag_,
        comm_
    );
}


// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

// General blocking/non-blocking MPI send
bool Foam::UPstream::mpi_send
(
    const UPstream::commsTypes commsType,
    const void* buf,                       // Type checking done by caller
    std::streamsize count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const int toProcNo,
    const int tag,
    const int communicator,
    UPstream::Request* req,
    const UPstream::sendModes sendMode
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    PstreamGlobals::reset_request(req);

    // Could check if nonBlocking and request are consistently specified...


    // TODO: some corrective action, at least when not nonBlocking
    #if 0
    if (count > std::streamsize(INT_MAX))
    {
        Perr<< "[mpi_send] : to rank " << toProcNo
            << " type:" << int(dataTypeId)
            << " exceeds INT_MAX values"
            << Foam::endl;

        error::printStack(Perr);
    }
    #endif

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        Perr<< "[mpi_send] : starting send to:" << toProcNo
            << " type:" << int(dataTypeId)
            << " count:" << label(count)
            << " tag:" << tag << " comm:" << communicator
            << " commType:" << UPstream::commsTypeNames[commsType]
            << " warnComm:" << UPstream::warnComm
            << Foam::endl;
        error::printStack(Perr);
    }
    else if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_send] : starting send to:" << toProcNo
            << " type:" << int(dataTypeId)
            << " count:" << label(count)
            << " tag:" << tag << " comm:" << communicator
            << " commType:" << UPstream::commsTypeNames[commsType]
            << Foam::endl;
    }

    PstreamGlobals::checkCommunicator(communicator, toProcNo);

    int returnCode = MPI_ERR_UNKNOWN;

    profilingPstream::beginTiming();

    if (commsType == UPstream::commsTypes::buffered)
    {
        returnCode = MPI_Bsend
        (
            buf,
            count,
            datatype,
            toProcNo,
            tag,
            PstreamGlobals::MPICommunicators_[communicator]
        );

        // Assume these are from scatters ...
        profilingPstream::addScatterTime();

        if (FOAM_UNLIKELY(UPstream::debug))
        {
            Perr<< "[mpi_send] : finished buffered send to:"
                << toProcNo
                << " count:" << label(count) << " tag:" << tag
                << Foam::endl;
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::sendModes::sync == sendMode)
        {
            returnCode = MPI_Ssend
            (
                buf,
                count,
                datatype,
                toProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator]
            );
        }
        else
        {
            returnCode = MPI_Send
            (
                buf,
                count,
                datatype,
                toProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator]
            );
        }

        // Assume these are from scatters ...
        profilingPstream::addScatterTime();

        if (FOAM_UNLIKELY(UPstream::debug))
        {
            Perr<< "[mpi_send] : finished send to:"
                << toProcNo
                << " type:" << int(dataTypeId)
                << " count:" << label(count) << " tag:" << tag
                << Foam::endl;
        }
    }
    else if (commsType == UPstream::commsTypes::nonBlocking)
    {
        MPI_Request request;

        if (UPstream::sendModes::sync == sendMode)
        {
            returnCode = MPI_Issend
            (
                buf,
                count,
                datatype,
                toProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );
        }
        else
        {
            returnCode = MPI_Isend
            (
                buf,
                count,
                datatype,
                toProcNo,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );
        }

        if (FOAM_UNLIKELY(UPstream::debug))
        {
            Perr<< "[mpi_send] : started non-blocking send to:"
                << toProcNo
                << " type:" << int(dataTypeId)
                << " count:" << label(count) << " tag:" << tag
                << " request:" <<
                (req ? label(-1) : PstreamGlobals::outstandingRequests_.size())
                << Foam::endl;
        }

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType)
            << Foam::abort(FatalError);
        return false;
    }

    return (returnCode == MPI_SUCCESS);
}


// ************************************************************************* //
