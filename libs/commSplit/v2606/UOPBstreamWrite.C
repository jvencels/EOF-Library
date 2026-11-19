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

#include "UOPstream.H"
#include "PstreamGlobals.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::UOPBstream::bufferIPCsend()
{
    // Uses double broadcast
    // 1. for the data size
    // 2. for the data itself
    // With this information, can determine and resize receive buffer

    PstreamGlobals::checkCommunicator(comm_, toProcNo_);

    // Broadcast #1 - data size
    // Same data type must be used in UIPBstream::bufferIPCrecv()

    int64_t count(sendBuf_.size());
    if (!PstreamGlobals::broadcast_int64(count, comm_))
    {
        FatalErrorInFunction
            << "Broadcast failure sending buffer size:"
            << label(count) << " comm:" << comm_ << nl
            << Foam::abort(FatalError);
        return false;
    }

    // Broadcast #2 - data content
    // - skip if there is no data to send
    if
    (
        (count > 0)  // ie, not empty
     && !UPstream::mpi_broadcast
        (
            sendBuf_.data(),
            sendBuf_.size(),  // same as count
            UPstream::dataTypes::type_byte,
            comm_
        )
    )
    {
        FatalErrorInFunction
            << "Broadcast failure sending buffer data:"
            << sendBuf_.size() << " comm:" << comm_ << nl
            << Foam::abort(FatalError);
        return false;
    }

    return true;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::UOPBstream::send(Foam::zero, int communicator, int root)
{
    if (UPstream::is_parallel(communicator))
    {
        int64_t count(0);
        PstreamGlobals::broadcast_int64(count, communicator, root);
    }
    else
    {
        // FatalError ?
    }
}


// ************************************************************************* //
