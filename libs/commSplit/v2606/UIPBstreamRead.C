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

#include "UIPstream.H"
#include "IOstreams.H"
#include "PstreamGlobals.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UIPBstream::bufferIPCrecv()
{
    // Uses double broadcast. Symmetric with UOPBstream::bufferIPCsend()
    // 1. for the data size
    // 2. for the data itself

    // Broadcast #1 - data size
    // Same data type must be used in UOPBstream::bufferIPCsend()

    int64_t count(0);
    if (!PstreamGlobals::broadcast_int64(count, comm_))
    {
        FatalErrorInFunction
            << "Broadcast failure receiving buffer size" << nl
            << " comm:" << comm_ << nl
            << Foam::abort(FatalError);
    }

    // This is not actually possible - sender uses List::size()
    //
    // if (FOAM_UNLIKELY(count > int64_t(UList<char>::max_size())))
    // {
    //     FatalErrorInFunction
    //         << "Broadcast list size larger than UList<char>::max_size()"
    //         << Foam::abort(FatalError);
    // }

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "UIPBstream IPC read buffer :"
            << " comm:" << comm_
            << " probed size:" << label(count)
            << " wanted size:" << recvBuf_.capacity()
            << Foam::endl;
    }


    // Set buffer size, avoiding any copying and resize doubling etc.
    recvBuf_.clear();
    if (recvBuf_.capacity() < label(count))
    {
        recvBuf_.setCapacity_nocopy(label(count));
    }
    recvBuf_.resize_nocopy(label(count));

    // This is the only real information we can trust
    messageSize_ = label(count);


    // Broadcast #2 - data content
    // - skip if there is no data to receive
    if
    (
        (count > 0)  // ie, not empty
     && !UPstream::mpi_broadcast
        (
            recvBuf_.data(),
            recvBuf_.size(),  // same as count
            UPstream::dataTypes::type_byte,
            comm_
        )
    )
    {
        FatalErrorInFunction
            << "Broadcast failure receiving buffer data:"
            << recvBuf_.size() << " comm:" << comm_ << nl
            << Foam::abort(FatalError);
    }

    if (recvBuf_.empty())
    {
        setEof();
    }
}


// ************************************************************************* //
