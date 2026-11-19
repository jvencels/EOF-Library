/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2025 OpenCFD Ltd.
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

#include "UPstream.H"
#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::Communicator::Communicator() noexcept
:
    UPstream::Communicator(MPI_COMM_NULL)
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::UPstream::Communicator
Foam::UPstream::Communicator::lookup(int communicator)
{
    MPI_Comm mpiComm = MPI_COMM_NULL;

    if (communicator < 0)
    {
        // Default/short-cut : world-comm
        communicator = UPstream::worldComm;
    }

    if (communicator < PstreamGlobals::MPICommunicators_.size())
    {
        mpiComm = PstreamGlobals::MPICommunicators_[communicator];
    }
    else
    {
        WarningInFunction
            << "Illegal communicator " << communicator << nl
            << "Should be within range [-1,"
            << PstreamGlobals::MPICommunicators_.size()
            << ')' << endl;
    }

    return UPstream::Communicator(mpiComm);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Communicator::good() const noexcept
{
    return MPI_COMM_NULL != PstreamUtils::Cast::to_mpi(*this);
}


void Foam::UPstream::Communicator::reset() noexcept
{
    *this = UPstream::Communicator(MPI_COMM_NULL);
}


int Foam::UPstream::Communicator::size() const
{
    int val = 0;

    MPI_Comm comm = PstreamUtils::Cast::to_mpi(*this);

    if (MPI_COMM_SELF == comm)
    {
        return 1;
    }
    else if
    (
        (MPI_COMM_NULL == comm)
     || (MPI_SUCCESS != MPI_Comm_size(comm, &val))
    )
    {
        val = 0;
    }

    return val;
}


// ************************************************************************* //
