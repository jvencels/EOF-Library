/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

MPI_Comm Foam::PstreamGlobals::MPI_COMM_FOAM;

Foam::DynamicList<MPI_Request> Foam::PstreamGlobals::outstandingRequests_;

int Foam::PstreamGlobals::nTags_ = 0;

Foam::DynamicList<int> Foam::PstreamGlobals::freedTags_;

Foam::DynamicList<MPI_Comm> Foam::PstreamGlobals::MPICommunicators_;
Foam::DynamicList<MPI_Group> Foam::PstreamGlobals::MPIGroups_;


void Foam::PstreamGlobals::checkCommunicator
(
    const label comm,
    const label toProcNo
)
{
    if
    (
        comm < 0
     || comm >= PstreamGlobals::MPICommunicators_.size()
    )
    {
        FatalErrorInFunction
            << "toProcNo:" << toProcNo << " : illegal communicator "
            << comm << nl
            << "Communicator should be within range [0,"
            << PstreamGlobals::MPICommunicators_.size()
            << ')' << abort(FatalError);
    }
}


// ************************************************************************* //
