/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "UPstreamWrapping.H"
#include "profilingPstream.H"
#include "PstreamGlobals.H"
#include "Map.H"
#include <cstring>  // memmove

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::PstreamDetail::broadcast
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    const int communicator,
    int root
)
{
    if (!UPstream::is_parallel(communicator))
    {
        return true;
    }

    int returnCode(MPI_SUCCESS);

    if (root < 0)
    {
        // Use last rank as the root
        root = static_cast<int>(UPstream::nProcs(communicator)-1);
    }

    profilingPstream::beginTiming();

    {
        returnCode =
            MPI_Bcast
            (
                values,
                count,
                datatype,
                root, // The broadcast root (normally 0 == UPstream::masterNo())
                PstreamGlobals::MPICommunicators_[communicator]
            );
    }

    profilingPstream::addBroadcastTime();

    return (returnCode == MPI_SUCCESS);
}


template<class Type>
void Foam::PstreamDetail::reduce
(
    const Type* sendData,
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const int communicator,

    UPstream::Request* req
)
{
    PstreamGlobals::reset_request(req);

    const bool immediate = (req);

    if (!UPstream::is_parallel(communicator))
    {
        return;
    }

    const void* send_buffer = sendData;
    if
    (
        UPstream::master(communicator)
     && (!sendData || (sendData == values))
    )
    {
        // Appears to be an in-place request.
        // - this setting only relevant (or usable) on the root rank
        send_buffer = MPI_IN_PLACE;
    }

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        if (immediate)
        {
            Perr<< "** MPI_Ireduce (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Reduce (blocking):";
        }
        if (UPstream::master(communicator) && (send_buffer == MPI_IN_PLACE))
        {
            Perr<< " [inplace]";
        }
        Perr<< " count:" << count
            << " comm:" << communicator
            << " warnComm:" << UPstream::warnComm << endl;
        error::printStack(Perr);
    }


    int returnCode(MPI_ERR_UNKNOWN);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        returnCode =
            MPI_Ireduce
            (
                send_buffer,
                values,
                count,
                datatype,
                optype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        returnCode =
            MPI_Reduce
            (
                send_buffer,
                values,
                count,
                datatype,
                optype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[communicator]
            );

        profilingPstream::addReduceTime();
    }

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalErrorInFunction<< "MPI Reduce ";
        if (immediate) FatalError<< "(non-blocking) ";

        FatalError
            << "failed for count:" << count
            << Foam::abort(FatalError);
    }
}


template<class Type>
void Foam::PstreamDetail::allReduce
(
    Type* values,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const int communicator,

    UPstream::Request* req
)
{
    PstreamGlobals::reset_request(req);

    const bool immediate = (req);

    if (!UPstream::is_parallel(communicator))
    {
        return;
    }

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        if (immediate)
        {
            Perr<< "** MPI_Iallreduce (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Allreduce (blocking):";
        }
        if constexpr (std::is_void_v<Type>)
        {
            Perr<< " count:" << count;
        }
        else
        {
            if (count == 1)
            {
                Perr<< (*values);
            }
            else
            {
                Perr<< UList<Type>(values, count);
            }
        }
        Perr<< " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm << endl;
        error::printStack(Perr);
    }


    int returnCode(MPI_ERR_UNKNOWN);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        returnCode =
            MPI_Iallreduce
            (
                MPI_IN_PLACE,  // recv is also send
                values,
                count,
                datatype,
                optype,
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        returnCode =
            MPI_Allreduce
            (
                MPI_IN_PLACE,  // recv is also send
                values,
                count,
                datatype,
                optype,
                PstreamGlobals::MPICommunicators_[communicator]
            );

        profilingPstream::addReduceTime();
    }

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalErrorInFunction<< "MPI Allreduce ";
        if (immediate) FatalError<< "(non-blocking) ";

        FatalError<< "failed for count:" << count
            << Foam::abort(FatalError);
    }
}


template<class Type>
void Foam::PstreamDetail::scanReduce
(
    const Type* sendData,
    Type* recvData,
    int count,
    MPI_Datatype datatype,
    MPI_Op optype,
    const int communicator,

    const int exclusive
)
{
    if (!UPstream::is_parallel(communicator))
    {
        return;
    }

    const void* send_buffer = sendData;
    if (!sendData || (sendData == recvData))
    {
        // Appears to be an in-place request.
        // - this setting only relevant (or usable) on the root rank
        send_buffer = MPI_IN_PLACE;
    }

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        if (exclusive)
        {
            Perr<< "** MPI_Exscan (blocking):";
        }
        else
        {
            Perr<< "** MPI_Scan (blocking):";
        }
        if (UPstream::master(communicator) && (send_buffer == MPI_IN_PLACE))
        {
            Perr<< " [inplace]";
        }
        Perr<< " count:" << count
            << " comm:" << communicator
            << " warnComm:" << UPstream::warnComm << endl;
        error::printStack(Perr);
    }


    int returnCode(MPI_ERR_UNKNOWN);

    profilingPstream::beginTiming();

    if (exclusive)
    {
        returnCode =
            MPI_Exscan
            (
                send_buffer,
                recvData,
                count,
                datatype,
                optype,
                PstreamGlobals::MPICommunicators_[communicator]
            );
    }
    else
    {
        returnCode =
            MPI_Scan
            (
                send_buffer,
                recvData,
                count,
                datatype,
                optype,
                PstreamGlobals::MPICommunicators_[communicator]
            );
    }
    profilingPstream::addReduceTime();

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalErrorInFunction<< "MPI Scan ";
        if (exclusive) FatalError<< "(exclusive) ";

        FatalError
            << "failed for count:" << count
            << Foam::abort(FatalError);
    }
}


template<class Type>
void Foam::PstreamDetail::allToAll
(
    const UList<Type>& sendData,
    UList<Type>& recvData,
    MPI_Datatype datatype,
    const int communicator,

    UPstream::Request* req
)
{
    static_assert(!std::is_void_v<Type>, "Does not handle void types");

    PstreamGlobals::reset_request(req);

    const bool immediate = (req);

    if (!UPstream::is_rank(communicator))
    {
        return;
    }

    const label numProc = UPstream::nProcs(communicator);

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        if (immediate)
        {
            Perr<< "** MPI_Ialltoall (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Alltoall (blocking):";
        }
        Perr<< " numProc:" << numProc
            << " sendData:" << sendData.size()
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    if
    (
        (sendData.size() != numProc || recvData.size() != numProc)
    )
    {
        FatalErrorInFunction
            << "Have " << numProc << " ranks, but size of sendData:"
            << sendData.size() << " or recvData:" << recvData.size()
            << " is different!"
            << Foam::abort(FatalError);
    }

    if (!UPstream::is_parallel(communicator))
    {
        recvData.deepCopy(sendData);
        return;
    }


    int returnCode(MPI_ERR_UNKNOWN);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        returnCode =
            MPI_Ialltoall
            (
                // NOTE: const_cast is a temporary hack for
                // backward-compatibility with versions of OpenMPI < 1.7.4
                const_cast<Type*>(sendData.cdata()),
                1,                      // one element per rank
                datatype,
                recvData.data(),
                1,                      // one element per rank
                datatype,
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        returnCode =
            MPI_Alltoall
            (
                // NOTE: const_cast is a temporary hack for
                // backward-compatibility with versions of OpenMPI < 1.7.4
                const_cast<Type*>(sendData.cdata()),
                1,                      // one element per rank
                datatype,
                recvData.data(),
                1,                      // one element per rank
                datatype,
                PstreamGlobals::MPICommunicators_[communicator]
            );

        profilingPstream::addAllToAllTime();
    }

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalErrorInFunction<< "MPI Alltoall ";
        if (immediate) FatalError<< "(non-blocking) ";

        FatalError
            << "[comm: " << communicator << "] failed for "
            << sendData << endl
            << Foam::abort(FatalError);
    }
}


template<class Type>
void Foam::PstreamDetail::allToAllv
(
    const Type* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,

    Type* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,

    MPI_Datatype datatype,
    const int communicator,

    UPstream::Request* req
)
{
    static_assert(!std::is_void_v<Type>, "Does not handle void types");

    PstreamGlobals::reset_request(req);

    const bool immediate = (req);

    if (!UPstream::is_rank(communicator))
    {
        return;
    }

    const label np = UPstream::nProcs(communicator);

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        if (immediate)
        {
            Perr<< "** MPI_Ialltoallv (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Alltoallv (blocking):";
        }
        Perr<< " sendCounts:" << sendCounts
            << " sendOffsets:" << sendOffsets
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    if
    (
        (sendCounts.size() != np || sendOffsets.size() < np)
     || (recvCounts.size() != np || recvOffsets.size() < np)
    )
    {
        FatalErrorInFunction
            << "Have " << np << " ranks, but sendCounts:" << sendCounts.size()
            << ", sendOffsets:" << sendOffsets.size()
            << ", recvCounts:" << recvCounts.size()
            << " or recvOffsets:" << recvOffsets.size()
            << " is different!"
            << Foam::abort(FatalError);
    }

    if (!UPstream::is_parallel(communicator))
    {
        if (FOAM_UNLIKELY(recvCounts[0] != sendCounts[0]))
        {
            FatalErrorInFunction
                << "Bytes to send " << sendCounts[0]
                << " does not equal bytes to receive " << recvCounts[0]
                << Foam::abort(FatalError);
        }
        std::memmove
        (
            recvData,
            (sendData + sendOffsets[0]),
            recvCounts[0]*sizeof(Type)
        );

        return;
    }


    int returnCode(MPI_ERR_UNKNOWN);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        returnCode =
            MPI_Ialltoallv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        returnCode =
            MPI_Alltoallv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                PstreamGlobals::MPICommunicators_[communicator]
            );

        profilingPstream::addAllToAllTime();
    }

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalErrorInFunction<< "MPI Alltoallv ";
        if (immediate) FatalError<< "(non-blocking) ";

        FatalError
            << "[comm: " << communicator << "] failed for "
            << " For sendCounts " << sendCounts
            << " recvCounts " << recvCounts << endl
            << Foam::abort(FatalError);
    }
}


template<class Type>
void Foam::PstreamDetail::allToAllConsensus
(
    const UList<Type>& sendData,
    UList<Type>& recvData,
    MPI_Datatype datatype,
    const int tag,
    const int communicator
)
{
    static_assert(!std::is_void_v<Type>, "Does not handle void types");

    const bool initialBarrier = (UPstream::tuning_NBX_ > 0);

    if (!UPstream::is_rank(communicator))
    {
        return;  // Process not in communicator
    }

    const label myProci = UPstream::myProcNo(communicator);
    const label numProc = UPstream::nProcs(communicator);

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        Perr<< "** non-blocking consensus Alltoall (list):";
        Perr<< " numProc:" << numProc
            << " sendData:" << sendData.size()
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    if (sendData.size() != numProc || recvData.size() != numProc)
    {
        FatalErrorInFunction
            << "Have " << numProc << " ranks, but size of sendData:"
            << sendData.size() << " or recvData:" << recvData.size()
            << " is different!"
            << Foam::abort(FatalError);
    }

    // Initial: assign zero everywhere. Values of zero are never transmitted
    const Type zeroValue = pTraits<Type>::zero;
    recvData = zeroValue;

    if (!UPstream::is_parallel(communicator))
    {
        // Non-parallel : deep copy
        recvData.deepCopy(sendData);
        return;
    }

    // Fake send/recv for myself
    {
        recvData[myProci] = sendData[myProci];
    }

    // Implementation description
    // --------------------------
    // "Scalable Communication Protocols for Dynamic Sparse Data Exchange",
    // Hoeffler, Siebert, Lumsdaine
    // May 2010 ACM SIGPLAN Notices 45(5):159-168
    // https://doi.org/10.1145/1837853.1693476
    //
    // - http://unixer.de/publications/img/hoefler-dsde-protocols.pdf
    //
    // Algorithm NBX: Nonblocking consensus

    // This specific specialization is largely just for integer data
    // so we initialise the receiving data with zero and then
    // do not send/recv them.
    // This is because we are dealing with a flat list of entries to
    // send and not a sparse Map etc.


    // ------------------------------------------------------------------------
    // Setup sends
    // ------------------------------------------------------------------------

    profilingPstream::beginTiming();

    // An initial barrier may help to avoid synchronisation problems
    // caused elsewhere
    if (initialBarrier)
    {
        MPI_Barrier(PstreamGlobals::MPICommunicators_[communicator]);
    }

    DynamicList<MPI_Request> sendRequests(sendData.size());

    // Start nonblocking synchronous send to destination rank
    for (label proci = 0; proci < numProc; ++proci)
    {
        if (sendData[proci] == zeroValue)
        {
            // Do not send/recv empty data
        }
        else if (proci != myProci)
        {
            // Has data to send

            MPI_Issend
            (
               &sendData[proci],
                1,              // one element per rank
                datatype,
                proci,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
               &sendRequests.emplace_back()
            );
        }
    }


    // ------------------------------------------------------------------------
    // Probe and receive
    // ------------------------------------------------------------------------

    MPI_Request barrierRequest;

    for (bool barrier_active = false, done = false; !done; /*nil*/)
    {
        int flag = 0;
        MPI_Status status;

        MPI_Iprobe
        (
            MPI_ANY_SOURCE,
            tag,
            PstreamGlobals::MPICommunicators_[communicator],
           &flag,
           &status
        );

        if (flag)
        {
            // Message found, receive into dest buffer location
            const int proci = status.MPI_SOURCE;

            // Only send/recv a single (fundamental) data type
            int count(0);
            MPI_Get_count(&status, datatype, &count);

            if (FOAM_UNLIKELY(count != 1))
            {
                FatalErrorInFunction
                    << "Incorrect message size from proc=" << proci
                    << ". Expected 1 but had " << count << nl
                    << exit(FatalError);
            }

            // Regular blocking receive [the data are small]

            MPI_Recv
            (
               &recvData[proci],
                count,          // count=1 (see above)
                datatype,
                proci,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                MPI_STATUS_IGNORE
            );
        }

        if (barrier_active)
        {
            // Test barrier for completion
            // - all received, or nothing to receive
            MPI_Test(&barrierRequest, &flag, MPI_STATUS_IGNORE);

            if (flag)
            {
                done = true;
            }
        }
        else
        {
            // Check if all sends have arrived
            MPI_Testall
            (
                sendRequests.size(),
                sendRequests.data(),
               &flag, MPI_STATUSES_IGNORE
            );

            if (flag)
            {
                MPI_Ibarrier
                (
                    PstreamGlobals::MPICommunicators_[communicator],
                   &barrierRequest
                );
                barrier_active = true;
            }
        }
    }

    profilingPstream::addAllToAllTime();
}


template<class Type>
void Foam::PstreamDetail::allToAllConsensus
(
    const Map<Type>& sendBufs,
    Map<Type>& recvBufs,
    MPI_Datatype datatype,
    const int tag,
    const int communicator
)
{
    static_assert(!std::is_void_v<Type>, "Does not handle void types");

    const bool initialBarrier = (UPstream::tuning_NBX_ > 0);

    const label myProci = UPstream::myProcNo(communicator);
    const label numProc = UPstream::nProcs(communicator);

    if (!UPstream::is_rank(communicator))
    {
        return;
    }

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        Perr<< "** non-blocking consensus Alltoall (map):";
        Perr<< " numProc:" << numProc
            << " sendData:" << sendBufs.size()
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    // Initial: clear out everything
    recvBufs.clear();

    // Fake send/recv for myself - parallel or non-parallel
    {
        const auto iter = sendBufs.find(myProci);
        if (iter.good())
        {
            // Do myself: insert_or_assign
            recvBufs(iter.key()) = iter.val();
        }
    }

    if (!UPstream::is_parallel(communicator))
    {
        // Nothing left to do
        return;
    }


    // ------------------------------------------------------------------------
    // Setup sends
    // ------------------------------------------------------------------------

    // Algorithm NBX: Nonblocking consensus
    // Implementation like above, but sending map data.

    DynamicList<MPI_Request> sendRequests(sendBufs.size());

    profilingPstream::beginTiming();

    // An initial barrier may help to avoid synchronisation problems
    // caused elsewhere
    if (initialBarrier)
    {
        MPI_Barrier(PstreamGlobals::MPICommunicators_[communicator]);
    }


    // Start nonblocking synchronous send to destination ranks

    // Same as forAllConstIters()
    const auto endIter = sendBufs.cend();
    for (auto iter = sendBufs.cbegin(); iter != endIter; ++iter)
    {
        const label proci = iter.key();
        const auto& sendData = iter.val();

        if (proci != myProci && proci >= 0 && proci < numProc)
        {
            MPI_Issend
            (
               &sendData,
                1,              // one element per rank
                datatype,
                proci,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
               &sendRequests.emplace_back()
            );
        }
    }


    // ------------------------------------------------------------------------
    // Probe and receive
    // ------------------------------------------------------------------------

    MPI_Request barrierRequest;

    for (bool barrier_active = false, done = false; !done; /*nil*/)
    {
        int flag = 0;
        MPI_Status status;

        MPI_Iprobe
        (
            MPI_ANY_SOURCE,
            tag,
            PstreamGlobals::MPICommunicators_[communicator],
           &flag,
           &status
        );

        if (flag)
        {
            // Message found, receive into dest buffer location
            const int proci = status.MPI_SOURCE;

            // Only send/recv a single (fundamental) data type
            int count(0);
            MPI_Get_count(&status, datatype, &count);

            if (FOAM_UNLIKELY(count != 1))
            {
                FatalErrorInFunction
                    << "Incorrect message size from proc=" << proci
                    << ". Expected 1 but had " << count << nl
                    << exit(FatalError);
            }

            auto& recvData = recvBufs(proci);

            // Regular blocking receive [the data are small]

            MPI_Recv
            (
               &recvData,
                count,          // count=1 (see above)
                datatype,
                proci,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                MPI_STATUS_IGNORE
            );
        }

        if (barrier_active)
        {
            // Test barrier for completion
            // - all received, or nothing to receive
            MPI_Test(&barrierRequest, &flag, MPI_STATUS_IGNORE);

            if (flag)
            {
                done = true;
            }
        }
        else
        {
            // Check if all sends have arrived
            MPI_Testall
            (
                sendRequests.size(),
                sendRequests.data(),
               &flag, MPI_STATUSES_IGNORE
            );

            if (flag)
            {
                MPI_Ibarrier
                (
                    PstreamGlobals::MPICommunicators_[communicator],
                   &barrierRequest
                );
                barrier_active = true;
            }
        }
    }

    profilingPstream::addAllToAllTime();
}


template<class Type>
void Foam::PstreamDetail::gather
(
    const Type* sendData,
    Type* recvData,

    int count,
    MPI_Datatype datatype,

    const int communicator,
    UPstream::Request* req
)
{
    PstreamGlobals::reset_request(req);

    const bool immediate = (req);

    if (!count || !UPstream::is_rank(communicator))
    {
        return;
    }
    else if (!UPstream::is_parallel(communicator))
    {
        if constexpr (std::is_void_v<Type>)
        {
            // Cannot copy data here since we don't know the number of bytes
            // - must be done by the caller.
        }
        else if (sendData && recvData && (sendData != recvData))
        {
            std::memmove(recvData, sendData, count*sizeof(Type));
        }
        return;
    }

    const void* send_buffer = sendData;
    if
    (
        UPstream::master(communicator)
     && (!sendData || (sendData == recvData))
    )
    {
        // Appears to be an in-place request.
        // - this setting only relevant (or usable) on the root rank
        send_buffer = MPI_IN_PLACE;
    }

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        if (immediate)
        {
            Perr<< "** MPI_Igather (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Gather (blocking):";
        }
        if (UPstream::master(communicator) && (send_buffer == MPI_IN_PLACE))
        {
            Perr<< " [inplace]";
        }
        Perr<< " numProc:" << UPstream::nProcs(communicator)
            << " count:" << count
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    int returnCode(MPI_ERR_UNKNOWN);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        returnCode =
            MPI_Igather
            (
                send_buffer, count, datatype,
                recvData, count, datatype,
                0,  // root: UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        returnCode =
            MPI_Gather
            (
                send_buffer, count, datatype,
                recvData, count, datatype,
                0,  // root: UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[communicator]
            );

        profilingPstream::addGatherTime();
    }

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalErrorInFunction << "MPI Gather ";
        if (immediate) FatalError<< "(non-blocking) ";

        FatalError
            << "[comm: " << communicator << "] failed."
            << " count:" << count << nl
            << Foam::abort(FatalError);
    }
}


template<class Type>
void Foam::PstreamDetail::scatter
(
    const Type* sendData,
    Type* recvData,

    int count,
    MPI_Datatype datatype,

    const int communicator,
    UPstream::Request* req
)
{
    PstreamGlobals::reset_request(req);

    const bool immediate = (req);

    if (!count || !UPstream::is_rank(communicator))
    {
        return;
    }
    else if (!UPstream::is_parallel(communicator))
    {
        if constexpr (std::is_void_v<Type>)
        {
            // Cannot copy data here since we don't know the number of bytes
            // - must be done by the caller.
        }
        else if (sendData && recvData && (sendData != recvData))
        {
            std::memmove(recvData, sendData, count*sizeof(Type));
        }
        return;
    }

    void* recv_buffer = recvData;

    if
    (
        UPstream::master(communicator)
     && (!recvData || (sendData == recvData)))
    {
        // Appears to be an in-place request.
        // - this setting only relevant (or usable) on the root rank
        recv_buffer = MPI_IN_PLACE;
    }

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        if (immediate)
        {
            Perr<< "** MPI_Iscatter (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Scatter (blocking):";
        }
        if (UPstream::master(communicator) && (recv_buffer == MPI_IN_PLACE))
        {
            Perr<< " [inplace]";
        }
        Perr<< " numProc:" << UPstream::nProcs(communicator)
            << " count:" << count
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }
    int returnCode(MPI_ERR_UNKNOWN);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        returnCode =
            MPI_Iscatter
            (
                sendData, count, datatype,
                recv_buffer, count, datatype,
                0,  // root: UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        returnCode =
            MPI_Scatter
            (
                sendData, count, datatype,
                recv_buffer, count, datatype,
                0,  // root: UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[communicator]
            );

        profilingPstream::addScatterTime();
    }

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalErrorInFunction << "MPI Scatter ";
        if (immediate) FatalError<< "(non-blocking) ";

        FatalError
            << "[comm: " << communicator << "] failed."
            << " count:" << count << nl
            << Foam::abort(FatalError);
    }
}


template<class Type>
void Foam::PstreamDetail::gatherv
(
    const Type* sendData,
    int sendCount,

    Type* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,

    MPI_Datatype datatype,
    const int communicator,

    UPstream::Request* req
)
{
    PstreamGlobals::reset_request(req);

    const bool immediate = (req);

    if (!UPstream::is_rank(communicator))
    {
        return;
    }
    else if (!UPstream::is_parallel(communicator))
    {
        if constexpr (std::is_void_v<Type>)
        {
            // Cannot copy data here since we don't know the number of bytes
            // - must be done by the caller.
        }
        else if (sendData && recvData)
        {
            // recvCounts[0] may be invalid - use sendCount instead
            std::memmove(recvData, sendData, sendCount*sizeof(Type));
        }
        return;
    }

    const label np = UPstream::nProcs(communicator);

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        if (immediate)
        {
            Perr<< "** MPI_Igatherv (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Gatherv (blocking):";
        }
        Perr<< " np:" << np
            << " recvCounts:" << recvCounts
            << " recvOffsets:" << recvOffsets
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    if
    (
        UPstream::master(communicator)
     && (recvCounts.size() != np || recvOffsets.size() < np)
    )
    {
        // Note: allow offsets to be e.g. 1 larger than nProc so we
        // can easily loop over the result

        FatalErrorInFunction
            << "Have " << np << " ranks, but recvCounts:" << recvCounts.size()
            << " or recvOffsets:" << recvOffsets.size()
            << " is too small!"
            << Foam::abort(FatalError);
    }

    // Ensure send/recv consistency on master
    if (UPstream::master(communicator) && !recvCounts[0])
    {
        sendCount = 0;
    }


    int returnCode(MPI_ERR_UNKNOWN);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        returnCode =
            MPI_Igatherv
            (
                const_cast<Type*>(sendData),
                sendCount,
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        returnCode =
            MPI_Gatherv
            (
                const_cast<Type*>(sendData),
                sendCount,
                datatype,
                recvData,
                const_cast<int*>(recvCounts.cdata()),
                const_cast<int*>(recvOffsets.cdata()),
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[communicator]
            );

        profilingPstream::addGatherTime();
    }

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalErrorInFunction << "MPI Gatherv ";
        if (immediate) FatalError<< "(non-blocking) ";

        FatalError
            << "[comm: " << communicator << "] failed."
            << " sendCount " << sendCount
            << " recvCounts " << recvCounts
            << Foam::abort(FatalError);
    }
}


template<class Type>
void Foam::PstreamDetail::scatterv
(
    const Type* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,

    Type* recvData,
    int recvCount,

    MPI_Datatype datatype,
    const int communicator,

    UPstream::Request* req
)
{
    PstreamGlobals::reset_request(req);

    const bool immediate = (req);

    if (!UPstream::is_rank(communicator))
    {
        return;
    }
    else if (!UPstream::is_parallel(communicator))
    {
        if constexpr (std::is_void_v<Type>)
        {
            // Cannot copy data here since we don't know the number of bytes
            // - must be done by the caller.
        }
        else if (sendData && recvData)
        {
            std::memmove(recvData, sendData, recvCount*sizeof(Type));
        }
        return;
    }

    const label np = UPstream::nProcs(communicator);

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        if (immediate)
        {
            Perr<< "** MPI_Iscatterv (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Scatterv (blocking):";
        }
        Perr<< " np:" << np
            << " sendCounts:" << sendCounts
            << " sendOffsets:" << sendOffsets
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }

    if
    (
        UPstream::master(communicator)
     && (sendCounts.size() != np || sendOffsets.size() < np)
    )
    {
        // Note: allow offsets to be e.g. 1 larger than nProc so we
        // can easily loop over the result

        FatalErrorInFunction
            << "Have " << np << " ranks, but sendCounts:" << sendCounts.size()
            << " or sendOffsets:" << sendOffsets.size()
            << " is too small!"
            << Foam::abort(FatalError);
    }


    int returnCode(MPI_ERR_UNKNOWN);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        returnCode =
            MPI_Iscatterv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                recvCount,
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        returnCode =
            MPI_Scatterv
            (
                const_cast<Type*>(sendData),
                const_cast<int*>(sendCounts.cdata()),
                const_cast<int*>(sendOffsets.cdata()),
                datatype,
                recvData,
                recvCount,
                datatype,
                0,  // (root rank) == UPstream::masterNo()
                PstreamGlobals::MPICommunicators_[communicator]
            );

        profilingPstream::addScatterTime();
    }

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalErrorInFunction << "MPI Scatterv ";
        if (immediate) FatalError<< "(non-blocking) ";

        FatalError
            << "[comm: " << communicator << "] failed."
            << " sendCounts " << sendCounts
            << " sendOffsets " << sendOffsets
            << Foam::abort(FatalError);
    }
}


template<class Type>
void Foam::PstreamDetail::allGather
(
    Type* allData,
    int count,

    MPI_Datatype datatype,
    const int communicator,

    UPstream::Request* req
)
{
    PstreamGlobals::reset_request(req);

    const bool immediate = (req);

    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do - ignore
        return;
    }

    if (FOAM_UNLIKELY(PstreamGlobals::warnCommunicator(communicator)))
    {
        if (immediate)
        {
            Perr<< "** MPI_Iallgather (non-blocking):";
        }
        else
        {
            Perr<< "** MPI_Allgather (blocking):";
        }
        Perr<< " numProc:" << UPstream::nProcs(communicator)
            << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Perr);
    }


    int returnCode(MPI_ERR_UNKNOWN);

#if defined(MPI_VERSION) && (MPI_VERSION >= 3)
    if (immediate)
    {
        // MPI-3 : eg, openmpi-1.7 (2013) and later
        profilingPstream::beginTiming();
        MPI_Request request;

        returnCode =
            MPI_Iallgather
            (
                MPI_IN_PLACE, count, datatype,
                allData, count, datatype,
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            );

        PstreamGlobals::push_request(request, req);
        profilingPstream::addRequestTime();
    }
    else
#endif
    {
        profilingPstream::beginTiming();

        returnCode =
            MPI_Allgather
            (
                MPI_IN_PLACE, count, datatype,
                allData, count, datatype,
                PstreamGlobals::MPICommunicators_[communicator]
            );

        // Is actually gather/scatter but we can't split it apart
        profilingPstream::addGatherTime();
    }

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalErrorInFunction << "MPI Allgather ";
        if (immediate) FatalError<< "(non-blocking) ";

        FatalError
            << "[comm: " << communicator << "] failed."
            << Foam::abort(FatalError);
    }
}


// ************************************************************************* //
