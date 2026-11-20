/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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

#include "UPstream.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"
#include "UPstreamWrapping.H"
#include "collatedFileOperation.H"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <numeric>
#include <string>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// The min value and default for MPI buffer length
constexpr int minBufLen = 20000000;

// Track size of attached MPI buffer
static int attachedBufLen = 0;

// Track if we initialized MPI
static bool ourMpi = false;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// Attach user-defined send buffer
static void attachOurBuffers()
{
#ifndef SGIMPI
    if (attachedBufLen)
    {
        return;  // Already attached
    }

    // Use UPstream::mpiBufferSize (optimisationSwitch),
    // but allow override with MPI_BUFFER_SIZE env variable (int value)

    int len = 0;

    const std::string str(Foam::getEnv("MPI_BUFFER_SIZE"));
    if (str.empty() || !Foam::read(str, len) || len <= 0)
    {
        len = Foam::UPstream::mpiBufferSize;
    }

    if (len < minBufLen)
    {
        len = minBufLen;
    }

    char* buf = new char[len];

    if (MPI_SUCCESS == MPI_Buffer_attach(buf, len))
    {
        // Properly attached
        attachedBufLen = len;

        if (Foam::UPstream::debug)
        {
            Foam::Perr<< "UPstream::init : buffer-size " << len << '\n';
        }
    }
    else
    {
        delete[] buf;
        Foam::Perr<< "UPstream::init : could not attach buffer\n";
    }
#endif
}


// Remove an existing user-defined send buffer
// IMPORTANT:
//     This operation will block until all messages currently in the
//     buffer have been transmitted.
static void detachOurBuffers()
{
#ifndef SGIMPI
    if (!attachedBufLen)
    {
        return;  // Nothing to detach
    }

    // Some MPI notes suggest that the return code is MPI_SUCCESS when
    // no buffer is attached.
    // Be extra careful and require a non-zero size as well.

    char* buf = nullptr;
    int len = 0;

    if (MPI_SUCCESS == MPI_Buffer_detach(&buf, &len) && len)
    {
        // This was presumably the buffer that we attached
        // and not someone else.
        delete[] buf;
    }

    // Nothing attached
    attachedBufLen = 0;
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NOTE:
// valid parallel options vary between implementations, but flag common ones.
// if they are not removed by MPI_Init(), the subsequent argument processing
// will notice that they are wrong
void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
    validParOptions.insert("machinefile", "machine file");
}


bool Foam::UPstream::initNull()
{
    int flag = 0;

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init\n"
            << Foam::abort(FatalError);

        return false;
    }

    MPI_Initialized(&flag);
    if (flag)
    {
        if (UPstream::debug)
        {
            Perr<< "UPstream::initNull : was already initialized\n";
        }
    }
    else
    {
        // Not already initialized

        MPI_Init_thread
        (
            nullptr,    // argc
            nullptr,    // argv
            MPI_THREAD_SINGLE,
            &flag       // provided_thread_support
        );

        ourMpi = true;
    }

    // Could also attach buffers etc.

    return true;
}


bool Foam::UPstream::init(int& argc, char**& argv, const bool needsThread)
{
    int flag = 0;
    int provided_thread_support = 0;

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init" << endl
            << Foam::abort(FatalError);

        return false;
    }

    MPI_Initialized(&flag);
    if (flag)
    {
        // Already initialized.
        // Warn if we've called twice, but skip if initialized externally

        if (ourMpi)
        {
            WarningInFunction
                << "MPI was already initialized - cannot perform MPI_Init" << nl
                << "This could indicate an application programming error!"
                << endl;

            return true;
        }
        else if (UPstream::debug)
        {
            Perr<< "UPstream::init : was already initialized\n";
        }

        MPI_Query_thread(&provided_thread_support);
    }
    else
    {
        // (SINGLE | FUNNELED | SERIALIZED | MULTIPLE)
        int required_thread_support =
        (
            needsThread
          ? MPI_THREAD_MULTIPLE
          : MPI_THREAD_SINGLE
        );

        MPI_Init_thread
        (
            &argc,
            &argv,
            required_thread_support,
           &provided_thread_support
        );

        ourMpi = true;
    }

    // Define data type mappings and user data types. Defined now so that
    // any OpenFOAM Pstream operations may make immediate use of them.
    PstreamGlobals::initDataTypes();
    PstreamGlobals::initOpCodes();

    if (UPstream::debug)
    {
        PstreamGlobals::printDataTypes();
    }


    // Check argument list for any of the following:
    // - local world
    //    -> Extract world name and filter out '-world <name>' from argv list
    // - mpi-no-comm-dup option
    //    -> disable initial comm_dup and filter out the option
    // - mpi-split-by-appnum option
    //   -> disable initial comm_dup, select split-by-appnum
    // and filter out the option

    // Default handling of initial MPI_Comm_dup(MPI_COMM_WORLD,...)
    UPstream::noInitialCommDup_ = false;
    bool split_by_appnum = false;

    // Local world name
    word worldName;

    for (int argi = 1; argi < argc; ++argi)
    {
        const char *optName = argv[argi];
        if (optName[0] != '-')
        {
            continue;
        }
        ++optName;  // Looks like an option, skip leading '-'

        if (strcmp(optName, "world") == 0)
        {
            if (argi+1 >= argc)
            {
                FatalErrorInFunction
                    << "Missing world name for option '-world'" << nl
                    << Foam::abort(FatalError);
            }
            worldName = argv[argi+1];

            // Remove two arguments (-world name)
            for (int i = argi+2; i < argc; ++i)
            {
                argv[i-2] = argv[i];
            }
            argc -= 2;
            --argi;  // re-examine
        }
        else if (strcmp(optName, "mpi-no-comm-dup") == 0)
        {
            UPstream::noInitialCommDup_ = true;

            // Remove one argument
            for (int i = argi+1; i < argc; ++i)
            {
                argv[i-1] = argv[i];
            }
            --argc;
            --argi;  // re-examine
        }
        else if (strcmp(optName, "mpi-split-by-appnum") == 0)
        {
            split_by_appnum = true;
            UPstream::noInitialCommDup_ = true;

            // Remove one argument
            for (int i = argi+1; i < argc; ++i)
            {
                argv[i-1] = argv[i];
            }
            --argc;
            --argi;  // re-examine
        }
    }

    const bool hasLocalWorld(!worldName.empty());

    // commSplit override: the OpenFOAM ranks always leave MPI_COMM_WORLD by
    // an MPI_Comm_split with a fixed colour (1), so that another program in
    // the same MPMD job (eg, ElmerSolver, colour 0) is never taken for an
    // OpenFOAM rank. Explicit -world or -mpi-split-by-appnum keep their
    // meaning; otherwise the colour split replaces the initial comm_dup.
    bool split_by_colour = false;
    if (!hasLocalWorld && !split_by_appnum)
    {
        split_by_colour = true;
        split_by_appnum = true;
        UPstream::noInitialCommDup_ = true;
    }

    if (hasLocalWorld && split_by_appnum)
    {
        FatalErrorInFunction
            << "Cannot specify both -world and -mpi-split-by-appnum" << nl
            << Foam::abort(FatalError);
    }

    int numProcs = 0, globalRanki = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &globalRanki);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    if (UPstream::debug)
    {
        Perr<< "UPstream::init :"
            << " thread-support : requested:" << needsThread
            << " provided:"
            << (
                   (provided_thread_support == MPI_THREAD_SINGLE)
                 ? "SINGLE"
                 : (provided_thread_support == MPI_THREAD_SERIALIZED)
                 ? "SERIALIZED"
                 : (provided_thread_support == MPI_THREAD_MULTIPLE)
                 ? "MULTIPLE"
                 : "other"
               )
            << " procs:" << numProcs
            << " rank:" << globalRanki
            << " world:" << worldName << endl;
    }

    if (numProcs <= 1 && !(hasLocalWorld || split_by_appnum))
    {
        FatalErrorInFunction
            << "attempt to run parallel on 1 processor"
            << Foam::abort(FatalError);
    }

    // Initialise parallel structure
    setParRun(numProcs, provided_thread_support == MPI_THREAD_MULTIPLE);

    if (hasLocalWorld)
    {
        // Using local worlds.
        // During startup, so commWorld() == commGlobal()
        const auto mpiGlobalComm =
            PstreamGlobals::MPICommunicators_[UPstream::commGlobal()];

        // Gather the names of all worlds and determine unique names/indices.
        //
        // Minimize communication and use low-level MPI to avoid relying on any
        // OpenFOAM structures which not yet have been created

        {
            // Include a trailing nul character in the lengths
            int stride = int(worldName.size()) + 1;

            // Use identical size on all ranks (avoids MPI_Allgatherv)
            MPI_Allreduce
            (
                MPI_IN_PLACE,
               &stride,
                1,
                MPI_INT,
                MPI_MAX,
                mpiGlobalComm
            );

            // Gather as an extended C-string with embedded nul characters
            auto buffer_storage = std::make_unique<char[]>(numProcs*stride);
            char* allStrings = buffer_storage.get();

            // Fill in local value, slot starts at (rank*stride)
            {
                char* slot = (allStrings + (globalRanki*stride));
                std::fill_n(slot, stride, '\0');
                std::copy_n(worldName.data(), worldName.size(), slot);
            }

            // Gather everything into the extended C-string
            MPI_Allgather
            (
                MPI_IN_PLACE, 0, MPI_CHAR,
                allStrings, stride, MPI_CHAR,
                mpiGlobalComm
            );

            worldIDs_.resize_nocopy(numProcs);

            // Transcribe and compact (unique world names)
            DynamicList<word> uniqWorlds(numProcs);

            for (label proci = 0; proci < numProcs; ++proci)
            {
                // Create from C-string at slot=(rank*stride),
                // relying on the embedded nul chars
                word world(allStrings + (proci*stride));

                worldIDs_[proci] = uniqWorlds.find(world);

                if (worldIDs_[proci] == -1)
                {
                    worldIDs_[proci] = uniqWorlds.size();
                    uniqWorlds.push_back(std::move(world));
                }
            }

            allWorlds_ = std::move(uniqWorlds);
        }

        const label myWorldId = worldIDs_[globalRanki];

        DynamicList<label> subRanks;
        forAll(worldIDs_, proci)
        {
            if (worldIDs_[proci] == myWorldId)
            {
                subRanks.push_back(proci);
            }
        }

        // New local-world communicator with comm-global as its parent.
        // - the updated (const) world comm does not change after this.

        UPstream::constWorldComm_ =
            UPstream::newCommunicator(UPstream::commGlobal(), subRanks);

        UPstream::worldComm = UPstream::constWorldComm_;
        UPstream::warnComm = UPstream::constWorldComm_;

        const int worldRanki = UPstream::myProcNo(UPstream::constWorldComm_);

        // MPI_COMM_SELF : the processor number wrt the new world communicator
        if (procIDs_[UPstream::commSelf()].size())
        {
            procIDs_[UPstream::commSelf()].front() = worldRanki;
        }

        // Name the old world communicator as '<openfoam:global>'
        // - it is the inter-world communicator
        if (MPI_COMM_NULL != mpiGlobalComm)
        {
            MPI_Comm_set_name(mpiGlobalComm, "<openfoam:global>");
        }

        const auto mpiWorldComm =
            PstreamGlobals::MPICommunicators_[UPstream::constWorldComm_];

        if (MPI_COMM_NULL != mpiWorldComm)
        {
            MPI_Comm_set_name(mpiWorldComm, ("world=" + worldName).data());
        }

        PstreamGlobals::MPI_COMM_FOAM = mpiWorldComm;

        if (UPstream::debug)
        {
            // Check
            int newRanki, newSize;
            MPI_Comm_rank(mpiWorldComm, &newRanki);
            MPI_Comm_size(mpiWorldComm, &newSize);

            Perr<< "UPstream::init : in world:" << worldName
                << " using local communicator:" << constWorldComm_
                << " rank " << newRanki << " of " << newSize << endl;
        }

        // Override Pout prefix (move to setParRun?)
        Pout.prefix() = '[' + worldName + '/' + Foam::name(worldRanki) + "] ";
        Perr.prefix() = Pout.prefix();
    }
    else if (split_by_appnum)
    {
        // Splitting by APPNUM.
        //
        // During startup, so commWorld() == commGlobal() and both are
        // guaranteed to be MPI_COMM_WORLD since the logic automatically
        // sets UPstream::noInitialCommDup_ = true (ie, no MPI_Comm_dup)

        const auto mpiGlobalComm =
            PstreamGlobals::MPICommunicators_[UPstream::commGlobal()];

        int appNum(0);

        if (split_by_colour)
        {
            // commSplit override: fixed OpenFOAM colour
            appNum = 1;
        }
        else
        {
            void* val;
            int flag;

            MPI_Comm_get_attr(mpiGlobalComm, MPI_APPNUM, &val, &flag);
            if (flag)
            {
                appNum = *static_cast<int*>(val);
            }
            else
            {
                appNum = 0;
                Perr<< "UPstream::init : used -mpi-split-by-appnum"
                       " with a single application??" << endl;
            }
        }

        // New world communicator with comm-global as its parent.
        // - the updated (const) world comm does not change after this.

        // Using MPI_APPNUM as the colour for splitting with MPI_Comm_split.
        // Do **NOT** use Allgather+Comm_create_group two-step process here
        // since other applications will not expect that (ie, deadlock)

        UPstream::constWorldComm_ =
            UPstream::splitCommunicator(UPstream::commGlobal(), appNum, false);

        UPstream::worldComm = UPstream::constWorldComm_;
        UPstream::warnComm = UPstream::constWorldComm_;

        const int worldRanki = UPstream::myProcNo(UPstream::constWorldComm_);

        // MPI_COMM_SELF : the processor number wrt the new world communicator
        if (procIDs_[UPstream::commSelf()].size())
        {
            procIDs_[UPstream::commSelf()].front() = worldRanki;
        }

        // Name the old world communicator as '<openfoam:global>'
        // - it is the inter-world communicator
        if (MPI_COMM_NULL != mpiGlobalComm)
        {
            MPI_Comm_set_name(mpiGlobalComm, "<openfoam:global>");
        }

        const auto mpiWorldComm =
            PstreamGlobals::MPICommunicators_[UPstream::constWorldComm_];

        const word commName
        (
            split_by_colour ? word("<openfoam:commSplit>")
          : word("app=" + Foam::name(appNum))
        );

        if (MPI_COMM_NULL != mpiWorldComm)
        {
            MPI_Comm_set_name(mpiWorldComm, commName.data());
        }

        PstreamGlobals::MPI_COMM_FOAM = mpiWorldComm;

        if (split_by_colour && UPstream::nProcs(UPstream::constWorldComm_) <= 1)
        {
            FatalErrorInFunction
                << "attempt to run parallel on 1 processor"
                << Foam::abort(FatalError);
        }

        if (UPstream::debug)
        {
            // Check
            int newRanki, newSize;
            MPI_Comm_rank(mpiWorldComm, &newRanki);
            MPI_Comm_size(mpiWorldComm, &newSize);

            Perr<< "UPstream::init : app:" << appNum
                << " using local communicator:" << constWorldComm_
                << " rank " << newRanki << " of " << newSize << endl;
        }

        // Override Pout prefix (move to setParRun?)
        Pout.prefix() = '[' + commName + '/' + Foam::name(worldRanki) + "] ";
        Perr.prefix() = Pout.prefix();
    }
    else
    {
        // All processors use world 0
        worldIDs_.resize_nocopy(numProcs);
        worldIDs_ = 0;

        const auto mpiWorldComm =
            PstreamGlobals::MPICommunicators_[UPstream::constWorldComm_];

        // Name the world communicator as '<openfoam:world>'
        if (MPI_COMM_NULL != mpiWorldComm)
        {
            MPI_Comm_set_name(mpiWorldComm, "<openfoam:world>");
        }

        PstreamGlobals::MPI_COMM_FOAM = mpiWorldComm;
    }


    // Define inter-node and intra-node communicators
    if (UPstream::nodeCommsControl_ >= 4)
    {
        // Debugging: split with given number per node
        setHostCommunicators(UPstream::nodeCommsControl_);
    }
    #ifndef MSMPI_VER  /* Uncertain if this would work with MSMPI */
    else if (UPstream::nodeCommsControl_ == 2)
    {
        // Defined based on shared-memory hardware information
        setSharedMemoryCommunicators();
    }
    #endif
    else
    {
        // Defined based on hostname, even if nominally disabled
        setHostCommunicators();
    }


    // Provide some names for these communicators
    if (MPI_COMM_NULL != PstreamGlobals::MPICommunicators_[commInterNode_])
    {
        MPI_Comm_set_name
        (
            PstreamGlobals::MPICommunicators_[commInterNode_],
            "<openfoam:inter-node>"
        );
    }
    if (MPI_COMM_NULL != PstreamGlobals::MPICommunicators_[commLocalNode_])
    {
        MPI_Comm_set_name
        (
            PstreamGlobals::MPICommunicators_[commLocalNode_],
            "<openfoam:local-node>"
        );
    }

    attachOurBuffers();

    return true;
}


void Foam::UPstream::shutdown(int errNo)
{
    int flag = 0;

    MPI_Initialized(&flag);
    if (!flag)
    {
        // MPI not initialized - we have nothing to do
        return;
    }

    MPI_Finalized(&flag);
    if (flag)
    {
        // MPI already finalized - we have nothing to do
        if (ourMpi)
        {
            WarningInFunction
                << "MPI was already finalized (by a connected program?)\n";
        }
        else if (UPstream::debug && errNo == 0)
        {
            Perr<< "UPstream::shutdown : was already finalized\n";
        }
        ourMpi = false;
        return;
    }

    if (!ourMpi)
    {
        WarningInFunction
            << "Finalizing MPI, but was initialized elsewhere\n";
    }
    ourMpi = false;


    // Abort - stop now, without any final synchonization steps!
    // -----

    if (errNo != 0)
    {
        UPstream::abort(errNo);
        return;
    }


    // Regular cleanup
    // ---------------

    if (UPstream::debug)
    {
        Perr<< "UPstream::shutdown\n";
    }

    // Check for any outstanding requests
    {
        label nOutstanding = 0;

        for (MPI_Request request : PstreamGlobals::outstandingRequests_)
        {
            if (MPI_REQUEST_NULL != request)
            {
                // TBD: MPI_Cancel(&request); MPI_Request_free(&request);
                ++nOutstanding;
            }
        }

        if (nOutstanding)
        {
            WarningInFunction
                << "Still have " << nOutstanding
                << " outstanding MPI requests."
                << " Should not happen for a normal code exit."
                << endl;
        }

        PstreamGlobals::outstandingRequests_.clear();
    }


    {
        detachOurBuffers();

        forAllReverse(myProcNo_, communicator)
        {
            freeCommunicatorComponents(communicator);
        }
    }

    // Free any user data types
    PstreamGlobals::deinitDataTypes();
    PstreamGlobals::deinitOpCodes();

    MPI_Finalize();
}


void Foam::UPstream::exit(int errNo)
{
    UPstream::shutdown(errNo);
    std::exit(errNo);
}


void Foam::UPstream::abort(int errNo)
{
    MPI_Comm abortComm = MPI_COMM_WORLD;

    // TBD: only abort on our own communicator?
    #if 0
    const label index = UPstream::commGlobal();

    if (index > 0 && index < PstreamGlobals::MPICommunicators_.size())
    {
        abortComm = PstreamGlobals::MPICommunicators_[index];
        if (MPI_COMM_NULL == abortComm)
        {
            abortComm = MPI_COMM_WORLD;
        }
    }
    #endif

    MPI_Abort(abortComm, errNo);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::allocateCommunicatorComponents
(
    const label parentIndex,
    const label index
)
{
    PstreamGlobals::initCommunicator(index);

    int returnCode = MPI_SUCCESS;

    if (parentIndex == -1)
    {
        // Global communicator. Same as world communicator for single-world

        if (index != UPstream::commGlobal())
        {
            FatalErrorInFunction
                << "base world communicator should always be index "
                << UPstream::commGlobal()
                << Foam::exit(FatalError);
        }
        auto& mpiNewComm = PstreamGlobals::MPICommunicators_[index];

        if (UPstream::noInitialCommDup_)
        {
            PstreamGlobals::pendingMPIFree_[index] = false;
            PstreamGlobals::MPICommunicators_[index] = MPI_COMM_WORLD;
        }
        else
        {
            PstreamGlobals::pendingMPIFree_[index] = true;
            MPI_Comm_dup(MPI_COMM_WORLD, &mpiNewComm);
        }

        MPI_Comm_rank(mpiNewComm, &myProcNo_[index]);

        // Set the number of ranks to the actual number
        int numProcs = 0;
        MPI_Comm_size(mpiNewComm, &numProcs);

        // identity [0-numProcs], as 'int'
        procIDs_[index].resize_nocopy(numProcs);
        std::iota(procIDs_[index].begin(), procIDs_[index].end(), 0);
    }
    else if (parentIndex == -2)
    {
        // MPI_COMM_SELF

        PstreamGlobals::pendingMPIFree_[index] = false;
        PstreamGlobals::MPICommunicators_[index] = MPI_COMM_SELF;

        MPI_Comm_rank(MPI_COMM_SELF, &myProcNo_[index]);

        // For MPI_COMM_SELF : the process IDs within the world communicator.
        // Uses MPI_COMM_WORLD in case called before UPstream::commGlobal()
        // was initialized

        procIDs_[index].resize_nocopy(1);
        MPI_Comm_rank(MPI_COMM_WORLD, &procIDs_[index].front());
    }
    else
    {
        // General sub-communicator.
        // Create based on the groupings predefined by procIDs_

        const auto mpiParentComm =
            PstreamGlobals::MPICommunicators_[parentIndex];

        auto& mpiNewComm =
            PstreamGlobals::MPICommunicators_[index];

        PstreamGlobals::pendingMPIFree_[index] = true;

        // Starting from parent
        MPI_Group parent_group;
        MPI_Comm_group(mpiParentComm, &parent_group);

        MPI_Group active_group;
        MPI_Group_incl
        (
            parent_group,
            procIDs_[index].size(),
            procIDs_[index].cdata(),
           &active_group
        );

        #if defined(MSMPI_VER)
        // ms-mpi (10.0 and others?) does not have MPI_Comm_create_group
        MPI_Comm_create
        (
            mpiParentComm,
            active_group,
           &mpiNewComm
        );
        #else
        // Create new communicator for this group
        MPI_Comm_create_group
        (
            mpiParentComm,
            active_group,
            UPstream::msgType(),
           &mpiNewComm
        );
        #endif

        // Groups not needed after this...
        MPI_Group_free(&parent_group);
        MPI_Group_free(&active_group);

        if (MPI_COMM_NULL == mpiNewComm)
        {
            // This process is not involved in the new communication pattern
            myProcNo_[index] = -1;
            PstreamGlobals::pendingMPIFree_[index] = false;

            // ~~~~~~~~~
            // IMPORTANT
            // ~~~~~~~~~
            // Always retain knowledge of the inter-node leaders,
            // even if this process is not on that communicator.
            // This will help when constructing topology-aware communication.

            if (index != commInterNode_)
            {
                procIDs_[index].clear();
            }
        }
        else
        {
            returnCode = MPI_Comm_rank(mpiNewComm, &myProcNo_[index]);

            if (FOAM_UNLIKELY(MPI_SUCCESS != returnCode))
            {
                FatalErrorInFunction
                    << "Problem :"
                    << " when allocating communicator at " << index
                    << " from ranks " << flatOutput(procIDs_[index])
                    << " of parent " << parentIndex
                    << " cannot find my own rank"
                    << Foam::exit(FatalError);
            }
        }
    }
}


void Foam::UPstream::dupCommunicatorComponents
(
    const label parentIndex,
    const label index
)
{
    PstreamGlobals::initCommunicator(index);

    PstreamGlobals::pendingMPIFree_[index] = true;
    MPI_Comm_dup
    (
        PstreamGlobals::MPICommunicators_[parentIndex],
       &PstreamGlobals::MPICommunicators_[index]
    );

    myProcNo_[index] = myProcNo_[parentIndex];
    procIDs_[index] = procIDs_[parentIndex];
}


void Foam::UPstream::splitCommunicatorComponents
(
    const label parentIndex,
    const label index,
    int colour,
    const bool two_step
)
{
    PstreamGlobals::initCommunicator(index);

    // ------------------------------------------------------------------------
    // Create sub-communicator according to its colouring
    //     => MPI_Comm_split().
    // Since other parts of OpenFOAM may still need a view of the siblings:
    //     => MPI_Group_translate_ranks().
    //
    // The MPI_Group_translate_ranks() step can be replaced with an
    // MPI_Allgather() of the involved parent ranks (since we alway maintain
    // the relative rank order when splitting).
    //
    // Since MPI_Comm_split() already does an MPI_Allgather() internally
    // to pick out the colours (and do any sorting), we can simply
    // do the same thing:
    //
    // Do the Allgather first and pickout identical colours to define the
    // group and create a communicator based on that.
    //
    // This is no worse than the Allgather communication overhead of using
    // MPI_Comm_split() and saves the extra translate_ranks step.
    // ------------------------------------------------------------------------

    const auto mpiParentComm = PstreamGlobals::MPICommunicators_[parentIndex];

    int parentRank = 0;
    int parentSize = 0;
    MPI_Comm_rank(mpiParentComm, &parentRank);
    MPI_Comm_size(mpiParentComm, &parentSize);

    auto& procIds = procIDs_[index];
    myProcNo_[index] = -1;

    if (two_step)
    {
        // First gather the colours
        procIds.resize_nocopy(parentSize);
        procIds[parentRank] = colour;

        MPI_Allgather
        (
            MPI_IN_PLACE, 0, MPI_INT,
            procIds.data(), 1, MPI_INT,
            mpiParentComm
        );

        if (colour < 0)
        {
            // Not involved
            procIds.clear();
        }
        else
        {
            // Select ranks based on the matching colour
            int nranks = 0;
            for (int i = 0; i < parentSize; ++i)
            {
                if (procIds[i] == colour)
                {
                    procIds[nranks++] = i;
                }
            }
            procIds.resize(nranks);
        }

        allocateCommunicatorComponents(parentIndex, index);
    }
    else
    {
        auto& mpiNewComm = PstreamGlobals::MPICommunicators_[index];

        MPI_Comm_split
        (
            mpiParentComm,
            (colour >= 0 ? colour : MPI_UNDEFINED),
            0,  // maintain relative ordering
           &mpiNewComm
        );

        if (MPI_COMM_NULL == mpiNewComm)
        {
            // Not involved
            PstreamGlobals::pendingMPIFree_[index] = false;
            procIds.clear();
        }
        else
        {
            PstreamGlobals::pendingMPIFree_[index] = true;
            MPI_Comm_rank(mpiNewComm, &myProcNo_[index]);

            // Starting from parent
            MPI_Group parent_group;
            MPI_Comm_group(mpiParentComm, &parent_group);

            MPI_Group new_group;
            MPI_Comm_group(mpiNewComm, &new_group);

            // Parent ranks: identity map
            List<int> parentIds(parentSize);
            std::iota(parentIds.begin(), parentIds.end(), 0);

            // New ranks:
            procIds.resize_nocopy(parentSize);
            procIds = -1;  // Some extra safety...

            MPI_Group_translate_ranks
            (
                parent_group, parentSize, parentIds.data(),
                new_group, procIds.data()
            );

            // Groups not needed after this...
            MPI_Group_free(&parent_group);
            MPI_Group_free(&new_group);

            // The corresponding ranks.
            // - since old ranks are an identity map, can just use position

            int nranks = 0;
            for (int i = 0; i < parentSize; ++i)
            {
                // Exclude MPI_UNDEFINED and MPI_PROC_NULL etc...
                if (procIds[i] >= 0 && procIds[i] < parentSize)
                {
                    procIds[nranks++] = i;
                }
            }
            procIds.resize(nranks);
        }
    }
}


void Foam::UPstream::freeCommunicatorComponents(const label index)
{
    if (UPstream::debug)
    {
        Perr<< "freeCommunicatorComponents: " << index
            << " from " << PstreamGlobals::MPICommunicators_.size() << endl;
    }

    // Only free communicators that we have specifically allocated ourselves
    //
    // Bounds checking needed since there are no UPstream communicator indices
    // when MPI is initialized outside of OpenFOAM

    if
    (
        (index >= 0 && index < PstreamGlobals::MPICommunicators_.size())
     && PstreamGlobals::pendingMPIFree_[index]
    )
    {
        PstreamGlobals::pendingMPIFree_[index] = false;

        // Free communicator. Sets communicator to MPI_COMM_NULL
        if (MPI_COMM_NULL != PstreamGlobals::MPICommunicators_[index])
        {
            MPI_Comm_free(&PstreamGlobals::MPICommunicators_[index]);
        }
    }
}


bool Foam::UPstream::setSharedMemoryCommunicators()
{
    // Uses the world communicator (not global communicator)

    // Skip if non-parallel
    if (!UPstream::parRun())
    {
        numNodes_ = 1;
        return false;
    }

    if (FOAM_UNLIKELY(commInterNode_ >= 0 || commLocalNode_ >= 0))
    {
        // Failed sanity check
        FatalErrorInFunction
            << "Node communicator(s) already created!" << endl
            << Foam::abort(FatalError);
        return false;
    }

    commInterNode_ = getAvailableCommIndex(constWorldComm_);
    commLocalNode_ = getAvailableCommIndex(constWorldComm_);

    PstreamGlobals::initCommunicator(commInterNode_);
    PstreamGlobals::initCommunicator(commLocalNode_);

    // Overwritten later
    myProcNo_[commInterNode_] = UPstream::masterNo();
    myProcNo_[commLocalNode_] = UPstream::masterNo();

    // Sorted order, purely cosmetic
    if (commLocalNode_ < commInterNode_)
    {
        std::swap(commLocalNode_, commInterNode_);
    }

    if (debug)
    {
        Perr<< "Allocating node communicators "
            << commInterNode_ << ", " << commLocalNode_ << nl
            << "    parent : " << constWorldComm_ << nl
            << endl;
    }


    const auto mpiParentComm =
        PstreamGlobals::MPICommunicators_[constWorldComm_];

    auto& mpiLocalNode =
        PstreamGlobals::MPICommunicators_[commLocalNode_];

    int parentRank = 0;
    int parentSize = 0;
    MPI_Comm_rank(mpiParentComm, &parentRank);
    MPI_Comm_size(mpiParentComm, &parentSize);

    List<int> nodeLeaders(parentSize);
    nodeLeaders = -1;

    MPI_Comm_split_type
    (
        mpiParentComm,
        MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
       &mpiLocalNode
    );

    if (FOAM_UNLIKELY(MPI_COMM_NULL == mpiLocalNode))
    {
        // This process is not involved in an intra-host communication?
        // - should never happen!

        const label index = commLocalNode_;
        PstreamGlobals::pendingMPIFree_[index] = false;

        myProcNo_[index] = -1;
        procIDs_[index].clear();

        FatalErrorInFunction
            << "Comm_split_type(shared) failed\n"
            << Foam::abort(FatalError);
    }
    else
    {
        // This process is involved in intra-host communication
        const label index = commLocalNode_;
        auto& procIds = procIDs_[index];

        PstreamGlobals::pendingMPIFree_[index] = true;

        int localRank = 0;
        int localSize = 0;
        MPI_Comm_rank(mpiLocalNode, &localRank);
        MPI_Comm_size(mpiLocalNode, &localSize);

        if (localRank == 0)
        {
            // This process is a host leader - mark its position
            nodeLeaders[parentRank] = parentRank;
        }

        procIds.resize_nocopy(localSize);
        procIds[localRank] = UPstream::myProcNo(UPstream::constWorldComm_);
        // OR: procIds[localRank] = parentRank;

        // Get all of the siblings (within the node)
        MPI_Allgather
        (
            MPI_IN_PLACE, 0, MPI_INT,
            procIds.data(), 1, MPI_INT,
            mpiLocalNode
        );
    }


    // Get all of the host-leader information and find who they are.
    {
        auto& procIds = procIDs_[commInterNode_];

        MPI_Allgather
        (
            MPI_IN_PLACE, 0, MPI_INT,
            nodeLeaders.data(), 1, MPI_INT,
            mpiParentComm
        );

        // Capture the size (number of nodes) before doing anything further
        numNodes_ = std::count_if
        (
            nodeLeaders.cbegin(),
            nodeLeaders.cend(),
            [](int rank){ return (rank >= 0); }
        );

        // ~~~~~~~~~
        // IMPORTANT
        // ~~~~~~~~~
        // Always retain knowledge of the inter-node leaders,
        // even if this process is not on that communicator.
        // This will help when constructing topology-aware communication.

        procIds.resize_nocopy(numNodes_);

        std::copy_if
        (
            nodeLeaders.cbegin(),
            nodeLeaders.cend(),
            procIds.begin(),
            [](int rank){ return (rank >= 0); }
        );
    }

    // From master to host-leader. Ranks between hosts.
    allocateCommunicatorComponents(UPstream::worldComm, commInterNode_);

    return true;
}


void Foam::UPstream::barrier(const int communicator, UPstream::Request* req)
{
    // No-op for non-parallel or not on communicator
    if (!UPstream::is_parallel(communicator))
    {
        PstreamGlobals::reset_request(req);
        return;
    }

    if (req)
    {
        MPI_Request request;

        // Non-blocking
        if
        (
            MPI_Ibarrier
            (
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Ibarrier returned with error"
                << Foam::abort(FatalError);
        }

        *req = UPstream::Request(request);
    }
    else
    {
        // Blocking
        if
        (
            MPI_Barrier
            (
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Barrier returned with error"
                << Foam::abort(FatalError);
        }
    }
}


void Foam::UPstream::send_done
(
    const int toProcNo,
    const int communicator,
    const int tag  // Message tag (must match on receiving side)
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }

    {
        MPI_Send
        (
            nullptr, 0, MPI_BYTE, toProcNo, tag,
            PstreamGlobals::MPICommunicators_[communicator]
        );
    }
}


int Foam::UPstream::wait_done
(
    const int fromProcNo,
    const int communicator,
    const int tag  // Message tag (must match on sending side)
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return -1;
    }
    else if (fromProcNo < 0)
    {
        MPI_Status status;
        MPI_Recv
        (
            nullptr, 0, MPI_BYTE, MPI_ANY_SOURCE, tag,
            PstreamGlobals::MPICommunicators_[communicator],
           &status
        );
        return status.MPI_SOURCE;
    }
    else
    {
        MPI_Recv
        (
            nullptr, 0, MPI_BYTE, fromProcNo, tag,
            PstreamGlobals::MPICommunicators_[communicator],
            MPI_STATUS_IGNORE
        );
        return fromProcNo;
    }
}


std::pair<int,int64_t>
Foam::UPstream::probeMessage
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
    const int tag,  // Message tag
    const int communicator
)
{
    std::pair<int,int64_t> result(-1, 0);

    // No-op for non-parallel or not on communicator
    if (!UPstream::is_parallel(communicator))
    {
        return result;
    }

    const int source = (fromProcNo < 0) ? MPI_ANY_SOURCE : fromProcNo;
    // Supporting MPI_ANY_TAG is not particularly useful...

    int flag(1);  // Flag only modified by MPI_Iprobe() version
    int returnCode = MPI_SUCCESS;
    MPI_Status status;

    profilingPstream::beginTiming();
    if (UPstream::commsTypes::nonBlocking == commsType)
    {
        // Non-blocking
        returnCode = MPI_Iprobe
        (
            source,
            tag,
            PstreamGlobals::MPICommunicators_[communicator],
           &flag,
           &status
        );
    }
    else
    {
        // Standard
        returnCode = MPI_Probe
        (
            source,
            tag,
            PstreamGlobals::MPICommunicators_[communicator],
           &status
        );
    }
    profilingPstream::addProbeTime();

    if (FOAM_UNLIKELY(MPI_SUCCESS != returnCode))
    {
        FatalErrorInFunction
            << "[mpi_probe] : returned with error"
            << Foam::abort(FatalError);
    }
    else if (flag)
    {
        // Unlikely to be used with large amounts of data,
        // but use MPI_Get_elements_x() instead of MPI_Count() anyhow

        MPI_Count num_recv(0);
        MPI_Get_elements_x(&status, MPI_BYTE, &num_recv);

        // Errors
        if (FOAM_UNLIKELY(num_recv == MPI_UNDEFINED || int64_t(num_recv) < 0))
        {
            FatalErrorInFunction
                << "MPI_Get_elements_x() : "
                   "returned undefined or negative value"
                << Foam::abort(FatalError);
        }
        else if (FOAM_UNLIKELY(int64_t(num_recv) > int64_t(INT_MAX)))
        {
            FatalErrorInFunction
                << "MPI_Get_elements_x() : "
                   "count is larger than INT_MAX bytes"
                << Foam::abort(FatalError);
        }


        result.first = status.MPI_SOURCE;
        result.second = int64_t(num_recv);
    }

    return result;
}


// Probe for messages coming from multiple sources.
//
// MPI_Probe will locally block (but no communication) until the
// corresponding message information arrives, which means that
// sequential probing also sequentially blocks. However, since blocking is
// just until message envelopes arrive and it does not transfer the
// contents, it should be reasonable. Using MPI_Iprobe isn't really much
// of an alternative, since that would move all of the looping, bookkeeping
// and polling to the caller!
//
void Foam::UPstream::probeMessages
(
    DynamicList<int64_t>& messageSizes,
    const UList<int>& fromProcs,
    const int tag,  // Message tag: assumed to be identical for all sources
    const int communicator
)
{
    messageSizes.clear();

    // No-op for non-parallel, not on communicator or no sources to probe
    if (!UPstream::is_parallel(communicator) || fromProcs.empty())
    {
        return;
    }

    // Reserve one larger: the caller may wish to create an offsets
    // array from the size information
    messageSizes.reserve_exact(fromProcs.size()+1);

    MPI_Status status;

    profilingPstream::beginTiming();
    for (const int source : fromProcs)
    {
        #ifdef FULLDEBUG
        if (FOAM_UNLIKELY(source < 0))
        {
            // Cannot use MPI_ANY_SOURCE at all here.
            FatalErrorInFunction
                << "Cannot use MPI_ANY_SOURCE when probing multiple sources!\n"
                << Foam::abort(FatalError);
        }
        #endif

        int returnCode = MPI_Probe
        (
            source,
            tag,
            PstreamGlobals::MPICommunicators_[communicator],
           &status
        );

        if (FOAM_UNLIKELY(MPI_SUCCESS != returnCode))
        {
            FatalErrorInFunction
                << "[mpi_probe] : returned with error"
                << Foam::abort(FatalError);
        }
        else
        {
            // Unlikely to be used with large amounts of data,
            // but use MPI_Get_elements_x() instead of MPI_Count() anyhow

            MPI_Count num_recv(0);
            MPI_Get_elements_x(&status, MPI_BYTE, &num_recv);

            // Errors
            if
            (
                FOAM_UNLIKELY
                (num_recv == MPI_UNDEFINED || int64_t(num_recv) < 0)
            )
            {
                FatalErrorInFunction
                    << "MPI_Get_elements_x() : "
                       "returned undefined or negative value"
                    << Foam::abort(FatalError);
            }
            else if (FOAM_UNLIKELY(int64_t(num_recv) > int64_t(INT_MAX)))
            {
                FatalErrorInFunction
                    << "MPI_Get_elements_x() : "
                       "count is larger than INT_MAX bytes"
                    << Foam::abort(FatalError);
            }

            messageSizes.push_back(int64_t(num_recv));
        }
    }
    profilingPstream::addProbeTime();
}


// ************************************************************************* //
