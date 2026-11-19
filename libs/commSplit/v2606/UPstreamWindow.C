/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "UPstreamWindow.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::Window::Window() noexcept
:
    UPstream::Window(MPI_WIN_NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Window::good() const noexcept
{
    return MPI_WIN_NULL != PstreamUtils::Cast::to_mpi(*this);
}


void Foam::UPstream::Window::reset() noexcept
{
    *this = UPstream::Window(MPI_WIN_NULL);
}


int Foam::UPstream::Window::size() const
{
    int val = 0;

    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);
    MPI_Group group;

    // Get num of ranks from the group information
    if
    (
        (MPI_WIN_NULL != win)
     && (MPI_SUCCESS == MPI_Win_get_group(win, &group))
    )
    {
        if (MPI_SUCCESS != MPI_Group_size(group, &val))
        {
            val = 0;
        }
        MPI_Group_free(&group);
    }

    return val;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//
// Allocate a local or shared memory window.
// Uses MPI_Win_allocate() or MPI_Win_allocate_shared(), respectively.
//
static std::pair<void*,int64_t>
call_window_allocate
(
    Foam::UPstream::Window* self,
    MPI_Comm communicator,

    // [in] number of elements (not bytes)
    std::streamsize num_elements,
    // [in] size of each element == sizeof(Type)
    const int disp_unit,
    const bool shared
)
{
    using namespace Foam;

    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        *self = UPstream::Window(MPI_WIN_NULL);
        return {nullptr, 0};
    }

    // if (FOAM_UNLIKELY(MPI_COMM_NULL == communicator))
    // {
    //     FatalErrorInFunction
    //         << "Attempt to use NULL communicator"
    //         << Foam::abort(FatalError);
    //         return false;
    // }

    MPI_Win win = PstreamUtils::Cast::to_mpi(*self);

    // Stringent handling of existing windows
    if (FOAM_UNLIKELY(MPI_WIN_NULL != win))
    {
        FatalErrorInFunction
            << "Window already exists. Use close() first"
            << Foam::abort(FatalError);
        return {nullptr, 0};
    }

    int returnCode(MPI_SUCCESS);
    void *baseptr = nullptr;

    if (shared)
    {
        returnCode = MPI_Win_allocate_shared
        (
            // From num elements -> num of bytes
            std::streamsize(num_elements * disp_unit),
            disp_unit,
            MPI_INFO_NULL,
            communicator,
           &baseptr,
           &win
        );
    }
    else
    {
        returnCode = MPI_Win_allocate
        (
            // From num elements -> num of bytes
            std::streamsize(num_elements * disp_unit),
            disp_unit,
            MPI_INFO_NULL,
            communicator,
           &baseptr,
           &win
        );
    }

    if (FOAM_UNLIKELY((MPI_SUCCESS != returnCode) || (MPI_WIN_NULL == win)))
    {
        if (shared)
        {
            FatalError("MPI_Win_allocate_shared()")
                << Foam::abort(FatalError);
        }
        else
        {
            FatalError("MPI_Win_allocate()")
                << Foam::abort(FatalError);
        }

        return {nullptr, 0};
    }

    // Now have a window
    *self = UPstream::Window(win);

    // The address and the type-specific count
    return {baseptr, num_elements};
}


// ------------------------------------------------------------------------- //

std::pair<void*,int64_t>
Foam::UPstream::Window::mpi_win_allocate
(
    std::streamsize num_elements,
    int disp_unit,
    UPstream::Communicator communicator,
    const bool shared
)
{
    return call_window_allocate
    (
        this,
        PstreamUtils::Cast::to_mpi(communicator),

        num_elements,
        disp_unit,
        shared
    );
}


std::pair<void*,int64_t>
Foam::UPstream::Window::mpi_win_allocate
(
    std::streamsize num_elements,
    int disp_unit,
    int communicator,   // Index into MPICommunicators_
    const bool shared
)
{
    return call_window_allocate
    (
        this,
        PstreamGlobals::MPICommunicators_[communicator],

        num_elements,
        disp_unit,
        shared
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// NOTE: Currently no
//   - MPI_Win_create_dynamic()
//   - MPI_Win_attach()
//   - MPI_Win_detach()
// since working with their addresses (and broadcasting them)
// is fairly painful and probably not particularly efficient either

//
// Create a window to existing memory with MPI_Win_create().
//
static bool call_window_create
(
    Foam::UPstream::Window* self,
    MPI_Comm communicator,

    // [in] base address
    void *baseptr,
    // [in] number of elements (not bytes)
    std::streamsize num_elements,
    // [in] size of each element == sizeof(Type)
    const int disp_unit
)
{
    using namespace Foam;

    // No-op for non-parallel
    if (!UPstream::parRun())
    {
        *self = UPstream::Window(MPI_WIN_NULL);
        return false;
    }

    // if (FOAM_UNLIKELY(MPI_COMM_NULL == communicator))
    // {
    //     using namespace Foam;
    //     FatalErrorInFunction
    //         << "Attempt to use NULL communicator"
    //         << Foam::abort(FatalError);
    //         return false;
    // }

    MPI_Win win = PstreamUtils::Cast::to_mpi(*self);

    // Stringent handling of existing windows
    if (FOAM_UNLIKELY(MPI_WIN_NULL != win))
    {
        FatalErrorInFunction
            << "Window already exists. Use close() first"
            << Foam::abort(FatalError);
        return false;
    }

    // Leave nothing to chance
    if (!baseptr || !num_elements)
    {
        baseptr = nullptr;
        num_elements = 0;
    }

    int returnCode = MPI_Win_create
    (
        baseptr,
        // From num elements -> num of bytes
        std::streamsize(num_elements * disp_unit),
        disp_unit,
        MPI_INFO_NULL,
        communicator,
       &win
    );

    if (FOAM_UNLIKELY((MPI_SUCCESS != returnCode) || (MPI_WIN_NULL == win)))
    {
        FatalError("MPI_Win_create()")
            << Foam::abort(FatalError);
        return false;
    }

    // Now have a window
    *self = UPstream::Window(win);

    return (MPI_SUCCESS == returnCode);
}


bool Foam::UPstream::Window::mpi_win_create
(
    void *baseptr,
    std::streamsize num_elements,
    const int disp_unit,
    UPstream::Communicator communicator
)
{
    return call_window_create
    (
        this,
        PstreamUtils::Cast::to_mpi(communicator),

        baseptr,
        num_elements,
        disp_unit
    );
}


bool Foam::UPstream::Window::mpi_win_create
(
    void *baseptr,
    std::streamsize num_elements,
    const int disp_unit,
    int communicator    // Index into MPICommunicators_
)
{
    return call_window_create
    (
        this,
        PstreamGlobals::MPICommunicators_[communicator],

        baseptr,
        num_elements,
        disp_unit
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

void Foam::UPstream::Window::close()
{
    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);

    if (UPstream::parRun() && (MPI_WIN_NULL != win))
    {
        MPI_Win_free(&win);
        *this = UPstream::Window(MPI_WIN_NULL);
    }
}


// * * * * * * * * * * * * * * * Synchronization * * * * * * * * * * * * * * //

void Foam::UPstream::Window::mpi_win_flushing(int rank, bool local)
{
    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);

    if (UPstream::parRun() && (MPI_WIN_NULL != win))
    {
        if (rank < 0)
        {
            if (local) MPI_Win_flush_local_all(win);
            else /* */ MPI_Win_flush_all(win);
        }
        else
        {
            if (local) MPI_Win_flush_local(rank, win);
            else /* */ MPI_Win_flush(rank, win);
        }
    }
}


void Foam::UPstream::Window::sync()
{
    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);

    if (UPstream::parRun() && (MPI_WIN_NULL != win))
    {
        MPI_Win_sync(win);
    }
}


void Foam::UPstream::Window::mpi_win_locking(int rank, bool exclusive)
{
    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);

    if (UPstream::parRun() && (MPI_WIN_NULL != win))
    {
        if (rank < 0)
        {
            MPI_Win_lock_all
            (
                (exclusive ? MPI_MODE_NOCHECK : 0),
                win
            );
        }
        else
        {
            MPI_Win_lock
            (
                (exclusive ? MPI_LOCK_EXCLUSIVE : MPI_LOCK_SHARED),
                rank,
                0,  // No assertion
                win
            );
        }
    }
}


void Foam::UPstream::Window::mpi_win_unlocking(int rank)
{
    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);

    if (UPstream::parRun() && (MPI_WIN_NULL != win))
    {
        if (rank < 0)
        {
            MPI_Win_unlock_all(win);
        }
        else
        {
            MPI_Win_unlock(rank, win);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Window::get_data
(
    void* origin,                   // Type checking done by caller
    std::streamsize count,
    const UPstream::dataTypes dataTypeId,
    int target_rank,
    int target_disp
) const
{
    if (!UPstream::parRun() || !origin || !count)
    {
        // Nothing to do
        return true;
    }

    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);
    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);

    if (FOAM_UNLIKELY(MPI_WIN_NULL == win))
    {
        FatalError("MPI_Get()")
            << "Called with MPI_WIN_NULL."
            << Foam::abort(FatalError);
        return false;
    }

    int returnCode = MPI_Get
    (
        // origin
        origin, count, datatype,
        // target
        target_rank, target_disp, count, datatype,
        // window
        win
    );

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalError("MPI_Get()")
            << Foam::abort(FatalError);
        return false;
    }

    return (MPI_SUCCESS == returnCode);
}


bool Foam::UPstream::Window::put_data
(
    const void* origin,             // Type checking done by caller
    std::streamsize count,
    const UPstream::dataTypes dataTypeId,
    int target_rank,
    int target_disp
) const
{
    if (!UPstream::parRun() || !origin || !count)
    {
        // Nothing to do
        return true;
    }

    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);
    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);

    if (FOAM_UNLIKELY(MPI_WIN_NULL == win))
    {
        FatalError("MPI_Put()")
            << "Called with MPI_WIN_NULL."
            << Foam::abort(FatalError);
        return false;
    }

    int returnCode = MPI_Put
    (
        // origin
        origin, count, datatype,
        // target
        target_rank, target_disp, count, datatype,
        // window
        win
    );

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalError("MPI_Put()")
            << Foam::abort(FatalError);
        return false;
    }

    return (MPI_SUCCESS == returnCode);
}


bool Foam::UPstream::Window::put_data
(
    const UPstream::opCodes opCodeId,
    const void* origin,             // Type checking done by caller
    std::streamsize count,
    const UPstream::dataTypes dataTypeId,
    int target_rank,
    int target_disp
) const
{
    if (UPstream::opCodes::invalid == opCodeId)
    {
        // Regular data put - doesn't use/need an op-type!
        return this->put_data
        (
            origin,
            count,
            dataTypeId,
            target_rank,
            target_disp
        );
    }

    if (!UPstream::parRun() || !origin || !count)
    {
        // Nothing to do
        return true;
    }

    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);
    MPI_Op optype = PstreamGlobals::getOpCode(opCodeId);
    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);

    if (FOAM_UNLIKELY(MPI_WIN_NULL == win))
    {
        FatalError("MPI_Accumulate()")
            << "Called with MPI_WIN_NULL."
            << Foam::abort(FatalError);
        return false;
    }
    if (FOAM_UNLIKELY(MPI_OP_NULL == optype))
    {
        FatalError("MPI_Accumulate()")
            << "Invalid opcode:" << int(opCodeId)
            << " type:" << int(dataTypeId) << " count:" << label(count) << nl
            << Foam::abort(FatalError);
        return false;
    }

    int returnCode = MPI_Accumulate
    (
        // origin
        origin, count, datatype,
        // target
        target_rank, target_disp, count, datatype,
        // operation
        optype,
        // window
        win
    );

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalError("MPI_Accumulate()")
            << Foam::abort(FatalError);
        return false;
    }

    return (MPI_SUCCESS == returnCode);
}


bool Foam::UPstream::Window::mpi_fetch_and_op
(
    const UPstream::opCodes opCodeId,
    const void* origin,             // Type checking done by caller
    void* result,                   // Type checking done by caller
    const UPstream::dataTypes dataTypeId,
    int target_rank,
    int target_disp
) const
{
    if (!UPstream::parRun())
    {
        // Fails in non-parallel
        return false;
    }

    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);
    MPI_Op optype = PstreamGlobals::getOpCode(opCodeId);
    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);

    if (FOAM_UNLIKELY(MPI_WIN_NULL == win))
    {
        FatalError("MPI_Fetch_and_op()")
            << "Called with MPI_WIN_NULL."
            << Foam::abort(FatalError);
        return false;
    }
    if (FOAM_UNLIKELY(MPI_OP_NULL == optype))
    {
        FatalError("MPI_Fetch_and_op()")
            << "Invalid opcode:" << int(opCodeId)
            << " type:" << int(dataTypeId) << nl
            << Foam::abort(FatalError);
        return false;
    }

    int returnCode = MPI_Fetch_and_op
    (
        origin, result, datatype,
        // target
        target_rank, target_disp,
        // operation
        optype,
        // window
        win
    );

    // Error handling
    if (FOAM_UNLIKELY(returnCode != MPI_SUCCESS))
    {
        FatalError("MPI_Fetch_and_op()")
            << Foam::abort(FatalError);
        return false;
    }

    return (MPI_SUCCESS == returnCode);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Check for failure of MPI_Win_get_attr
#undef  CheckFail_Win_get_attr
#define CheckFail_Win_get_attr(returnCode, flag, attribute)       \
{                                                                 \
    if (FOAM_UNLIKELY((MPI_SUCCESS != returnCode) || !flag))      \
    {                                                             \
        FatalError("MPI_Win_get_attr()")                          \
            << "Failed getting attribute " << attribute << endl   \
            << Foam::abort(FatalError);                           \
    }                                                             \
}


bool Foam::UPstream::Window::is_shared(const bool failNonShared) const
{
    if (!UPstream::parRun())
    {
        // Nothing to do
        return false;
    }

    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);

    if (FOAM_UNLIKELY(MPI_WIN_NULL == win))
    {
        return false;
    }

    // Error handling flags
    int returnCode(MPI_ERR_UNKNOWN);
    int flag(1);
    int flavour(0);

    // MPI_WIN_CREATE_FLAVOR : Type (int *)
    {
        // const auto key = MPI_WIN_CREATE_FLAVOR;
        typedef int value_type;
        void* val(nullptr);

        returnCode = MPI_Win_get_attr(win, MPI_WIN_CREATE_FLAVOR, &val, &flag);
        CheckFail_Win_get_attr(returnCode, flag, "MPI_WIN_CREATE_FLAVOR");

        flavour = int
        (
            *static_cast<value_type*>(val)
        );
    }

    if (failNonShared && (MPI_WIN_FLAVOR_SHARED != flavour))
    {
        FatalErrorInFunction
            << "Expecting a shared window but had ("
            << flavour << ") flavour instead" << endl
            << Foam::abort(FatalError);
    }

    return (MPI_WIN_FLAVOR_SHARED == flavour);
}


std::pair<void*,int64_t>
Foam::UPstream::Window::mpi_win_query
(
    UPstream::Window window,
    const int expected_disp_unit
)
{
    if (!UPstream::parRun())
    {
        // Nothing to do
        return {nullptr, 0};
    }

    MPI_Win win = PstreamUtils::Cast::to_mpi(window);

    if (FOAM_UNLIKELY(MPI_WIN_NULL == win))
    {
        FatalError("MPI_Win_get_attr()")
            << "Called with MPI_WIN_NULL."
            << Foam::abort(FatalError);
        return {nullptr, 0};
    }


    // Error handling flags
    int returnCode(MPI_ERR_UNKNOWN);
    int flag(1);

    // Debugging
    // MPI_WIN_CREATE_FLAVOR : Type (int *)
    // if (FOAM_UNLIKELY(UPstream::debug & 2))
    // {
    //     // const auto key = MPI_WIN_CREATE_FLAVOR;
    //     typedef int value_type;
    //     void* val(nullptr);
    //
    //     returnCode =
    //         MPI_Win_get_attr(win, MPI_WIN_CREATE_FLAVOR, &val, &flag);
    //     CheckFail_Win_get_attr(returnCode, flag, "MPI_WIN_CREATE_FLAVOR");
    //
    //     int flavour = *static_cast<value_type*>(val);
    //     Perr<< "Window created with flavour (" << flavour << ')' << endl;
    // }

    std::pair<void*,int64_t> result(nullptr, 0);

    // The window size
    // MPI_WIN_SIZE : Type (MPI_Aint *)
    {
        // const auto key = MPI_WIN_SIZE;
        typedef MPI_Aint value_type;
        void* val(nullptr);

        returnCode = MPI_Win_get_attr(win, MPI_WIN_SIZE, &val, &flag);
        CheckFail_Win_get_attr(returnCode, flag, "MPI_WIN_SIZE");

        result.second = *static_cast<value_type*>(val);
    }

    // Early exit
    if (result.second == 0)
    {
        return {nullptr, 0};
    }

    // The base address
    // MPI_WIN_BASE : Type (void *)
    {
        // const auto key = MPI_WIN_BASE;
        void* value(nullptr);

        returnCode = MPI_Win_get_attr(win, MPI_WIN_BASE, &value, &flag);
        CheckFail_Win_get_attr(returnCode, flag, "MPI_WIN_BASE");

        result.first = value;
    }

    // Early exit - this probably can never happen
    // (ie, nullptr but non-zero size)
    if (result.first == nullptr)
    {
        return {nullptr, 0};
    }

    // Scale count by the expected displacement unit
    if (expected_disp_unit)
    {
        result.second /= expected_disp_unit;

        int disp_unit = 1;

        // The displacement units
        // MPI_WIN_DISP_UNIT : Type (int *)
        {
            // const auto key = MPI_WIN_DISP_UNIT;
            typedef int value_type;
            void* val(nullptr);

            returnCode = MPI_Win_get_attr(win, MPI_WIN_DISP_UNIT, &val, &flag);
            CheckFail_Win_get_attr(returnCode, flag, "MPI_WIN_DISP_UNIT");

            disp_unit = *static_cast<value_type*>(val);
        }

        // Error if the expected disp_unit is incorrect
        // - ignore this check if the window is empty

        if (expected_disp_unit != disp_unit)
        {
            FatalErrorInFunction
                << "Window [size=" << result.second
                << "] created with Type size=" << disp_unit
                << " but expecting Type size=" << expected_disp_unit << endl
                << Foam::abort(FatalError);
        }
    }

    return result;
}


std::pair<void*,int64_t>
Foam::UPstream::Window::mpi_win_query_shared
(
    UPstream::Window window,
    int target_rank,
    const int expected_disp_unit
)
{
    if (!UPstream::parRun())
    {
        // Nothing to do
        return {nullptr, 0};
    }

    MPI_Win win = PstreamUtils::Cast::to_mpi(window);

    if (FOAM_UNLIKELY(MPI_WIN_NULL == win))
    {
        FatalError("MPI_Win_shared_query()")
            << "Called with MPI_WIN_NULL."
            << Foam::abort(FatalError);
        return {nullptr, 0};
    }

    // Fail if window is not shared
    const bool shared = window.is_shared(true);

    if (!shared)
    {
        return {nullptr, 0};
    }

    // Initial values and fallback
    MPI_Aint num_bytes = 0;
    void *baseptr = nullptr;
    int disp_unit = 1;

    int returnCode = MPI_Win_shared_query
    (
        win,
        target_rank,
       &num_bytes,
       &disp_unit,
       &baseptr
    );

    if (FOAM_UNLIKELY(MPI_SUCCESS != returnCode))
    {
        FatalError("MPI_Win_shared_query()")
            << Foam::abort(FatalError);
        return {nullptr, 0};
    }

    std::pair<void*,int64_t> result(baseptr, num_bytes);

    // Scale count by the expected displacement unit
    // - probably Fatal not to supply this value
    //
    // Note that with share the baseptr will be non-null even if the
    // local window has zero bytes. This maintains the contiguous
    // addressing across all ranks

    if (result.second && expected_disp_unit)
    {
        result.second /= expected_disp_unit;

        // Error if the expected disp_unit is incorrect
        // - ignore this check if the window is empty

        if (expected_disp_unit != disp_unit)
        {
            FatalErrorInFunction
                << "Window on rank(" << target_rank
                << ") [size=" << result.second
                << "] created with Type size=" << disp_unit
                << " but expecting Type size=" << expected_disp_unit << endl
                << Foam::abort(FatalError);
        }
    }

    return result;
}


#undef CheckFail_Win_get_attr


// ************************************************************************* //
