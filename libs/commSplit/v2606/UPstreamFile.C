/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Mark Olesen
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

#include "fileName.H"
#include "UPstreamFile.H"
#include "PstreamGlobals.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Has _c() version?
#undef Foam_UPstream_largeCountFile

#if (MPI_VERSION >= 4)
#define Foam_UPstream_largeCountFile
#endif

// Macros for calling versions with or without '_c'
#ifdef Foam_UPstream_largeCountFile
#define Foam_mpiCall(Function)  Function##_c
#else
#define Foam_mpiCall(Function)  Function
#endif


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

inline bool checkCount(std::streamsize count, const char* what)
{
    #ifndef Foam_UPstream_largeCountFile
    if (FOAM_UNLIKELY(count > std::streamsize(INT_MAX)))
    {
        using namespace Foam;
        FatalErrorInFunction
            << "Write size " << label(count)
            << " exceeds INT_MAX bytes for '" << what << "'\n"
            << Foam::abort(Foam::FatalError);
        return false;
    }
    #endif

    return true;
}

}  // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class UPstream::File::Impl Declaration
\*---------------------------------------------------------------------------*/

class UPstream::File::Impl
{
    //- The file-handle
    MPI_File fh_;

    //- Path of the open file
    fileName name_;

    //- The current state (open|read|write|closed etc)
    int state_;

    //- The associated rank when openned
    int rank_;

public:

    //- The file states
    enum states : int { CLOSED = 0, READ, WRITE, ATOMIC_WRITE };

    // Constructors

        //- Default construct
        Impl()
        :
            fh_(MPI_FILE_NULL),
            state_(CLOSED),
            rank_(0)
        {}


    // Member Functions

        // The file handle
        const MPI_File& handle() const noexcept { return fh_; }
        MPI_File& handle() noexcept { return fh_; }

        // Path to the open file
        const fileName& name() const noexcept { return name_; }
        fileName& name() noexcept { return name_; }

        // Change the file state, return the old value
        int state(states val) noexcept
        {
            int old(state_);
            state_ = val;
            return old;
        }

        //- Is rank 0 ? (master rank)
        bool master() const noexcept { return (rank_ == 0); }

        //- Get the associated rank
        int rank() const noexcept { return rank_; }

        //- Set the associated rank
        void rank(int val) noexcept { rank_ = val; }


    // Checks

        // The file state
        bool is_open() const noexcept { return state_; }

        // The file read state
        bool is_read() const noexcept
        {
            return (states::READ == state_);
        }

        // The file write atomic state
        bool is_atomic() const noexcept
        {
            return (states::ATOMIC_WRITE == state_);
        }

        // The file write state (atomic or non-atomic)
        bool is_write() const noexcept
        {
            return (states::ATOMIC_WRITE == state_ || states::WRITE == state_);
        }

        //- Assert is_read() or FatalError
        inline bool checkReadable(const char* what) const
        {
            if (FOAM_UNLIKELY(!is_read()))
            {
                FatalErrorInFunction
                    << "File handler not open for reading '" << what << "'\n"
                    << "name: " << name() << nl
                    << Foam::exit(Foam::FatalError);
                return false;
            }
            return true;
        }

        //- Assert is_write() or FatalError
        inline bool checkWritable(const char* what) const
        {
            if (FOAM_UNLIKELY(!is_write()))
            {
                FatalErrorInFunction
                    << "File handler not open for writing'" << what << "'\n"
                    << "name: " << name() << nl
                    << Foam::exit(Foam::FatalError);
                return false;
            }
            return true;
        }
};

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::UPstream::File::supported()
{
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::File::File()
:
    file_(new UPstream::File::Impl)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UPstream::File::~File()
{
    if (FOAM_UNLIKELY(file_ && file_->is_open()))
    {
        WarningInFunction
            << "Exited scope without close()" << nl
            << "    FIX YOUR CODE!!" << endl;
        // Do not call close() since we don't know where that collective
        // should have been called
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fileName& Foam::UPstream::File::name() const
{
    return (file_ ? file_->name() : fileName::null);
}


bool Foam::UPstream::File::is_open() const
{
    return bool(file_ && file_->is_open());
}


bool Foam::UPstream::File::close()
{
    if (FOAM_UNLIKELY(!file_->is_open()))
    {
        WarningInFunction
            << "Called without an open file handler !" << endl;
        return false;
    }

    MPI_File_close(&(file_->handle()));

    // Atomic rename of file (master only)
    const fileName& pathname = file_->name();

    if (file_->master() && file_->is_atomic() && !pathname.empty())
    {
        std::rename
        (
            (pathname + "~tmp~").c_str(),
            pathname.c_str()
        );
    }

    file_->state(Impl::CLOSED);
    file_->name() = "";
    file_->rank(0);

    return true;
}


// * * * * * * * * * * * * Member Functions (Reading)  * * * * * * * * * * * //

#if 0
bool Foam::UPstream::File::open_read
(
    const int communicator,
    const fileName& pathname
)
{
    //Needed?  PstreamGlobals::checkCommunicator(communicator, 0);

    if (FOAM_UNLIKELY(file_->is_open()))
    {
        WarningInFunction
            << "Previous use of file handler did not call close()" << nl
            << "    FIX YOUR CODE!!" << endl;
        // Do not call close() since we don't know where that collective
        // should have been called
    }
    file_->state(Impl::CLOSED);
    file_->name() = pathname;  // <- set now for external error messages
    file_->rank(0);

    int returnCode = MPI_File_open
    (
        PstreamGlobals::MPICommunicators_[communicator],
        pathname.c_str(),
        (MPI_MODE_RDONLY),
        MPI_INFO_NULL,
        &(file_->handle())
    );

    if (FOAM_UNLIKELY(MPI_SUCCESS != returnCode))
    {
        FatalErrorInFunction
            << "Error encounted in MPI_File_open() : "
            << pathname << nl
            << Foam::exit(Foam::FatalError);

        return false;
    }

    file_->state(Impl::READ);
    file_->name() = pathname;
    file_->rank(UPstream::myProcNo(communicator));

    return true;  // ie, is_read()
}
#endif


// * * * * * * * * * * * *  Non-Collective Reading * * * * * * * * * * * * * //

#if 0
bool Foam::UPstream::File::get_header(DynamicList<char>& content)
{
    std::streamsize headerSize(4096);

    // constexpr const char* const func = "MPI_File_read_at";
    file_->checkReadable("MPI_File_read_at");

    if (off_t fileLen = Foam::fileSize(this->name()); fileLen >= 0)
    {
        std::streamsize size = std::streamsize(fileLen);
        if (headerSize > size)
        {
            headerSize = size;
        }
    }
    else
    {
        content.clear();
        return false;
    }

    // Get the first header content:
    content.resize_nocopy(headerSize);

    int returnCode = Foam_mpiCall(MPI_File_read_at)
    (
        file_->handle(),
        0,  // offset
        content.data(),
        content.size(),
        MPI_BYTE,
        MPI_STATUS_IGNORE
    );

    // Wrap as ISpanStream headerStream(content);
    if (MPI_SUCCESS == returnCode)
    {
        ISpanStream is(content);
        dictionary headerDict;

        // Read the regular "FoamFile" header
        bool ok = io.readHeader(headerDict, is);

        // Probably collated - extract class from "data.class"
        if
        (
            decomposedBlockData::isCollatedType(io)
         && headerDict.readIfPresent("data.class", io.headerClassName())
        )
        {
            return ok;
        }
    }

    return (MPI_SUCCESS == returnCode);
}
#endif


// * * * * * * * * * * * * Member Functions (Writing)  * * * * * * * * * * * //

bool Foam::UPstream::File::open_write
(
    const int communicator,
    const fileName& pathname,
    IOstreamOption::atomicType atomicType
)
{
    //Needed?  PstreamGlobals::checkCommunicator(communicator, 0);

    if (FOAM_UNLIKELY(file_->is_open()))
    {
        WarningInFunction
            << "Previous use of file handler did not call close()" << nl
            << "    FIX YOUR CODE!!" << endl;
        // Do not call close() since we don't know where that collective
        // should have been called
    }
    file_->state(Impl::CLOSED);
    file_->name() = pathname;  // <- set now for external error messages
    file_->rank(0);

    const bool atomic = (IOstreamOption::atomicType::ATOMIC == atomicType);

    // When opening new files, remove file variants out of the way.
    // Eg, opening "file1"
    // - remove old "file1.gz" (compressed)
    // - also remove old "file1" if it is a symlink

    const fileName pathname_gz(pathname + ".gz");
    const fileName pathname_tmp(pathname + "~tmp~");

    // File to open with MPI_File_open
    const auto& target = (atomic ? pathname_tmp : pathname);

    // Remove old compressed version (if any)
    if
    (
        auto fType = Foam::type(pathname_gz, false);
        (fType == fileName::SYMLINK || fType == fileName::FILE)
    )
    {
        Foam::rm(pathname_gz);
    }

    // Avoid writing into symlinked files (non-append mode)
    if
    (
        auto fType = Foam::type(target, false);
        fType == fileName::SYMLINK
    )
    {
        Foam::rm(target);
    }

    int returnCode = MPI_File_open
    (
        PstreamGlobals::MPICommunicators_[communicator],
        target.c_str(),
        (MPI_MODE_CREATE | MPI_MODE_WRONLY),
        MPI_INFO_NULL,
        &(file_->handle())
    );

    if (FOAM_UNLIKELY(MPI_SUCCESS != returnCode))
    {
        FatalErrorInFunction
            << "Error encounted in MPI_File_open() : "
            << target << nl
            << Foam::exit(Foam::FatalError);

        return false;
    }

    file_->state(atomic ? Impl::ATOMIC_WRITE : Impl::WRITE);
    file_->name() = pathname;
    file_->rank(UPstream::myProcNo(communicator));

    return true;  // ie, is_write()
}


bool Foam::UPstream::File::set_size(std::streamsize num_bytes)
{
    if (FOAM_UNLIKELY(!file_->is_open()))
    {
        WarningInFunction
            << "Called without an open file handler !" << endl;
        return false;
    }

    int returnCode = MPI_File_set_size(file_->handle(), num_bytes);

    return (MPI_SUCCESS == returnCode);
}


// * * * * * * * * * * * *  Non-Collective Writing * * * * * * * * * * * * * //

bool Foam::UPstream::File::write_data
(
    const void* data,
    std::streamsize count,
    const UPstream::dataTypes dataTypeId
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    // constexpr const char* const func = "MPI_File_write";
    file_->checkWritable("MPI_File_write");
    checkCount(count, "MPI_File_write");

    int returnCode = Foam_mpiCall(MPI_File_write)
    (
        file_->handle(),
        data,
        count,
        datatype,
        MPI_STATUS_IGNORE
    );

    return (MPI_SUCCESS == returnCode);
}


bool Foam::UPstream::File::write_data_at
(
    std::streamsize offset,
    const void* data,
    std::streamsize count,
    const UPstream::dataTypes dataTypeId
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    // constexpr const char* const func = "MPI_File_write_at";
    file_->checkWritable("MPI_File_write_at");
    checkCount(count, "MPI_File_write_at");

    int returnCode = Foam_mpiCall(MPI_File_write_at)
    (
        file_->handle(),
        offset,
        data,
        count,
        datatype,
        MPI_STATUS_IGNORE
    );

    return (MPI_SUCCESS == returnCode);
}


// * * * * * * * * * * * * *  Collective Writing * * * * * * * * * * * * * * //

bool Foam::UPstream::File::write_data_all
(
    const void* data,
    std::streamsize count,
    const UPstream::dataTypes dataTypeId
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    // constexpr const char* const func = "MPI_File_write_all";
    file_->checkWritable("MPI_File_write_all");
    checkCount(count, "MPI_File_write_all");

    int returnCode = Foam_mpiCall(MPI_File_write_all)
    (
        file_->handle(),
        data,
        count,
        datatype,
        MPI_STATUS_IGNORE
    );

    return (MPI_SUCCESS == returnCode);
}


bool Foam::UPstream::File::write_data_at_all
(
    std::streamsize offset,
    const void* data,
    std::streamsize count,
    const UPstream::dataTypes dataTypeId
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    // constexpr const char* const func = "MPI_File_write_at_all";
    file_->checkWritable("MPI_File_write_at_all");
    checkCount(count, "MPI_File_write_at_all");

    int returnCode = Foam_mpiCall(MPI_File_write_at_all)
    (
        file_->handle(),
        offset,
        data,
        count,
        datatype,
        MPI_STATUS_IGNORE
    );

    return (MPI_SUCCESS == returnCode);
}


// bool Foam::UPstream::File::write_data_all_begin
// (
//     const void* data,
//     std::streamsize count,
//     const UPstream::dataTypes dataTypeId
// )
// {
//     MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);
//
//     // constexpr const char* const func = "MPI_File_write_all_begin";
//     file_->checkWritable("MPI_File_write_all_begin");
//     checkCount(count, "MPI_File_write_all_begin");
//
//     int returnCode = Foam_mpiCall(MPI_File_write_all_begin)
//     (
//         file_->handle(),
//         data,
//         count,
//         datatype,
//         MPI_STATUS_IGNORE
//     );
//
//     return (MPI_SUCCESS == returnCode);
// }


// bool Foam::UPstream::File::write_data_all_end
// (
//     const void* data
// )
// {
//     file_->checkWritable("MPI_File_write_all_end");
//     int returnCode = Foam_mpiCall(MPI_File_write_all_end)
//     (
//         file_->handle(),
//         data
//         MPI_STATUS_IGNORE
//     );
//
//     return (MPI_SUCCESS == returnCode);
// }


// bool Foam::UPstream::File::write_data_at_all_begin
// (
//     std::streamsize offset,
//     const void* data,
//     std::streamsize count,
//     const UPstream::dataTypes dataTypeId
// )
// {
//     MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);
//
//     // constexpr const char* const func = "MPI_File_write_at_all_begin";
//     file_->checkWritable("MPI_File_write_at_all_begin");
//     checkCount(count, "MPI_File_write_at_all_begin");
//
//     int returnCode = Foam_mpiCall(MPI_File_write_at_all_begin)
//     (
//         file_->handle(),
//         offset,
//         data,
//         count,
//         datatype,
//         MPI_STATUS_IGNORE
//     );
//
//     return (MPI_SUCCESS == returnCode);
// }


// bool Foam::UPstream::File::write_data_at_all_end
// (
//     const void* data
// )
// {
//     file_->checkWritable("MPI_File_write_at_all_end");
//
//     int returnCode = Foam_mpiCall(MPI_File_write_at_all_end)
//     (
//         file_->handle(),
//         data
//         MPI_STATUS_IGNORE
//     );
//
//     return (MPI_SUCCESS == returnCode);
// }


// ************************************************************************* //
