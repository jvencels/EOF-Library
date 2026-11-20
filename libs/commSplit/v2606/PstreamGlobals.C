/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2023-2025 OpenCFD Ltd.
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

MPI_Comm Foam::PstreamGlobals::MPI_COMM_FOAM = MPI_COMM_WORLD;

Foam::DynamicList<bool> Foam::PstreamGlobals::pendingMPIFree_;
Foam::DynamicList<MPI_Comm> Foam::PstreamGlobals::MPICommunicators_;
Foam::DynamicList<MPI_Request> Foam::PstreamGlobals::outstandingRequests_;

Foam::PstreamGlobals::DataTypeCountLookupTable
Foam::PstreamGlobals::dataTypesCount_(1);

Foam::PstreamGlobals::DataTypeLookupTable
Foam::PstreamGlobals::MPIdataTypes_(MPI_DATATYPE_NULL);

Foam::PstreamGlobals::OpCodesLookupTable
Foam::PstreamGlobals::MPIopCodes_(MPI_OP_NULL);


// * * * * * * * * * * * * * * * Communicators * * * * * * * * * * * * * * * //

void Foam::PstreamGlobals::initCommunicator(const label index)
{
    if (FOAM_UNLIKELY(index < 0 || index > MPICommunicators_.size()))
    {
        FatalErrorInFunction
            << "PstreamGlobals out of sync with UPstream data. Problem."
            << Foam::abort(FatalError);
    }
    else if (index == MPICommunicators_.size())
    {
        // Extend storage with null values
        pendingMPIFree_.emplace_back(false);
        MPICommunicators_.emplace_back(MPI_COMM_NULL);
    }
    else
    {
        // Init with null values
        pendingMPIFree_[index] = false;
        MPICommunicators_[index] = MPI_COMM_NULL;
    }
}


// * * * * * * * * * * * * * * * * Data Types  * * * * * * * * * * * * * * * //

void Foam::PstreamGlobals::initDataTypes()
{
    static_assert
    (
        PstreamGlobals::DataTypeCountLookupTable::max_size()
     == (int(UPstream::dataTypes::DataTypes_end)+1),
        "Data count lookup table size != number of dataTypes enumerations"
    );
    static_assert
    (
        PstreamGlobals::DataTypeLookupTable::max_size()
     == (int(UPstream::dataTypes::DataTypes_end)+1),
        "Lookup table size != number of dataTypes enumerations"
    );

    // From enumeration to MPI datatype for fundamental types
    // (count is always 1)
    #undef  defineType
    #define defineType(Idx, BaseType)                                         \
    {                                                                         \
        dataTypesCount_[int(UPstream::dataTypes::Idx)] = 1;                   \
        MPIdataTypes_[int(UPstream::dataTypes::Idx)] = BaseType;              \
    }

    // Fundamental Types [10]:
    defineType(type_byte,   MPI_BYTE);
    defineType(type_int16,  MPI_INT16_T);
    defineType(type_int32,  MPI_INT32_T);
    defineType(type_int64,  MPI_INT64_T);
    defineType(type_uint16, MPI_UINT16_T);
    defineType(type_uint32, MPI_UINT32_T);
    defineType(type_uint64, MPI_UINT64_T);
    defineType(type_float,  MPI_FLOAT);
    defineType(type_double, MPI_DOUBLE);
    defineType(type_long_double, MPI_LONG_DOUBLE);

    #undef defineType

    // User-define types
    #undef  defineUserType
    #define defineUserType(Idx, Count, BaseType, Name)                        \
    {                                                                         \
        dataTypesCount_[int(UPstream::dataTypes::Idx)] = Count;               \
        auto& dt = MPIdataTypes_[int(UPstream::dataTypes::Idx)];              \
        MPI_Type_contiguous(Count, BaseType, &dt);                            \
        MPI_Type_set_name(dt, Name);                                          \
        MPI_Type_commit(&dt);                                                 \
    }

    // User Types [6]:
    defineUserType(type_3float,  3, MPI_FLOAT,  "float[3]");
    defineUserType(type_3double, 3, MPI_DOUBLE, "double[3]");
    defineUserType(type_6float,  6, MPI_FLOAT,  "float[6]");
    defineUserType(type_6double, 6, MPI_DOUBLE, "double[6]");
    defineUserType(type_9float,  9, MPI_FLOAT,  "float[9]");
    defineUserType(type_9double, 9, MPI_DOUBLE, "double[9]");

    #undef defineUserType
}


void Foam::PstreamGlobals::deinitDataTypes()
{
    // User types only
    auto first =
    (
        MPIdataTypes_.begin() + int(UPstream::dataTypes::User_begin)
    );
    const auto last =
    (
        MPIdataTypes_.begin() + int(UPstream::dataTypes::User_end)
    );

    for (; first != last; ++first)
    {
        if (MPI_DATATYPE_NULL != *first)
        {
            MPI_Type_free(&(*first));
        }
    }
}


// Debugging
bool Foam::PstreamGlobals::checkDataTypes()
{
    // Check all types, not just user types
    auto first =
    (
        MPIdataTypes_.begin()
    );
    const auto last =
    (
        MPIdataTypes_.begin() + int(UPstream::dataTypes::DataTypes_end)
    );

    for (; (first != last); ++first)
    {
        if (MPI_DATATYPE_NULL == *first)
        {
            return false;
        }
    }

    return true;
}


// Debugging
void Foam::PstreamGlobals::printDataTypes(bool all)
{
    int rank = -1;
    if
    (
        (MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &rank))
     || (rank != 0)
    )
    {
        return;
    }

    const auto print = [&](auto firstIndex, auto lastIndex)
    {
        auto first =
        (
            PstreamGlobals::MPIdataTypes_.begin() + int(firstIndex)
        );
        const auto last =
        (
            PstreamGlobals::MPIdataTypes_.begin() + int(lastIndex)
        );

        for (; (first != last); ++first)
        {
            std::cerr
                << "  name = "
                << PstreamGlobals::dataType_name(*first) << '\n';
        }
    };

    if (all)
    {
        std::cerr << "enumerated data types:\n";
        print
        (
            UPstream::dataTypes::Basic_begin,
            UPstream::dataTypes::Basic_end
        );
    }
    else
    {
        // User types only.
        std::cerr << "enumerated user-defined data types:\n";
    }
    print
    (
        UPstream::dataTypes::User_begin,
        UPstream::dataTypes::User_end
    );
}


std::string Foam::PstreamGlobals::dataType_name(MPI_Datatype datatype)
{
    if (MPI_DATATYPE_NULL == datatype)
    {
        return std::string("(null)");
    }

    char buf[MPI_MAX_OBJECT_NAME];
    int len;

    if (MPI_SUCCESS == MPI_Type_get_name(datatype, buf, &len))
    {
        if (len > 0)
        {
            return std::string(buf, len);
        }
        else
        {
            return std::string("(anon)");
        }
    }

    return std::string("???");
}


// * * * * * * * * * * * * * * * * Op Codes  * * * * * * * * * * * * * * * * //

void Foam::PstreamGlobals::initOpCodes()
{
    static_assert
    (
        PstreamGlobals::OpCodesLookupTable::max_size()
     == (int(UPstream::opCodes::OpCodes_end)+1),
        "Lookup table size != number of opCodes enumerations"
    );

    // From enumeration to MPI datatype
    #undef  defineCode
    #define defineCode(Idx, CodeType) \
    MPIopCodes_[int(UPstream::opCodes::Idx)] = CodeType;

    defineCode(op_min,  MPI_MIN);
    defineCode(op_max,  MPI_MAX);
    defineCode(op_sum,  MPI_SUM);
    defineCode(op_prod, MPI_PROD);

    // TBD: still need to sort out if they are MPI_C_BOOL or MPI_CXX_BOOL
    // ...
    defineCode(op_bool_and, MPI_LAND);
    defineCode(op_bool_or,  MPI_LOR);
    defineCode(op_bool_xor, MPI_LXOR);

    defineCode(op_bit_and,  MPI_BAND);
    defineCode(op_bit_or,   MPI_BOR);
    defineCode(op_bit_xor,  MPI_BXOR);

    // Do not include MPI_MINLOC, MPI_MAXLOC since they are tied to
    // float_int, double_int and larger or other types

    // window-only
    defineCode(op_replace, MPI_REPLACE);
    defineCode(op_no_op, MPI_NO_OP);

    #undef defineCode
}


void Foam::PstreamGlobals::deinitOpCodes()
{}


bool Foam::PstreamGlobals::checkOpCodes()
{
    auto first = MPIopCodes_.begin();
    const auto last =
    (
        MPIopCodes_.begin() + int(UPstream::opCodes::OpCodes_end)
    );

    for (; (first != last); ++first)
    {
        if (MPI_OP_NULL == *first)
        {
            return false;
        }
    }

    return true;
}


// ************************************************************************* //
