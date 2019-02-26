// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

//! [exercise]
#include <iostream>   // std::cerr, std::endl
#include <seqan3/alphabet/all.hpp>

// Definition of our alphabet.
enum class dna2
{
    S,   // strong (GC)
    W    // weak (AT)
};

// Overload alphabet_size.
template <>
struct seqan3::alphabet_size<dna2> : std::integral_constant<size_t, 2> {};

// Overload underlying_rank.
template <>
struct seqan3::underlying_rank<dna2>
{
    using type = uint8_t;
};

// Define to_rank.
seqan3::underlying_rank<dna2>::type to_rank(dna2 const alph) noexcept
{
    return (uint8_t)(alph == dna2::W ? 1u : 0u);
}

// Define assign_rank.
dna2 & assign_rank (dna2 & alph, seqan3::underlying_rank<dna2>::type const rank) noexcept
{
    return alph = (rank == 1 ? dna2::W : dna2::S);
}

dna2 assign_rank (dna2 && alph, seqan3::underlying_rank<dna2>::type const rank) noexcept
{
    alph = (rank == 1 ? dna2::W : dna2::S);
    return std::move(alph);
}

// Alphabet requirements.

// Overload underlying_char.
template <>
struct seqan3::underlying_char<dna2>
{
    using type = char;
};

// Define to_char.
seqan3::underlying_char<dna2>::type to_char(dna2 const alph) noexcept
{
    return alph == dna2::W ? 'W' : 'S';
}

// Define assign_char.
dna2 & assign_char (dna2 & alph, seqan3::underlying_char<dna2>::type const ch) noexcept
{
    alph = (ch == 'W' ? dna2::W : dna2::S);
    return alph;
}

dna2 assign_char (dna2 && alph, seqan3::underlying_char<dna2>::type const ch) noexcept
{
    alph = (ch == 'W' ? dna2::W : dna2::S);
    return std::move(alph);
}

// Define char_is_valid_for.
template <typename T>
    requires std::Same<seqan3::remove_cvref_t<T>, dna2>
bool char_is_valid_for(typename seqan3::underlying_char<T>::type const ch) noexcept
{
    return to_char(assign_char(dna2{}, ch)) == ch;
}

// Define assign_char_strict.
dna2 & assign_char_strict (dna2 & alph, seqan3::underlying_char<dna2>::type const ch)
{
    return alph = (ch == 'W' ? dna2::W : dna2::S);
}

dna2 assign_char_strict (dna2 && alph, seqan3::underlying_char<dna2>::type const ch)
{
    alph = (ch == 'W' ? dna2::W : dna2::S);
    return std::move(alph);
}




// Concept checks.
static_assert(std::Regular<dna2>);
static_assert(std::StrictTotallyOrdered<dna2>);
static_assert(seqan3::Semialphabet<dna2>);
//static_assert(seqan3::Alphabet<dna2>);

// Constrained function that works only for seqan3::Alphabet types.
template <typename alph_type>
    requires seqan3::Semialphabet<alph_type>
void test_function(alph_type)
{
    std::cerr << "You're good!" << std::endl;
    std::cerr << "The alphabet size is " << (seqan3::alphabet_size<dna2>::value) << "." << std::endl;
}

int main ()
{
    // Let's test our new alphabet class here. The compilation fails, if members are missing.
    test_function(dna2{});
    return 0;
}
//! [exercise]
