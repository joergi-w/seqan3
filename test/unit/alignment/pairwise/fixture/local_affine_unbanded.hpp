// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "alignment_fixture.hpp"

namespace seqan3::test::alignment::fixture::local::affine::unbanded
{

inline constexpr auto align_config = align_cfg::mode{local_alignment} |
                                     align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}};

static auto dna4_01 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        // score: 11 (4 matches, 1 mismatch)
        // alignment:
        // GTTTA
        // || ||
        // GTCTA
        "AACCGGTTTAACCGGTT"_dna4,
        "ACGTCTACGTA"_dna4,
        align_config | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        11,
        "GTTTA",
        "GTCTA",
        alignment_coordinate{column_index_type{5u}, row_index_type{2u}},
        alignment_coordinate{column_index_type{10u}, row_index_type{7u}}
    };
}();

static auto dna4_02 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "ACGTCTACGTA"_dna4,
        "AACCGGTTTAACCGGTT"_dna4,
        align_config | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        11,
        "GTCTA",
        "GTTTA",
        alignment_coordinate{column_index_type{2u}, row_index_type{5u}},
        alignment_coordinate{column_index_type{7u}, row_index_type{10u}}
    };
}();

} // namespace seqan3::test::alignment::fixture::local::affine::unbanded
