// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Tests for the alignment graph.
 */

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/align/alignment_graph.hpp>

using namespace seqan3;
using namespace literal;

template <typename T>
class align_graph : public ::testing::Test
{};

// ------------------------------------------------------------------
// basics
// ------------------------------------------------------------------

std::vector<rna5_vector> sequences{"ACG"_rna5, "UCG"_rna5, "CCG"_rna5};

TEST(align_graph, constructors)
{
    using TRandomAccessContainer = std::vector<rna5_vector>;
    EXPECT_FALSE((std::is_default_constructible_v<alignment_graph<TRandomAccessContainer>>));
    EXPECT_TRUE((std::is_copy_constructible_v<alignment_graph<TRandomAccessContainer>>));
    EXPECT_TRUE((std::is_move_constructible_v<alignment_graph<TRandomAccessContainer>>));
    EXPECT_TRUE((std::is_copy_assignable_v<alignment_graph<TRandomAccessContainer>>));
    EXPECT_TRUE((std::is_move_assignable_v<alignment_graph<TRandomAccessContainer>>));
}

TEST(align_graph, input_container)
{
    // std::vector<std::vector<rna5>>
    alignment_graph graph_vv{sequences};

    // std::array of std::vector
    std::array<dna15_vector, 3> seq_array{"ACG"_dna15, "TCG"_dna15, "CCG"_dna15};
    alignment_graph graph_av{seq_array};

    // std::array of std::array
    std::array<std::array<aa27, 10>, 10> seq_array2{};
    alignment_graph graph_aa{seq_array2};

    // std::deque of std::deque
    std::deque<std::deque<wuss51>> seq_deque{};
    alignment_graph graph_dd{seq_deque};
}

TEST(align_graph, size)
{
    alignment_graph graph{sequences};
    EXPECT_EQ(graph.size(), 9u);
}
