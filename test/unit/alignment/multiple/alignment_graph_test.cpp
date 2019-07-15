// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Tests for the alignment graph.
 */

#include <gtest/gtest.h>

#include <array>
#include <deque>
#include <type_traits>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alignment/multiple/alignment_graph.hpp>

using namespace seqan3;

template <typename T>
class align_graph : public ::testing::Test
{};

TEST(align_graph, constructors)
{
    EXPECT_TRUE((std::is_default_constructible_v<detail::alignment_graph<>>));
    EXPECT_TRUE((std::is_copy_constructible_v<detail::alignment_graph<>>));
    EXPECT_TRUE((std::is_move_constructible_v<detail::alignment_graph<>>));
    EXPECT_TRUE((std::is_copy_assignable_v<detail::alignment_graph<>>));
    EXPECT_TRUE((std::is_move_assignable_v<detail::alignment_graph<>>));
    EXPECT_TRUE((std::is_destructible_v<detail::alignment_graph<>>));
}

TEST(align_graph, default_constructor)
{
    detail::alignment_graph ag{};
    EXPECT_EQ(ag.size(), 0ul);
}

TEST(align_graph, input_container)
{
    // std::vector<std::vector<rna5>>
    std::vector<rna5_vector> seq_vec{"ACG"_rna5, "UCG"_rna5, "CCG"_rna5};
    detail::alignment_graph graph_vv{seq_vec};

    // std::array of std::vector
    std::array<rna5_vector, 3> seq_array{"ACG"_rna5, "UCG"_rna5, "CCG"_rna5};
    detail::alignment_graph graph_av{seq_array};

    // std::array of std::array
    std::array<std::array<rna5, 10>, 10> seq_array2{};
    detail::alignment_graph graph_aa{seq_array2};

    // std::deque of std::deque
    std::deque<std::deque<rna5>> seq_deque{};
    detail::alignment_graph graph_dd{seq_deque};
}

TEST(align_graph, size)
{
    std::vector<rna5_vector> sequences{"ACG"_rna5, "UCG"_rna5, "CCG"_rna5};
    detail::alignment_graph graph{sequences};
    EXPECT_EQ(graph.size(), 9ul);
    EXPECT_EQ(graph.num_sequences(), 3ul);
}
