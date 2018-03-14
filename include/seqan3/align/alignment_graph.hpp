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
 * \brief Contains the alignment graph for (multiple) sequence alignments.
 */

#pragma once

#include <seqan3/core/platform.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/alphabet/concept.hpp>

#if SEQAN3_WITH_LEMON

#include <lemon/list_graph.h>

namespace seqan3
{

template <typename TRandomAccessContainer>
requires random_access_range_concept<TRandomAccessContainer>
      && random_access_range_concept<typename TRandomAccessContainer::value_type>
      && sized_range_concept<typename TRandomAccessContainer::value_type>
      && alphabet_concept<typename TRandomAccessContainer::value_type::value_type>
class alignment_graph
{
protected:
    using fragment_type = std::tuple<unsigned, unsigned, unsigned>;
    //                               seq_id    offset    length

    lemon::ListDigraph graph;
    lemon::ListDigraph::ArcMap<double> weights;
    TRandomAccessContainer & sequences;
    std::vector<fragment_type> fragments;
    unsigned num_nodes;

public:
    /*!
     * \brief Constructor for an alignment graph.
     * \param seq A sequence container that the graph is linked to.
     */
    alignment_graph(TRandomAccessContainer & seq) : weights(graph), sequences(seq), num_nodes(0u)
    {
        // count the number of characters in all the sequences
        for (auto sequence : sequences)
            num_nodes += sequence.size();

        // create fragments and nodes
        fragments.resize(num_nodes);
        graph.reserveNode(num_nodes);
        for (unsigned cnt = 0u; cnt < num_nodes; ++cnt)
            graph.addNode();
    }

    //! delete default constructor
    alignment_graph() = delete;

    //! copy constructor
    alignment_graph(alignment_graph const & origin) : sequences(origin.sequences), fragments(origin.fragments)
    {
        lemon::digraphCopy(origin.graph, graph).arcMap(origin.weights, weights).run();
    }

    //! copy assignment
    alignment_graph & operator=(alignment_graph const & origin)
    {
        return alignment_graph(origin);
    }

    //! default move constructor / assigment
    alignment_graph(alignment_graph &&) = default;
    alignment_graph & operator=(alignment_graph &&) = default;

    /*!
     * \brief Return the number of nodes in the graph (constant time access).
     * \return number of nodes
     */
    unsigned size()
    {
        return num_nodes;
    }
};

} // namespace seqan3

#else
    #error "This module is only available if the LEMON library is installed. See README for more information."
#endif // SEQAN3_WITH_LEMON
