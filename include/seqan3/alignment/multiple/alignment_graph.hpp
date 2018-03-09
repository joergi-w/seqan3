// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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

template <std::ranges::RandomAccessRange sequence_container_type, typename weight_type = int>
requires std::ranges::RandomAccessRange<typename sequence_container_type::value_type>
      && std::ranges::SizedRange<typename sequence_container_type::value_type>
      && Alphabet<typename sequence_container_type::value_type::value_type>
class alignment_graph
{
protected:
    //                               seq_id    offset    length
    using fragment_type = std::tuple<unsigned, unsigned, unsigned>;

    //! \brief The k-partite graph that stores the connections between fragments.
    lemon::ListDigraph graph;

    //! \brief The edges of the graph contain weights.
    lemon::ListDigraph::ArcMap<weight_type> weights;

    //! \brief The underlying sequences.
    sequence_container_type sequences;


    std::vector<fragment_type> fragments;
    unsigned num_nodes;

public:
    /*!
     * \brief Constructor for an alignment graph.
     * \param seq A sequence container that the graph is linked to.
     */
    alignment_graph(sequence_container_type & seq) : weights(graph), sequences(seq), num_nodes(0u)
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
