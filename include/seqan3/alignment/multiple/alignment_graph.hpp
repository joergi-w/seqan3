// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \brief Contains the alignment graph for (multiple) sequence alignments.
 */

#pragma once

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/std/ranges>

#if SEQAN3_WITH_LEMON

#include <lemon/list_graph.h>

namespace seqan3::detail
{

template <Arithmetic weight_type = int>
class alignment_graph
{
private:
    //! \brief The k-partite graph that stores the connections between fragments.
    lemon::ListDigraph graph;

    //! \brief The edges of the graph contain weights.
    lemon::ListDigraph::ArcMap<weight_type> weights;

    //! \brief The number of nodes in the graph.
    size_t num_nodes;

    //! \brief The number of sequences in the graph.
    size_t num_seqs;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */

    //! Default constructor.
    constexpr alignment_graph() :
        graph{},
        weights{graph},
        num_nodes{0},
        num_seqs{0}
    {}

    /*!
     * \brief Constructor for an alignment graph.
     * \param seq A sequence container that the graph is linked to.
     */
    template <std::ranges::ForwardRange seqs_t>
        requires std::ranges::SizedRange<typename seqs_t::value_type>
    constexpr alignment_graph(seqs_t const & sequences) : alignment_graph()
    {
        // count the number of characters in all the sequences
        for (auto const & sequence : sequences)
        {
            num_nodes += std::ranges::size(sequence);
            ++num_seqs;
        }

        // create fragments and nodes
        graph.reserveNode(num_nodes);
        for (size_t cnt = 0ul; cnt < num_nodes; ++cnt)
            graph.addNode();
    }

    //! Copy constructor.
    constexpr alignment_graph(alignment_graph const & origin) :
        graph{},
        weights{graph},
        num_nodes{origin.num_nodes},
        num_seqs{origin.num_seqs}
    {
        lemon::digraphCopy(origin.graph, graph).arcMap(origin.weights, weights).run();
    }

    //! Copy assignment.
    constexpr alignment_graph & operator=(alignment_graph const & origin)
    {
        return alignment_graph{origin};
    }

    //! Default move constructor.
    constexpr alignment_graph(alignment_graph &&) noexcept = default;

    //! Default move assignment.
    constexpr alignment_graph & operator=(alignment_graph &&) noexcept = default;

    //! Default destructor.
    ~alignment_graph() noexcept = default;
    //!\}

    /*!
     * \brief Return the number of nodes in the graph (constant time access).
     * \return number of nodes
     */
    constexpr size_t size() const noexcept
    {
        return num_nodes;
    }

    /*!
     * \brief Return the number of sequences in the graph (constant time access).
     * \return number of sequences
     */
    constexpr size_t num_sequences() const noexcept
    {
        return num_seqs;
    }
};

} // namespace seqan3::detail

#else
    #error "This module is only available if the LEMON graph library is installed. See README for more information."
#endif // SEQAN3_WITH_LEMON
