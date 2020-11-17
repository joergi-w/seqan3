// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::format_clustal.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

//#include <iterator>
//#include <string>
//#include <vector>
//
//#include <seqan3/core/char_operations/predicate.hpp>
//#include <seqan3/core/concept/core_language.hpp>
//#include <seqan3/core/concept/tuple.hpp>
//#include <seqan3/core/detail/to_string.hpp>
//#include <seqan3/core/type_traits/range.hpp>
//#include <seqan3/core/type_traits/template_inspection.hpp>
//#include <seqan3/io/alignment_file/detail.hpp>
//#include <seqan3/io/alignment_file/format_sam_base.hpp>
#include <seqan3/io/alignment_file/header.hpp>
#include <seqan3/io/alignment_file/input_format_concept.hpp>
#include <seqan3/io/alignment_file/input_options.hpp>
//#include <seqan3/io/alignment_file/misc.hpp>
#include <seqan3/io/alignment_file/output_format_concept.hpp>
#include <seqan3/io/alignment_file/output_options.hpp>
//#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
//#include <seqan3/io/detail/ignore_output_iterator.hpp>
//#include <seqan3/io/detail/misc.hpp>
//#include <seqan3/io/stream/iterator.hpp>
//#include <seqan3/io/sequence_file/input_format_concept.hpp>
//#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/range/detail/misc.hpp>
//#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/istreambuf.hpp>
//#include <seqan3/range/views/slice.hpp>
//#include <seqan3/range/views/take_until.hpp>
//#include <seqan3/range/views/to_char.hpp>
//#include <seqan3/range/views/to.hpp>
//#include <seqan3/std/algorithm>
//#include <seqan3/std/concepts>
//#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief       The SAM format (tag).
 * \implements  AlignmentFileFormat
 * \ingroup     alignment_file
 *
 * \details
 *
 * ### Introduction
 *
 * SAM is often used for storing alignments of several read sequences against one
 * or more reference sequences. See the
 * [article on wikipedia](https://en.wikipedia.org/wiki/SAM_(file_format)) for an
 * introduction of the format or look into the official
 * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
 * **SeqAn implements version 1.6 of the SAM specification**.
 *
 * Take a look at our tutorial \ref tutorial_alignment_file for a walk through of how to read alignment files.
 *
 * ### fields_specialisation
 *
 * The SAM format provides the following fields:
 * seqan3::field::alignment, seqan3::field::seq, seqan3::field::qual,
 * seqan3::field::id, seqan3::field::ref_seq, seqan3::field::ref_id
 * seqan3::field::ref_ossfet, seqan3::field::offset, seqan3::field::flag,
 * seqan3::field::mapq and seqan3::field::mate.
 * In addition there is the seqan3::field::header_ptr, which is usually only used internally
 * to provide the range-based functionality of the file.
 *
 * **None of the fields are required** when writing but will be defaulted
 * to '0' for numeric fields and '*' for other fields.
 *
 * ### SAM format columns -> fields
 *
 * Since many users will be accustomed to the columns of the SAM format, here is a
 * mapping of the common SAM format columns to the SeqAn record fields:
 *
 * | #  | SAM Column ID |  FIELD name                                       |
 * |:--:|:--------------|:--------------------------------------------------|
 * | 1  | QNAME         | seqan3::field::id                                 |
 * | 2  | FLAG          | seqan3::field::flag                               |
 * | 3  | RNAME         | seqan3::field::ref_id                             |
 * | 4  | POS           | seqan3::field::ref_offset                         |
 * | 5  | MAPQ          | seqan3::field::mapq                               |
 * | 6  | CIGAR         | implicilty stored in seqan3::field::alignment     |
 * | 7  | RNEXT         | seqan3::field::mate (tuple pos 0)                 |
 * | 8  | PNEXT         | seqan3::field::mate (tuple pos 1)                 |
 * | 9  | TLEN          | seqan3::field::mate (tuple pos 2)                 |
 * | 10 | SEQ           | seqan3::field::seq                                |
 * | 11 | QUAL          | seqan3::field::qual                               |
 *
 * The (read sequence/query) **OFFSET** will be required to store the soft
 * clipping information at the read start (end clipping will be automatically
 * deduced by how much the read sequence length + offset is larger than the
 * alignment length).
 *
 * Note: SeqAn currently does not support hard clipping. When reading SAM,
 * hard-clipping is discarded; but the resulting alignment/sequence combination
 * is still valid.
 *
 * ### Format Check
 *
 * The format checks are implemented according to the official
 * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf)
 * in order to ensure correct SAM file output.
 *
 * If a non-recoverable format violation is encountered on reading, or you specify
 * invalid values/combinations when writing, seqan3::format_error is thrown.
 *
 * ### Header implementation
 *
 * The SAM header (if present) is read/written once in the beginning before the
 * first record is read/written.
 */
class format_clustal
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    // construction cannot be noexcept because this class has a std::string variable as a quality string buffer.
    format_clustal() = default; //!< Defaulted.
    format_clustal(format_clustal const &) = default; //!< Defaulted.
    format_clustal & operator=(format_clustal const &) = default; //!< Defaulted.
    format_clustal(format_clustal &&) = default; //!< Defaulted.
    format_clustal & operator=(format_clustal &&) = default; //!< Defaulted.
    ~format_clustal() = default; //!< Defaulted.
    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "aln" }
    };

protected:
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type,
              typename ref_seqs_type,
              typename ref_ids_type,
              typename seq_type,
              typename id_type,
              typename offset_type,
              typename ref_seq_type,
              typename ref_id_type,
              typename ref_offset_type,
              typename align_type,
              typename cigar_type,
              typename flag_type,
              typename mapq_type,
              typename qual_type,
              typename mate_type,
              typename tag_dict_type,
              typename e_value_type,
              typename bit_score_type>
    void read_alignment_record(stream_type & stream,
                               alignment_file_input_options<seq_legal_alph_type> const & SEQAN3_DOXYGEN_ONLY(options),
                               ref_seqs_type & SEQAN3_DOXYGEN_ONLY(ref_seqs),
                               alignment_file_header<ref_ids_type> & SEQAN3_DOXYGEN_ONLY(header),
                               seq_type & SEQAN3_DOXYGEN_ONLY(seq),
                               qual_type & SEQAN3_DOXYGEN_ONLY(qual),
                               id_type & id,
                               offset_type & SEQAN3_DOXYGEN_ONLY(offset),
                               ref_seq_type & SEQAN3_DOXYGEN_ONLY(ref_seq),
                               ref_id_type & SEQAN3_DOXYGEN_ONLY(ref_id),
                               ref_offset_type & SEQAN3_DOXYGEN_ONLY(ref_offset),
                               align_type & SEQAN3_DOXYGEN_ONLY(align),
                               cigar_type & SEQAN3_DOXYGEN_ONLY(cigar_vector),
                               flag_type & SEQAN3_DOXYGEN_ONLY(flag),
                               mapq_type & SEQAN3_DOXYGEN_ONLY(mapq),
                               mate_type & SEQAN3_DOXYGEN_ONLY(mate),
                               tag_dict_type & SEQAN3_DOXYGEN_ONLY(tag_dict),
                               e_value_type & SEQAN3_DOXYGEN_ONLY(e_value),
                               bit_score_type & SEQAN3_DOXYGEN_ONLY(bit_score))
    {
        auto stream_view = views::istreambuf(stream);
        auto it = stream_view.begin();
        auto e = stream_view.end();

        // skip initial whitespace
        while (it != e && is_space(*it))
            ++it;

        // check if file starts with "CLUSTAL"
        for (char chr : std::string{"CLUSTAL"})
        {
            if (*it != chr)
                throw parse_error{std::string{"Expected to read '"} + chr + std::string{"', but found "} +
                                  detail::make_printable(*it) + std::string{" in the CLUSTAL header"}};

            if (++it == e)
                throw unexpected_end_of_input{"CLUSTAL header does not end in newline."};
        }

        // skip until space
        while (it != e && (!is_char<'\n'>)(*it))
            ++it;

        // skip until id starts
        while (it != e && is_space(*it))
            ++it;

        // read id
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            for (; it != e; ++it)
            {
                id.push_back(assign_char_to(*it, std::ranges::range_value_t<id_type>{}));
            }
        }
    }

    template <typename stream_type,
              typename header_type,
              typename seq_type,
              typename id_type,
              typename ref_seq_type,
              typename ref_id_type,
              typename align_type,
              typename qual_type,
              typename mate_type,
              typename tag_dict_type,
              typename e_value_type,
              typename bit_score_type>
    void write_alignment_record(stream_type & stream,
                                alignment_file_output_options const & options,
                                header_type && header,
                                seq_type && seq,
                                qual_type && qual,
                                id_type && id,
                                int32_t const offset,
                                ref_seq_type && SEQAN3_DOXYGEN_ONLY(ref_seq),
                                ref_id_type && ref_id,
                                std::optional<int32_t> ref_offset,
                                align_type && align,
                                std::vector<cigar> const & cigar_vector,
                                sam_flag const flag,
                                uint8_t const mapq,
                                mate_type && mate,
                                tag_dict_type && tag_dict,
                                e_value_type && SEQAN3_DOXYGEN_ONLY(e_value),
                                bit_score_type && SEQAN3_DOXYGEN_ONLY(bit_score))
    {}

} // namespace seqan3
