// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/io/alignment_file/format_clustal.hpp>
#include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/io/alignment_file/output.hpp>

TEST(clustal_format, read_aln)
{
    std::vector<std::string> ids
    {
        {"M83762.1-1031_1093"},
        {"AC008670.6-83725_83795"},
        {"Z82044.1-16031_16103"},
        {"AE004843.1-4972_4900"},
        {"AB042432.1-14140_14072"}
    };

    std::string input
    {
        "CLUSTAL FORMAT\n"
        "\n"
        "M83762.1-1031_1093      gcuuuaaaagc-uuu---gcugaagcaacggcc----uuguaagucguag\n"
        "AC008670.6-83725_83795  acuuuuaaagg-aua-acagccauccguugguc----uuaggccccaaaa\n"
        "Z82044.1-16031_16103    gcgguuguggcgaag-ugguuaacgcaccagauuguggcucuggcacuc-\n"
        "AE004843.1-4972_4900    gcucauguagc-ucaguugguagagcacacccu----ugguaagggugag\n"
        "AB042432.1-14140_14072  guuucuguagu-ugaau---uacaacgaugauu----uuucaugucauug\n"
        "                                 *               *                        \n"
        "\n"
        "M83762.1-1031_1093      aa-aacu--a-ua---cguuuuaaagcu\n"
        "AC008670.6-83725_83795  au-uuuggugcaacuccaaauaaaagua\n"
        "Z82044.1-16031_16103    ---guggguucgauucccaucaaucgcc\n"
        "AE004843.1-4972_4900    gucagcgguucaaauccgcucaugagcu\n"
        "AB042432.1-14140_14072  gu-cgcaguugaaugcuguguagaaaua\n"
        "                                    *               "
    };

    std::stringstream istream{input};

    seqan3::alignment_file_input fin{istream, seqan3::format_clustal{}, seqan3::fields<seqan3::field::id>{}};
//    fin.options = options;

    auto it = fin.begin();
//    std::cerr << seqan3::get<seqan3::field::id>(*it) << std::endl;
//    EXPECT_EQ(seqan3::get<seqan3::field::id>(*it), ids[0]);
//    for (unsigned i = 0u; i < 1u; ++i, ++it)
//    {
//        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::seq>(*it), seqs[i])));
//        EXPECT_EQ(seqan3::get<seqan3::field::id>(*it), ids[i]);
//    }
}
