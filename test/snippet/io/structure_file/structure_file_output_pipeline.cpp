#include <sstream>

#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/io/structure_file/output.hpp>
#include <seqan3/range/views/take.hpp>

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
> example
UUGGAGUACACAACCUGUACACUCUUUC
..(((((..(((...)))..)))))... (-3.71))";

int main()
{
#if !SEQAN3_WORKAROUND_GCC_93983
    seqan3::structure_file_input my_in{std::istringstream{input}, seqan3::format_vienna{}};
    my_in | seqan3::views::take(2) | seqan3::structure_file_output{std::ostringstream{}, seqan3::format_vienna{}};
#endif // !SEQAN3_WORKAROUND_GCC_93983
}
