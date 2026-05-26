#include <catch2/catch_test_macros.hpp>

#include <igor/Core/Aligner.h>

#include <algorithm>
#include <cstdio>
#include <forward_list>
#include <string>
#include <vector>

static std::vector<int> sorted_list(std::forward_list<int> xs)
{
    std::vector<int> out(xs.begin(), xs.end());
    std::sort(out.begin(), out.end());
    return out;
}

TEST_CASE("parse extended CIGAR", "[cigar]")
{
    auto ops = parse_cigar("3=1X2I4D5M");
    REQUIRE(ops.size() == 5);
    REQUIRE(ops[0] == std::make_pair(3, '='));
    REQUIRE(ops[4] == std::make_pair(5, 'M'));
    REQUIRE_THROWS_AS(parse_cigar("0="), std::invalid_argument);
    REQUIRE_THROWS_AS(parse_cigar("=10"), std::invalid_argument);
    REQUIRE_THROWS_AS(parse_cigar("10Z"), std::invalid_argument);
}

TEST_CASE("alignment data to CIGAR basic cases", "[cigar]")
{
    REQUIRE(alignment_data_to_cigar(Alignment_data("g", 0, 0, 9, 10, {}, {}, {}, 10)) == "10=");
    REQUIRE(alignment_data_to_cigar(Alignment_data("g", 0, 0, 9, 10, {}, {}, {0,1,2,3,4,5,6,7,8,9}, 10)) == "10X");
    REQUIRE(alignment_data_to_cigar(Alignment_data("g", 0, 0, 12, 13, {7,6,5}, {}, {}, 10)) == "5=3I5=");
    REQUIRE(alignment_data_to_cigar(Alignment_data("g", 0, 0, 9, 12, {}, {5,4}, {}, 10)) == "4=2D6=");
}

TEST_CASE("alignment data CIGAR mixed round trip", "[cigar]")
{
    Alignment_data aln("gene", -2, 0, 9, 12, {6}, {3, 4}, {2, 8}, 42.0);
    std::string cigar = alignment_data_to_cigar(aln);
    REQUIRE(cigar == "1=2D1=1X3=1I1=1X1=");

    Alignment_data round = alignment_data_from_cigar(
        aln.gene_name, cigar,
        alignment_data_sequence_start(aln), alignment_data_sequence_end(aln),
        alignment_data_germline_start(aln), alignment_data_germline_end(aln), aln.score);
    REQUIRE(round.gene_name == aln.gene_name);
    REQUIRE(round.offset == aln.offset);
    REQUIRE(round.five_p_offset == aln.five_p_offset);
    REQUIRE(round.three_p_offset == aln.three_p_offset);
    REQUIRE(round.align_length == aln.align_length);
    REQUIRE(sorted_list(round.insertions) == sorted_list(aln.insertions));
    REQUIRE(sorted_list(round.deletions) == sorted_list(aln.deletions));
    std::sort(round.mismatches.begin(), round.mismatches.end());
    auto expected_mismatches = aln.mismatches;
    std::sort(expected_mismatches.begin(), expected_mismatches.end());
    REQUIRE(round.mismatches == expected_mismatches);
}

TEST_CASE("CIGAR M fallback loses mismatch detail by design", "[cigar]")
{
    Alignment_data aln = alignment_data_from_cigar("gene", "5M", 1, 5, 1, 5, 10);
    REQUIRE(aln.align_length == 5);
    REQUIRE(aln.mismatches.empty());
    REQUIRE(aln.five_p_offset == 0);
    REQUIRE(aln.three_p_offset == 4);
}

TEST_CASE("CSV reader preserves alignment length and offsets", "[cigar]")
{
    std::string path = "alignment_cigar_tmp.csv";
    Aligner aligner;
    std::unordered_map<int, std::forward_list<Alignment_data>> alignments;
    Alignment_data aln("gene", 1, 1, 10, 11, {6}, {}, {4}, 12.0);
    alignments[7].push_front(aln);
    aligner.write_alignments_seq_csv(path, alignments);

    auto loaded = aligner.read_alignments_seq_csv(path, 0.0, true);
    std::remove(path.c_str());
    REQUIRE(loaded.count(7) == 1);
    REQUIRE(!loaded[7].empty());
    const Alignment_data& got = loaded[7].front();
    REQUIRE(got.align_length == aln.align_length);
    REQUIRE(got.five_p_offset == aln.five_p_offset);
    REQUIRE(got.three_p_offset == aln.three_p_offset);
    REQUIRE(alignment_data_to_cigar(got) == alignment_data_to_cigar(aln));
}
