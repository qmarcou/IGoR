#include <catch2/catch_test_macros.hpp>

#include <igor/Core/SeqAn2Aligner.h>

#ifdef IGOR_WITH_SEQAN2
static Matrix<double> simple_matrix()
{
    std::vector<double> v(15 * 15, -1.0);
    for (int i = 0; i < 15; ++i) v[i + 15 * i] = 2.0;
    return Matrix<double>(15, 15, v);
}

TEST_CASE("SeqAn2Aligner basic global alignment", "[seqan2]")
{
    SeqAn2Aligner aligner(simple_matrix(), 2, V_gene, AlignmentPreset::with_band(20));
    aligner.set_genomic_sequences({{"g1", "ACGTACGT"}});
    auto alns = aligner.align_seq("ACGTACGT", 1.0, true, true);
    REQUIRE(!alns.empty());
    auto aln = alns.front();
    REQUIRE(aln.gene_name == "g1");
    REQUIRE(aln.five_p_offset == 0);
    REQUIRE(aln.three_p_offset == 7);
    REQUIRE(alignment_data_to_cigar(aln) == "8=");
}

TEST_CASE("SeqAn2Aligner allow_in_dels filters gaps", "[seqan2]")
{
    SeqAn2Aligner aligner(simple_matrix(), 2, V_gene, AlignmentPreset::with_band(20));
    aligner.set_genomic_sequences({{"g1", "ACGTACGT"}});
    auto with_gaps = aligner.align_seq("ACGTTACGT", -100.0, true, false);
    auto no_gaps = aligner.align_seq("ACGTTACGT", -100.0, false, false);
    REQUIRE(!with_gaps.empty());
    REQUIRE(no_gaps.empty());
}

TEST_CASE("known diagonal calibration updates band", "[seqan2][calibration]")
{
    auto preset = AlignmentPreset::with_band(1);
    BandCalibrationParams params;
    params.min_band_slack = 5;
    auto calibrated = calibrate_preset(preset, std::vector<int>{60, 70, 80}, params);
    REQUIRE(calibrated.band_lower_diag <= 55);
    REQUIRE(calibrated.band_upper_diag >= 85);
}
#endif
