#include <catch2/catch_test_macros.hpp>

#include <igor/Core/AlignmentPreset.h>

TEST_CASE("AlignmentPreset band factories", "[preset]")
{
    auto v_race = AlignmentPreset::v_gene_5p_race(80, 15);
    REQUIRE(v_race.band_lower_diag == 65);
    REQUIRE(v_race.band_upper_diag == 95);

    auto v_pcr = AlignmentPreset::v_gene_multiplex_pcr(50, 15);
    REQUIRE(v_pcr.band_lower_diag == -65);
    REQUIRE(v_pcr.band_upper_diag == -35);

    auto j_c = AlignmentPreset::j_gene_c_primer(60, 5);
    REQUIRE(j_c.band_lower_diag == 0);
    REQUIRE(j_c.band_upper_diag == 65);

    auto j_j = AlignmentPreset::j_gene_j_primer(60, 5);
    REQUIRE(j_j.band_lower_diag == 0);
    REQUIRE(j_j.band_upper_diag == 65);

    auto exact = AlignmentPreset::v_gene_5p_race(80, 0);
    REQUIRE(exact.band_lower_diag == 80);
    REQUIRE(exact.band_upper_diag == 80);

    auto zero_primer = AlignmentPreset::v_gene_multiplex_pcr(0, 10);
    REQUIRE(zero_primer.band_lower_diag == -10);
    REQUIRE(zero_primer.band_upper_diag == 10);
}

TEST_CASE("AlignmentPreset free-end booleans", "[preset]")
{
    auto v_race = AlignmentPreset::v_gene_5p_race();
    REQUIRE_FALSE(v_race.seq2_leading_free);
    REQUIRE(v_race.seq1_leading_free);
    REQUIRE(v_race.seq1_trailing_free);
    REQUIRE_FALSE(v_race.seq2_trailing_free);
    REQUIRE_FALSE(v_race.local_alignment);

    auto v_pcr = AlignmentPreset::v_gene_multiplex_pcr();
    REQUIRE(v_pcr.seq2_leading_free);
    REQUIRE_FALSE(v_pcr.seq1_leading_free);
    REQUIRE(v_pcr.seq1_trailing_free);
    REQUIRE_FALSE(v_pcr.seq2_trailing_free);
    REQUIRE_FALSE(v_pcr.local_alignment);

    auto j_c = AlignmentPreset::j_gene_c_primer();
    REQUIRE_FALSE(j_c.seq2_leading_free);
    REQUIRE(j_c.seq1_leading_free);
    REQUIRE(j_c.seq1_trailing_free);
    REQUIRE_FALSE(j_c.seq2_trailing_free);
    REQUIRE_FALSE(j_c.local_alignment);

    auto j_j = AlignmentPreset::j_gene_j_primer();
    REQUIRE_FALSE(j_j.seq2_leading_free);
    REQUIRE(j_j.seq1_leading_free);
    REQUIRE_FALSE(j_j.seq1_trailing_free);
    REQUIRE(j_j.seq2_trailing_free);
    REQUIRE_FALSE(j_j.local_alignment);

    auto d = AlignmentPreset::d_gene(3);
    REQUIRE(d.band_lower_diag == -3);
    REQUIRE(d.band_upper_diag == 3);
    REQUIRE(d.local_alignment);

    auto band = AlignmentPreset::with_band(20);
    REQUIRE(band.band_lower_diag == -20);
    REQUIRE(band.band_upper_diag == 20);
    REQUIRE_FALSE(band.seq2_leading_free);
    REQUIRE_FALSE(band.seq1_leading_free);
    REQUIRE_FALSE(band.seq1_trailing_free);
    REQUIRE_FALSE(band.seq2_trailing_free);
    REQUIRE_FALSE(band.local_alignment);
}
