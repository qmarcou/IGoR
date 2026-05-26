#include <catch2/catch_test_macros.hpp>

#include <igor/Core/SeqAn2Aligner.h>

#ifdef IGOR_WITH_SEQAN2
TEST_CASE("calibration is reproducible from known diagonals", "[calibration]")
{
    BandCalibrationParams params;
    params.min_band_slack = 15;
    params.seed = 123;
    std::vector<int> diagonals{60, 70, 80, 90, 100};
    auto a = calibrate_preset(AlignmentPreset::v_gene_5p_race(), diagonals, params);
    auto b = calibrate_preset(AlignmentPreset::v_gene_5p_race(), diagonals, params);
    REQUIRE(a.band_lower_diag == b.band_lower_diag);
    REQUIRE(a.band_upper_diag == b.band_upper_diag);
    REQUIRE(a.band_lower_diag <= 45);
    REQUIRE(a.band_upper_diag >= 115);
}

TEST_CASE("empty known diagonal calibration keeps preset", "[calibration]")
{
    auto preset = AlignmentPreset::with_band(20);
    auto calibrated = calibrate_preset(preset, std::vector<int>{});
    REQUIRE(calibrated.band_lower_diag == -20);
    REQUIRE(calibrated.band_upper_diag == 20);
}
#endif
