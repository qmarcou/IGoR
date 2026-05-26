#pragma once

#include <igorCoreExport.h>

struct CORE_EXPORT BandCalibrationParams {
    int n_seed_seqs = 50;
    int n_verify_seqs = 100;
    int min_band_slack = 15;
    int band_expand_step = 10;
    int max_iterations = 3;
    int score_tolerance = 0;
    int seed = 42;
};

struct CORE_EXPORT AlignmentPreset {
    int band_lower_diag;
    int band_upper_diag;
    bool seq2_leading_free;
    bool seq1_leading_free;
    bool seq1_trailing_free;
    bool seq2_trailing_free;
    bool local_alignment;

    static AlignmentPreset v_gene_5p_race(int max_leader_nt = 80, int max_indel = 15)
    {
        return {max_leader_nt - max_indel, max_leader_nt + max_indel, false, true, true, false, false};
    }

    static AlignmentPreset v_gene_multiplex_pcr(int primer_pos_in_gene = 0, int max_indel = 15)
    {
        return {-(primer_pos_in_gene + max_indel), -(primer_pos_in_gene - max_indel), true, false, true, false,
                false};
    }

    static AlignmentPreset j_gene_c_primer(int max_cdr3_nt = 60, int max_indel = 5)
    {
        return {0, max_cdr3_nt + max_indel, false, true, true, false, false};
    }

    static AlignmentPreset j_gene_j_primer(int max_cdr3_nt = 60, int max_indel = 5)
    {
        return {0, max_cdr3_nt + max_indel, false, true, false, true, false};
    }

    static AlignmentPreset d_gene(int max_indel = 3)
    {
        return {-max_indel, max_indel, false, false, false, false, true};
    }

    static AlignmentPreset with_band(int half_width)
    {
        return {-half_width, half_width, false, false, false, false, false};
    }
};
