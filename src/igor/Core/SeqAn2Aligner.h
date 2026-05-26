#pragma once

#include <igor/Core/Aligner.h>
#include <igor/Core/AlignmentPreset.h>
#include <igor/Core/Utils.h>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef IGOR_WITH_SEQAN2
#include <seqan/score.h>
#endif

struct SeqAn2AlignConfig {
    int band_lower_diag = -20;
    int band_upper_diag = 20;
    bool seq2_leading_free = false;
    bool seq1_leading_free = false;
    bool seq1_trailing_free = true;
    bool seq2_trailing_free = true;
    bool local_alignment = false;
};

#ifdef IGOR_WITH_SEQAN2
class CORE_EXPORT SeqAn2Aligner {
public:
    SeqAn2Aligner(Matrix<double> subst, int gap_penalty, Gene_class gene, SeqAn2AlignConfig config = {});
    SeqAn2Aligner(Matrix<double> subst, int gap_penalty, Gene_class gene, AlignmentPreset preset,
                  std::vector<std::string> calibration_reads = {}, BandCalibrationParams calib_params = {});
    ~SeqAn2Aligner();

    std::forward_list<Alignment_data> align_seq(const std::string& seq, double threshold,
                                                bool allow_in_dels, bool best_only = false) const;
    std::unordered_map<int, std::forward_list<Alignment_data>>
    align_seqs(const std::vector<std::pair<int, std::string>>& seqs, double threshold, bool allow_in_dels) const;
    void set_genomic_sequences(std::vector<std::pair<std::string, std::string>>);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

CORE_EXPORT AlignmentPreset calibrate_preset(AlignmentPreset preset, const std::vector<std::string>& reads,
                                             const std::vector<std::string>& germlines,
                                             const seqan2::Score<int, seqan2::Simple>& scoring,
                                             BandCalibrationParams params = {});
CORE_EXPORT AlignmentPreset calibrate_preset(AlignmentPreset preset, const std::vector<int>& known_diagonals,
                                             BandCalibrationParams params = {});
#endif
