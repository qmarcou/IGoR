#include <igor/Core/SeqAn2Aligner.h>

#ifdef IGOR_WITH_SEQAN2
#include <seqan/align.h>
#include <seqan/sequence.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

using TSeq = seqan2::String<seqan2::Dna5>;
using TAlign = seqan2::Align<TSeq, seqan2::ArrayGaps>;
using TScore = seqan2::Score<int, seqan2::Simple>;

static SeqAn2AlignConfig to_config(const AlignmentPreset& p)
{
    return {p.band_lower_diag, p.band_upper_diag, p.seq2_leading_free, p.seq1_leading_free,
            p.seq1_trailing_free, p.seq2_trailing_free, p.local_alignment};
}

static TSeq to_dna5(const std::string& s)
{
    TSeq out;
    seqan2::resize(out, s.size());
    for (size_t i = 0; i < s.size(); ++i) out[i] = seqan2::Dna5(s[i]);
    return out;
}

static bool has_indels(const TAlign& ali)
{
    auto const& r0 = seqan2::row(ali, 0);
    auto const& r1 = seqan2::row(ali, 1);
    auto begin = std::min(seqan2::clippedBeginPosition(r0), seqan2::clippedBeginPosition(r1));
    auto end = std::max(seqan2::clippedEndPosition(r0), seqan2::clippedEndPosition(r1));
    for (decltype(end) col = begin; col < end; ++col) {
        if (seqan2::isGap(r0, col) || seqan2::isGap(r1, col)) return true;
    }
    return false;
}

static Alignment_data from_seqan2_align(const TAlign& ali, const std::string& gene_name, double score)
{
    auto const& read_row = seqan2::row(ali, 0);
    auto const& germ_row = seqan2::row(ali, 1);
    int t = static_cast<int>(seqan2::clippedBeginPosition(read_row));
    int g = static_cast<int>(seqan2::clippedBeginPosition(germ_row));
    int seq_start = t + 1;
    int ref_start = g + 1;
    int last_t = t;
    int last_g = g;
    std::string cigar;
    char last_op = 0;
    int run = 0;
    auto emit = [&](char op) {
        if (run && op == last_op) { ++run; return; }
        if (run) { cigar += std::to_string(run); cigar += last_op; }
        last_op = op; run = 1;
    };

    auto begin = std::min(seqan2::clippedBeginPosition(read_row), seqan2::clippedBeginPosition(germ_row));
    auto end = std::max(seqan2::clippedEndPosition(read_row), seqan2::clippedEndPosition(germ_row));
    for (decltype(end) col = begin; col < end; ++col) {
        if (seqan2::isGap(read_row, col)) { emit('D'); last_g = g++; }
        else if (seqan2::isGap(germ_row, col)) { emit('I'); last_t = t++; }
        else { emit(read_row[col] == germ_row[col] ? '=' : 'X'); last_t = t++; last_g = g++; }
    }
    if (run) { cigar += std::to_string(run); cigar += last_op; }
    return alignment_data_from_cigar(gene_name, cigar, seq_start, last_t + 1, ref_start, last_g + 1, score);
}

template <bool A, bool B, bool C, bool D>
static int run_global(TAlign& ali, const TScore& score, int lower, int upper)
{
    return seqan2::globalAlignment(ali, score, seqan2::AlignConfig<A, B, C, D>(), lower, upper);
}

static int global_dispatch(TAlign& ali, const TScore& score, const SeqAn2AlignConfig& cfg)
{
#define CASE(a,b,c,d) if (cfg.seq2_leading_free==a && cfg.seq1_leading_free==b && cfg.seq1_trailing_free==c && cfg.seq2_trailing_free==d) return run_global<a,b,c,d>(ali, score, cfg.band_lower_diag, cfg.band_upper_diag)
    CASE(false,false,false,false); CASE(false,false,false,true); CASE(false,false,true,false); CASE(false,false,true,true);
    CASE(false,true,false,false); CASE(false,true,false,true); CASE(false,true,true,false); CASE(false,true,true,true);
    CASE(true,false,false,false); CASE(true,false,false,true); CASE(true,false,true,false); CASE(true,false,true,true);
    CASE(true,true,false,false); CASE(true,true,false,true); CASE(true,true,true,false); CASE(true,true,true,true);
#undef CASE
    return 0;
}

struct SeqAn2Aligner::Impl {
    TScore scoring;
    Gene_class gene;
    SeqAn2AlignConfig config;
    std::vector<std::pair<std::string, TSeq>> germlines;
    std::vector<std::string> germline_strings;

    Impl(Matrix<double> subst, int gap_penalty, Gene_class g, SeqAn2AlignConfig c)
        : scoring(static_cast<int>(subst(0,0)), -std::abs(gap_penalty), -std::abs(gap_penalty)), gene(g), config(c)
    {
        // SeqAn Simple scoring cannot hold IUPAC matrix scores; use match/mismatch from matrix as a robust baseline.
        int mismatch = static_cast<int>(subst(0,1));
        scoring = TScore(static_cast<int>(subst(0,0)), mismatch, -std::abs(gap_penalty));
    }
};

SeqAn2Aligner::SeqAn2Aligner(Matrix<double> subst, int gap_penalty, Gene_class gene, SeqAn2AlignConfig config)
    : impl_(new Impl(subst, gap_penalty, gene, config)) {}

SeqAn2Aligner::SeqAn2Aligner(Matrix<double> subst, int gap_penalty, Gene_class gene, AlignmentPreset preset,
                             std::vector<std::string> calibration_reads, BandCalibrationParams calib_params)
    : impl_(new Impl(subst, gap_penalty, gene, to_config(preset)))
{
    if (!calibration_reads.empty()) {
        AlignmentPreset calibrated = calibrate_preset(preset, calibration_reads, impl_->germline_strings, impl_->scoring, calib_params);
        impl_->config = to_config(calibrated);
    }
}

SeqAn2Aligner::~SeqAn2Aligner() = default;

void SeqAn2Aligner::set_genomic_sequences(std::vector<std::pair<std::string, std::string>> gs)
{
    impl_->germlines.clear();
    impl_->germline_strings.clear();
    for (auto& g : gs) {
        impl_->germlines.push_back({g.first, to_dna5(g.second)});
        impl_->germline_strings.push_back(g.second);
    }
}

std::forward_list<Alignment_data> SeqAn2Aligner::align_seq(const std::string& seq, double threshold,
                                                           bool allow_in_dels, bool best_only) const
{
    std::forward_list<Alignment_data> out;
    TSeq read = to_dna5(seq);
    double best = -1e300;
    std::forward_list<Alignment_data> all;
    for (const auto& germ : impl_->germlines) {
        if (impl_->config.local_alignment) {
            TAlign ali; seqan2::resize(seqan2::rows(ali), 2);
            seqan2::assignSource(seqan2::row(ali, 0), read);
            seqan2::assignSource(seqan2::row(ali, 1), germ.second);
            int sc = seqan2::localAlignment(ali, impl_->scoring);
            if (sc >= threshold && (allow_in_dels || !has_indels(ali))) all.push_front(from_seqan2_align(ali, germ.first, sc));
        } else {
            TAlign ali; seqan2::resize(seqan2::rows(ali), 2);
            seqan2::assignSource(seqan2::row(ali, 0), read);
            seqan2::assignSource(seqan2::row(ali, 1), germ.second);
            int sc = global_dispatch(ali, impl_->scoring, impl_->config);
            if (sc >= threshold && (allow_in_dels || !has_indels(ali))) all.push_front(from_seqan2_align(ali, germ.first, sc));
        }
    }
    if (!best_only) return all;
    for (const auto& a : all) best = std::max(best, a.score);
    for (const auto& a : all) if (a.score == best) out.push_front(a);
    return out;
}

std::unordered_map<int, std::forward_list<Alignment_data>>
SeqAn2Aligner::align_seqs(const std::vector<std::pair<int, std::string>>& seqs, double threshold, bool allow_in_dels) const
{
    std::unordered_map<int, std::forward_list<Alignment_data>> out;
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < seqs.size(); ++i) {
        auto alns = align_seq(seqs[i].second, threshold, allow_in_dels);
#pragma omp critical
        out[seqs[i].first] = alns;
    }
    return out;
}

static int percentile(std::vector<int> values, double p)
{
    if (values.empty()) return 0;
    std::sort(values.begin(), values.end());
    size_t idx = static_cast<size_t>(std::round((values.size() - 1) * p));
    return values[idx];
}

AlignmentPreset calibrate_preset(AlignmentPreset preset, const std::vector<int>& known_diagonals, BandCalibrationParams params)
{
    if (known_diagonals.empty()) return preset;
    int med = percentile(known_diagonals, 0.5);
    int p5 = percentile(known_diagonals, 0.05);
    int p95 = percentile(known_diagonals, 0.95);
    preset.band_lower_diag = p5 - params.min_band_slack;
    preset.band_upper_diag = p95 + params.min_band_slack;
    (void)med;
    return preset;
}

AlignmentPreset calibrate_preset(AlignmentPreset preset, const std::vector<std::string>& reads,
                                 const std::vector<std::string>& germlines, const TScore& scoring,
                                 BandCalibrationParams params)
{
    std::vector<int> diagonals;
    int n = std::min<int>(params.n_seed_seqs, reads.size());
    for (int i = 0; i < n; ++i) {
        TSeq read = to_dna5(reads[i]);
        for (const auto& gstr : germlines) {
            TSeq germ = to_dna5(gstr);
            TAlign ali; seqan2::resize(seqan2::rows(ali), 2);
            seqan2::assignSource(seqan2::row(ali, 0), read);
            seqan2::assignSource(seqan2::row(ali, 1), germ);
            if (preset.local_alignment) seqan2::localAlignment(ali, scoring);
            else seqan2::globalAlignment(ali, scoring);
            diagonals.push_back(static_cast<int>(seqan2::clippedBeginPosition(seqan2::row(ali, 0))) -
                                static_cast<int>(seqan2::clippedBeginPosition(seqan2::row(ali, 1))));
        }
    }
    return calibrate_preset(preset, diagonals, params);
}
#endif
