#include <igor/Core/Aligner.h>
#include <igor/Core/AlignmentPreset.h>
#ifdef IGOR_WITH_SEQAN2
#include <igor/Core/SeqAn2Aligner.h>
#endif

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

static Matrix<double> simple_matrix()
{
    std::vector<double> v(15 * 15, -1.0);
    for (int i = 0; i < 15; ++i) v[i + 15 * i] = 2.0;
    return Matrix<double>(15, 15, v);
}

static std::vector<std::pair<const int, const std::string>> load_reads(const char* path)
{
    try { return read_txt(path); }
    catch (...) { return {{0, "ACGTACGTACGT"}, {1, "TGCATGCATGCA"}}; }
}

int main(int argc, char** argv)
{
    std::string seq_path = argc > 1 ? argv[1] : std::string(IGOR_SOURCE_DIR) + "/demo/murugan_naive1_noncoding_demo_seqs.txt";
    auto reads_const = load_reads(seq_path.c_str());
    std::vector<std::pair<std::string, std::string>> germlines{{"synthetic", "ACGTACGTACGT"}};
    Matrix<double> subst = simple_matrix();

    Aligner legacy(subst, 2, V_gene);
    legacy.set_genomic_sequences(germlines);
    auto t0 = std::chrono::high_resolution_clock::now();
    auto legacy_res = legacy.align_seqs(reads_const, -100.0, true);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "{\n";
    std::cout << "  \"legacy_ms\": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << ",\n";
    std::cout << "  \"legacy_sequences\": " << legacy_res.size();

#ifdef IGOR_WITH_SEQAN2
    std::vector<std::pair<int, std::string>> reads;
    for (const auto& r : reads_const) reads.push_back({r.first, r.second});
    SeqAn2Aligner seqan2(subst, 2, V_gene, AlignmentPreset::with_band(50));
    seqan2.set_genomic_sequences(germlines);
    auto s0 = std::chrono::high_resolution_clock::now();
    auto seqan_res = seqan2.align_seqs(reads, -100.0, true);
    auto s1 = std::chrono::high_resolution_clock::now();
    std::cout << ",\n  \"seqan2_ms\": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(s1 - s0).count() << ",\n";
    std::cout << "  \"seqan2_sequences\": " << seqan_res.size();
#endif
    std::cout << "\n}\n";
    return 0;
}
