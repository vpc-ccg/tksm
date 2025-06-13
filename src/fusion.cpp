#include "fusion.h"

#include <cgranges/IITree.h>

#include <algorithm>
#include <cxxopts.hpp>
#include <fstream>
#include <optional>
#include <random>
#include <ranges>
#include <set>
#include <string>
#include <variant>
#include <vector>

#include "gtf.h"
#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

using std::set;
using std::string;
using std::unordered_map;
using std::vector;
namespace sr = std::ranges;


template <class B>
inline void
print_tsv(std::ostream &ost, B b) {
    ost << b << "\n";
}

template <class B, class... A>
inline void
print_tsv(std::ostream &ost, B b, A... a) {
    ost << b << "\t";
    print_tsv(ost, a...);
}

string fusion_separator = "::";
class locus {
    string chr;
    int position;
    bool plus_strand;

public:
    locus() : chr{}, position{0}, plus_strand{true} {}
    locus(const string &chr, int position, bool plus_strand) : chr{chr}, position{position}, plus_strand{plus_strand} {}
    bool operator<(const locus &other) const {
        if (chr == other.chr) {
            return position < other.position;
        }
        return chr < other.chr;
    }

    locus operator+(int offset) const { return locus{chr, position + offset, plus_strand}; }

    locus operator-(int offset) const { return locus{chr, position - offset, plus_strand}; }

    bool operator==(const locus &other) const {
        return chr == other.chr && position == other.position && plus_strand == other.plus_strand;
    }

    bool operator!=(const locus &other) const { return !(*this == other); }

    auto operator<=>(const locus &other) const {
        if (chr == other.chr) {
            return position <=> other.position;
        }
        return chr <=> other.chr;
    }

    auto get_chr() const -> const string & { return chr; }

    auto get_position() const -> int { return position; }

    auto is_plus_strand() const -> bool { return plus_strand; }
    auto to_string() const -> string { return chr + ":" + std::to_string(position) + (plus_strand ? "+" : "-"); }

    friend std::ostream &operator<<(std::ostream &os, const locus &loc) {
        os << loc.to_string();
        return os;
    }

    friend std::istream &operator>>(std::istream &is, locus &loc) {
        string chr;
        int position;
        string strand;
        is >> chr >> position >> strand;
        loc = locus{chr, position, strand == "+"};
        return is;
    }

    operator std::pair<string, int>() const { return {chr, position}; }
};
template <>
struct fmt::formatter<locus> : ostream_formatter {};

enum class EventType { DELETION, INVERSION, TRANSLOCATION, DUPLICATION, INSERTION, NONE };
inline auto
event_type_to_string(EventType type) -> string {
    switch (type) {
        case EventType::DELETION:
            return "deletion";
        case EventType::INVERSION:
            return "inversion";
        case EventType::TRANSLOCATION:
            return "translocation";
        case EventType::DUPLICATION:
            return "duplication";
        case EventType::INSERTION:
            return "insertion";
        case EventType::NONE:
            return "none";
        default:
            throw runtime_error("Invalid event type " + std::to_string(static_cast<int>(type)));
    }
}

inline auto
string_to_event_type(const string &st) -> EventType {
    string str{st};
    sr::transform(str.begin(), str.end(), str.begin(), ::tolower);
    if (str == "INVERSION") {
        return EventType::INVERSION;
    }
    else if (str == "DELETION") {
        return EventType::DELETION;
    }
    else if (str == "TRANSLOCATION") {
        return EventType::TRANSLOCATION;
    }
    else if (str == "DUPLICATION") {
        return EventType::DUPLICATION;
    }
    else if (str == "INSERTION") {
        return EventType::INSERTION;
    }
    else if (str == "NONE") {
        return EventType::NONE;
    }
    else {
        throw runtime_error("Invalid event type " + str);
    }
}

class chimeric_event : public ginterval {
public:
    enum class CUT { HEAD, TAIL };
    EventType event_type;
    double event_ratio;

    string event_name;
    string chr2;
    string orientation2;
    int count;
    // Event
    chimeric_event(EventType event_type) : ginterval{}, event_type{event_type}, event_ratio{0.5}, event_name{"::"} {}

    chimeric_event(const string &chr, int start, int end, const string &orientation, const string &orientation2,
                   const string &chr2, EventType event_type, const string &event_name, int count)
        : ginterval{chr, start, end, orientation},
          event_type{event_type},
          event_ratio{0.5},
          event_name{event_name},
          chr2{chr2},
          orientation2{orientation2},
          count{count} {}
    chimeric_event(const string &chr, int start, int end, const string &orientation, const string &orientation2,
                   const string &chr2, EventType event_type, const string &event_name)
        : ginterval{chr, start, end, orientation},
          event_type{event_type},
          event_ratio{0.5},
          event_name{event_name},
          chr2{chr2},
          orientation2{orientation2},
          count{0} {}

    chimeric_event(const chimeric_event &other) = default;
    auto cut_transcript(const transcript &t, int cut_position, CUT cut) const
        -> std::pair<transcript, std::optional<gtf>> {
        if (t.start > cut_position || t.end < cut_position) {
            logw("Cut position {} is not within transcript {}={}", cut_position, t.info.at("transcript_id"),
                 static_cast<ginterval>(t));
        }
        ginterval i = cut == CUT::TAIL ? ginterval{t.chr, cut_position, t.end, t.plus_strand}
                                       : ginterval{t.chr, t.start, cut_position, t.plus_strand};

        transcript newt{gtf{i, gtf::entry_type::transcript}, 0};
        bool exon_has_been_cut = false;
        gtf cpy_exon;
        for (const auto &exon : t.cget_exons()) {
            int overlap = i.overlap(exon);
            if (overlap == 0) {
                continue;
            }
            else if (overlap == exon.size()) {
                newt.add_exon(exon);
            }
            else {
                cpy_exon = exon;
                if (cut == CUT::TAIL) {
                    cpy_exon.start = cut_position;
                }
                else {
                    cpy_exon.end = cut_position;
                }
                cpy_exon.info["exon_id"] = cpy_exon.info["exon_id"];
                exon_has_been_cut        = true;
                newt.add_exon(cpy_exon);
            }
        }
        if (exon_has_been_cut) {
            return {newt, cpy_exon};
        }
        else {
            return {newt, std::nullopt};
        }
    }

    auto fuse_transcripts(const transcript &t1, const transcript &t2, auto rand_gen) const -> transcript {
        logd("{} {}", gtf::type_to_string(t1.type), gtf::type_to_string(t2.type));

        auto fusion_transcript_name =
            fmt::format("{}{}{}", t1.info.at("transcript_name"), fusion_separator, t2.info.at("transcript_name"));

        auto fusion_gene_name =
            fmt::format("{}{}{}", t1.info.at("gene_name"), fusion_separator, t2.info.at("gene_name"));
        auto fusion_transcript_id =
            fmt::format("{}{}{}", t1.info.at("transcript_id"), fusion_separator, t2.info.at("transcript_id"));
        auto fusion_gene_id = fmt::format("{}{}{}", t1.info.at("gene_id"), fusion_separator, t2.info.at("gene_id"));

        auto head_cut_orientation = [&]() -> CUT {
            if (event_type == EventType::DELETION) {
                return CUT::HEAD;
            }
            else if (event_type == EventType::INVERSION) {
                std::uniform_real_distribution<double> dist(0, 1);
                // Pick a random orientation, we might want to generate both, then we have to modify upstream code
                return dist(rand_gen) > 0.5 ? CUT::HEAD : CUT::TAIL;
            }
            else {
                throw runtime_error("Invalid event type " + event_type_to_string(event_type));
            }
        }();

        auto tail_cut_orientation = [&]() -> CUT {
            if (event_type == EventType::DELETION) {
                return CUT::TAIL;
            }
            else if (event_type == EventType::INVERSION) {
                return head_cut_orientation;  // Inversion cuts the same orientation as the head
            }
            else {
                throw runtime_error("Invalid event type " + event_type_to_string(event_type));
            }
        }();
        auto [fusion_transcript, head_cut_exon]      = cut_transcript(t1, this->start, head_cut_orientation);
        auto [pseudo_tail_transcript, tail_cut_exon] = cut_transcript(t2, this->end, tail_cut_orientation);
        if (head_cut_exon) {
            if (!tail_cut_exon) {
                head_cut_exon.value().end = this->start;  // Retain exon
            }
            fusion_transcript.add_exon(head_cut_exon.value());
            if (tail_cut_exon) {
                fusion_transcript.add_exon(tail_cut_exon.value());
            }
        }
        if (head_cut_orientation == CUT::TAIL) {
            sr::swap(fusion_transcript, pseudo_tail_transcript);
            // TODO: Reverse the exons if required, I might merge deletion and inversion after this
            // preprocessing
        }
        for (const auto &exon : pseudo_tail_transcript.get_exons()) {
            fusion_transcript.add_exon(exon);
        }
        int number = 1;

        // Fix GTF info
        for (auto &exon : fusion_transcript.get_exons()) {
            string exon_id = exon.info.at("exon_id");
            exon.info.clear();
            exon.info["exon_number"]       = std::to_string(number++);
            exon.info["transcript_name"]   = fusion_transcript_name;
            exon.info["gene_name"]         = fusion_gene_name;
            exon.info["transcript_id"]     = fusion_transcript_id;
            exon.info["gene_id"]           = fusion_gene_id;
            exon.info["transcript_source"] = "TKSM_" + event_type_to_string(event_type);
            exon.info["tag"]               = "TKSM_fusion";
        }

        fusion_transcript.info["transcript_name"] = fusion_transcript_name;
        fusion_transcript.info["gene_name"]       = fusion_gene_name;
        fusion_transcript.info["transcript_id"]   = fusion_transcript_id;
        fusion_transcript.info["gene_id"]         = fusion_gene_id;
        fusion_transcript.end                     = pseudo_tail_transcript.end;
        return fusion_transcript;
    }

    auto crunch_fusions(vector<transcript> &fusion_events) const -> void {
        std::unordered_map<string, transcript> id_to_transcript;
        for (auto &t : fusion_events) {
            auto [it, inserted] = id_to_transcript.try_emplace(t.info.at("transcript_id"), t);
            if (!inserted) {
                it->second.set_abundance(it->second.get_abundance() + t.get_abundance());
            }
        }
        fusion_events.clear();
        for (auto &[id, t] : id_to_transcript) {
            fusion_events.push_back(t);
        }
        sr::sort(fusion_events,
                 [](const auto &t1, const auto &t2) { return t1.info.at("gene_id") < t2.info.at("gene_id"); });
    }

    auto execute_event(const map<locus, vector<string>> &relavant_transcripts_by_loci, auto &rand_gen,
                       const auto &tid_to_tpcm, const auto &transcript_templates) const -> vector<transcript> {
        set<locus> loci;
        for (const auto &[loc, transcripts] : relavant_transcripts_by_loci) {
            loci.insert(loc);
        }
        auto it = loci.lower_bound(get_start());
        if (it == loci.end()) {
            throw runtime_error("Could not find start locus");
        }
        auto start_locus = *it;
        it               = loci.lower_bound(get_end());
        if (it == loci.end()) {
            for (const auto &l : loci) {
                logd("{}", l);
            }
            throw runtime_error("Could not find end locus");
        }
        auto end_locus = *it;

        auto &start_molecules = relavant_transcripts_by_loci.at(start_locus);
        auto &end_molecules   = relavant_transcripts_by_loci.at(end_locus);

        std::uniform_real_distribution<double> dist(0, 1);

        vector<transcript> fused_transcripts;

        double total_head_abundance = count;
        for (const auto &t : start_molecules) {
            auto iter = tid_to_tpcm.find(t);
            if (iter == tid_to_tpcm.end()) {
                continue;
            }
            for (const auto &cell_and_abundance : iter->second) {
                total_head_abundance += cell_and_abundance.second;
            }
        }

        double total_tail_abundance = 0;

        map<string, double> tail_abundance_flat;
        for (const auto &t : end_molecules) {
            auto iter = tid_to_tpcm.find(t);
            if (iter == tid_to_tpcm.end()) {
                continue;
            }
            for (const auto &cell_and_abundance : iter->second) {
                tail_abundance_flat[t] += cell_and_abundance.second;
                total_tail_abundance += cell_and_abundance.second;
            }
        }

        if (total_tail_abundance <= 0) {  // If tail abundance is 0 we will generate a uniform distribution
            logd("Total tail abundance is 0");
            logd("Generating uniform tail abundance");
            total_tail_abundance = end_molecules.size();
            for (auto &t : tail_abundance_flat) {
                t.second = 1;
            }
        }
        for (const auto &start_molecule : start_molecules) {
            auto sm_iter = tid_to_tpcm.find(start_molecule);
            if (sm_iter == tid_to_tpcm.end()) {
                logd("Could not find start molecule {} abundance", start_molecule);
                continue;
            }
            const auto sm_per_cell = sm_iter->second;
            for (const auto &cell_and_abundance : sm_per_cell) {
                double sma = cell_and_abundance.second;
                if (sma <= 0) {
                    logd("Head abundance is 0 skipping");
                    continue;
                }
                for (const auto &end_molecule : end_molecules) {
                    double ema = tail_abundance_flat[end_molecule];
                    logd("Fusing {} and {}", start_molecule, end_molecule);
                    fmtlog::poll(true);

                    if (ema <= 0) {
                        logd("Tail abundance is 0 skipping");
                        continue;
                    }
                    double abundance =
                        total_head_abundance * (sma / total_head_abundance) * (ema / total_tail_abundance);

                    fused_transcripts.push_back(fuse_transcripts(transcript_templates.at(start_molecule),
                                                                 transcript_templates.at(end_molecule), rand_gen));
                    logd(
                        "Total head_abundance: {}, head_abundance: {}, total_tail_abundance: {}, tail_abundance: {}, "
                        "abundance: {}",
                        total_head_abundance, sma, total_tail_abundance, ema, abundance);
                    fused_transcripts.back().set_abundance(abundance);
                    fused_transcripts.back().info["CB"] = cell_and_abundance.first;
                }
            }
        }
        crunch_fusions(fused_transcripts);
        return fused_transcripts;
    }

    locus get_start(bool strand = true) const { return locus{chr, start, strand}; }
    locus get_end(bool strand = true) const { return locus{chr, end, strand}; }
    friend ostream &operator<<(ostream &os, const chimeric_event &event) {
        os << event.chr << "\t" << event.start << "\t" << event.end << "\t" << event_type_to_string(event.event_type)
           << "\t" << event.chr2 << "\t" << event.event_name;
        return os;
    }
};
template <>
struct fmt::formatter<chimeric_event> : ostream_formatter {};

inline vector<chimeric_event>
read_fusions(std::istream &fusion_file) {
    vector<chimeric_event> fusions;
    string line;
    while (std::getline(fusion_file, line)) {
        std::istringstream iss{line};
        string chr1, chr2, orientation1, orientation2, event_name;
        int start1, end1;
        double count;
        iss >> chr1 >> start1 >> end1 >> orientation1 >> orientation2 >> chr2 >> event_name >> count;
        std::stringstream oss;
        oss << chr1 << ":" << start1 << "-" << end1 << orientation1;
        EventType et = [&]() -> EventType {
            if (chr1 == chr2) {
                if (orientation1 == orientation2) {
                    return EventType::DELETION;
                }
                else {
                    return EventType::DUPLICATION;
                }
            }
            else {
                return EventType::TRANSLOCATION;
            }
        }();
        chimeric_event fusion{chr1, start1, end1, orientation1, orientation2, chr2, et, event_name, (int)count};
        fusions.push_back(fusion);
    }
    return fusions;
}

class Fusion_submodule : public tksm_submodule {
    void add_options(cxxopts::Options &options) override {
        // clang-format off
        options.add_options("submodule:fusion")
            (
                "fusion-gtf",
                "Path to output GTF annotation file",
                cxxopts::value<string>()
            )(
                "fusion-file",
                "Path to tab separated fusion file",
                cxxopts::value<string>()
            )(
                "fusion-output",
                "Description of generated fusions",
                cxxopts::value<string>()
            )(
                "fusion-count",
                "Number of random fusions to generate",
                cxxopts::value<int>()->default_value("0")
            )(
                "disable-deletions",
                "Disables deletions (from fusions) that removes expression on the overlapping genes",
                cxxopts::value<bool>()->default_value("false")->implicit_value("true")
            )(
                "translocation-ratio",
                "Ratio of translocated fusions",
                cxxopts::value<double>()->default_value("0")
            )(
                "expression-fallback",
                "Fallback expression distribution for transcripts that are not expressed. "
                "If not provided only expressed transcripts will be fused. "
                "Comma separated list of distribution type and parameters. "
                "Available distributions: uniform,s,e; normal,μ,σ. ",
                cxxopts::value<vector<string>>()
            )
            ;
        // clang-format on
    }

    auto generate_random_breakpoint(const gtf &gene) -> locus {
        std::uniform_int_distribution<int> breakpoint_selector(gene.start, gene.end);
        return locus{gene.chr, breakpoint_selector(rand_gen), gene.plus_strand};
    }

    auto process_expression_parameters(const std::vector<string> &vec) const
        -> std::variant<std::uniform_real_distribution<>, std::normal_distribution<>> {
        if (vec.size() != 3) {
            throw runtime_error("Invalid number of parameters for expression distribution");
        }
        if (vec[0] == "uniform") {
            return std::uniform_real_distribution<>{std::stod(vec[1]), std::stod(vec[2])};
        }
        else if (vec[0] == "normal") {
            return std::normal_distribution<>{std::stod(vec[1]), std::stod(vec[2])};
        }
        else {
            throw runtime_error("Invalid distribution type " + vec[0]);
        }
    }

    auto generate_fusions(vector<chimeric_event> &fusions_so_far, uint64_t total_fusion_count,
                          const map<string, std::pair<gtf, vector<gtf>>> &genes, double translocation_ratio)
        -> vector<chimeric_event> {
        // Convert genes vector to a map of chr -> vector<gtf>

        std::map<string, vector<gtf>> genes_by_chr;
        for (auto &[gene_id, gene_pair] : genes) {
            const auto &gene = gene_pair.first;
            genes_by_chr[gene.chr].push_back(gene);
        }
        uint64_t fusion_count = total_fusion_count - fusions_so_far.size();
        // Assuming distance between the first and last gene of chromosome estimates the size of the chromosome
        std::map<string, unsigned long> chr_size;
        std::map<string, uint64_t> fusion_count_per_chr;
        unsigned long total_size = 0;
        for (auto &gp : genes_by_chr) {
            auto &chr        = gp.first;
            auto &genes_copy = gp.second;
            sr::sort(genes_copy, [](const auto &a, const auto &b) { return a.start < b.start; });
            chr_size[chr] = genes_copy.back().end - genes_copy.front().start;
            total_size += chr_size[chr];
            logd("Chr {} size is {}, {}:{}-{}", chr, chr_size[chr], chr, genes_copy.front().start,
                 genes_copy.back().end);
        }

        for (const auto &gp : genes_by_chr) {
            auto &chr = gp.first;

            fusion_count_per_chr[chr] = std::round(static_cast<double>(fusion_count) * chr_size[chr] / total_size);
            logd("Fusion count for chr {} is {}", chr, fusion_count_per_chr[chr]);
        }
        uint64_t calculated_fusions = std::accumulate(fusion_count_per_chr.begin(), fusion_count_per_chr.end(), 0UL,
                                                      [](uint64_t a, const auto &b) { return a + b.second; });
        if (fusion_count > calculated_fusions) {
            logd("Adding {} fusions to random chromosomes", fusion_count - calculated_fusions);

            for (uint64_t i = 0; i < fusion_count - calculated_fusions; i++) {
                std::vector<std::pair<string, int>> v;
                sr::sample(fusion_count_per_chr, std::back_inserter(v), 1, rand_gen);
                fusion_count_per_chr[v.begin()->first]++;
            }
        }
        // Convert fusion vector to a map of chr -> vector<chimeric_event>

        std::map<string, vector<chimeric_event>> fusions_by_chr;
        for (auto &fusion : fusions_so_far) {
            fusions_by_chr[fusion.chr].push_back(fusion);
        }

        for (const auto &fp : fusion_count_per_chr) {
            auto &chr                         = fp.first;
            auto &fusions                     = fusions_by_chr[chr];
            auto &genes_copy                  = genes_by_chr[chr];
            uint64_t fusion_count_of_this_chr = fp.second;
            for (auto &fusion : fusions) {
                auto it = std::remove_if(genes_copy.begin(), genes_copy.end(), [&fusion](const gtf &gene) {
                    return gene.chr == fusion.chr && gene.start >= fusion.start && gene.end <= fusion.end;
                });
                genes_copy.erase(it, genes_copy.end());
            }

            sr::shuffle(genes_copy, rand_gen);
            // Take the first fusion_count_per_chr[chr] * 2 genes and sort it back
            if (fusion_count_of_this_chr * 2 > genes_copy.size()) {
                logw("Not enough genes to generate {} fusions on chr {}", fusion_count_of_this_chr, chr);
                logw("Reducing fusion count to {}", genes_copy.size() / 2);
                fusion_count_of_this_chr = genes_copy.size() / 2;
            }
            sr::sort(genes_copy.begin(), genes_copy.begin() + fusion_count_of_this_chr * 2,
                     [](const gtf &a, const gtf &b) { return a.start < b.start; });

            // Generate chimeric events for each adjacent gene pair in the chromosome as many as
            // fusion_count_per_chr[chr]
            for (uint64_t i = 0; i < fusion_count_of_this_chr; i += 2) {
                auto &gene1 = genes_copy[i];
                auto &gene2 = genes_copy[i + 1];
                if (gene1.overlap(gene2) > 0) {
                    logd("Skipping fusion between overlapping genes {} and {}", gene1.info["gene_name"],
                         gene2.info["gene_name"]);
                    continue;
                }
                // if genes are on the same strand generate deletion else inversion

                EventType event_type =
                    gene1.plus_strand != gene2.plus_strand ? EventType::INVERSION : EventType::DELETION;

                locus g1 = generate_random_breakpoint(gene1);
                locus g2 = generate_random_breakpoint(gene2);
                chimeric_event fusion{g1.get_chr(),
                                      g1.get_position(),
                                      g2.get_position(),
                                      g1.is_plus_strand() ? "+" : "-",
                                      g2.is_plus_strand() ? "+" : "-",
                                      g2.get_chr(),
                                      event_type,
                                      gene1.info.at("gene_name") + "::" + gene2.info.at("gene_name")};
                logd("Generated fusion: {}, on genes {} and {}", fusion, gene1.info["gene_name"],
                     gene2.info["gene_name"]);
                fusions_so_far.push_back(fusion);
            }
        }
        return fusions_so_far;
    }

    auto index_fusions(const vector<chimeric_event> &fusions)
        -> std::tuple<IITree<locus, chimeric_event>, IITree<locus, chimeric_event>> const {
        IITree<locus, chimeric_event> fusion_tree;
        IITree<locus, chimeric_event> deletion_tree;
        int RANGE = 0;
        logi("Indexing fusions");
        for (const auto &fusion : fusions) {
            fusion_tree.add(fusion.get_start() - RANGE, fusion.get_start() + RANGE, fusion);
            fusion_tree.add(fusion.get_end() - RANGE, fusion.get_end() + RANGE, fusion);
            fusion_tree.add(fusion.get_start(false) - RANGE, fusion.get_start(false) + RANGE, fusion);
            fusion_tree.add(fusion.get_end(false) - RANGE, fusion.get_end(false) + RANGE, fusion);

            if (fusion.event_type == EventType::DELETION) {
                deletion_tree.add(fusion.get_start(), fusion.get_end(), fusion);
                deletion_tree.add(fusion.get_start(false), fusion.get_end(false), fusion);
            }
        }
        fusion_tree.index();
        deletion_tree.index();
        return {fusion_tree, deletion_tree};
    }

    auto get_fusions(const auto &genes,           // string->gtf map
                     const auto &expression_map,  // string->double map
                     const cxxopts::ParseResult &args) {
        unsigned fusion_count      = args["fusion-count"].as<int>();
        double translocation_ratio = args["translocation-ratio"].as<double>();

        vector<chimeric_event> fusions;
        if (args["fusion-file"].count() > 0) {
            string fusion_file_path{args["fusion-file"].as<string>()};
            std::ifstream fusion_file{fusion_file_path};
            fusions = read_fusions(fusion_file);
        }

        logi("Generating fusions");
        if (args["expression-fallback"].count() > 0) {
            while (fusions.size() < fusion_count) {
                logd("Generating fusions");
                fmtlog::poll(true);
                generate_fusions(fusions, fusion_count, genes, translocation_ratio);
                logd("Generated {}/{} fusions", fusions.size(), fusion_count);
                fmtlog::poll(true);
            }
            logi("Using expression fallback");
        }
        else {  // Only simulate fusions from expressed transcripts
            while (fusions.size() < fusion_count) {
                logd("Generating fusions");
                fmtlog::poll(true);

                map<string, std::pair<gtf, vector<gtf>>> expressed_genes;
                for (const auto &[name, gene] : genes) {
                    if (gene.first.info.at("gene_biotype") != "protein_coding") continue;
                    auto iter = expression_map.find(name);
                    if (iter != expression_map.end() && iter->second > 0) {
                        expressed_genes[name] = gene;
                    }
                }
                generate_fusions(fusions, fusion_count, expressed_genes, translocation_ratio);
                logd("Generated {}/{} fusions", fusions.size(), fusion_count);
                fmtlog::poll(true);
            }
        }
        return fusions;
    }
    gtf combine_gene_entries(const gtf &gene1, const gtf &gene2) {
        gtf gene;
        gene.chr                  = gene1.chr;
        gene.source               = "TKSM";
        gene.start                = std::min(gene1.start, gene2.start);
        gene.end                  = std::max(gene1.end, gene2.end);
        gene.plus_strand          = gene1.plus_strand;
        gene.info                 = gene1.info;
        gene.info["gene_name"]    = gene1.info.at("gene_name") + "::" + gene2.info.at("gene_name");
        gene.info["gene_id"]      = gene1.info.at("gene_id") + "::" + gene2.info.at("gene_id");
        gene.info["gene_source"]  = "TKSM_fusion";
        gene.info["gene_version"] = "1";
        gene.type                 = gtf::entry_type::gene;

        return gene;
    }

    auto read_genes(const vector<string> &gtf_files) -> map<string, std::pair<gtf, vector<gtf>>> {
        map<string, std::pair<gtf, vector<gtf>>> genes;
        for (const string &gtf_file : gtf_files) {
            for (const auto &[g, trans_vec] : read_gtf_genes(gtf_file, true, false)) {
                genes[g.info.at("gene_id")] = {g, trans_vec};
            }
        }
        return genes;
    }

    auto count_expressions(auto &abundances, auto &transcript_templates, const auto &tid_to_gene)
        -> std::pair<map<string, double>, unordered_map<string, unordered_map<string, double>>> {
        map<string, double> expression_map;
        unordered_map<string, unordered_map<string, double>> transcript_cell_to_abundance;
        for (const auto &[tid, tpm, comment] : abundances) {
            auto iter = transcript_templates.find(tid);
            if (iter != transcript_templates.end()) {
                //                iter->second->set_abundance(tpm);
                //                transcript cpy            = *iter->second;
                //                cpy.info["transcript_id"] = tid;
                //                cpy.info["CB"]            = comment;
                //                transcripts.emplace_back(cpy);
                const auto &gene                           = tid_to_gene.at(tid);
                transcript_cell_to_abundance[tid][comment] = tpm;
                expression_map[gene.info.at("gene_id")] += tpm;
            }
        }
        return {expression_map, transcript_cell_to_abundance};
    }

    auto find_relevant_molecules(const auto &transcripts, const auto &fusion_tree) {
        //        map<chimeric_event, map<locus, vector<transcript>>> relevant_molecules;
        map<chimeric_event, map<locus, vector<string>>> relevant_molecules_indices;
        for (const auto &[name, t] : transcripts) {
            string chr       = t.chr;
            bool plus_strand = t.plus_strand;
            locus start      = {chr, t.start, plus_strand};
            locus end        = {chr, t.end, plus_strand};

            vector<size_t> overlaps;
            fusion_tree.overlap(start, end, overlaps);

            for (auto &overlap : overlaps) {
                const chimeric_event &event       = fusion_tree.data(overlap);
                const locus &start                = fusion_tree.start(overlap);
                [[maybe_unused]] const locus &end = fusion_tree.end(overlap);
                relevant_molecules_indices[event][start].push_back(name);
                //                relevant_molecules[event][start].push_back(t);
            }
            if (overlaps.size() > 0) {
                logd("Found {} overlaps, for {}:{}-{}", overlaps.size(), chr, start, end);
            }
            overlaps.clear();
        }

        return relevant_molecules_indices;
    }

    auto update_expression_of_affected_transcripts(const auto &relevant_molecules, auto &abundances) {
        // Build a tid to modification ratio map by iterating relevant molecules

        unordered_map<string, double> tid_to_modification_ratio;

        for (const auto &[event, locus_to_tids] : relevant_molecules) {
            for (const auto &[locus, tids] : locus_to_tids) {
                for (const auto &tid : tids) {
                    tid_to_modification_ratio[tid] = 1 - event.event_ratio;
                }
            }
        }

        // Iterate over all transcripts and update their expression

        for (auto &[tid, tpm, comment] : abundances) {
            auto iter = tid_to_modification_ratio.find(tid);
            if (iter != tid_to_modification_ratio.end()) {
                tpm *= iter->second;
            }
        }
    }

    auto compute_fusion_transcripts(map<chimeric_event, map<locus, vector<string>>> &relevant_molecules,
                                    const auto &genes, auto &gene_expression_map,
                                    const unordered_map<string, unordered_map<string, double>> &tid_to_tpm,
                                    const auto &tid_to_gene, const auto &transcript_templates, const auto &args)
        -> generator<std::tuple<chimeric_event, transcript>> {
        bool do_fall_back_to_expression = args.count("expression-fallback") > 0;
        auto tpm_dist                   = [&]() {
            if (do_fall_back_to_expression) {
                return process_expression_parameters(args["expression-fallback"].template as<vector<string>>());
            }
            else {
                return std::variant<std::uniform_real_distribution<>, std::normal_distribution<>>{
                    std::normal_distribution<>{0, 0}};
            }
        }();
        map<gtf, vector<transcript>> gene2transcripts_of_fusions;
        string current_gene = "";
        gtf current_gene_obj;
        for (auto &[event, loci] : relevant_molecules) {
            if (do_fall_back_to_expression) {
                for (auto &l : loci) {
                    for (const string &tid : l.second) {
                        const string &gene_id = tid_to_gene.at(tid).info.at("gene_id");
                        const auto &it        = gene_expression_map.find(gene_id);
                        if (it == gene_expression_map.end()) {
                            logd("Could not find tpm for {}", tid);
                            double tpm = std::visit([&](auto &&arg) { return arg(rand_gen); }, tpm_dist);
                            logd("Using fallback tpm={}", tpm);
                            gene_expression_map[gene_id] = tpm;
                        }
                    }
                }
            }
            auto transcript_vec = event.execute_event(loci, rand_gen, tid_to_tpm, transcript_templates);
            logd("Creating {} fusion transcripts on event {}", transcript_vec.size(), event);

            for (auto &t : transcript_vec) {
                if (t.info["gene_id"] != current_gene) {
                    current_gene              = t.info["gene_id"];
                    vector<string> gene_names = rsplit(current_gene, "::");
                    current_gene_obj =
                        combine_gene_entries(genes.at(gene_names[0]).first, genes.at(gene_names[1]).first);
                }

                if (t.cget_exons().empty()) {
                    logd("Skipping empty fusion transcript {}", t.info.at("transcript_id"));
                    continue;
                }
                co_yield {event, t};
            }
        }
    }

public:
    Fusion_submodule(cxxopts::Options &opt, std::mt19937 &rand_gen)
        : tksm_submodule{"fusion", "Fusion module", rand_gen} {
        add_options(opt);
    }
    ~Fusion_submodule() = default;
    submodule_status receive_arguments(const cxxopts::ParseResult &args) override {
        if (args["fusion-count"].as<int>() == 0 && args["fusion-file"].count() == 0) {
            return update_status(submodule_status::DONT_RUN);
        }

        return update_status(submodule_status::RUN);
    }

    template <class TopModule>
    auto run(TopModule *top_module, vector<std::tuple<string, double, string>> &abundances,
             unordered_map<string, transcript> &transcript_templates) -> int {
        cxxopts::ParseResult &args = top_module->args;

        vector<string> gtf_files                    = args["gtf"].as<vector<string>>();
        [[maybe_unused]] size_t fusion_count        = args["fusion-count"].as<int>();
        [[maybe_unused]] double translocation_ratio = args["translocation-ratio"].as<double>();
        [[maybe_unused]] bool disable_deletions     = args["disable-deletions"].as<bool>();

        if (args.count("fusion-output") == 0) {
            loge("No fusion output file specified");
            return 1;
        }

        auto genes = read_genes(gtf_files);
        unordered_map<string, gtf> transcript_id_to_genes;

        for (auto &[gene, trans_vec] : genes) {
            for (auto &transcript : trans_vec.second) {
                transcript_id_to_genes[transcript.info.at("transcript_id")] = trans_vec.first;
            }
        }

        auto [gene_expression_map, transcript_to_tpcm] =
            count_expressions(abundances, transcript_templates, transcript_id_to_genes);
        vector<chimeric_event> fusions    = get_fusions(genes, gene_expression_map, args);
        auto [fusion_tree, deletion_tree] = index_fusions(fusions);

        std::map<chimeric_event, std::map<locus, vector<string>>> relevant_molecules =
            find_relevant_molecules(transcript_templates, fusion_tree);

        std::ofstream fusion_out{args["fusion-output"].as<string>()};

        logi("Creating fusion transcripts");
        for (const auto &[event, t] :
             compute_fusion_transcripts(relevant_molecules, genes, gene_expression_map, transcript_to_tpcm,
                                        transcript_id_to_genes, transcript_templates, args)) {
            abundances.emplace_back(t.info.at("transcript_id"), t.get_abundance(), t.info.at("CB"));
            transcript_templates.emplace(t.info.at("transcript_id"), t);
            print_tsv(fusion_out, event, t.info.at("gene_id"), t.info.at("gene_name"), t.info.at("transcript_id"),
                      t.info.at("transcript_name"), t.get_abundance());
        }
        update_expression_of_affected_transcripts(relevant_molecules, abundances);
        return 0;
    }

    void describe_program(const cxxopts::ParseResult &args) override {
        if (status != submodule_status::RUN) {
            return;
        }
        logi("Fusion submodule");

        if (args["fusion-file"].count() > 0) {
            logi("Fusion file: {}", args["fusion-file"].as<string>());
        }
        logi("Fusion count: {}", args["fusion-count"].as<int>());
        if (args["disable-deletions"].as<bool>()) {
            logi("Deletions are disabled");
        }
        logi("Translocation ratio: {}", args["translocation-ratio"].as<double>());

        if (args["expression-fallback"].count() > 0) {
            logi("Expression fallback: {}", fmt::join(args["expression-fallback"].as<vector<string>>(), ","));
        }
        else {
            logi("Expression fallback: Disabled");
        }
    }
};
