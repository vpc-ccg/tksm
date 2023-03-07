#include "fusion.h"

#include <cgranges/IITree.h>

#include <cxxopts.hpp>
#include <optional>
#include <random>
#include <ranges>
#include <set>
#include <string>
#include <vector>

#include "gtf.h"
#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "util.h"

using std::set;
using std::string;
using std::vector;

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
            throw std::runtime_error("Invalid event type " + std::to_string(static_cast<int>(type)));
    }
}

class chimeric_event : public ginterval {
public:
    enum class CUT { HEAD, TAIL };
    EventType event_type;
    double event_ratio;
    ginterval supplementary_interval;  // To be used for translocations

    // Event
    chimeric_event(EventType event_type) : ginterval{}, event_type{event_type}, event_ratio{0.5} {}
    chimeric_event(const string &chr, int start, int end, const string &orientation, EventType event_type)
        : ginterval{chr, start, end, orientation}, event_type{event_type}, event_ratio{0.5} {}

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
        for (const auto &exon : t.get_exons()) {
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

    auto fuse_transcripts(const transcript &t1, const transcript &t2) const -> transcript {
        logd("{} {}", gtf::type_to_string(t1.type), gtf::type_to_string(t2.type));

        auto fusion_transcript_name =
            fmt::format("{}{}{}", t1.info.at("transcript_name"), fusion_separator, t2.info.at("transcript_name"));

        auto fusion_gene_name =
            fmt::format("{}{}{}", t1.info.at("gene_name"), fusion_separator, t2.info.at("gene_name"));
        auto fusion_transcript_id =
            fmt::format("{}{}{}", t1.info.at("transcript_id"), fusion_separator, t2.info.at("transcript_id"));
        auto fusion_gene_id = fmt::format("{}{}{}", t1.info.at("gene_id"), fusion_separator, t2.info.at("gene_id"));

        switch (event_type) {
            case EventType::DELETION:
                if (t1.plus_strand != t2.plus_strand) {
                    // Cannot fuse
                }
                else {
                    auto [fusion_transcript, head_cut_exon] =
                        cut_transcript(t1, this->start, CUT::HEAD);  // Start from the head
                    auto [pseudo_tail_transcript, tail_cut_exon] = cut_transcript(t2, this->end, CUT::TAIL);

                    if (head_cut_exon) {
                        if (!tail_cut_exon) {
                            head_cut_exon.value().end = this->start;  // Retain exon
                        }
                        fusion_transcript.add_exon(head_cut_exon.value());
                        fusion_transcript.add_exon(tail_cut_exon.value());
                    }

                    for (const auto &exon : pseudo_tail_transcript.get_exons()) {
                        fusion_transcript.add_exon(exon);
                    }

                    int number = 1;
                    // Fix GTF info
                    for (auto &exon : fusion_transcript.get_exons()) {
                        exon.info["exon_number"]       = std::to_string(number++);
                        exon.info["transcript_name"]   = fusion_transcript_name;
                        exon.info["gene_name"]         = fusion_gene_name;
                        exon.info["transcript_id"]     = fusion_transcript_id;
                        exon.info["gene_id"]           = fusion_gene_id;
                        exon.info["transcript_source"] = "TKSM";
                        exon.info["tag"]               = "TKSM-fusion";
                    }

                    fusion_transcript.info["transcript_name"] = fusion_transcript_name;
                    fusion_transcript.info["gene_name"]       = fusion_gene_name;
                    fusion_transcript.info["transcript_id"]   = fusion_transcript_id;
                    fusion_transcript.info["gene_id"]         = fusion_gene_id;
                }
                break;
            case EventType::INVERSION:
                break;
            case EventType::TRANSLOCATION:
                break;
            case EventType::DUPLICATION:
                break;
            case EventType::INSERTION:
                break;
            case EventType::NONE:
                break;
            default:
                throw std::runtime_error("Invalid event type " + std::to_string(static_cast<int>(event_type)));
        }

        return t1;
    }

    auto execute_event(const map<locus, vector<transcript>> &molecules, auto &rand_gen) const -> vector<transcript> {
        set<locus> loci;
        for (const auto &[loc, transcripts] : molecules) {
            loci.insert(loc);
        }
        auto it = loci.lower_bound(get_start());
        if (it == loci.end()) {
            throw std::runtime_error("Could not find start locus");
        }
        auto start_locus = *it;
        it               = loci.lower_bound(get_end());
        if (it == loci.end()) {
            throw std::runtime_error("Could not find end locus");
        }
        auto end_locus = *it;

        auto &start_molecules = molecules.at(start_locus);
        auto &end_molecules   = molecules.at(end_locus);

        std::uniform_real_distribution<double> dist(0, 1);

        vector<transcript> fused_transcripts;
        for (const auto &start_molecule : start_molecules) {
            for (const auto &end_molecule : end_molecules) {
                logd("Fusing {} and {}", start_molecule, end_molecule);
                fmtlog::poll(true);
                fused_transcripts.push_back(fuse_transcripts(start_molecule, end_molecule));
            }
        }
        return fused_transcripts;
    }

    locus get_start(bool strand = true) const { return locus{chr, start, strand}; }
    locus get_end(bool strand = true) const { return locus{chr, end, strand}; }
    friend ostream &operator<<(ostream &os, const chimeric_event &event) {
        os << event.chr << "\t" << event.start << "\t" << event.end << "\t" << event_type_to_string(event.event_type);
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
        string chr1, chr2, orientation1, orientation2;
        int start1, end1, start2, end2;
        iss >> chr1 >> start1 >> end1 >> orientation1 >> chr2 >> start2 >> end2 >> orientation2;
        chimeric_event fusion{chr1, start1, end1, orientation1, EventType::NONE};
        fusions.push_back(fusion);
    }
    return fusions;
}

class Fusion_module::impl : public tksm_module {
    cxxopts::ParseResult parse(int argc, char **argv) {
        // clang-format off
        options.add_options("main")
            (
                "a,abundance",
                "input mdf file",
                cxxopts::value<string>()
            )(
                "o,output",
                "output mdf file",
                cxxopts::value<string>()
            )(
                "gtfi",
                "Path to input GTF annotation file",
                cxxopts::value<string>()
            )(
                "gtfo",
                "Path to output GTF annotation file",
                cxxopts::value<string>()
            )(
                "fusion-file",
                "Path to tab separated fusion file",
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
                "use-whole-id",
                "Use whole gene ids",
                cxxopts::value<bool>()->default_value("false")->implicit_value("true")
            )(
                "translocation-ratio",
                "Ratio of translocated fusions",
                cxxopts::value<double>()->default_value("0")
            )
            ;
        // clang-format on

        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

    auto generate_fusions(vector<chimeric_event> &fusions_so_far, int fusion_count, const vector<gtf> &genes,
                          double translocation_ratio) -> vector<chimeric_event> {
        // Convert genes vector to a map of chr -> vector<gtf>

        std::map<string, vector<gtf>> genes_by_chr;
        for (auto &gene : genes) {
            genes_by_chr[gene.chr].push_back(gene);
        }

        // Assuming distance between the first and last gene of chromosome estimates the size of the chromosome
        std::map<string, unsigned long> chr_size;
        std::map<string, int> fusion_count_per_chr;
        unsigned long total_size = 0;
        for (auto &gp : genes_by_chr) {
            auto &chr   = gp.first;
            auto &genes = gp.second;
            std::ranges::sort(genes, [](const auto &a, const auto &b) { return a.start < b.start; });
            chr_size[chr] = genes.back().end - genes.front().start;
            total_size += chr_size[chr];
            logd("Chr {} size is {}, {}:{}-{}", chr, chr_size[chr], chr, genes.front().start, genes.back().end);
        }

        for (const auto &gp : genes_by_chr) {
            auto &chr = gp.first;

            fusion_count_per_chr[chr] = std::round(static_cast<double>(fusion_count) * chr_size[chr] / total_size);
            logi("Fusion count for chr {} is {}", chr, fusion_count_per_chr[chr]);
        }
        // Convert fusion vector to a map of chr -> vector<chimeric_event>

        std::map<string, vector<chimeric_event>> fusions_by_chr;
        for (auto &fusion : fusions_so_far) {
            fusions_by_chr[fusion.chr].push_back(fusion);
        }

        for (const auto &fp : fusion_count_per_chr) {
            auto &chr        = fp.first;
            auto &fusions    = fusions_by_chr[chr];
            auto &genes_copy = genes_by_chr[chr];
            int fusion_count = fp.second;
            for (auto &fusion : fusions) {
                auto it = std::remove_if(genes_copy.begin(), genes_copy.end(), [&fusion](const gtf &gene) {
                    return gene.chr == fusion.chr && gene.start >= fusion.start && gene.end <= fusion.end;
                });
                genes_copy.erase(it, genes_copy.end());
            }
            std::shuffle(genes_copy.begin(), genes_copy.end(), rand_gen);
            // Take the first fusion_count_per_chr[chr] * 2 genes and sort it back

            std::sort(genes_copy.begin(), genes_copy.begin() + fusion_count_per_chr[chr] * 2,
                      [](const gtf &a, const gtf &b) { return a.start < b.start; });

            // Generate chimeric events for each adjacent gene pair in the chromosome as many as
            // fusion_count_per_chr[chr]
            for (int i = 0; i < fusion_count; i += 2) {
                auto &gene1 = genes_copy[i];
                auto &gene2 = genes_copy[i + 1];

                // if genes are on the same strand generate deletion else inversion

                EventType event_type =
                    gene1.plus_strand != gene2.plus_strand ? EventType::INVERSION : EventType::DELETION;

                chimeric_event fusion{gene1.chr, gene1.end, gene2.start, gene1.plus_strand ? "+" : "-", event_type};
                logi("Generated fusion: {}", fusion);
                fusions_so_far.push_back(fusion);
            }
        }
        return fusions_so_far;
    }

public:
    impl(int argc, char **argv) : tksm_module{"fusion", "Fusion module"}, args(parse(argc, argv)) {}
    ~impl() = default;
    int validate_arguments() {
        std::vector<string> mandatory = {"abundance", "output", "gtfi", "gtfo"};
        int missing_parameters        = 0;
        for (string &param : mandatory) {
            if (args.count(param) == 0) {
                report_missing_parameter(param);
                ++missing_parameters;
            }
        }

        if (missing_parameters > 0) {
            std::cerr << options.help() << std::endl;
            return 1;
        }

        if (args["fusion-count"].as<int>() == 0 && args["fusion-file"].as<string>().empty()) {
            loge("Either fusion-file or fusion-count must be specified");
            return 1;
        }

        return 0;
    }

    int run() {
        if (process_utility_arguments(args)) {
            return 0;
        }

        if (validate_arguments()) {
            return 1;
        }
        describe_program();

        std::string gtf_file       = args["gtfi"].as<string>();
        std::string abundance_file = args["abundance"].as<string>();
        std::string output_file    = args["output"].as<string>();

        int fusion_count                        = args["fusion-count"].as<int>();
        double translocation_ratio              = args["translocation-ratio"].as<double>();
        [[maybe_unused]] bool disable_deletions = args["disable-deletions"].as<bool>();

        vector<chimeric_event> fusions;
        if (args["fusion-file"].count() > 0) {
            string fusion_file_path{args["fusion-file"].as<string>()};
            std::ifstream fusion_file{fusion_file_path};
            fusions = read_fusions(fusion_file);
        }

        vector<transcript> transcripts = read_gtf_transcripts_deep(gtf_file);
        vector<gtf> genes              = read_gtf_genes(gtf_file);
        map<string, gtf> transcripts_by_id;
        for (auto &transcript : transcripts) {
            transcripts_by_id[transcript.info["transcript_id"]] = transcript;
        }

        logi("Generating fusions");
        generate_fusions(fusions, fusion_count, genes, translocation_ratio);

        IITree<locus, chimeric_event> fusion_tree;
        IITree<locus, chimeric_event> deletion_tree;

        logi("Indexing fusions");
        for (auto &fusion : fusions) {
            fusion_tree.add(fusion.get_start(), fusion.get_start() + 1, fusion);
            fusion_tree.add(fusion.get_end(), fusion.get_end() + 1, fusion);
            fusion_tree.add(fusion.get_start(false), fusion.get_start(false) + 1, fusion);
            fusion_tree.add(fusion.get_end(false), fusion.get_end(false) + 1, fusion);

            if (fusion.event_type == EventType::DELETION) {
                deletion_tree.add(fusion.get_start(), fusion.get_end(), fusion);
                deletion_tree.add(fusion.get_start(false), fusion.get_end(false), fusion);
            }
        }
        fusion_tree.index();
        deletion_tree.index();

        std::ifstream abundance_file_stream{abundance_file};
        string buffer;
        std::getline(abundance_file_stream, buffer);  // Header
        std::map<chimeric_event, std::map<locus, vector<transcript>>> relevant_molecules;
        std::set<string> processed_transcripts;
        std::ofstream output_file_stream{output_file};

        logi("Processing abundance file");
        int total_count = 0;
        int skipped     = 0;

        for (const auto &t : transcripts) {
            locus start{t.chr, t.start, t.plus_strand};
            locus end{t.chr, t.end, t.plus_strand};
            vector<size_t> overlaps;
            fusion_tree.overlap(start, end, overlaps);
            logi("Found {} overlaps, for {}:{}-{}", overlaps.size(), t.chr, t.start, t.end);
            for (auto &overlap : overlaps) {
                const chimeric_event &event       = fusion_tree.data(overlap);
                const locus &start                = fusion_tree.start(overlap);
                [[maybe_unused]] const locus &end = fusion_tree.end(overlap);
                relevant_molecules[event][start].push_back(t);
            }
            if (overlaps.size() > 0) {
                processed_transcripts.insert(t.info.at("transcript_id"));
            }
            overlaps.clear();
        }

        map<string, double> fusion_tpm;
        while (std::getline(abundance_file_stream, buffer)) {
            ++total_count;
            string tid;
            double tpm;
            string comment;
            std::istringstream(buffer) >> tid >> tpm >> comment;

            format_annot_id(tid, !args["use-whole-id"].as<bool>());
            auto iter = transcripts_by_id.find(tid);
            if (iter == transcripts_by_id.end()) {
                ++skipped;
                continue;
            }
            if (processed_transcripts.find(tid) != processed_transcripts.end()) {
                output_file_stream << buffer << "\n";
                continue;
            }
            else {
                fusion_tpm[tid] = tpm;
            }
        }

        logi("Skipped {} transcripts out of {}", skipped, total_count);
        logi("Creating fusion transcripts");
        std::ofstream gtfo_file{args["gtfo"].as<string>()};
        for (auto &[event, loci] : relevant_molecules) {
            for (auto &l : loci) {
                for (auto &t : l.second) {
                    const auto &it = fusion_tpm.find(t.info.at("transcript_id"));
                    if (it != fusion_tpm.end()) {
                        t.set_abundance(it->second);
                    }
                }
            }
            auto transcript_vec = event.execute_event(loci, rand_gen);
            logi("Creating {} fusion transcripts on event {}", transcript_vec.size(), event);
            for (auto &t : transcript_vec) {
                output_file_stream << t.to_abundance_str() << "\n";
                gtfo_file << t << "\n";
            }
        }
        return 0;
    }

    void describe_program() {
        logi("Running Fusion module");
        logi("Abundance file: {}", args["abundance"].as<string>());
        logi("Output file: {}", args["output"].as<string>());
        logi("GTF Input file: {}", args["gtfi"].as<string>());
        logi("GTF Output file: {}", args["gtfo"].as<string>());
        if (args["fusion-file"].count() > 0) {
            logi("Fusion file: {}", args["fusion-file"].as<string>());
        }
        logi("Fusion count: {}", args["fusion-count"].as<int>());
        if (args["disable-deletions"].as<bool>()) {
            logi("Deletions are disabled");
        }
        logi("Translocation ratio: {}", args["translocation-ratio"].as<double>());

        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Fusion_module);
