#ifndef _FUSION_H_
#define _FUSION_H_
#include <cgranges/IITree.h>

#include <cxxopts.hpp>
#include <random>
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




class locus {
    string chr;
    int position;

public:
    locus(const string &chr, int position) : chr{chr}, position{position} {}

    bool operator<(const locus &other) const {
        if (chr == other.chr) {
            return position < other.position;
        }
        return chr < other.chr;
    }

    locus operator+(int offset) const { return locus{chr, position + offset}; }

    locus operator-(int offset) const { return locus{chr, position - offset}; }

    bool operator==(const locus &other) const { return chr == other.chr && position == other.position; }

    bool operator!=(const locus &other) const { return !(*this == other); }

    auto get_chr() const -> const string & { return chr; }

    auto get_position() const -> int { return position; }

    auto to_string() const -> string { return chr + ":" + std::to_string(position); }

    friend std::ostream &operator<<(std::ostream &os, const locus &loc) {
        os << loc.to_string();
        return os;
    }

    friend std::istream &operator>>(std::istream &is, locus &loc) {
        string chr;
        int position;
        is >> chr >> position;
        loc = locus{chr, position};
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
    EventType event_type;
    double event_ratio;
    ginterval supplementary_interval;  // To be used for translocations

    auto fuse_molecules(const molecule_descriptor &md1, const molecule_descriptor &md2) const -> molecule_descriptor {
        switch (event_type) {
            case EventType::DELETION:

                break;
            case EventType::INVERSION:

                break;
            case EventType::TRANSLOCATION:
            case EventType::DUPLICATION:
            case EventType::INSERTION:
            case EventType::NONE:
                throw std::runtime_error("Invalid event type " + event_type_to_string(event_type));
        }
        return {};
    }

public:
    // Event
    chimeric_event(EventType event_type) : ginterval{}, event_type{event_type}, event_ratio{0.5} {}
    chimeric_event(const string &chr, int start, int end, const string &orientation, EventType event_type)
        : ginterval{chr, start, end, orientation}, event_type{event_type}, event_ratio{0.5} {}

    auto execute_event(const map<locus, vector<molecule_descriptor>> &molecules, auto &rand_gen) const
        -> vector<molecule_descriptor> {
        std::unordered_set<locus> keys;
        for (auto &[key, md] : molecules) {
            keys.insert(key);
        }

        assert(keys.size() == 2);

        locus start = *keys.begin();
        locus end   = *(std::next(keys.begin()));

        // shuffle the molecules
        std::shuffle(molecules.at(start).begin(), molecules.at(start).end(), rand_gen);
        std::shuffle(molecules.at(end).begin(), molecules.at(end).end(), rand_gen);
        auto random_position_picker = std::uniform_int_distribution<size_t>(0, molecules.at(end).size() - 1);

        // Fix the order if needed considering the position and strands

        std::uniform_real_distribution<double> dist(0, 1);
        vector<molecule_descriptor> fused_molecules;
        set<molecule_descriptor> used_at_the_end;

        for (const auto &md : molecules.at(start)) {
            if (dist(rand_gen) < event_ratio) {
                // Fuse the md with a random molecule from the end
                auto &end_molecules = molecules.at(end);
                auto &random_md     = end_molecules[random_position_picker(rand_gen)];
                used_at_the_end.insert(random_md);

                // Create a new molecule
                molecule_descriptor fused_md = fuse_molecules(md, random_md);
                fused_molecules.push_back(fused_md);
            }
            else {
                fused_molecules.push_back(md);
            }
        }

        // Add the molecules that were not used at the end
        for (const auto &md : molecules.at(end)) {
            if (used_at_the_end.find(md) == used_at_the_end.end()) {
                fused_molecules.push_back(md);
            }
        }

        return fused_molecules;
    }

    locus get_start() const { return locus{chr, start}; }
    locus get_end() const { return locus{chr, end}; }
};

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

class Fusion_module : public tksm_module {
    cxxopts::ParseResult parse(int argc, char **argv) {
        // clang-format off
        options.add_options("main")
            (
                "i,input",
                "input mdf file",
                cxxopts::value<string>()
            )(
                "o,output",
                "output mdf file",
                cxxopts::value<string>()
            )(
                "g,gtf",
                "Path to GTF annotation file",
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
        unsigned long total_size = 0;
        for (const auto &gp : genes_by_chr) {
            auto &chr     = gp.first;
            auto &genes   = gp.second;
            chr_size[chr] = genes.back().end - genes.front().start;
            total_size += chr_size[chr];
        }

        std::map<string, int> fusion_count_per_chr;
        for (const auto &fp : fusions_so_far) {
            fusion_count_per_chr[fp.chr] =
                std::round(static_cast<double>(fusion_count) * chr_size[fp.chr] / total_size);
        }

        // Convert fusion vector to a map of chr -> vector<chimeric_event>

        std::map<string, vector<chimeric_event>> fusions_by_chr;
        for (auto &fusion : fusions_so_far) {
            fusions_by_chr[fusion.chr].push_back(fusion);
        }

        for (const auto &fp : fusions_by_chr) {
            auto &chr        = fp.first;
            auto &fusions    = fp.second;
            auto &genes_copy = genes_by_chr[chr];
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
            for (int i = 0; i < fusion_count_per_chr[chr]; i += 2) {
                auto &gene1 = genes_copy[i];
                auto &gene2 = genes_copy[i + 1];

                // if genes are on the same strand generate deletion else inversion

                EventType event_type =
                    gene1.plus_strand != gene2.plus_strand ? EventType::INVERSION : EventType::DELETION;

                chimeric_event fusion{gene1.chr, gene1.end, gene2.start, gene1.plus_strand ? "+" : "-", event_type};
                fusions_so_far.push_back(fusion);
            }
        }
        return fusions_so_far;
    }

public:
    Fusion_module(int argc, char **argv) : tksm_module{"fusion", "Fusion module"}, args(parse(argc, argv)) {}

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output", "gtf"};
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

        string mdf_file_path{args["input"].as<string>()};
        string output_file_path{args["output"].as<string>()};
        string gtf_file_path{args["gtf"].as<string>()};

        std::ifstream mdf_file{mdf_file_path};
        auto streamer = stream_mdf(mdf_file);

        vector<gtf> transcripts = read_gtf(gtf_file_path);

        int fusion_count                        = args["fusion-count"].as<int>();
        double translocation_ratio              = args["translocation-ratio"].as<double>();
        [[maybe_unused]] bool disable_deletions = args["disable-deletions"].as<bool>();

        vector<chimeric_event> fusions;
        if (args["fusion-file"].count() > 0) {
            string fusion_file_path{args["fusion-file"].as<string>()};
            std::ifstream fusion_file{fusion_file_path};
            fusions = read_fusions(fusion_file);
        }
        generate_fusions(fusions, fusion_count, transcripts, translocation_ratio);

        // Build an interval tree of the fusion events
        IITree<locus, chimeric_event> fusion_tree;
        for (auto &fusion : fusions) {
            fusion_tree.add(fusion.get_start(), fusion.get_start() + 1, fusion);
            fusion_tree.add(fusion.get_end(), fusion.get_end() + 1, fusion);
        }

        fusion_tree.index();

        std::ofstream output_file{output_file_path};

        std::uniform_real_distribution<double> dist(0.0, 1.0);

        map<chimeric_event, map<locus, vector<molecule_descriptor>>> relevant_molecules;
        // Find relavent molecule descriptions save them in a vector and print the rest
        while (streamer) {
            auto md = streamer();
            // Combine the segments of the molecule description to a single ginterval

            ginterval g = md.cget_segments().front();
            for (const ginterval &gi : md.cget_segments()) {
                // Check if the ginterval overlaps with any of the fusion events
                if (g.start > gi.start) {
                    g.start = gi.start;
                }
                if (g.end < gi.end) {
                    g.end = gi.end;
                }
            }

            // Check if the ginterval overlaps with any of the fusion events

            locus start{g.chr, g.start};
            locus end{g.chr, g.end};
            vector<size_t> overlapping_positions;
            fusion_tree.overlap(start, end, overlapping_positions);

            if (overlapping_positions.empty()) {
                output_file << md;
            }
            else {
                const auto &target_event = fusion_tree.data(overlapping_positions[0]);
                const auto &start        = fusion_tree.start(overlapping_positions[0]);
                relevant_molecules[target_event][start].push_back(md);
            }
        }

        // Generate md for the fusion transcripts

        for (auto &[fusion, md_map] : relevant_molecules) {
            // Generate md for the fusion transcripts
            auto fused_molecules = fusion.execute_event(md_map, rand_gen);  // Also contains the unfused molecules
            for (auto &md : fused_molecules) {
                output_file << md;
            }
        }

        return 0;
    }

    void describe_program() {
        logi("Running Fusion module");
        logi("Input file: {}", args["input"].as<string>());
        logi("Output file: {}", args["output"].as<string>());
        logi("GTF file: {}", args["gtf"].as<string>());
        logi("Fusion file: {}", args["fusion-file"].as<string>());
        logi("Fusion count: {}", args["fusion-count"].as<int>());
        if (args["disable-deletions"].as<bool>()) {
            logi("Deletions are disabled");
        }
        logi("Translocation ratio: {}", args["translocation-ratio"].as<double>());

        fmtlog::poll(true);
    }
};

#endif
