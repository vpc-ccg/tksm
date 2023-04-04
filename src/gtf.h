#pragma once
#include <unordered_map>
#ifndef GTF_H
#define GTF_H

#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "graph.h"
#include "interval.h"
#include "util.h"
/*
inline graph<gene, double>
build_gene_graph(std::string path_to_gtf, const std::map<gene, int> &gene2count, bool coding_only = true) {
    graph<gene, double> gene_graph;

    std::ifstream file(path_to_gtf);

    std::string str;
    while (std::getline(file, str)) {
        if (str[0] == '#') {
            continue;
        }
        std::vector<std::string> fields = rsplit(str, "\t");

        if (fields[2] != "gene") {
            continue;
        }
        std::string chr{fields[0]};
        int s = stoi(std::string{fields[3]});
        int e = stoi(std::string{fields[4]});

        std::string st{fields[6]};

        std::string info_str{fields[8]};
        std::vector<std::string> info = rsplit(info_str, ";");
        strip_for_each(info, " ");

        std::string gene_id   = "";
        std::string gene_name = "";
        std::string biotype   = "";
        int count             = 0;
        for (auto iter = info.begin(); iter != info.end(); ++iter) {
            if ((*iter).size() <= 1) {
                continue;
            }
            std::string f{*iter};

            std::vector<std::string> fs = rsplit(f, " ");
            strip_for_each(fs, "\"");

            if (fs[0] == "gene_id") {
                gene_id = fs[1];
                ++count;
            }
            if (fs[0] == "gene_name") {
                gene_name = fs[1];
                ++count;
            }
            if (fs[0] == "gene_biotype") {
                biotype = fs[1];

                ++count;
            }
        }
        if (coding_only && biotype != "protein_coding") {
            continue;
        }
        gene g{chr, s, e, st, gene_id, gene_name};
        if (gene2count.find(g) != gene2count.end()) {
            gene_graph.add(g);
        }
    }
    return gene_graph;
}

inline auto
read_gtf_exons(std::string path_to_gtf, bool coding_only = true) {
    std::vector<exon> annots;
    std::ifstream file(path_to_gtf);
    std::string str;

    std::vector<gene *> gptrs;
    gene *current_gene;

    std::map<std::string, gene> t2g;

    while (std::getline(file, str)) {
        if (str[0] == '#') {
            continue;
        }
        std::vector<std::string> fields = rsplit(str, "\t");
        std::string type                = fields[2];
        std::string chr{fields[0]};
        int s = stoi(std::string{fields[3]}) - 1;  // GTF is 1 based
        int e = stoi(std::string{fields[4]});      // 1 based but inclusive so we add 1 to make it exclusive - 1 + 1
        std::string st{fields[6]};

        if (type == "gene") {
            std::string info_str{fields[8]};
            std::vector<std::string> info = rsplit(info_str, ";");
            strip_for_each(info, " \n\t");

            std::string gene_id   = "";
            std::string gene_name = "";
            std::string biotype   = "";
            int count             = 0;
            for (auto iter = info.begin(); iter != info.end(); ++iter) {
                if ((*iter).size() <= 1) {
                    continue;
                }
                std::string f{*iter};
                std::vector<std::string> fs = rsplit(f, " ");
                strip_for_each(fs, "\"");
                if (fs[0] == "gene_id") {
                    gene_id = fs[1];
                    ++count;
                }
                if (fs[0] == "gene_name") {
                    gene_name = fs[1];
                    ++count;
                }
                if (fs[0] == "gene_biotype") {
                    biotype = fs[1];
                    ++count;
                }
            }
            if (coding_only && (biotype != "protein_coding")) {
                continue;
            }
            gptrs.push_back(new gene(chr, s, e, st, gene_id, gene_name));
            current_gene = gptrs.back();
        }
        else if (type == "transcript") {
            std::string info_str{fields[8]};
            std::vector<std::string> info = rsplit(info_str, ";");
            strip_for_each(info, " \n\t");

            std::string transcript_id = "";
            std::string biotype       = "";
            int count                 = 0;
            for (auto iter = info.begin(); iter != info.end(); ++iter) {
                if ((*iter).size() <= 1) {
                    continue;
                }
                std::string f{*iter};
                std::vector<std::string> fs = rsplit(f, " ");
                strip_for_each(fs, "\"");
                if (fs[0] == "transcript_id") {
                    transcript_id = fs[1];
                    ++count;
                }
                if (fs[0] == "transcript_biotype") {
                    biotype = fs[1];
                    ++count;
                }
            }
            if (coding_only && (biotype != "protein_coding")) {
                continue;
            }
            // current_transcript = std::make_shared<transcript>(chr,s,e,st,transcript_id,current_gene);
            t2g[transcript_id.substr(0, 15)] = *current_gene;
        }
        else if (type == "exon") {
            std::string info_str{fields[8]};

            std::vector<std::string> info = rsplit(info_str, ";");
            strip_for_each(info, " \n\t");

            std::string biotype  = "";
            std::string tbiotype = "";

            std::string exon_id       = "";
            std::string transcript_id = "";
            int count                 = 0;
            for (auto iter = info.begin(); iter != info.end(); ++iter) {
                if ((*iter).size() <= 1) {
                    continue;
                }
                std::string f{*iter};
                std::vector<std::string> fs = rsplit(f, " ");
                strip_for_each(fs, "\"");

                if (fs[0] == "exon_id") {
                    exon_id = fs[1];
                    ++count;
                }

                if (fs[0] == "transcript_id") {
                    transcript_id = fs[1];
                    ++count;
                }
                if (fs[0] == "gene_biotype") {
                    biotype = fs[1];
                    ++count;
                }
                if (fs[0] == "transcript_biotype") {
                    tbiotype = fs[1];
                    ++count;
                }
            }
            if (coding_only && ((biotype != "protein_coding") || (tbiotype != "protein_coding"))) {
                continue;
            }
            annots.emplace_back(chr, s, e, st, exon_id, transcript_id, current_gene);
        }
    }

    return make_tuple(annots, gptrs, t2g);
}
*/

inline auto
read_gtf(const std::string &path2gtf) -> std::vector<gtf> {
    std::ifstream gtfile(path2gtf);
    std::string buffer;
    std::vector<gtf> gtf_entries;

    while (std::getline(gtfile, buffer)) {
        if (buffer[0] == '#') {
            continue;
        }
        gtf_entries.push_back(gtf{buffer});
    }
    return gtf_entries;
}

inline void
fillnames(gtf &entry) {
    auto g_name_iter = entry.info.find("gene_name");
    if (g_name_iter == entry.info.end()) {
        entry.info["gene_name"] = entry.info["gene_id"];
    }
    if (entry.type == gtf::entry_type::gene) {
        return;
    }
    auto t_name_iter = entry.info.find("transcript_name");
    if (t_name_iter == entry.info.end()) {
        entry.info["transcript_name"] = entry.info["transcript_id"];
    }
}

inline auto
read_gtf_genes(const std::string &path2gtf, bool fill_names = true, bool skip_lnc = true)
    -> std::vector<std::pair<gtf, vector<gtf>>> {
    std::ifstream gtfile(path2gtf);
    std::string buffer;
    std::vector<std::pair<gtf, vector<gtf>>> genes;

    while (std::getline(gtfile, buffer)) {
        if (buffer[0] == '#') {
            continue;
        }
        gtf entry{buffer};
        if (skip_lnc && entry.info["gene_biotype"] != "protein_coding") {
            continue;
        }
        if (fill_names) {
            fillnames(entry);
        }
        if (entry.type == gtf::entry_type::gene) {
            genes.emplace_back(entry, vector<gtf>{});
        }
        else if (entry.type == gtf::entry_type::transcript) {
            genes.back().second.push_back(entry);
        }
    }
    return genes;
}

inline auto
read_gtf_transcripts_deep(const std::string &path2gtf, bool skip_lnc = true, bool fill_names = true)
    -> std::unordered_map<std::string, transcript> {
    std::ifstream gtfile(path2gtf);
    std::string buffer;
    std::unordered_map<string, transcript> transcripts;
    std::string current_transcript = "";
    while (std::getline(gtfile, buffer)) {
        if (buffer[0] == '#') {
            continue;
        }
        gtf entry{buffer};
        if (entry.info["gene_biotype"] != "protein_coding" && skip_lnc) {
            continue;
        }
        if (fill_names) {
            fillnames(entry);
        }
        //        if( entry.info["transcript_biotype"] != "protein_coding" && skip_lnc){
        //            continue;
        //        }
        if (entry.type == gtf::entry_type::transcript) {
            transcripts.emplace(entry.info["transcript_id"], entry);
            current_transcript = entry.info["transcript_id"];
        }
        if (entry.type == gtf::entry_type::exon) {
            transcripts.at(current_transcript).add_exon(entry);
        }
    }
    return transcripts;
}

inline auto
read_gtf_transcripts(const std::string &path2gtf, int default_depth = 1) {
    std::map<std::string, molecule_descriptor> isoforms;
    std::ifstream gtfile(path2gtf);
    std::string buffer;

    gtf *current_gene       = nullptr;
    gtf *current_transcript = nullptr;

    while (std::getline(gtfile, buffer)) {
        if (buffer[0] == '#') {
            continue;
        }
        gtf *entry = new gtf(buffer);

        switch (entry->type) {
            case gtf::entry_type::gene:
                delete current_gene;
                current_gene = entry;
                break;
            case gtf::entry_type::transcript: {
                delete current_transcript;
                current_transcript = entry;
                auto iter          = isoforms.emplace(
                    current_transcript->info["transcript_id"],
                    molecule_descriptor{current_transcript->info["transcript_id"], !entry->plus_strand});
                iter.first->second.depth(default_depth);
                break;
            }
            case gtf::entry_type::exon:
                isoforms[current_transcript->info["transcript_id"]].append_segment(*entry);
                delete entry;
                break;
            default:
                break;
        }
    }

    return isoforms;
}

#endif
