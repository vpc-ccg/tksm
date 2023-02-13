#pragma once
#ifndef GTF_H
#define GTF_H

#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "graph.h"
#include "interval.h"
#include "util.h"

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
