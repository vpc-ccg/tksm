
#pragma once
#ifndef MDF_H
#define MDF_H
#include <cstdio>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "generator.h"
#include "interval.h"
#include "util.h"

using std::istream;
using std::ostream;
using std::string;
using std::vector;

/*
inline Generator<molecule_descriptor>
stream_mdf(FILE *file, bool unroll = false){
    char id[1024];
    int depth;
    char comment[10000];
    while(std::scanf("%s\t%d\t%s", id, &depth, comment) == 3){
        char chr[1024];
        int start;
        int end;
        char strand;
        char mutations[10000];
        molecule_descriptor md{id, strand == '-'};
        md.depth(depth)->comment(comment);

        vector<std::pair<int, char>> errors_so_far;
        while(std::scanf("%s\t%d\t%d\t%c\t%s", chr, &start, &end, &strand, mutations) == 5){

            char *token = strtok(mutations, ",");
            while(token != NULL){
                int mutation_pos;
                char target;
                sscanf(token, "%d%c", &mutation_pos, &target);
                errors_so_far.push_back(std::make_pair(mutation_pos, target));
                token = strtok(NULL, ",");
            }

            md.append_segment(ginterval{chr, start, end, std::to_string(strand)});

        }
        md.update_errors(errors_so_far);
        if(unroll){
            molecule_descriptor mdc = md;
            mdc.depth(1);
            for(int i = 0; i < md.get_depth(); ++i){
                mdc.id(md.get_id() + "_" + std::to_string(i));
                co_yield mdc;
            }
        }
        else{
            co_yield md;
        }
    }
}
*/
inline generator<molecule_descriptor>
stream_mdf(std::istream &ist, bool unroll = false) {
    string buffer;
    buffer.reserve(1000);
    std::getline(ist, buffer);
    while (ist) {
        auto fields = rsplit(buffer, "\t");
        string id{fields[0].substr(1)};
        int depth{stoi(fields[1])};
        //        int exon_count{stoi(fields[2])};
        string comment{fields[2]};
        vector<einterval> segments;
        vector<std::pair<int, char>> errors_so_far;

        std::getline(ist, buffer);
        while (ist && buffer[0] != '+') {
            //        for(int i = 0; i< exon_count; ++i){

            auto fields = rsplit(buffer, "\t");
            string chr{fields[0]};
            int start{stoi(fields[1])};
            int end{stoi(fields[2])};
            string strand{fields[3]};
            string error_str;
            if (fields.size() > 4) {
                error_str = fields[4];
            }

            segments.emplace_back(chr, start, end, strand).parse_and_add_errors(error_str);
            std::getline(ist, buffer);
        }
        molecule_descriptor md{id, !segments[0].plus_strand};
        md.depth(depth)->assign_segments(segments)->comment(comment);
        if (unroll && md.get_depth() > 1) {
            molecule_descriptor mdc = md;
            mdc.depth(1);
            for (int i = 0; i < md.get_depth(); ++i) {
                molecule_descriptor mdcx = mdc;
                mdcx.id(md.get_id() + "_" + std::to_string(i));
                co_yield mdcx;
            }
        }
        else {
            co_yield md;
        }
    }
}

inline generator<molecule_descriptor>
stream_mdf(const string &filename, bool unroll = false) {
    std::basic_ifstream<char> ist{filename};
    if (!ist) {
        throw std::runtime_error("Could not open file " + filename);
    }

    string buffer;
    buffer.reserve(1000);
    std::getline(ist, buffer);
    while (ist) {
        auto fields = rsplit(buffer, "\t");
        string id{fields[0].substr(1)};
        int depth{stoi(fields[1])};
        //        int exon_count{stoi(fields[2])};
        string comment{fields[2]};
        vector<einterval> segments;
        vector<std::pair<int, char>> errors_so_far;

        std::getline(ist, buffer);
        while (ist && buffer[0] != '+') {
            //        for(int i = 0; i< exon_count; ++i){

            auto fields = rsplit(buffer, "\t");
            string chr{fields[0]};
            int start{stoi(fields[1])};
            int end{stoi(fields[2])};
            string strand{fields[3]};
            string error_str;
            if (fields.size() > 4) {
                error_str = fields[4];
            }

            segments.emplace_back(chr, start, end, strand).parse_and_add_errors(error_str);
            std::getline(ist, buffer);
        }
        molecule_descriptor md{id, !segments[0].plus_strand};
        md.depth(depth)->assign_segments(segments)->comment(comment);
        if (unroll && md.get_depth() > 1) {
            molecule_descriptor mdc = md;
            mdc.depth(1);
            for (int i = 0; i < md.get_depth(); ++i) {
                molecule_descriptor mdcx = mdc;
                mdcx.id(md.get_id() + "_" + std::to_string(i));
                co_yield mdcx;
            }
        }
        else {
            co_yield md;
        }
    }
}

inline vector<molecule_descriptor>
parse_mdf(std::istream &ist, bool unroll = false) {
    vector<molecule_descriptor> mdfs;
    for (const molecule_descriptor &md : stream_mdf(ist, unroll)) {
        mdfs.push_back(md);
    }
    return mdfs;
}

#endif
