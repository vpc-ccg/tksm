
#pragma once
#ifndef MDF_H
#define MDF_H
#include <fstream>
#include <string>
#include <vector>
#include <utility>

#include "generator.h"
#include "interval.h"

using std::vector;
using std::string;
using std::ostream;
using std::ifstream;



inline Generator<molecule_descriptor> stream_mdf(ifstream &ist, bool unroll = false){
    string buffer;
    buffer.reserve(1000);
    std::getline(ist, buffer);
    while(ist){
        auto fields = rsplit(buffer, "\t");
        string id {fields[0].substr(1)};
        int depth{stoi(fields[1])};
//        int exon_count{stoi(fields[2])};
        string comment{fields[2]};
        vector<ginterval> segments;
        vector<std::pair<int, char>> errors_so_far;
        size_t size_so_far = 0;
        std::getline(ist, buffer);
        while(ist && buffer[0] != '+'){
        //        for(int i = 0; i< exon_count; ++i){

            auto fields = rsplit(buffer, "\t");
            string chr{fields[0]};
            int start {stoi(fields[1])};
            int end   {stoi(fields[2])};
            string strand{fields[3]};
            if(fields.size() > 4){
                auto mutations = rsplit(fields[4], ",");
                for(string mutation : mutations){
                    if(mutation == ""){
                        continue;
                    }
                    char target = mutation[mutation.size()-1];
                    int mutation_pos = stoi(mutation.substr(0,mutation.size()-1));
                    errors_so_far.push_back(std::make_pair( size_so_far + mutation_pos, target));
                }
            }
            size_so_far += (end-start);
            segments.emplace_back(chr, start, end, strand);
            std::getline(ist, buffer);
        }
        molecule_descriptor md{id, !segments[0].plus_strand};
        md.depth(depth)->update_errors(errors_so_far)->assign_segments(segments)->comment(comment);
        if (unroll){
            molecule_descriptor mdc = md;
            mdc.depth(1);
            for(int i =0;i<md.get_depth();++i){
                mdc.id(md.get_id() + "_" + std::to_string(i));
                co_yield mdc;
            }
        }
        else{
            co_yield md;
        }
    }
}
inline vector<molecule_descriptor> parse_mdf(ifstream &ist, bool unroll = false){
    auto streamer = stream_mdf(ist, unroll);
    vector<molecule_descriptor> mdfs;
    while(streamer){
        mdfs.push_back(streamer());
    }
    return mdfs;
}
/*
inline vector<molecule_descriptor> parse_mdf(ifstream &ist){
    vector<molecule_descriptor> molecules;
    string buffer;
    buffer.reserve(1000);
    std::getline(ist, buffer);
    while(ist){
        auto fields = rsplit(buffer, "\t");
        string id {fields[0].substr(1)};
        int depth{stoi(fields[1])};
//        int exon_count{stoi(fields[2])};
        string comment{fields[2]};
        vector<ginterval> segments;
        vector<std::pair<int, char>> errors_so_far;
        size_t size_so_far = 0;
        std::getline(ist, buffer);
        while(ist && buffer[0] != '+'){
        //        for(int i = 0; i< exon_count; ++i){

            auto fields = rsplit(buffer, "\t");
            string chr{fields[0]};
            int start {stoi(fields[1])};
            int end   {stoi(fields[2])};
            string strand{fields[3]};
            if(fields.size() > 4){
                auto mutations = rsplit(fields[4], ",");
                for(string mutation : mutations){
                    if(mutation == ""){
                        continue;
                    }
                    char target = mutation[mutation.size()-1];
                    int mutation_pos = stoi(mutation.substr(0,mutation.size()-1));
                    errors_so_far.push_back(std::make_pair( size_so_far + mutation_pos, target));
                }
            }
            size_so_far += (end-start);
            segments.emplace_back(chr, start, end, strand);
            std::getline(ist, buffer);
        }
//        molecules.emplace_back(id, segments, errors_so_far, depth);
        molecules.emplace_back(id, !segments[0].plus_strand).depth(depth)->update_errors(errors_so_far)->assign_segments(segments)->comment(comment);
    }
    return molecules;
}
*/
inline void print_mdf(ostream &ost, const pcr_copy& molecule){
    //For now use the depth of first pcr_copy
    print_tsv(ost, "+"+molecule.id, molecule.depth/*depth*/, molecule.comment);

    const vector<std::pair<int, char>> &errors = molecule.errors_so_far;
    int size_so_far = 0;
    for( const ginterval &ival : molecule.segments){
        string error_str = "";


        for(std::pair<int, char> error : errors){
            if( error.first > size_so_far && error.first < size_so_far + ival.end - ival.start){
                error_str += (std::to_string(error.first-size_so_far) + error.second + ",");
            }
        }
        if(error_str != ""){
            error_str.pop_back();
        }

        print_tsv(ost, ival.chr, ival.start, ival.end, (ival.plus_strand?"+":"-"), error_str);
        size_so_far += (ival.end-ival.start);
    }
}

/*
inline void print_mdf(ostream &ost, const string &id, const pcr_molecule &molecule, const vector< vector<std::pair<int, char>>> &errors_per_segment){


    //For now use the depth of first pcr_copy
    print_tsv(ost, "+"+id, molecule.paired[0].depth, molecule.paired[0].comment);

    int interval_counter = 0;
    for( const pcr_copy &pcp : molecule.paired){
        int size_so_far = 0;
        for( const ginterval &ival : pcp.segments){
            string error_str = "";

            const vector<std::pair<int, char>> &errors = errors_per_segment[interval_counter];

            for(std::pair<int, char> error : errors){
                if( error.first > size_so_far && error.first < size_so_far + ival.end - ival.start){
                    error_str += (std::to_string(error.first-size_so_far) + error.second + ",");
                }
            }
            if(error_str != ""){
                error_str.pop_back();
            }

            print_tsv(ost, ival.chr, ival.start, ival.end, (ival.plus_strand?"+":"-"), error_str);
            size_so_far += (ival.end-ival.start);
            ++interval_counter;
        }
    }
}
*/
inline void print_all_mdf(ostream &ost, const vector<pcr_copy> &molecules){

    for( const pcr_copy &molecule : molecules){
        //For now use the depth of first pcr_copy
        print_tsv(ost, "+"+molecule.id, molecule.depth/*depth*/, molecule.comment);

        int size_so_far = 0;
        for( const ginterval &ival : molecule.segments){
            string error_str = "";

            const vector<std::pair<int, char>> &errors = molecule.errors_so_far;

            for(std::pair<int, char> error : errors){
                if( error.first > size_so_far && error.first < size_so_far + ival.end - ival.start){
                    error_str += (std::to_string(error.first-size_so_far) + error.second + ",");
                }
            }
            if(error_str != ""){
                error_str.pop_back();
            }

            print_tsv(ost, ival.chr, ival.start, ival.end, (ival.plus_strand?"+":"-"), error_str);
            size_so_far += (ival.end-ival.start);

        }
    }
}

#endif

