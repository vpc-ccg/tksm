#ifndef MDF_H
#define MDF_H
#include <fstream>
#include <string>
#include <vector>
#include <utility>

#include "interval.h"

using std::vector;
using std::string;
using std::ostream;
using std::ifstream;

inline vector<pcr_copy> parse_mdf(ifstream &ist){

    vector<pcr_copy> molecules;
    string buffer;
    buffer.reserve(1000);
    while(std::getline(ist, buffer)){
        auto fields = rsplit(buffer, "\t");
        string id {fields[0].substr(1)};
        int depth{stoi(fields[1])};
        int exon_count{stoi(fields[2])};
        string comment{fields[3]};
        vector<ginterval> segments;
        vector<std::pair<int, char>> errors_so_far;
        size_t size_so_far = 0;
        for(int i = 0; i< exon_count; ++i){
            std::getline(ist, buffer);
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
        }
        molecules.emplace_back(id,segments, errors_so_far, depth);
    }
    return molecules;
}

inline void print_mdf(ostream &ost, const string &id, const pcr_molecule &molecule, const vector< vector<std::pair<int, char>>> &errors_per_segment){

    size_t interval_count = 0;//molecule.paired.size();
    for(const pcr_copy &pcp : molecule.paired){
        interval_count += pcp.segments.size();
    }
    //For now use the depth of first pcr_copy
    print_tsv(ost, "+"+id, molecule.paired[0].depth/*depth*/, interval_count, "comment");

    int interval_counter = 0;
    for( const pcr_copy &pcp : molecule.paired){
        size_t size_so_far = 0;
        for( const ginterval &ival : pcp.segments){
            string error_str = "";

            const vector<std::pair<int, char>> &errors = errors_per_segment[interval_counter];

            for(std::pair<int, char> error : errors){
                error_str += (std::to_string(error.first-size_so_far) + error.second + ",");
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


#endif

