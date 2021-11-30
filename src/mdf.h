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
void print_mdf(ostream &ost, const string &id, const pcr_molecule &molecule, const vector< vector<std::pair<int, char>>> &errors_per_segment){

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
            size_t error_cnt = 0;

            const vector<std::pair<int, char>> &errors = errors_per_segment[interval_counter];
            error_cnt = errors.size();
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

