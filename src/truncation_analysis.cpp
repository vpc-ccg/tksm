
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "util.h"
#include "interval.h"
#include "gtf.h"

using std::map;
using std::string;
using std::ifstream;
using std::vector;
using std::stoi;
/*
 * Picks up the single isoform gene mappings to calculate rate of truncation
 * 
 */
/*
 * 1       havana  gene    11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
1       havana  transcript      11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "lncRNA"; tag "basic"; transcript_support_level "1";
1       havana  exon    11869   12227   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE00002234944"; exon_version "1"; tag "basic"; transcript_support_level "1";
1       havana  exon    12613   12721   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE00003582793"; exon_version "1"; tag "basic"; transcript_support_level "1";
1       havana  exon    13221   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";
*/
int main(int argc, char **argv){

    string  gtf_file_path       {argv[1]};
    string  cdna_paf_file_path  {argv[2]};

    map<string, isoform>   single_isoform_genes;
    string buffer;

    ifstream gtf_file {gtf_file_path};
    gene *current_gene = nullptr;

    vector<gene *> gptrs;
    vector<transcript> transcripts_so_far;
    vector<exon> exons_so_far;
    transcript current_transcript;
    exon current_exon;
    //Assumes gtf is properly ordered gene->transcript->exons

    while(std::getline(gtf_file, buffer)){
        if(buffer[0] == '#'){
            continue;
        }
        gtf entry{buffer};

        switch(entry.type){
            case gtf::entry_type::gene:
                if(transcripts_so_far.size() == 1){ //Meaning single isoform gene
                    gptrs.push_back(current_gene);
                    single_isoform_genes.emplace(current_gene->gene_id, isoform{exons_so_far});
                }
                else{
                    delete current_gene;
                }
                current_gene = new gene{entry};
                transcripts_so_far.clear();
                exons_so_far.clear();
                break;
            case gtf::entry_type::transcript:
                exons_so_far.clear();
                current_transcript = transcript{entry, current_gene};
                transcripts_so_far.push_back(current_transcript);
                break;
            case gtf::entry_type::exon:
                current_exon = exon{entry, current_gene};
                exons_so_far.push_back(current_exon);
                break;
            default: ;// Skip others
        }
    }

    //ENSG00000000419_0       644     0       370     +       ENSG00000000419.12_ENST00000371582.8    1161    791     1161    370     370     0       NM:i:0  ms:i:370        AS:i:370        nn:i:0  tp:A:P  cm:i:117
    //       s1:i:365        s2:i:365        de:f:0  rl:i:0  cg:Z:370M
    int largest = 0;
    std::vector<int> start_position_counts(50000,0);
    std::vector<int> end_position_counts(50000,0);
    ifstream paf_file {cdna_paf_file_path};
    while(std::getline(paf_file, buffer)){
        vector<string> fields = rsplit(buffer, "\t");
        mapping m{buffer, 10, 50};
        if( !m.primary){ //Skip supplementary aligs
            continue;
        }
        if( m.segments.empty()){//Empty alignment check just in case
            std::cerr << "Cannot find any alignment on the following paf line\n";
            std::cerr << buffer << "\n";
            continue;
        }
        string mapping_transcript_id = (m.segments.begin()->tmplt.chr); //chr here refers to contig id, I might change this

        string gene_id = rsplit(mapping_transcript_id, ".")[0];
        
        if( single_isoform_genes.find(gene_id) == single_isoform_genes.end()) { //gene is not single isoform
            continue;
        }

        start_position_counts[stoi(fields[7])] +=1;
        end_position_counts[stoi(fields[8])] +=1;
        largest = std::max(largest,stoi(fields[8]));
    }
    for(int i =0; i< largest; ++i){
        std::cout << i << "\t" << start_position_counts[i] << "\t" << end_position_counts[i] << "\n";
    }
    //Cleanup
    for( gene *g :gptrs){
        delete g;
    }
    return 0;
}
