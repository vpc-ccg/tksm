#include <exception>
#include <string>
#include <random>
#include <vector>
#include <map>
#include <fstream>
#include <set>

#include <climits>

#include "extern/IITree.h"
#include "extern/cxxopts.hpp"
#include "interval.h"
#include "tree.h"
#include "gtf.h"
#include "mdf.h"
#define FMT_HEADER_ONLY 1
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fmt/format.h>

using std::ofstream;
using std::set;
using std::vector;
using std::string;
using std::map;


double variant_ratio(const string &genotype){
    size_t pos = 0;
    double sum = 0;
    int cc = 0;
    while( pos < genotype.size()){
        size_t p = genotype.find_first_of("|/", pos );
        if(p==std::string::npos){
            break;
        }
        string f = genotype.substr(pos,p-pos + 1);
        if(f == "."){
            f = "0";
        }
        int variant = stoi(f);
        sum += ((variant == 0) ? 0 : 1);
        cc+=1;
        pos = p + 1;
    }
    string f = genotype.substr(pos);
    if(f == "."){
        f = "0";
    }
    int variant = stoi(f);
    sum += ((variant == 0) ? 0 : 1);
    cc+=1; 
    return sum / cc;
}

std::mt19937 rand_gen{std::random_device{}()};

vector<exon> fuse_isoforms(const vector<exon> &g1, const vector<exon> &g2, std::pair<int,int> bps){
    vector<exon> fusiso;
    fusiso.reserve(g1.size() + g2.size());
    auto cut_isoform = [ &fusiso] (const vector<exon> &iso, int b1, bool first) mutable -> void{
        auto iter1 = iso.begin();
        
        bool strand = iso.front().plus_strand;
        bool before_or_after = !strand ^ first;
        if(!before_or_after){
            while(iter1->end < b1 && iter1 != iso.end()){
                fusiso.push_back(*iter1);
                ++iter1;
            }
            //Create a break exon ending at breakpoint if breakpoint is on an exon
            if(iter1 != iso.end() && iter1->start < b1){
                fusiso.push_back( exon{iter1->chr, iter1->start, b1, iter1->strand, iter1->exon_id + std::to_string(b1), iter1->transcript_id + std::to_string(b1), iter1->gene_ref});
            }
        }
        else{
            while(iter1 != iso.end() && iter1->end < b1){
                ++iter1;
            }
            //Create a break exon ending at breakpoint if breakpoint is on an exon
            if(iter1 != iso.end() && iter1->start < b1){
                fusiso.push_back( exon{iter1->chr, iter1->start, b1, iter1->strand, iter1->exon_id + std::to_string(b1), iter1->transcript_id + std::to_string(b1), iter1->gene_ref});
                ++iter1;
            }
            while(iter1 != iso.end()){
                fusiso.push_back(*iter1);
                ++iter1;
            }
        }
    };

    cut_isoform(g1, bps.first, true);
    cut_isoform(g2, bps.second, false);

    return fusiso;

}



auto generate_fusion_isoforms(
        const vector<std::pair<gene, gene>> &fusions, 
        const tree<exon, int> &exon_tree,
        int gene_index){

    //Returns tuple<
    map<string, isoform> isoforms;
    vector<intra_event> event_intervals;
    //        >

    map<gene, vector<exon>> gene2firstexons; //reverse_index
    map<gene, int>          gene2count;
    map<gene, std::pair<int, int>> gene2largestrange;
    
    auto find_largest_range = [] (const tree<exon, int> &t) -> std::pair<int,int>{
        int start = INT_MAX;
        int end = 0;
        auto finder = [&start, &end] (int depth, const tree<exon, int> *current) mutable -> void {
            if(current->identity.start < start){
                start = current->identity.start;
            }
            if(current->identity.end > end){
                end = current->identity.end;
            }
        };
        t.df_execute2(finder);
        return std::make_pair(start, end);
    };

    for( const auto &pair : exon_tree.children){
        const exon &e = pair.first;
        gene2count[e.gene_ref]+=  (exon_tree.value(e));
        gene2firstexons[e.gene_ref].push_back(e);
        gene2largestrange[e.gene_ref] = find_largest_range(pair.second);
    }
    auto find_breakpoint = [] (const std::pair<int,int> &range) -> int{
        int start = range.first;
        int end = range.second;
        std::uniform_int_distribution<> dist(start, end);
        return dist(rand_gen);
    };

    map<gene, int> gene2bp;
    for( auto p : gene2largestrange){
        gene2bp[p.first] = find_breakpoint(p.second);
    }
/*
            .   .   .   .   .
    _______/                 \________
            .   .   .   .   .
    _______/          ______/
    .   .   .   .   .   .
    \______              \______

   */
    auto find_isoform = [&gene2count, &gene2firstexons, &exon_tree] (const gene &g) -> vector<exon>{
        vector<exon> isoform;
//        bool is_forward = g.plus_strand;
        const tree<exon, int> *current = nullptr;
        int count;
        try{
            count = gene2count.at(g);
        }
        catch( std::exception &e){
            std::cerr << g.gene_name << " not in gene2count!\n";
            count  = 1;
        }
        std::uniform_int_distribution<> dist(0, count);
        int sum = 0;
        int rand_val = dist(rand_gen);
        for(const exon &e : gene2firstexons[g]){
            const auto &p = exon_tree.children.at(e);
            int count = p.data;
            current = &p;
            sum += count;
            if(sum > rand_val){
                break;
            }
        }
        if(current != nullptr){
            while((current->children.size() > 0 ))
            {
                std::uniform_int_distribution<> dist(0, current->data);
                int rand_val = dist(rand_gen);


                int sum = 0;
                for(const auto &pc : current->children){
                    int cc = pc.second.data;
                    sum+=cc;
                    if(sum >= rand_val){
                        isoform.push_back(current->identity);
                        current = &pc.second;

                        break;
                    }   //Decide
                }
            }
            isoform.push_back(current->identity);
        }
        return isoform;
    };


    for( auto fusion : fusions){
        const gene &g1 = fusion.first;
        const gene &g2 = fusion.second;
        std::uniform_int_distribution<> dist(0, 1);
        int direction = dist(rand_gen); // Coin flip 
        int expression = gene2count[g1];
        if(direction > 0){
            expression = gene2count[g2];
        }
        if(g1.chr == g2.chr){
            char pm[2] = {'+', '-'};
            gene _g1 = g1;
            gene _g2 = g2;
            int start1 = gene2bp[g1];
            int end2 = gene2bp[g2];
            interval ix (start1, end2);
            if( start1 > end2){
                _g1 = g2;
                _g2 = g1;
                ix = interval{end2, start1};
            }
            event_intervals.emplace_back(g1.chr, ix.start, ix.end, fmt::format("{}{}",pm[_g1.plus_strand],pm[_g2.plus_strand]));
        }
        map<isoform, string> transcript_ids;

        int transcript_index = 0;
        for(int i =0; i < expression; ++i){
            auto iso1 = find_isoform(g1);
            auto iso2 = find_isoform(g2);
            
            isoform fusiso;
            if( g1.plus_strand && g2.plus_strand) { // ++ deletion
                fusiso = fuse_isoforms(iso1, iso2, std::make_pair(gene2bp[g1],gene2bp[g2]));
            }
            else if( (!g1.plus_strand) && (!g2.plus_strand)){ // -- deletion
                fusiso = fuse_isoforms(iso2, iso1, std::make_pair(gene2bp[g2],gene2bp[g1]));
            }
            else if( (g1.plus_strand) && (!g2.plus_strand)){ // +- inversion
                if(direction){
                    fusiso = fuse_isoforms(iso1, iso2, std::make_pair(gene2bp[g1],gene2bp[g2]));
                }
                else{
                    fusiso = fuse_isoforms(iso2, iso1, std::make_pair(gene2bp[g2],gene2bp[g1]));
                }
            }
            else if( (!g1.plus_strand) && (g2.plus_strand)){ // -+ inversion
                if(direction){
                    fusiso = fuse_isoforms(iso2, iso1, std::make_pair(gene2bp[g2],gene2bp[g1]));
                }
                else{
                    fusiso = fuse_isoforms(iso1, iso2, std::make_pair(gene2bp[g1],gene2bp[g2]));
                }
            }
            else{
                throw std::exception();
            }
            auto fusiter = transcript_ids.find(fusiso);
            if( fusiter == transcript_ids.end()){
                std::string new_id = fmt::format("FUST{:05}{:06}", gene_index, transcript_index);
                transcript_ids[fusiso] =new_id;
                isoforms[new_id] = fusiso;
                auto &ff = isoforms[new_id];
                ff.gene =  fmt::format("{}-{}", g1.gene_name,g2.gene_name);
                ff.transcript_id = new_id;
                ++transcript_index;
            }
            else{
                isoforms[transcript_ids[fusiso]].depth++;
            }

        }
        ++gene_index;
    }

    return std::make_tuple(isoforms, event_intervals);
}

map<string, vector<isoform>> generate_normal_isoforms(const tree<exon, int> &exon_tree){
    map<gene, vector<exon>> gene2firstexons; //reverse_index
    map<gene, int>          gene2count;
    for( const auto &pair : exon_tree.children){
        const exon &e = pair.first;
//        int count = pair.second.first;
//        const tree<exon, int> &child = pair.second.second;
        gene2count[e.gene_ref]+=  (exon_tree.value(e));
        gene2firstexons[e.gene_ref].push_back(e);
    }
    
    auto find_isoform = [&gene2count, &gene2firstexons, &exon_tree] (const gene &g) -> vector<exon>{
        vector<exon> isoform;
//        bool is_forward = g.plus_strand;
        const tree<exon, int> *current = nullptr;
        int count = gene2count.at(g);
        std::uniform_int_distribution<> dist(0, count);
        int sum = 0;
        int rand_val = dist(rand_gen);

        for(const exon &e : gene2firstexons[g]){
            const auto &p = exon_tree.children.at(e);
            int count = p.data;
            current = &p;
            sum += count;
            if(sum > rand_val){
                break;
            }
        }

        while((current->children.size() > 0 )){
            std::uniform_int_distribution<> dist(0, current->data);
            int rand_val = dist(rand_gen);


            int sum = 0;
            for(const auto &pc : current->children){
                int cc = pc.second.data;
                sum+=cc;
                if(sum >= rand_val){
                    isoform.push_back(current->identity);
                    current = &pc.second;

                    break;
                }   //Decide
            }
        }
        isoform.push_back(current->identity);

        return isoform;
    };
    map<string, vector<isoform>> isoforms;
    for( const auto &p : gene2count){
        const gene &g = p.first;
        int cnt = p.second;

        for(int i =0; i < cnt; ++i){
            auto iso = find_isoform(g);

            isoforms[g.gene_id].push_back(isoform{iso});
        }

    }
    return isoforms;
}

map<string, std::pair<int,int>> find_gene_counts_per_contig(const vector<gene> &gptrs){
    size_t index = 0;
    string prev_chr = "-1";

    map<string, int> starts_at;
    map<string, int> ends_at;
    for( gene p:gptrs){
        if ( p.chr != prev_chr){
            starts_at[p.chr] = index;
            ends_at[prev_chr] = index - 1;
        }
        prev_chr = p.chr;
        ++index;
    }
    ends_at[prev_chr] = index;
    map<string, std::pair<int, int>> chr_range;
    for( auto &p : starts_at){
        int end = ends_at[p.first];
        int start = p.second;
        chr_range.emplace(p.first, std::make_pair(start, end));
    }
    return chr_range;
}

vector< std::pair<gene,gene>> generate_random_fusions( const vector<gene> &gptrs, map<string, int> fusion_count_per_chrX2, int translocation_count){

    vector< std::pair<gene,gene>> fusions;
    vector<gene> translocation_targets;


    map<string, std::pair<int,int>> gene_ranges_per_contig = find_gene_counts_per_contig(gptrs);
    for( auto &p : gene_ranges_per_contig){
        auto it = fusion_count_per_chrX2.find(p.first);
        if( it == fusion_count_per_chrX2.end()){
            continue;
        }
        vector<size_t> values(p.second.second-p.second.first);
        std::iota(values.begin(), values.end(), p.second.first);
        std::shuffle(values.begin(),values.end(), rand_gen);

        std::vector<size_t> picked_values{values.begin(), values.begin() + 3 * it->second};
        

        if(picked_values.begin() == picked_values.end()){
            continue;
        }
        sort(picked_values.begin(),picked_values.end());

        int odd = 0;
        int cnt = 0;
        for( auto iter = picked_values.begin(); iter != picked_values.end() && (iter+1) != picked_values.end(); ++iter){

            if( odd == 0){

                fusions.push_back(std::make_pair(gptrs[*iter], gptrs[*(iter+1)]));
            }
            else if( odd == 2){
                gene g = gptrs[*iter];
                translocation_targets.push_back(g); 
            }
            ++odd;
            if(odd>2){
                odd = 0;
            }
            ++cnt;
        }  
    }

    std::shuffle(translocation_targets.begin(), translocation_targets.end(), rand_gen);
    auto current = translocation_targets.begin();
    auto iter = std::next(current);
    int cnt = 0;

    while(iter != translocation_targets.end() && cnt < translocation_count){
        while( iter != translocation_targets.end() && current->chr == iter->chr){
            ++iter;
        }
        if(iter->chr != current->chr){
            fusions.push_back(std::make_pair(*current, *iter));
            ++current;
            ++iter;
            ++cnt;
        }
    }
    return fusions;
}

vector<std::pair<gene,gene>> generate_random_fusions( const vector<gene> &gptrs, int count, int tloc_count){
        
    map<string, std::pair<int,int>> gene_ranges_per_contig = find_gene_counts_per_contig(gptrs);
   
    map<string, int> fusion_count_per_chrX2;
    double sum = 0;

    for( auto &p : gene_ranges_per_contig){
        int cnt = p.second.second - p.second.first;
        sum+=cnt;
    }

    double total=0;
    for( auto &p : gene_ranges_per_contig){
        double ratio = (p.second.second - p.second.first) / sum;
        fusion_count_per_chrX2[p.first] = int( count * ratio );
        total+= int(count*ratio);
    }


    while(total < count){
        for( auto &p : gene_ranges_per_contig){
        
            if (total >= count){
                break;
            }
            auto iter = fusion_count_per_chrX2.find(p.first);
            if(iter == fusion_count_per_chrX2.end()){
                continue;
            }
            ++(iter->second);
            ++total;
        }
    }
        

    return generate_random_fusions(gptrs, fusion_count_per_chrX2, tloc_count);
}

auto count_isoforms_from_mdf(
        vector<pcr_copy> &molecules,
        const map<string, gene> &t2g
        ){

    //Returns
    tree<exon, int> exon_transitions{};
    //

    for( const pcr_copy &pcp : molecules){
        int cc =1;

        const gene &g = t2g.at(pcp.id);

        vector<exon> iso;
        auto it = pcp.segments.begin();
        const auto &seg = *it;
        ++it;
        exon e(seg.chr, seg.start, seg.end, (seg.plus_strand?"+":"-"),fmt::format("{}-{}",pcp.id, 0),pcp.id,g);  
        tree<exon, int> *prev = &exon_transitions.try_get(e, 0);
        if(prev->parent != nullptr){
            prev->parent->value(e)++;
        }
        exon preve = e;
        for(;it!=pcp.segments.end();++it){
            const auto &seg = *it;
            exon e(seg.chr, seg.start, seg.end, (seg.plus_strand?"+":"-"),fmt::format("{}-{}",pcp.id, cc),pcp.id,g);  
            
            tree<exon, int> *tre = &(prev->try_get(e,0));
            prev->value(e)++;
            prev = tre;
            preve = e;
            ++cc;
        }

    }

    return exon_transitions;
}

auto count_reads_on_tree(
        vector<exon> &exon_ref,
        map<string,IITree<int, size_t>> &exon_forest,
        vector<mapping> &reads,
        int min_dist,
        int max_dist){

    //Returns
    tree<exon, int> exon_transitions{};
    //

    for( const mapping &r : reads){
        std::map<string, gene> segmap;
        if( r.segments.begin() == r.segments.end() ){
            continue;   
        }

        //Find the gene by votes of segments
        auto counts_of_genes = [&exon_forest, &exon_ref] (const mapping &q) -> map<string, int> {
            vector<set<string>> gene_set_per_segment;
            vector<size_t> overlaps;
            for( auto iter = q.segments.begin(); iter != q.segments.end(); ++iter){
                gene_set_per_segment.emplace_back();
                //Find overlapping exons
                const ginterval &inter = iter->tmplt;
                exon_forest[inter.chr].overlap(inter.start, inter.end,overlaps);
                if(overlaps.empty()){
                    continue;
                }
                for(size_t i : overlaps){
                    size_t index = exon_forest[inter.chr].data(i);
                    if( index < 0){
                        continue;
                    }
                    gene_set_per_segment.back().insert(exon_ref[index].gene_ref.gene_id);
                }
            }
            map<string, int> gene_counts;
            for( auto &s : gene_set_per_segment){
                for( auto &g : s){
                    gene_counts[g] += 1;
                }
            }
            return gene_counts;
        }(r);
        
        auto find_exon = [&counts_of_genes,&exon_ref,&exon_forest] (const ginterval &g) -> exon {
            auto &exon_tree = exon_forest[g.chr];
            vector<size_t> overlaps; 
            exon_tree.overlap(g.start,g.end,overlaps);
            
            if( overlaps.empty()){
                return exon{};
            }
            std::sort(overlaps.begin(), overlaps.end(), [&counts_of_genes, &exon_tree, &exon_ref](size_t a, size_t b){
                string ga = exon_ref[exon_tree.data(a)].gene_ref.gene_id;
                string gb = exon_ref[exon_tree.data(b)].gene_ref.gene_id;
                return counts_of_genes[ga] > counts_of_genes[gb];
            });
            return exon_ref[exon_tree.data(overlaps[0])];
        };

        auto first = r.segments.begin();
        exon e = find_exon(first->tmplt);
        if( e == exon{}){
            continue;
        }
        
        tree<exon, int> *prev = &exon_transitions.try_get(e, 0);
        if(prev->parent != nullptr){
            prev->parent->value(e)++;
        }

        exon preve = find_exon(first->tmplt);//{};
        for( auto siter = std::next(first); siter != r.segments.end(); ++siter){
            exon e = find_exon(siter->tmplt);
            if( e == exon{}){
                continue;
            }
            if(preve.reciprocal(e) > 0.75){
                continue;
            }
            tree<exon, int> *tre = &(prev->try_get(e,0));
            prev->value(e)++;
            prev = tre;
            preve = e;
        }
    }
    return exon_transitions;
}

map<string, IITree<int, size_t>> make_exon_interval_tree( const vector<exon> &annots){

    map<string,IITree<int, size_t>> itree;
    size_t index = 0;

    //for(const exon &annot: annots){
    for(auto iter = annots.begin(); iter != annots.end(); ++iter){
        auto &annot = *iter;
        itree[annot.chr].add(annot.start, annot.end, index);

        ++index;
    }

    for( auto &p : itree){
        p.second.index();
    }
    return itree;
}
vector<mapping> convert_aligs_to_segments(const string &path_to_aligs, int max_skip = 25, int min_segment = 100){
    vector<mapping> mappings;

    std::ifstream file( path_to_aligs);
    string str;

    //0.180525.0      1636    1192    1592    +       18      80373285        3254010 3256236 56      401     60      tp:A:P  mm:i:2  gn:i:343        go:i:3  cg:Z:17M1D3M341I4M2167N20M1I14M
    while( std::getline(file, str)){ 
        mappings.emplace_back(str, max_skip, min_segment);
    }
    return mappings;
}



/*
 *  Read bedlike file from path and simulate fusions ([x] means x is an optional field)
 *  chr1 start1 end1 chr2 start2 end2 gene1:gene2 genotype|genotype [abundance] [isoform1:isoform2]
 *  where we take the segments between start and end
 *  or
 *  chr1 pos1 chr2 pos2 gene1:gene2 genotype|genotype [abundance] [isoform1:isoform2]
 *  where we take the breakpoint as input and compute the segments from the gene positions and strands
 */
auto simulate_given_fusions( const string &path,
                                        const tree<exon, int> &exon_tree,
                                        const vector<gene *> &gptrs){
    
    //Returns tuple<
    map<string,isoform> fusion_isoforms;
    vector<intra_event> event_intervals;
    int gene_index = 0;
    // >

    string str;
    map<gene, int> gene2count;
    map<gene,vector<exon>> gene2firstexons;
   
    map<string, gene> name2gene;

    for( const auto &pair : exon_tree.children){
        const exon &e = pair.first;
        gene2count[e.gene_ref]+=  (exon_tree.value(e));
        gene2firstexons[e.gene_ref].push_back(e);
    


    }

    for( const gene *gpt : gptrs){
        
        string id = gpt->gene_id;
        string name = gpt->gene_name;

        //Adding both symbol and id of the gene to a table
        if( name2gene.find(id) == name2gene.end()){
            name2gene[id] = *gpt;
        }
        if( name2gene.find(name) == name2gene.end()){
            name2gene[name] = *gpt;
        }
    }

    auto find_isoform = [&gene2count, &gene2firstexons, &exon_tree] (const gene &g, const string &isoname) -> vector<exon>{

        vector<exon> isoform;
        if(isoname == "find"){
            const tree<exon, int> *current = nullptr;
            int count = gene2count[g];
            std::uniform_int_distribution<> dist(0, count);
            int sum = 0;
            int rand_val = dist(rand_gen);
            for(const exon &e : gene2firstexons[g]){
                const auto &p = exon_tree.children.at(e);
                int count = p.data;
                current = &p;
                sum += count;
                if(sum > rand_val){
                    break;
                }
            }

            while((current->children.size() > 0 ))
            {
                std::uniform_int_distribution<> dist(0, current->data);
                int rand_val = dist(rand_gen);
                int sum = 0;
                for(const auto &pc : current->children){
                    int cc = pc.second.data;
                    sum+=cc;
                    if(sum >= rand_val){
                        isoform.push_back(current->identity);
                        current = &pc.second;
                        break;
                    }   //Decide
                }
            }
            isoform.push_back(current->identity);
            std::sort(isoform.begin(), isoform.end());
            return isoform;
        }
        else{
            return isoform;
        }
    };


    auto isoform_overlaps_bp = [] (const vector<exon> &v, int bp) -> bool{
        interval i {v.front().start, v.back().end};
        return (bp > i.start && bp < i.end);
    };
//chr1 start chr2 end gene1:gene2 genotype|genotype [info]
// info: ab:[abundance] is:[isoform1:isoform2] id:[id]
    std::ifstream file(path);

    while(std::getline(file, str)){
        //chr   start   chr2 end    gene1   gene2   count   comment
        //Transcription starts from the first locus
        vector<string> fields { rsplit(str, "\t")};
        string chr1 {fields[0]};
        int start1  {stoi(fields[1])};
        string chr2 {fields[2]};
        int end2    {stoi(fields[3])};
        
        vector<string> genes { rsplit(fields[4],":")};
        string genotype {fields[5]};

        map<string, string> info;
        for(const string &s : rsplit(fields[6]," ")){
            string key = s.substr(0,2);
            string value = s.substr(3);
            info[key] = value;
        }

        gene g1 = name2gene.at(genes[0]);
        gene g2 = name2gene.at(genes[1]);

        if(g1.chr == g2.chr){
            char pm[2] = {'+', '-'};
            gene _g1 = g1;
            gene _g2 = g2;
            interval ix (start1, end2);
            if( start1 > end2){
                _g1 = g2;
                _g2 = g1;
                ix = interval{end2, start1};
            }
            event_intervals.emplace_back(chr1, ix.start, ix.end, fmt::format("{}{}",pm[_g1.plus_strand],pm[_g2.plus_strand]));
        }
        auto abiter = info.find("ab");
        int count;
        if( abiter != info.end()){
            count = stoi(abiter->second);
        }
        else{ //Calculate count from head gene
            count = gene2count.at(g1);
            count = count * variant_ratio(genotype);
        }

        auto isiter = info.find("is");
        string isoname1{"find"};
        string isoname2{"find"};
        if( isiter != info.end()){//set flag to skip find_isoform
            auto ii = rsplit(isiter->second, ":");
            isoname1 = ii[0];
            isoname2 = ii[1];
        }

        int transcript_index = 0;
        map<isoform, string> transcript_ids;
        for( int i = 0; i < count; ++i){
            auto iso1 = find_isoform(g1, isoname1);
            int try_c = 0;
            while(!isoform_overlaps_bp(iso1, start1)){
                iso1 = find_isoform(g1, isoname1);
                ++try_c;
                if(try_c > 1000){
                    std::cerr << "Cannot find good isoform for " << genes[0] << " for " << genes[0] << ":" << genes[1] << " fusion\n";
                    std::cerr << fmt::format("{}, {}-{}\n", start1, iso1.front().start, iso1.back().end);
                    break;
                }
            }
            auto iso2 = find_isoform(g2, isoname2);
            try_c = 0;
            while(!isoform_overlaps_bp(iso2, end2)){
                iso2 = find_isoform(g2, isoname2);
                ++try_c;
                if(try_c > 1000){
                    std::cerr << "Cannot find good isoform for " << genes[1] << " for " << genes[0] << ":" << genes[1] << " fusion\n";
                    std::cerr << fmt::format("{}, {}-{}\n", start1, iso1.front().start, iso1.back().end);
                    break;
                }
            }


            auto fusion = fuse_isoforms(iso1, iso2, std::make_pair(start1,end2));

            auto tid_iter = transcript_ids.find(fusion);

            if(tid_iter == transcript_ids.end()){
            
                transcript_ids[fusion] = fmt::format("FUST{:05}{:06}", gene_index, transcript_index);
                transcript_index++;
                tid_iter = transcript_ids.find(fusion);

                fusion_isoforms[tid_iter->second] = (fusion);
                auto &ff =fusion_isoforms[tid_iter->second];
                ff.gene =  fmt::format("{}-{}", g1.gene_name,g2.gene_name);
                ff.transcript_id = tid_iter->second;


            }
            else{
                auto &ff =fusion_isoforms[tid_iter->second];
                ff.depth++;
            }
        }
        ++gene_index;
    }
    return std::make_tuple(fusion_isoforms,event_intervals, gene_index);
}


int main(int argc, char **argv){

    cxxopts::Options options("RNAInfuser Fusion", "Fusion module of RNAInfuser");

    options.add_options()
        //("p,paf",  "Path to whole genome mappings in paf format", cxxopts::value<string>())
        ("g,gtf",  "Path to gtf annotation file", cxxopts::value<string>())
        ("d,dna",  "Path to human genome reference file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("m,input", "Path to input mdf", cxxopts::value<string>())
        ("fusion-count", "Number of gene fusions to simulate", cxxopts::value<int>()->default_value("0"))
        ("translocation-ratio", "Percentage of fusions to be simulated as translocations", cxxopts::value<double>()->default_value("0.1"))
        ("f,fusion-file", "Tab separated file to describe fusions", cxxopts::value<string>())
        ("disable-deletions", "Disables deletions (from fusions) to remove expression on overlapping genes", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))       
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory {{"input", "gtf", "output"}};

    int missing_parameters = 0;
    for( string &param : mandatory){
        if(args.count(param) == 0){
            std::cerr << param << " is required!\n";
            ++missing_parameters;
        }
    }
    if(missing_parameters  > 0){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    string path_to_gtf   {args["gtf"].as<string>()};
    string path_to_dna   {args["dna"].as<string>()};

    string out_mdf_path {args["output"].as<string>()};
    string input_mdf_path {args["m"].as<string>()};

    int seed = args["seed"].as<int>();;
    rand_gen.seed(seed);

 
    int fusion_count = args["fusion-count"].as<int>();
    double translocation_ratio = args["translocation-ratio"].as<double>();

    bool delete_genes = !args["disable-deletions"].as<bool>();

    std::cerr << "Reading GTF!\n";

    auto [gtf_exons, gene_ptrs, t2g] = read_gtf_exons(path_to_gtf,false);
    map<string, IITree<int, size_t>> exon_interval_tree = make_exon_interval_tree(gtf_exons);

    ifstream mf(input_mdf_path);
    vector<pcr_copy> molecules = parse_mdf(mf);
    tree<exon, int> exon_trie = count_isoforms_from_mdf(molecules,t2g);

    set<gene> expressed_genes;
    for(const auto &tr : exon_trie.children){
        expressed_genes.insert(tr.first.gene_ref);
    }

    std::cerr << "Generating fusions\n";

    vector<gene> exp_gene_vec;
    exp_gene_vec.insert(exp_gene_vec.begin(),expressed_genes.begin(),expressed_genes.end());
    
    std::sort(exp_gene_vec.begin(), exp_gene_vec.end() ,[](const gene &g1, const gene &g2) -> bool {
        return static_cast<ginterval>(g1) < g2; 
    });

    vector<isoform> fusion_isoforms;
    set<string> deleted_genes;
    int event_count_so_far = 0;
    if(args["f"].count() > 0){
        auto [given_fus, event_positions, _event_count] = simulate_given_fusions(args["f"].as<string>(), exon_trie, gene_ptrs);
        event_count_so_far = _event_count;
        for( const auto &ff : given_fus){
            fusion_isoforms.push_back(ff.second);
        }
        std::sort(event_positions.begin(), event_positions.end());
        
        auto gene_iter = exp_gene_vec.begin();
        std::vector<gene> valid_genes;
        for(const intra_event &gi : event_positions){
            while( gene_iter != exp_gene_vec.end() && gene_iter->chr < gi.chr){
                ++gene_iter;

                valid_genes.push_back(*gene_iter);
            }
            while(gene_iter != exp_gene_vec.end()  && gene_iter->start < gi.end){
                if(gene_iter->end < gi.start){
                    valid_genes.push_back(*gene_iter);
                }
                else if( gi.layout[0] == gi.layout[1]){ //Deletion
                    deleted_genes.insert(gene_iter->gene_id);
                }
                ++gene_iter;
            }
        }
        exp_gene_vec = valid_genes;
    }


    vector< std::pair<gene,gene>> fusions = generate_random_fusions(
            exp_gene_vec, fusion_count * ( 1- translocation_ratio), fusion_count * translocation_ratio);


    auto [gen_fusion, event_positions] = generate_fusion_isoforms(
            fusions, 
            exon_trie,
            event_count_so_far);


    auto gene_iter = exp_gene_vec.begin();
    for(const intra_event &gi : event_positions){
        while( gene_iter != exp_gene_vec.end() && gene_iter->chr < gi.chr){
            ++gene_iter;
        }
        while(gene_iter != exp_gene_vec.end()  && gene_iter->start < gi.end){
            if(gene_iter->end < gi.start){
            }
            else if( gi.layout[0] == gi.layout[1]){ //Deletion
                deleted_genes.insert(gene_iter->gene_id);
            }
            ++gene_iter;
        }
    }

    std::cerr << "Simulating fusions\n";

    ofstream outfile {out_mdf_path};

    for(const auto &pcp : molecules){
        if(delete_genes && deleted_genes.find(t2g.at(pcp.id).gene_id) != deleted_genes.end()){ // Gene is deleted skip
            continue;
        }
        print_mdf(outfile, pcp); 
    }

    for(const auto &iso : fusion_isoforms){
        print_mdf(outfile, iso.transcript_id, pcr_molecule{iso.transcript_id, iso}, vector<vector<std::pair<int, char>>>{iso.segments.size()});
    }
    for(const auto &iso : gen_fusion){
        print_mdf(outfile, iso.second.transcript_id, pcr_molecule{iso.second.transcript_id, iso.second}, vector<vector<std::pair<int, char>>>{iso.second.segments.size()});
    }

    return 0;
}
