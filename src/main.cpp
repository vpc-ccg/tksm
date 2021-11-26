/*
 *  Left to do:
 *      Remove unused graph code
 *      Reimplement readthrough with tree
 *      Refactor 
 *      Breakpoint on exons case
 */

//Precompiled headers
#include "libs.h"

//Included again for clarity, include guard should take care of this
#include <filesystem>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <fstream>
#include <utility>
#include <sstream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <array>

#include <cctype>
#include <cassert>
#include <climits>

#include <zlib.h>
#include "IITree.h"
#include "kekseq.h"
#include "cigar.h"
#include "interval.h"
#include "graph.h"
#include "tree.h"
#include "reverse_complement.h"
#include "util.h"
#include "gtf.h"

#include "extern/cxxopts.hpp"

using std::string;
using std::vector;
using std::cout;
using std::map;
using std::ostream;
using std::set;

#ifdef DEBUG
    using std::cerr;
    std::cerr << "DEBUG MODE!\n";
#else

class mockstream{
    template<class K>
    friend mockstream& operator << (mockstream &mock, const K &k){
        return mock;
    }
};
mockstream cerr;
#endif

//Random number generator, seed is set in the main function
std::mt19937 rand_gen{std::random_device{}()};

map<string, IITree<int, size_t>> make_exon_interval_tree( const vector<exon> &annots){

    map<string,IITree<int, size_t>> itree;
    size_t index = 0;

    //for(const exon &annot: annots){
    for(auto iter = annots.begin(); iter != annots.end(); ++iter){
        auto &annot = *iter;
        itree[annot.chr].add(annot.start, annot.end, index);

        ++index;
    }
    int cnt = 0;
    for( auto &p : itree){
        p.second.index();
        cout << ++cnt <<"\n";
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

void convert_counts_to_probabilities( graph<gene, double> &gene_graph){

    for( std::pair<gene, double>& g : gene_graph){
        auto nei = gene_graph.neighbours(g.first);
        int sum = g.second;
        for( const gene &n : nei){
            sum+= gene_graph.arc(g.first, n);
        }
        for( const gene &n : nei){
            double val = gene_graph.arc(g.first, n);
            gene_graph.arc(g.first, n) = val / sum;
        }
    }
}

auto count_reads_on_tree(
        vector<exon> &exon_ref,
        map<string,IITree<int, size_t>> &exon_forest,
        vector<mapping> &reads,
        int min_dist,
        int max_dist){

    tree<exon, int> exon_transitions{};
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
                    gene_set_per_segment.back().insert(exon_ref[index].gene_ref->gene_id);
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
                string ga = exon_ref[exon_tree.data(a)].gene_ref->gene_id;
                string gb = exon_ref[exon_tree.data(b)].gene_ref->gene_id;
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

//Rewrite
void count_reads_on_graph(
        graph<gene, double> &gene_graph,
        vector<exon> &exon_ref,
        map<string,IITree<int, size_t>> &exon_forest,
        vector<mapping> &reads,
        int min_dist,
        int max_dist){

    tree<exon, int> exon_transitions{};
    for( const mapping &r : reads){
        std::map<string, gene> segmap;
        if( r.segments.begin() == r.segments.end() ){
            continue;   
        }

        //Find the gene by votes of segments
        auto counts_of_genes = [&exon_forest, &exon_ref] (const mapping &q){
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
                    gene_set_per_segment.back().insert(exon_ref[index].gene_ref->gene_id);
                }
            }
            map<string, int> gene_counts;
            for( auto &s : gene_set_per_segment){
                for( auto &g : s){
                    gene_counts[g] += 1;
                }
            }
            return gene_counts;
            //vector<std::pair<string,int>> ret{gene_counts.begin(),gene_counts.end()};
            //std::sort(ret.begin(),ret.end(),[] (auto &a, auto &b) -> bool {
            //    return a.second < b.second;
            //});
            //return ret;
        }(r);
        
        auto find_exon = [&counts_of_genes,&exon_ref,&exon_forest] (const ginterval &g) -> exon{
            auto &exon_tree = exon_forest[g.chr];
            vector<size_t> overlaps; 
            exon_tree.overlap(g.start,g.end,overlaps);
            std::sort(overlaps.begin(), overlaps.end(), [&counts_of_genes, &exon_tree, &exon_ref](size_t a, size_t b){
                string ga = exon_ref[exon_tree.data(a)].gene_ref->gene_id;
                string gb = exon_ref[exon_tree.data(b)].gene_ref->gene_id;
                return counts_of_genes[ga] < counts_of_genes[gb];
            });
            return exon_ref[exon_tree.data(overlaps[0])];
        };

        auto first = r.segments.begin();
        
        tree<exon, int> *prev = &exon_transitions.try_get(find_exon(first->tmplt), 0);
        for( auto siter = std::next(first); siter != r.segments.end(); ++siter){
            exon e = find_exon(siter->tmplt);
            tree<exon, int> *tre = &(prev->try_get(e,0));
            prev->value(e)++;
            prev = tre;
        }
/*
        long prev_index = -1;
        bool gene_switch_flag = false;
        gene *last_gene;
        bool lg_set = false;
        for( auto siter = r.segments.begin(); std::next(siter) != r.segments.end(); ++siter){
            auto nex_it = std::next(siter);
            long index = l_first_overlap(nex_it->tmplt);

            exon &e2 = exon_ref[index];
            if(index > 0){
                lg_set = true;
                last_gene = e2.gene_ref;
            }
            if( index > 0 && prev_index > 0){
                exon &e1 = exon_ref[prev_index];
                int distance = e1.gene_ref->distance(*e2.gene_ref);
                if( *e1.gene_ref != *e2.gene_ref && e1.chr == e2.chr && e1.plus_strand == e2.plus_strand && distance > min_dist && distance < max_dist ){
                    gene_switch_flag = true;
                    ++gene_graph.arc(*e1.gene_ref, *e2.gene_ref);
                    //cout << r.rid << "\t" << siter->tmplt << "\t" << * exon_ref[prev_index].gene_ref <<
                    //    " -> " << nex_it->tmplt << "\t" <<  * exon_ref[index].gene_ref << "\t" << index << "\n";
                }

            }
            else if(index > 0 ){
                //                cout << r.rid << "\tNULL -> " << exon_ref[index] << "\n";
            }
            else if(prev_index > 0){
                //                cout << r.rid << "\t" << exon_ref[prev_index] << " -> NULL\n";
            }
            else{
                //                cout << r.rid << "\tNULL -> NULL\n";
            }
            prev_index = index;
        }
        //Single segment case
        if( std::next(r.segments.begin()) == r.segments.end()){
            long index = l_first_overlap(r.segments.begin()->tmplt);
            if(index < 0 ){
                continue;
            }
            exon &e1 = exon_ref[index];
            last_gene = e1.gene_ref;
            lg_set=true;
        }
        if( lg_set && !gene_switch_flag && (r.segments.begin() != r.segments.end()) && gene_graph.in(*last_gene)){
            gene_graph.value(*last_gene) +=1;
        }
    
    */
    }
}

map<string, std::pair<int,int>> find_gene_counts_per_contig(graph<gene, double> &gene_graph){

    size_t index = 0;
    string prev_chr = "-1";

    map<string, int> starts_at;
    map<string, int> ends_at;
    for( std::pair<gene,double> &p: gene_graph.nodes){
        if ( p.first.chr != prev_chr){
            starts_at[p.first.chr] = index;
            ends_at[prev_chr] = index - 1;
        }
        prev_chr = p.first.chr;
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

vector< std::pair<gene,gene>> generate_random_fusions( graph<gene, double> &gene_graph, map<string, int> fusion_count_per_chrX2, int translocation_count){
    map<string, std::pair<int,int>> gene_ranges_per_contig = find_gene_counts_per_contig(gene_graph);
    vector< std::pair<gene,gene>> fusions;
    vector<gene> translocation_targets;
    for( auto &p : gene_ranges_per_contig){
        auto iter = fusion_count_per_chrX2.find(p.first);
        if( iter == fusion_count_per_chrX2.end()){
            continue;
        }
        vector<size_t> values(p.second.second-p.second.first);
        std::iota(values.begin(), values.end(), p.second.first);
        std::shuffle(values.begin(),values.end(), rand_gen);
        std::vector<size_t> picked_values{values.begin(), values.begin() + 2 * iter->second};
        
        if(picked_values.begin() == picked_values.end()){
            continue;
        }
        sort(picked_values.begin(),picked_values.end());

        int odd = 0;
        for( auto iter = picked_values.begin(); std::next(iter) != picked_values.end(); ++iter){
            if( odd == 0){
                //cout << gene_graph.nodes[*iter].first << "\t" << gene_graph.nodes[*std::next(iter)].first << "\n";
                fusions.push_back(std::make_pair(gene_graph.nodes[*iter].first, gene_graph.nodes[*std::next(iter)].first));
            }
            else if( odd == 2){
                gene g = gene_graph.nodes[*iter].first;
                translocation_targets.push_back(g); 
            }
            ++odd;
            if(odd>2){
                odd = 0;
            }
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

vector<std::pair<gene,gene>> generate_random_fusions( graph<gene, double> &gene_graph, int count, int tloc_count){
        
    map<string, std::pair<int,int>> gene_ranges_per_contig = find_gene_counts_per_contig(gene_graph);
   
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
        

    return generate_random_fusions(gene_graph, fusion_count_per_chrX2, tloc_count);
}


template< class B>
void print_tsv(B b){
    cout << b << "\n";
}

template <class B, class... A>
void print_tsv(B b, A... a){
    cout << b << "\t";
    print_tsv(a...);
}



template <size_t L>
class awful_fasta{
    const kekseq::kstring &ks;
    static constexpr size_t line_length = L;
    public:
    awful_fasta(kekseq::kstring &seq) : ks(seq){
        
    }
    friend ostream &operator<<(ostream &ost, const awful_fasta<L> &f){
        
        for(int i =0; i < f.ks.l; i+=f.line_length){
            for(int j=0; j< std::min(f.line_length,f.ks.l-i);++j){
                ost << f.ks.s[i+j];
            }
            ost << "\n";
        }
        return ost;
    }
};


enum class bpstrategy{
    uniform,
    exon_biased,
};

vector<std::pair<int,int>> generate_fusion_breakpoints( const vector<std::pair<gene,gene>> &gene_pairs, bpstrategy strat){

    vector<std::pair<int,int>> breakpoints;
    switch (strat){
        case bpstrategy::uniform:
            for(auto p : gene_pairs){
                std::uniform_int_distribution<> dist1(p.first.start, p.first.end);
                std::uniform_int_distribution<> dist2(p.second.start, p.second.end);
                breakpoints.push_back(std::make_pair(dist1(rand_gen),dist2(rand_gen)));
            }
            break;
        case bpstrategy::exon_biased:
            cerr << "Fusion strategy not implemented!\n";
            exit(-1);
            break;
    }

    return breakpoints;
}

map<gene, int> generate_transcript_expression( const tree<exon, int> &exon_tree){

    map<gene,int> expressions;
    for( auto p: exon_tree.children){
        expressions[*p.first.gene_ref]+= p.second.data;
    }
    return expressions;

}

vector<int>  generate_fusion_expression( const vector<std::pair<gene,gene>> &fusions, const tree<exon,int> &exon_tree){
    vector<int> expressions;
    for( auto p: fusions){
        expressions.push_back(10);
    }
    return expressions;
}

vector<exon> fuse_isoforms(const vector<exon> &g1, const vector<exon> &g2, std::pair<int,int> bps){
    vector<exon> fusiso;
    fusiso.reserve(g1.size() + g2.size());
    auto cut_isoform = [ &fusiso] (const vector<exon> &g, int b1) mutable -> void{
        auto iter1 = g.begin();

        if(g.front().plus_strand){
            while(iter1->end < b1 && iter1 != g.end()){
                fusiso.push_back(*iter1);
                ++iter1;
            }
            //Create a break exon ending at breakpoint if breakpoint is on an exon
            if(iter1 != g.end() && iter1->start < b1){
                fusiso.push_back( exon{iter1->chr, iter1->start, b1, iter1->strand, iter1->exon_id + std::to_string(b1), iter1->gene_ref});
            }
        }
        else{
            while(iter1 != g.end() && iter1->end < b1){
                ++iter1;
            }
            //Create a break exon ending at breakpoint if breakpoint is on an exon
            if(iter1 != g.end() && iter1->start < b1){
                fusiso.push_back( exon{iter1->chr, iter1->start, b1, iter1->strand, iter1->exon_id + std::to_string(b1), iter1->gene_ref});
                ++iter1;
            }
            while(iter1 != g.end()){
                fusiso.push_back(*iter1);
                ++iter1;
            }
        }
    };

    cut_isoform(g1, bps.first);
    cut_isoform(g2, bps.second);

    return fusiso;

}

map<string, vector<isoform>> generate_normal_isoforms(const tree<exon, int> &exon_tree){
    map<gene, vector<exon>> gene2firstexons; //reverse_index
    map<gene, int>          gene2count;
    for( const auto &pair : exon_tree.children){
        const exon &e = pair.first;
//        int count = pair.second.first;
//        const tree<exon, int> &child = pair.second.second;
        gene2count[*e.gene_ref]+=  (exon_tree.value(e));
        gene2firstexons[*e.gene_ref].push_back(e);
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
            bool flag = false;

            int sum = 0;
            for(const auto &pc : current->children){
                int cc = pc.second.data;
                sum+=cc;
                if(sum >= rand_val){
                    isoform.push_back(current->identity);
                    current = &pc.second;
                    flag = true;
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

map<string, vector<isoform>> generate_fusion_isoforms(
        const vector<std::pair<gene, gene>> &fusions, 
        const vector<int> &fusion_expression,
        const tree<exon, int> &exon_tree){

    map<string, vector<isoform>> isoforms;

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
        gene2count[*e.gene_ref]+=  (exon_tree.value(e));
        gene2firstexons[*e.gene_ref].push_back(e);
        gene2largestrange[*e.gene_ref] = find_largest_range(pair.second);
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

    cerr << "There are " << gene2firstexons.size()  << " genes!\n";
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

        while((current->children.size() > 0 ))
        {
            std::uniform_int_distribution<> dist(0, current->data);
            int rand_val = dist(rand_gen);
            bool flag = false;

            int sum = 0;
            for(const auto &pc : current->children){
                int cc = pc.second.data;
                sum+=cc;
                if(sum >= rand_val){
                    isoform.push_back(current->identity);
                    current = &pc.second;
                    flag = true;
                    break;
                }   //Decide
            }
        }
        isoform.push_back(current->identity);

        return isoform;
    };

    auto expiter = fusion_expression.begin();
    for( auto fusion : fusions){
        const gene &g1 = fusion.first;
        const gene &g2 = fusion.second;
        std::uniform_int_distribution<> dist(0, 1);
        int direction = dist(rand_gen); // Coin flip 
        for(int i =0; i < *expiter; ++i){
            auto iso1 = find_isoform(g1);
            auto iso2 = find_isoform(g2);

            if( g1.plus_strand && g2.plus_strand) { // ++ deletion
                auto fusiso = fuse_isoforms(iso1, iso2, std::make_pair(gene2bp[g1],gene2bp[g2]));
                isoforms[g1.gene_id + "_" + g2.gene_id].push_back(isoform{fusiso});
            }
            else if( (!g1.plus_strand) && (!g2.plus_strand)){ // -- deletion
                auto fusiso = fuse_isoforms(iso2, iso1, std::make_pair(gene2bp[g2],gene2bp[g1]));
                isoforms[g2.gene_id + "_" + g1.gene_id].push_back(isoform{fusiso});
            }
            else if( (g1.plus_strand) && (!g2.plus_strand)){ // +- inversion
                if(direction){
                    auto fusiso = fuse_isoforms(iso1, iso2, std::make_pair(gene2bp[g1],gene2bp[g2]));
                    isoforms[g1.gene_id + "_" + g2.gene_id].push_back(isoform{fusiso});
                }
                else{
                    auto fusiso = fuse_isoforms(iso2, iso1, std::make_pair(gene2bp[g2],gene2bp[g1]));
                    isoforms[g2.gene_id + "_" + g1.gene_id].push_back(isoform{fusiso});
                }
            }
            else if( (!g1.plus_strand) && (g2.plus_strand)){ // -+ inversion
                if(direction){
                    auto fusiso = fuse_isoforms(iso2, iso1, std::make_pair(gene2bp[g2],gene2bp[g1]));
                    isoforms[g2.gene_id + "_" + g1.gene_id].push_back(isoform{fusiso});
                }
                else{
                    auto fusiso = fuse_isoforms(iso1, iso2, std::make_pair(gene2bp[g1],gene2bp[g2]));
                    isoforms[g1.gene_id + "_" + g2.gene_id].push_back(isoform{fusiso});
                }
            }
        }
        ++expiter;
    }

    return isoforms;
}


template<typename K>
constexpr int roundup32(K x){
    --x;
    x|=x>>1;
    x|=x>>2;
    x|=x>>4;
    x|=x>>8;
    x|=x>>16;
    ++x;
    return x;
}

map<string, string> read_fasta_fast( const string &fasta_path){
    string fai_path = fasta_path + ".fai";
    map <string, string> contig2seq;

    string buf;
    if(std::filesystem::exists(fai_path)){ //This helps saves little time, but keeping it for reference
        std::cerr << "Index exists... Allocating memory ahead of time!\n";
        std::ifstream f(fai_path);

        while(std::getline(f, buf)){
            std::istringstream str(buf);
            string contig;
            int base_count;
            str >> contig >> base_count;
            contig2seq[contig].reserve( roundup32(base_count));
        }
    }
    
    std::cerr << "Reading fasta file!\n";
    string contig;
    std::ifstream f(fasta_path);
    while (std::getline(f, buf)){
        if( buf[0] == '>'){ //new contig
            std::istringstream str(buf.c_str()+1);
            str >> contig;
        }
        else{
            contig2seq[contig] += buf;   
        }
    }
    return contig2seq;
}

/*
 *  Read bedlike file from path and simulate fusions
 */
vector<isoform> simulate_given_fusions( const string &path,
                                        const tree<exon, int> &exon_tree){
    
    vector<isoform> fusion_isoforms;

    std::ifstream file(path, std::ios::in);

    string str;
    map<gene, int> gene2count;
    map<gene,vector<exon>> gene2firstexons;
    
    map<string, gene> name2gene;

    for( const auto &pair : exon_tree.children){
        const exon &e = pair.first;
        gene2count[*e.gene_ref]+=  (exon_tree.value(e));
        gene2firstexons[*e.gene_ref].push_back(e);
        
        if(e.gene_ref != nullptr){
            string id = e.gene_ref->gene_id;
            string name = e.gene_ref->gene_name;

            //Adding both symbol and id of the gene to a table
            if( name2gene.find(id) == name2gene.end()){
                name2gene[id] = *e.gene_ref;
            }
            if( name2gene.find(name) == name2gene.end()){
                name2gene[name] = *e.gene_ref;
            }
        }
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

        while((current->children.size() > 0 ))
        {
            std::uniform_int_distribution<> dist(0, current->data);
            int rand_val = dist(rand_gen);
            bool flag = false;

            int sum = 0;
            for(const auto &pc : current->children){
                int cc = pc.second.data;
                sum+=cc;
                if(sum >= rand_val){
                    isoform.push_back(current->identity);
                    current = &pc.second;
                    flag = true;
                    break;
                }   //Decide
            }
        }
        isoform.push_back(current->identity);

        return isoform;
    };

    while(std::getline(file, str)){
        //chr   start   chr2 end    gene1   gene2   count   comment
        //Transcription starts from the first locus
        vector<string> fields { rsplit(str, "\t")};
        string chr1 {fields[0]};
        string chr2 {fields[2]};

        int start  {stoi(fields[1])};
        int end    {stoi(fields[3])};

        string gs1  {fields[4]};
        string gs2  {fields[5]};
        
        int count {stoi(fields[6])};
        string comment {fields[7]};
        
        gene g1 = name2gene.at(gs1);
        gene g2 = name2gene.at(gs2);

        for( int i = 0; i < count; ++i){
            auto iso1 = find_isoform(g1);
            auto iso2 = find_isoform(g2);

            auto fusion = fuse_isoforms(iso1, iso2, std::make_pair(start,end));
            
            fusion_isoforms.push_back(fusion);
        }
    }
    return fusion_isoforms;
}

void generate_and_print_fasta(
        const map<string, vector<isoform>> &isoforms,
        const map<string,string> &sequences,
        std::ofstream &ost){

    for( const auto &pair : isoforms){
        string base = pair.first;
        int cnt=0;
        for( const isoform &i : pair.second){
            ost << ">" << base << "_" << cnt << "\n";
            for(const exon &e : i.segments){
                if( sequences.find(e.chr) == sequences.end()){
                    cerr << e << "!\n";
                    throw std::runtime_error("Chr not in fasta");
                }
                if( e.end > sequences.at(e.chr).size()){
                    cerr << e << "!\n";
                    throw std::runtime_error("Exon is outside of the chr");
                }
                string seq = sequences.at(e.chr).substr( e.start, e.end - e.start); // Maybe string view in the future?
                if(! e.plus_strand){
                    reverse_complement::complement_inplace(seq);
                }
#ifdef DEBUG
                ost << e << "\t";
#endif
                ost << seq << "\n";
            }
            ++cnt;
        }

    }
}

//pcr copy structure that tracks pcr errors introduced
struct pcr_copy{
    isoform iso;
    vector< std::pair< int, char>> errors_so_far;
    bool reversed;

    pcr_copy( const isoform &iso, auto errors_so_far) :iso(iso), errors_so_far(errors_so_far) {}
    pcr_copy( const isoform &iso) :iso(iso) {}
};

//pcr molecule structure that can model paired molecules
struct pcr_molecule{
    vector<pcr_copy> paired;

    pcr_molecule() {}

    pcr_molecule(const pcr_molecule &other) : paired(other.paired) {}
    pcr_molecule(const isoform &other) {
        paired.emplace_back(other);
    }
    pcr_molecule(const pcr_molecule &first, const pcr_molecule &second) : paired(first.paired) {
        paired.reserve(first.paired.size() + second.paired.size());
        paired.insert(paired.end(), second.paired.begin(), second.paired.end());
    }
};

void generate_and_print_fasta_with_pcr( 
        const map<string, vector<isoform>> &isoforms,
        const map<string, string> &sequences,
        std::ofstream &ost,
        int cycles,
        double pcr_duplication_rate,
        double error_rate,
        int umi_length,
        int number_of_target_reads,
        double random_pairing_rate_per_cycle){

    vector< vector< pcr_molecule>> saved_for_rp(cycles);

    using mut_tree = tree<int, vector<std::pair< int, char >>>;

//Calculate number of molecules to be produced with pcr
    int64_t molecule_count = std::accumulate(isoforms.begin(), isoforms.end(), 0, [] ( int64_t sum, const std::pair<string, vector<isoform>> &pair) -> int64_t{
            return sum + pair.second.size();
            });

//Calculate expected number of molecules after pcr
    int64_t expected_number_after_pcr = std::pow((1 + pcr_duplication_rate), cycles) * molecule_count;

    double drop_ratio = static_cast<double>(number_of_target_reads) / expected_number_after_pcr;
    
//Convert error_rate to random_error_rate where characters can mutate into same chars (A->A or T->T etc.)
    error_rate = error_rate * 4.0 / 3.0; 

    std::uniform_real_distribution<> z1dist(0,1);
    std::uniform_int_distribution<> basedist(0,3);
    char bases[] = {'A','C','T','G'};

    //PCR lambda function
    std::function<void(const pcr_molecule&, mut_tree &,int, int, vector<int>&)> do_pcr = 
        [&do_pcr, &z1dist, pcr_duplication_rate, &bases, &basedist, random_pairing_rate_per_cycle, &saved_for_rp, cycles]
        (const pcr_molecule &iso, mut_tree &mt, int current_cycle, int expected_mutation_count, vector<int> &positions) mutable -> void{
        for(int cyc = current_cycle + 1; cyc <= cycles; ++cyc){
            if(z1dist(rand_gen) < pcr_duplication_rate){ // Molecule is caught in this cycle of pcr
                                                         // GC bias can be introduced here
                std::shuffle(positions.begin(),positions.end(), rand_gen);
                vector<std::pair<int, char>> mutations;
                for(int j = 0; j < expected_mutation_count; ++j){ //TODO: We can make this a normal distribution    
                    mutations.push_back(std::make_pair(positions[j], bases[basedist(rand_gen)]));
                }
                if(z1dist(rand_gen) < random_pairing_rate_per_cycle){
                    vector<std::pair<int,char>> mutations_so_far;
                    mut_tree *current = &mt;
                    while(current !=nullptr){
                        for( auto pr : current->data){
                            mutations_so_far.push_back(pr);
                        }
                        current = current->parent;
                    }
                    pcr_molecule other(iso);
                    other.paired[0].errors_so_far = mutations_so_far;
                    saved_for_rp[current_cycle].push_back(other);
                }
                else{
                    mt.add_child(cyc, mutations);
                    do_pcr(iso, mt[cyc], cyc, expected_mutation_count, positions);
                }
            }
        }
    };

    //Makes a printer lambda function to be run on the tree
    auto make_pcr_printer = [&z1dist, drop_ratio] (ostream &ost, string *_base_id, int copy_number, const map<string, string> *sequences, const pcr_molecule *_pm){
        return  [&ost, &z1dist, drop_ratio, _base_id, copy_number, sequences, _pm] (int depth, const mut_tree *mt) mutable -> void{
            string base_id{*_base_id};
            const pcr_molecule &pm {*_pm};
            if( z1dist(rand_gen) < drop_ratio) { //Molecule is captured by sequencing
                ost << ">" << base_id << "_" << copy_number << " ";
                mut_tree *current = mt->parent;
                ost << "lineage:";
                ost << mt->identity;
                while( current != nullptr){ //Print pcr Lineage
                    ost << "_" << current->identity;
                    current = current->parent;
                }
                ost << "\n";

                for( const pcr_copy &pcp : pm.paired){
                    int seq_index = 0;
                    const isoform &iso = pcp.iso;
                    for(const exon &e : iso.segments){
                        if( sequences->find(e.chr) == sequences->end()){
                            cerr << e << "!\n";
                            throw std::runtime_error("Chr not in fasta");
                        }
                        if( e.end > sequences->at(e.chr).size()){
                            cerr << e << "!\n";
                            throw std::runtime_error("Exon is outside of the chr");
                        }
                        string seq = sequences->at(e.chr).substr( e.start, e.end - e.start); // Maybe string view in the future?
                       
                        for( const std::pair<int, char> &pr : pcp.errors_so_far){//Errors so far, only used for random_pairings
                            if(pr.first > seq_index && pr.first < seq_index + seq.size()){
                                seq[ pr.first - seq_index] = pr.second;
                            }
                        }
                        std::function<void(const mut_tree *, string&)> apply_mutation = [&apply_mutation, seq_index] ( const mut_tree *mt, string &seq) mutable -> void{
                            if(mt->parent != nullptr){
                                apply_mutation(mt->parent, seq);
                            }
                            for(const auto &pr: mt->data){
                                if(pr.first > seq_index && pr.first < seq_index + seq.size()){
                                    seq[ pr.first - seq_index] = pr.second;
                                }
                            }
                        };

                        apply_mutation(mt, seq);
                        if(! e.plus_strand){
                            reverse_complement::complement_inplace(seq);
                        }

#ifdef DEBUG
                        ost << e << "\t";
#endif
                        ost << seq << "\n";
                        seq_index += seq.size();
                    }
                }
            }
        };
    };
    for( const auto &pair : isoforms){
        string base_id = pair.first;
        int cnt=0;
        for( const isoform &iso : pair.second){
            int molecule_size = 0;
            for(const exon &e : iso.segments){
                molecule_size += (e.end - e.start);
            }
            pcr_molecule pcm{iso};
            mut_tree mutation_tree;
            vector<int> positions(molecule_size);
            std::iota(positions.begin(), positions.end(), 0);
            int expected_mutation_count = error_rate * molecule_size; 

            do_pcr(pcm, mutation_tree, 0, expected_mutation_count, positions);
            auto printer = make_pcr_printer(ost, &base_id, cnt, &sequences, &pcm);
            mutation_tree.df_execute2(printer);
            ++cnt;
        }

    }
    int cycle = 1;

    vector<pcr_molecule> random_pairings;
    for( auto vec_iter = saved_for_rp.begin(); vec_iter !=saved_for_rp.end(); ++vec_iter){
        vector<pcr_molecule> &pc_vec = *vec_iter; 
        //After shuffle, pair adjacent elements
        std::shuffle(pc_vec.begin(), pc_vec.end(), rand_gen);
        for(auto iter = pc_vec.begin(); iter + 1 != pc_vec.end() && iter !=pc_vec.end(); iter= std::next(iter,2)){
            auto next_iter = std::next(iter);
            pcr_molecule rp(*iter, *next_iter);
            int molecule_size = 0;
            for( const pcr_copy &cp :rp.paired){
                for(const exon &e : cp.iso.segments){
                    molecule_size += (e.end - e.start);
                }
            }
            mut_tree mutation_tree;
            vector<int> positions(molecule_size);
            std::iota(positions.begin(), positions.end(), 0);
            int expected_mutation_count = error_rate * molecule_size; 

            do_pcr(rp, mutation_tree, cycle, expected_mutation_count, positions);
            string base_id = "RP";
            for( const pcr_copy &pc : rp.paired){
                base_id += "_";
                base_id += pc.iso.segments[0].gene_ref->gene_id;
            }
            auto printer = make_pcr_printer(ost, &base_id, cycle, &sequences, &rp);
            mutation_tree.df_execute2(printer);
        }
        
        //Continue PCR
        ++cycle;
    }
}

int main(int argc, char **argv){
    cxxopts::Options options("RNAInfuser", "RNA simulator");

    options.add_options()
        ("p,paf",  "Path to whole genome mappings in paf format", cxxopts::value<string>())
        ("g,gtf",  "Path to gtf annotation file", cxxopts::value<string>())
        ("c,cdna", "Path to cdna reference file", cxxopts::value<string>())
        ("d,dna",  "Path to human genome reference file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())

        ("read-count", "Number of reads to simulate", cxxopts::value<int>()->default_value("100000"))
        ("fusion-count", "Number of gene fusions to simulate", cxxopts::value<int>()->default_value("0"))
        ("translocation-ratio", "Percentage of fusions to be simulated as translocations", cxxopts::value<double>()->default_value("0.1"))

        ("rt-minimum-distance", "Minimum distance between readthrough genes to simulate",cxxopts::value<int>()->default_value("1000"))
        ("rt-maximum-distance", "Maximum distance between readthrough genes to simulate",cxxopts::value<int>()->default_value("500000"))
    
        ("umi", "Insert umi sequences to the molecules before pcr. You can specify the length.", cxxopts::value<int>()->implicit_value("16"))
        ("pcr", "Number of pcr cycles to simulate", cxxopts::value<int>()->default_value("6"))
        ("pcr-efficiency", "Probability of a molecule being duplicated during each pcr cycle", cxxopts::value<double>()->default_value("0.75"))
        ("pcr-error-rate", "Probability of substition errors for each base during each pcr cycle", cxxopts::value<double>()->default_value("0.01"))
        ("pcr-random-pairing-rate", "Probability of molecules randomly attaching to each other during each pcr cycle", cxxopts::value<double>()->default_value("0.0001"))

        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("h,help", "Help screen")
    ;

    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory {{"paf", "gtf", "cdna", "dna", "output"}};

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


    //parameters will be moved to argument parser when done
    string path_to_aligs {args["paf"].as<string>()};
    string path_to_gtf   {args["gtf"].as<string>()};
    string path_to_cdna  {args["cdna"].as<string>()};
    string path_to_dna   {args["dna"].as<string>()};
    string out_cdna_path {args["out"].as<string>()};
    int seed = args["seed"].as<int>();;

    rand_gen.seed(seed);
    int min_dist = args["rt-minimum-distance"].as<int>();
    int max_dist = args["rt-maximum-distance"].as<int>();
    
    gzFile gzf = gzopen(path_to_cdna.c_str(), "r");
    kekseq::kseq<gzFile,gzread> fr(gzf); 

    int r;
    map<string, string> transcript2sequence;
    map<string, string> transcript2comment;
    while((r = fr.read()) >= 0){
        if(transcript2sequence.find(fr.name.s)!=transcript2sequence.end()){
            cerr << "Duplicate transcript " << fr.name.s << "\n";
        }
        transcript2sequence[fr.name.s] = string{fr.seq.s};
        transcript2comment[fr.name.s] = string{fr.comment.s};
    }
    gzclose(gzf);

    map<string, string> chr2contig = read_fasta_fast(path_to_dna);

    std::cerr << "Reading FASTQ!\n";
    vector<mapping> reads = convert_aligs_to_segments(path_to_aligs);   
    std::cerr << "Reading GTF!\n";
    auto [gtf_exons, gene_ptrs] = read_gtf_exons(path_to_gtf);

    std::cerr << "Flattening GTF!\n";
//    auto [merged_exons, merge_index] = flatten_gtf(gtf_exons,0.75); // vector<exon>, map<string, string>

//    std::cerr << "Making GTF interval tree!\n" << " exons: " << merged_exons.size() << "\n";
    map<string, IITree<int, size_t>> exon_interval_tree = make_exon_interval_tree(gtf_exons);

    std::cerr << "Counting reads on tree\n";
    tree<exon, int> ec = count_reads_on_tree(gtf_exons, exon_interval_tree, reads, min_dist,  max_dist);
    std::cerr << "Done.\n";
    
    map<gene, vector<exon>> gene2firstexons; //reverse_index
    map<gene, int>          gene2count;
    for( const auto &pair : ec.children){
        const exon &e = pair.first;
        gene2count[*e.gene_ref]+=  (ec.value(e));
        gene2firstexons[*e.gene_ref].push_back(e);
    }

    std::cerr << "Building graphs\n";
    graph<gene, double> gene_graph = build_gene_graph(path_to_gtf, gene2count); // I should replace this

    std::cerr << "Generating fusions\n";
    int fusion_count = args["fusion-count"].as<int>();
    double translocation_ratio = args["translocation-ratio"].as<double>();
    
    vector< std::pair<gene,gene>> fusions = generate_random_fusions(gene_graph, fusion_count * ( 1- translocation_ratio), fusion_count * translocation_ratio);
//    vector<std::pair<int,int>> fusion_breakpoints = generate_fusion_breakpoints( fusions, bpstrategy::uniform);

//    map<gene,int> transcript_expression = generate_transcript_expression(  ec);//Convert this to transcript 
    vector<int> fusion_expression = generate_fusion_expression( fusions, ec);

    map<string, vector<isoform>> normal_isoforms = generate_normal_isoforms(ec);

    std::cerr << "Simulating fusions\n";
    map<string, vector<isoform>> fusion_isoforms = generate_fusion_isoforms(fusions, fusion_expression, ec);
    
    std::cerr << "Generating molecule sequences\n";
    std::ofstream output_cdna_stream(out_cdna_path);

    int pcr_cycles = args["pcr"].as<int>();
    double pcr_efficiency = args["pcr-efficiency"].as<double>();
    double pcr_error_rate = args["pcr-error-rate"].as<double>();
    double pcr_random_pairing_rate = args["pcr-random-pairing-rate"].as<double>();
    int umi_length = args["umi"].as<int>();
    generate_and_print_fasta_with_pcr(normal_isoforms, chr2contig, output_cdna_stream, pcr_cycles, pcr_efficiency, pcr_error_rate, umi_length, 10000000, pcr_random_pairing_rate);
    generate_and_print_fasta(fusion_isoforms, chr2contig, output_cdna_stream);

    std::cerr << "Cleaning Up\n";
    //Cleanup
    for( auto &pp : gene_ptrs){
        delete pp;
    }
    return 0;
/*
    cerr << "Printing tree\n";
    ec.df_execute([](int depth, exon data, int arcval){
        for(int i=0;i<depth;++i){
            cout << "\t";
        }
        cout << data << "\t"  << arcval <<"\n";
    });
    return 0;
*/
}
