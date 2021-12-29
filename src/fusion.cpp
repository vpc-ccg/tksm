#ifndef FUSION_H
#define FUSION_H

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

using std::ofstream;
using std::set;
using std::vector;
using std::string;
using std::map;

std::mt19937 rand_gen{std::random_device{}()};

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
        //    std::cerr << e.what() << "\n";
        //    throw e;
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

vector<int>  generate_fusion_expression( const vector<std::pair<gene,gene>> &fusions, const tree<exon,int> &exon_tree){
    vector<int> expressions;
    for( auto p: fusions){
        expressions.push_back(10);
    }
    return expressions;
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
            std::cout << gptrs.size() << "\t" << *iter  << "\t" << *(iter+1) << "\t" << cnt << "\t" << picked_values.size()<< "\n";

            if( odd == 0){
                //cout << gene_graph.nodes[*iter].first << "\t" << gene_graph.nodes[*std::next(iter)].first << "\n";
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

int main(int argc, char **argv){

    cxxopts::Options options("RNAInfuser Fusion", "Fusion module of RNAInfuser");

    options.add_options()
        ("p,paf",  "Path to whole genome mappings in paf format", cxxopts::value<string>())
        ("g,gtf",  "Path to gtf annotation file", cxxopts::value<string>())
        ("d,dna",  "Path to human genome reference file", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())

        ("molecule-count", "Number of molecules to simulate", cxxopts::value<int>()->default_value("100000"))
        ("fusion-count", "Number of gene fusions to simulate", cxxopts::value<int>()->default_value("0"))
        ("translocation-ratio", "Percentage of fusions to be simulated as translocations", cxxopts::value<double>()->default_value("0.1"))

        ("rt-minimum-distance", "Minimum distance between readthrough genes to simulate",cxxopts::value<int>()->default_value("1000"))
        ("rt-maximum-distance", "Maximum distance between readthrough genes to simulate",cxxopts::value<int>()->default_value("500000"))

        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("h,help", "Help screen")
    ;
    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    std::vector<string> mandatory {{"paf", "gtf", "output"}};

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
//    string path_to_cdna  {args["cdna"].as<string>()};
    string path_to_dna   {args["dna"].as<string>()};
    string out_cdna_path {args["output"].as<string>()};
    int seed = args["seed"].as<int>();;

    rand_gen.seed(seed);
    int min_dist = args["rt-minimum-distance"].as<int>();
    int max_dist = args["rt-maximum-distance"].as<int>();
 

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



    set<gene> expressed_genes;
    for(const auto &tr : ec.children){
        expressed_genes.insert(*tr.first.gene_ref);
    }
    

    std::cerr << "Generating fusions\n";
    int fusion_count = args["fusion-count"].as<int>();
    double translocation_ratio = args["translocation-ratio"].as<double>();
    vector<gene> exp_gene_vec;
    exp_gene_vec.insert(exp_gene_vec.begin(),expressed_genes.begin(),expressed_genes.end());
    
    std::sort(exp_gene_vec.begin(), exp_gene_vec.end() ,[](const gene &g1, const gene &g2) -> bool {
        return g1.chr > g2.chr; 
    });

    std::cout << exp_gene_vec.size() << "\n";
    vector< std::pair<gene,gene>> fusions = generate_random_fusions(exp_gene_vec, fusion_count * ( 1- translocation_ratio), fusion_count * translocation_ratio);
//    vector<std::pair<int,int>> fusion_breakpoints = generate_fusion_breakpoints( fusions, bpstrategy::uniform);
    vector<int> fusion_expression = generate_fusion_expression( fusions, ec);

    map<string, vector<isoform>> normal_isoforms = generate_normal_isoforms(ec);

    std::cerr << "Simulating fusions\n";
    map<string, vector<isoform>> fusion_isoforms = generate_fusion_isoforms(fusions, fusion_expression, ec);

    ofstream outfile {out_cdna_path};
    for( const auto &pv : normal_isoforms){
        int cnt = 0;
        for(const isoform &iso : pv.second){
            print_mdf(outfile, pv.first+"_"+std::to_string(cnt), pcr_molecule{pv.first+"_"+std::to_string(cnt),iso}, vector<vector<std::pair<int, char>>>{iso.segments.size()});
            ++cnt;
        }
    }
    std::cerr << "LE: " << fusion_isoforms.size() << "\n";
    for( const auto &pv : fusion_isoforms){
        int cnt = 0;
        for(const isoform &iso : pv.second){
            print_mdf(outfile, pv.first+"_"+std::to_string(cnt), pcr_molecule{pv.first+"_"+std::to_string(cnt),iso},vector<vector<std::pair<int, char>>>{iso.segments.size()});
            ++cnt;
        }
    }

    return 0;
}


#endif
