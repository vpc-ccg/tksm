/*
 *  Left to do:
 *      Translocations
 *      Multiple isoforms (Just run the single isoform code #expression times)
 *      Sanity checks on various types of gene 2 gene orientations verify sequences are correct
 *      User entered gene list 2 fuse
 *      reverse_complement function (is empty now)
 *      Remove unused graph code
 *      Reimplement readthrough with tree
 *      Refactor 
 *      
 *
 *
 *
 */

#include <filesystem>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <fstream>
#include <string_view>
#include <utility>
#include <sstream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>

#include <cctype>
#include <cassert>

#include "IITree.h"

#include <zlib.h>
#include "kekseq.h"


#if DEBUG
    using std::cerr;
#else
class mockstream{
    template<class K>
    friend mockstream& operator << (mockstream &mock, const K &k){
        return mock;
    }
};
mockstream cerr;
#endif


using std::string;
using std::vector;
using std::cout;
using std::string_view;
using std::map;
using std::ostream;
using std::shared_ptr;
using std::set;


vector<string> rsplit(string str, string delim){
    std::vector<std::string> splits;
    size_t p1 = 0;
    size_t p2 = 0;
    while((p2= str.find(delim,p1)) != string::npos){
        splits.push_back(str.substr(p1,p2-p1));
        p1 = p2+delim.size();
    }
    splits.push_back(str.substr(p1));
    return splits;
}


enum class cigar_character_type{
    matched, onquery, ontemplate, hardclip, softclip, notcigar,
};

#define CIGAR_CHARACTERS "DHIMNPSX="
cigar_character_type what_is_this_cigar(char c){
    switch(c){
        case 'M':
        case '=':
        case 'X':
            return cigar_character_type::matched;
        case 'D':
        case 'N':
            return cigar_character_type::ontemplate;
        case 'I':
        case 'P':
            return cigar_character_type::onquery;
        case 'H':
            return cigar_character_type::hardclip;
        case 'S':
            return cigar_character_type::softclip;
        default:
            return cigar_character_type::notcigar;
    }
}

class cigar{
    std::vector<std::pair<int,char>> cigarray;
    public:
    cigar(const string &cgr){
        size_t prev = -1;
        size_t index = cgr.find_first_of(CIGAR_CHARACTERS);
        while(index != string::npos){
            int len = std::stoi(cgr.substr(1+prev,index));
            cigarray.push_back(std::make_pair(len,cgr[index]));
            prev = index;
            index = cgr.find_first_of(CIGAR_CHARACTERS,index + 1);
        }
    }
    decltype(cigar::cigarray.begin()) begin(){
        return cigarray.begin();
    }

    decltype(cigar::cigarray.end()) end(){
        return cigarray.end();
    }

    auto operator [](size_t index) const{
        return cigarray[index];
    }
};


class interval{
    public:
        int start;
        int end;

        interval():start(0),end(0){}
        interval(int start, int end): start(start), end(end){}
        interval(const interval &i1, const interval &i2): start(std::min(i1.start, i2.start)), end(std::max(i1.end, i2.end)){}
        int distance( const interval &other) const{
            if( other.start < this->start){
                return other.distance(*this);
            }
            return other.start - this->end;
        }
        int overlap(const interval &other)  const{
            if (other.end <= start){ //BEFORE
                return 0;
            }
            else if (other.start >= end){ //AFTER
                return 0;
            }
            else if (other.start >= start && other.end <= end){ // IN
                return other.end - other.start;
            }
            else if (other.start < start && other.end > end) { // AROUND
                return end - start;
            }
            else if (other.start < start && other.end < end && other.end > start){ // LEFT OVERLAP
                return other.end - start;
            }
            else if (other.start > start && other.start < end && other.end > end){ //RIGHT OVERLAP
                return end - other.start;
            }
            return 0;
        }
        friend int larger_interval( const interval &i1, const interval &i2){
            return std::max(i1.end - i1.start, i2.end - i2.start);
        }
        double reciprocal(const interval &other) const {
            return static_cast<double>(overlap(other)) / larger_interval(*this, other);
        }
};


class ginterval: public interval{
    public:
        string chr;
        bool plus_strand;
        ginterval(): interval(0,0), chr(""), plus_strand(true){}
        ginterval(string chr, int start, int end, const string &strand): 
            interval(start, end), chr(chr), plus_strand(strand=="+"){}
        ginterval( const ginterval &g1, const ginterval &g2) : interval(g1,g2), chr(g1.chr), plus_strand(g1.plus_strand) {} // Merge constructor
        int overlap( const ginterval &other) const{
            if( chr != other.chr){
                return 0;
            }
            return interval::overlap(other);
        } 
        double reciprocal(const ginterval &other) const {
            return static_cast<double>(overlap(other)) / larger_interval(*this, other);
        }
};


struct gene: public ginterval{
    string gene_id;
    string gene_name;
    gene() : ginterval(),gene_id("NAN"), gene_name("NAN") {}
    gene(const string &id) : gene_id(id) {} //Mock constructor for map access
    gene(string chr, int start, int end, const string &strand, const string &gene_id,
            const string &gene_name) : ginterval(chr, start, end, strand),
    gene_id(gene_id), gene_name(gene_name){}



    friend auto operator ==(const gene &a, const gene &b){
        return a.gene_id == b.gene_id;
    }
    friend auto operator !=(const gene &a, const gene &b){
        return ! (a==b);
    }
    friend auto operator <=> (const gene &a, const gene &b){
        return a.gene_id <=> b.gene_id;
    }
};

namespace std
{
    template<> struct hash<gene>
    {
        std::size_t operator()(gene const& s) const noexcept
        {
            return std::hash<std::string>{}(s.gene_id);
        }
    };
}

struct transcript: public ginterval{
    string transcript_id;
    gene* gene_ref;
    transcript(string chr, int start, int end, const string &strand, const string &transcript_id, gene *gref):
        ginterval(chr, start, end, strand), transcript_id(transcript_id), gene_ref(gref) {}
};

struct exon: public ginterval{
    string exon_id;
    string strand;

    virtual ~exon() {}
    gene * gene_ref;

    exon(string chr, int start, int end, const string &strand, const string &exon_id, gene* gref, shared_ptr<transcript> tref): 
        ginterval(chr, start, end, strand), exon_id(exon_id), gene_ref(gref){}
    exon() : exon_id("NULL") {}

    exon(const exon &g1, const exon &g2, const string &id) : ginterval(g1,g2), exon_id(id), strand(g1.strand), gene_ref(g1.gene_ref)  {
    }
    exon(const exon &g1, const exon &g2) : ginterval(g1,g2), exon_id(g1.exon_id), strand(g1.strand), gene_ref(g1.gene_ref){
    } // Merge constructor

    bool operator<( const exon &other) const{
        if( gene_ref-> chr != other.gene_ref->chr){
            return gene_ref->chr < other.gene_ref->chr;
        }
        if( gene_ref != other.gene_ref){
            return *gene_ref < *other.gene_ref;
        }
        if( start == other.start){
            return end < other.end;
        }
        return start < other.start;
        //}

    }
    bool operator==(const exon &other) const{
        return exon_id == other.exon_id;
    }
};


ostream &operator << ( ostream &ost, const ginterval &ex){
    ost << ex.chr << ":" << ex.start << "-" << ex.end;
    return ost;
}


ostream &operator << ( ostream &ost, const exon &ex){
    ost << dynamic_cast<const ginterval&>(ex) << " " << ex.strand << " " << ex.exon_id;
    if(ex.gene_ref != nullptr){
        ost << " " << ex.gene_ref->gene_id;
    }
    return ost;
}
ostream &operator << ( ostream &ost, const gene &ex){
    ost << (ginterval) ex <<  " " << ex.gene_name ;
    return ost;
}


struct segment{
    ginterval tmplt;
    interval query;
    segment(string chr, int start, int end, int s2, int e2, const string &strand): 
        tmplt(chr, s2, e2, strand), query(start,end){}
};

struct mapping{
    string rid;
    vector<segment> segments;
    bool primary;
    mapping(): rid("-1"){}
    mapping(const string &paf, int max_skip, int min_segment){

        std::istringstream ps(paf);


        ps >> rid;

        int rlen;
        int rstart;
        int rend;
        ps >> rlen;
        ps >> rstart;
        ps >> rend;

        std::string strand;
        ps >> strand;
        bool complemented = strand == "-";
        std::string chr;
        ps >> chr;

        if(chr.find("chr")!= std::string::npos){
            chr = chr.substr(3);
        }

        int template_len;
        int template_start;
        int template_end;
        int num_matches;
        int alig_block_len;
        int mapq;
        std::string field;
        ps >> template_len >> template_start >> template_end >> num_matches >> alig_block_len >> mapq;

        ps >> field;
        string cigar_str = "-1";

        while(!ps.eof()){ 
            ps >> field;
            if (field.substr(0,2) == "tp"){
                primary = (field == "tp:A:P");
            }
            if (field.substr(0,2) == "cg"){
                cigar_str = field.substr(5);
            }
        }

        if( cigar_str == "-1"){
            cerr << "Paf line doesn't include alignment cigar! Exiting!.\n";
            exit(-1);
        }
        cigar cig(cigar_str);

        std::vector<segment> aligs;

        int st = template_start;
        int sq = rstart;
        int et, eq;
        if(complemented){
            sq = rend - 1;
        }
        for( auto pair : cig){
            int length = pair.first;
            char c = pair.second;
            cigar_character_type type = what_is_this_cigar(c);

            if( type == cigar_character_type::matched){
                et = st + length;

                if(complemented){
                    eq = sq - length;
                    aligs.emplace_back(chr,eq,sq-1,st,et-1,strand);
                }
                else{
                    eq = sq + length;
                    aligs.emplace_back(chr,sq,eq-1,st,et-1,strand);
                }

                sq = eq;
                st = et;
            }
            else if( type == cigar_character_type::ontemplate){
                st = st + length;
            }
            else if( type == cigar_character_type::onquery){
                if(complemented){
                    sq = sq - length;
                }   
                else{
                    sq = sq + length;
                }
            }
            else{

            }
        }
        if( aligs.size() == 0){
            cerr << "NOALIG " << rid << "\n";
        }
        if( max_skip > 0 && aligs.size() > 0){

            st = aligs[0].tmplt.start;
            et = aligs[0].tmplt.end;
            sq = aligs[0].query.start;
            eq = aligs[0].query.end;

            for(auto iter = std::begin(aligs); std::next(iter) != std::end(aligs); ++iter){

                if(complemented){

                    if( std::next(iter)->tmplt.start - et < max_skip){
                        sq = std::next(iter)->query.start;
                        et = std::next(iter)->tmplt.end;

                    }
                    else{
                        if(eq-sq >= min_segment && et-st >= min_segment){
                            segments.emplace_back(chr,sq,eq,st,et,strand);
                        }
                        auto nxt = std::next(iter);
                        st = nxt->tmplt.start;
                        et = nxt->tmplt.end;
                        sq = nxt->query.start;
                        eq = nxt->query.end;
                    }
                }
                else{

                    if( std::next(iter)->tmplt.start - et < max_skip){
                        et = std::next(iter)->tmplt.end;
                        eq = std::next(iter)->query.end;
                    }
                    else{
                        if(eq-sq >= min_segment && et-st >= min_segment){
                            segments.emplace_back(chr,sq,eq,st,et,strand);
                        }
                        auto nxt = std::next(iter);
                        st = nxt->tmplt.start;
                        et = nxt->tmplt.end;
                        sq = nxt->query.start;
                        eq = nxt->query.end;
                    }
                }
            }
            if(eq-sq >= min_segment && et-st >= min_segment){
                segments.emplace_back(chr,sq,eq,st,et,strand);
            }
        }
        else{
            segments = aligs;
        }
    }
};


template <class N, class A>
class tree{

    public:
    map<N,std::pair<A,tree<N,A>>> children;
    tree<N,A> *parent;
    N data;


    tree() : parent(nullptr), data(N{}){}
    tree(const N &data) :parent(nullptr), data(data){}
    void add_child( const N &n, const A &a){
        children.emplace(n, std::make_pair(a,tree<N,A>{n}));
        children.at(n).second.parent = this;
    }

    tree<N,A> &operator[](const N &n){
        return children.at(n).second;
    }

    tree<N,A> &try_get(const N &n, const A &a){
        if(children.find(n) == children.end()){
            add_child(n,a);
        }
        return children.at(n).second;
    }
    A& value(const N &n) {
        return children.at(n).first;
    }
    A value(const N &n) const {
        return children.at(n).first;
    }


    template< class Func>
    void df_execute2(int depth, Func foo) const{
        foo(depth, this);
        for( auto &p : children){
            p.second.second.df_execute2(depth+1,foo);
        }
    }

    template< class Func>
    void df_execute2(Func foo)const {
        foo(0,this);
        for( auto &p : children){
            p.second.second.df_execute2(1,foo);
        }
    }
    template< class Func>
    void df_execute(int depth, A& arc_val, Func foo){
        foo(depth,data, arc_val);
        for( auto &p : children){
            p.second.second.df_execute(depth+1,p.second.first,foo);
        }
    }
    template< class Func>
    void df_execute(Func foo){
        foo(0,data,A{});
        for( auto &p : children){
            p.second.second.df_execute(1,p.second.first,foo);
        }
    }
};

template <class N, class A>
class graph{
    public:
    vector<std::pair<N,A>> nodes;
    vector< map<size_t, A> > arcs;
    map< N, size_t> reverse_index;

    graph(){}

    void add( N node){
        reverse_index[node] = nodes.size();
        nodes.push_back(std::make_pair(node,A{}));
        arcs.emplace_back();
    }

    A &arc(size_t i, size_t j){
        assert(i < nodes.size());
        assert(j < nodes.size());
        return arcs[i][j];
    }
    A &arc(const N &a, const N &b){
        size_t i = reverse_index[a];
        size_t j = reverse_index[b];
        return arc(i,j);
    }
    auto begin(){
        return nodes.begin();
    }
    auto end(){
        return nodes.end();
    }
    class neighbour{
        graph &owner;
        size_t target;
        public:
        neighbour(graph &owner, size_t target) : owner(owner), target(target) {}
        class nei_iter{
            graph &owner;
            size_t target;
            decltype(owner.arcs[0].begin()) index;
            public:

            nei_iter(graph &owner, size_t target) : owner(owner), target(target), index(owner.arcs[target].begin()) {}
            nei_iter(graph &owner, size_t target, decltype(owner.arcs[0].begin()) index) : owner(owner), target(target), index(index) {}
            auto &operator *(){
                return owner.nodes[index->first].first;
            }
            nei_iter &operator++(){
                ++index;
                return *this;
            }
            nei_iter operator++(int){
                nei_iter tmp = *this;
                ++(*this);
                return tmp;
            }
            friend bool operator== (const nei_iter& a, const nei_iter& b) { return a.target == b.target && a.index == b.index; };
            friend bool operator!= (const nei_iter& a, const nei_iter& b) { return !(a==b); };
        };
        nei_iter begin(){
            return nei_iter(owner,target);
        }

        nei_iter end(){
            return nei_iter(owner,target,owner.arcs[target].end());
        }
    };
    neighbour neighbours(size_t index){
        assert(index < nodes.size());
        assert(index >= 0);
        return neighbour(*this, index);
    }
    neighbour neighbours(const N &key){
        size_t index = reverse_index[key];
        return neighbour(*this, index);
    }

    bool in( const N &key){
        return reverse_index.find(key) != reverse_index.end();
    }
    A& value(const N &key){
        size_t index = reverse_index[key];
        assert(in(key));

        return nodes[index].second;
    }
};


string strip_str(const std::string &inpt, const string &chrs)
{

    auto frst = inpt.find_first_not_of(chrs);
    auto last = inpt.find_last_not_of(chrs);
    if(frst == string::npos || last == string::npos){
        return "";
    }
    return inpt.substr(frst, last-frst+1);
/*
    auto start_it = inpt.begin();
    auto end_it = inpt.rbegin();
    while (std::isspace(*start_it))
        ++start_it;
    while (std::isspace(*end_it))
        ++end_it;
    return std::string(start_it, end_it.base());
*/
}
void strip_for_each(vector<string> &vec, const string &chrs = " "){

    for( string &st : vec){
        st = strip_str(st, chrs);
    }
}
/*
class lazy_split{
    const string &source;
    string delim;

    public:
    vector<size_t> known_delimiters;
    string strip;


    lazy_split(string &source) : source(source), delim(" "), known_delimiters({0}), strip(""){}

    lazy_split(const string &source, const string &delim) : source(source), delim(delim), known_delimiters({0}), strip("") {}
    lazy_split(const string &source, const string &delim, const string &strip) : source(source), delim(delim), known_delimiters({0}), strip(strip) {}

    string_view operator[](size_t index){

        while(index + 1 >= known_delimiters.size()){
            auto it = std::search(source.begin()+known_delimiters.back(),source.end(), std::default_searcher(delim.begin(),delim.end()));

            if( it == source.end()){
                known_delimiters.push_back(source.size() + delim.size() - 1);
                assert( index + 1< known_delimiters.size());
                break;
            }
            known_delimiters.push_back(it - source.begin() + delim.size());  
        }


        string_view sv(source);

        sv.remove_prefix(known_delimiters[index] );
        sv.remove_suffix(source.size() - known_delimiters[index+1] + 1 );
        if(strip.size()>0){
            sv.remove_prefix(std::min(sv.find_first_not_of(strip),sv.size()));
            sv.remove_suffix( sv.size() - std::min(1 + sv.find_last_not_of(strip),sv.size()));
        }
        return sv;

    }
    class iterator{
        public:
            size_t index;
            lazy_split &owner;
            bool end;



            iterator(lazy_split &owner) : index(0), owner(owner), end(owner.source.size()==0) {
            }

            iterator(lazy_split &owner, size_t index) : index(index), owner(owner), end(owner.source.size()==0) {}
            iterator(lazy_split &owner, size_t index, bool end) : index(index), owner(owner), end(end) {}

            string_view operator *(){
                return owner[index];
            }

            iterator &operator++(){
                ++index;
                if(owner.known_delimiters.back() - owner.delim.size() + 1 == owner.source.size() && index + 2 > owner.known_delimiters.size()){
                    end = true;
                    return *this;
                }
                owner[index];
                return *this;
            }
            iterator operator++(int){
                iterator tmp = *this;
                ++(*this);

                if(owner.known_delimiters.back() - owner.delim.size() + 1== owner.source.size()  && index + 2 == owner.known_delimiters.size()){
                    end = true;
                    return *this;
                }
                owner[index];
                return tmp;
            }
            friend bool operator== (const iterator& a, const iterator& b) { 
                if(a.end && b.end){
                    return true;
                }
                return a.index == b.index;
            }
            friend bool operator!= (const iterator& a, const iterator& b) { return !(a==b); };

    };
    iterator begin(){
        return iterator(*this,0);
    }
    iterator end(){
        return iterator(*this, -1 ,true);
    }

};
*/



graph<gene, double> build_gene_graph( string path_to_gtf, const  map<gene, int> &gene2count,bool coding_only = true){


    graph<gene, double> gene_graph;

    std::ifstream file(path_to_gtf);

    string str;
    while( std::getline(file, str)){
        if(str[0] == '#'){
            continue;
        }
        vector<string> fields = rsplit(str, "\t");

        if(fields[2] != "gene"){
            continue;
        }
        string chr{ fields[0]};
        int s = stoi(string{fields[3]});
        int e =  stoi(string{fields[4]});

        string st{ fields[6]};

        string info_str{fields[8]};
        vector<string> info = rsplit(info_str, ";");
        strip_for_each(info , " ");

        string gene_id = ""; 
        string gene_name = "";
        string biotype = "";
        int count = 0;
        for( auto iter = info.begin(); iter != info.end(); ++iter){
            if((*iter).size() <= 1){ continue;}
            string f{*iter};

            vector<string> fs = rsplit(f , " ");
            strip_for_each(fs, "\"");
        
            if(fs[0] == "gene_id"){
                gene_id = fs[1];
                ++ count;
            }
            if(fs[0] == "gene_name"){
                gene_name = fs[1];
                ++ count;
            }
            if(fs[0] == "gene_biotype"){
                biotype = fs[1];

                ++ count;
            }
        }
        if( coding_only && biotype != "protein_coding"){
            continue;
        }
        gene g{chr, s, e, st, gene_id, gene_name};
        if( gene2count.find(g) != gene2count.end()){
            gene_graph.add(g);
        }
    }
    return gene_graph;
}
auto read_gtf_exons( string path_to_gtf, bool coding_only = true){

    vector<exon> annots;

    std::ifstream file(path_to_gtf);

    string str;

    vector<gene *> gptrs;
    gene* current_gene;
    shared_ptr<transcript> current_transcript;
    while( std::getline(file, str)){
        if(str[0] == '#'){
            continue;
        }
        vector<string> fields = rsplit(str, "\t");
        string_view type = fields[2];
        string chr{ fields[0]};
        int s = stoi(string{fields[3]});
        int e =  stoi(string{fields[4]});
        string st{ fields[6]};
        
        if( type == "gene"){
            string info_str{fields[8]};
            vector<string> info = rsplit(info_str, ";");
            strip_for_each(info, " \n\t");

            string gene_id = ""; 
            string gene_name = "";
            string biotype = "";
            int count = 0;
            for( auto iter = info.begin(); iter != info.end(); ++iter){
                if((*iter).size() <= 1){ continue;}
                string f{*iter};
                vector<string> fs = rsplit(f , " ");
                strip_for_each(fs, "\"");
                if(fs[0] == "gene_id"){
                    gene_id = fs[1];
                    ++ count;
                }
                if(fs[0] == "gene_name"){
                    gene_name = fs[1];
                    ++ count;
                }
                if(fs[0] == "gene_biotype"){
                    biotype = fs[1];
                    ++ count;
                }
            }
            if( coding_only && (biotype != "protein_coding")){
                continue;
            }
            gptrs.push_back(new gene(chr,s,e,st,gene_id,gene_name));
            current_gene = gptrs.back();
        }
        else if( type == "transcript"){
            string info_str{fields[8]};
            vector<string> info = rsplit(info_str, ";");
            strip_for_each(info, " \n\t");

            string transcript_id = "";
            string biotype = "";
            int count = 0;
            for( auto iter = info.begin(); iter != info.end(); ++iter){
                if((*iter).size() <= 1){ continue;}
                string f{*iter};
                vector<string> fs = rsplit(f , " ");
                strip_for_each( fs, "\"");
                if(fs[0] == "transcript_id"){
                    transcript_id = fs[1];
                    ++ count;
                }
                if(fs[0] == "transcript_biotype"){
                    biotype = fs[1];
                    ++ count;
                }
            }
            if( coding_only && (biotype != "protein_coding")){
                continue;
            }
            current_transcript = std::make_shared<transcript>(chr,s,e,st,transcript_id,current_gene);

        }
        else if( type == "exon"){

            string info_str{fields[8]};

            vector<string> info = rsplit(info_str, ";");
            strip_for_each(info, " \n\t");

            string biotype = "";
            string tbiotype = "";
            string exon_id = ""; 
            int count = 0;
            for( auto iter = info.begin(); iter != info.end(); ++iter){
                if((*iter).size() <= 1){ continue;}
                string f{*iter};
                vector<string> fs = rsplit(f , " ");
                strip_for_each( fs, "\"");

                if(fs[0] == "exon_id"){
                    exon_id = fs[1];
                    ++ count;
                }
                if(fs[0] == "gene_biotype"){
                    biotype = fs[1];
                    ++ count;
                }if(fs[0] == "transcript_biotype"){
                    tbiotype = fs[1];
                    ++ count;
                }
            }
            if (coding_only && ( (biotype!="protein_coding") || (tbiotype!="protein_coding"))) {
                continue;
            }
            annots.emplace_back(chr,s,e,st,exon_id, current_gene, current_transcript);
        }
    }
    /*
       while( std::getline(file, str)){
       if(str[0] == '#'){
       continue;
       }
       lazy_split fields{str, "\t"};

       if(fields[2] != "exon"){
       continue;
       }
       string chr{ fields[0]};
       int s = stoi(string{fields[3]});
       int e =  stoi(string{fields[4]});

       string st{ fields[6]};

       string info_str{fields[8]};
       lazy_split info{info_str, ";", " \n\t"};


       string gene_id = ""; 
       string transcript_id = ""; 
       string exon_id = ""; 
       int count = 0;
       for( auto iter = info.begin(); iter != info.end(); ++iter){
       if((*iter).size() <= 1){ continue;}
       string f{*iter};
       lazy_split fs{f , " ", "\""};

       if(fs[0] == "gene_id"){
       gene_id = fs[1];
       ++ count;
       }
       else if(fs[0] == "transcript_id"){
       transcript_id = fs[1];
       ++ count;
       }
       else if(fs[0] == "exon_id"){
       exon_id = fs[1];
       ++ count;
       }
       if(count == 3){
       break;
       }
       }

       annots.emplace_back(chr, s, e, st, exon_id);
       }
       */
    return make_pair(annots,gptrs);
}

map<string, IITree<int, size_t>> make_exon_tree( const vector<exon> &annots){

    map<string,IITree<int, size_t>> itree;

    size_t index = 0;

    //for(const exon &annot: annots){
    for(auto iter = annots.begin(); iter != annots.end(); ++iter){
        auto &annot = *iter;
        itree[annot.chr].add(annot.start, annot.end, index);
        //       cerr << annot.chr << "\t" << annot << "\t" << index << "\n";
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
        cerr << "1!\n";

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
        }(r);
        
        cerr << "2!\n";
        auto find_exon = [&counts_of_genes,&exon_ref,&exon_forest] (const ginterval &g) -> exon{
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
        cerr << "3!\n";

        auto first = r.segments.begin();
        exon e = find_exon(first->tmplt);
        if( e == exon{}){
            continue;
        }
        
        tree<exon, int> *prev = &exon_transitions.try_get(e, 0);
        if(prev->parent != nullptr){
            prev->parent->value(e)++;
        }
        cerr << "4!\n";

        exon preve{};
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
        cerr << "\nDONE\n";

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

void build_and_print_updated_cdna_ref( const string &in_path,
        const string &out_path, graph<gene, double> &gene_graph){

    std::ifstream in_file(in_path);
    std::ofstream out_file(out_path);

    string header;
    string seq;
    while( std::getline(in_file, header)){
        std::getline(in_file, seq);
        vector<string> head_fields = rsplit(header, " ");
        strip_for_each(head_fields, ">");
//        string_view transcript_id = head_fields[0];

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

vector< std::pair<gene,gene>> generate_random_fusions( graph<gene, double> &gene_graph, map<string, int> fusion_count_per_chrX2, int seed = 121){
    std::mt19937 rand_gen{std::random_device{}()};
    rand_gen.seed(seed);
    map<string, std::pair<int,int>> gene_ranges_per_contig = find_gene_counts_per_contig(gene_graph);
    vector< std::pair<gene,gene>> fusions;

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

        bool odd = true;
        for( auto iter = picked_values.begin(); std::next(iter) != picked_values.end(); ++iter){
            if( odd){
                //cout << gene_graph.nodes[*iter].first << "\t" << gene_graph.nodes[*std::next(iter)].first << "\n";
                fusions.push_back(std::make_pair(gene_graph.nodes[*iter].first, gene_graph.nodes[*std::next(iter)].first));
            }
            odd = ! odd;
        }  
    }
    return fusions;
}

vector<std::pair<gene,gene>> generate_random_fusions( graph<gene, double> &gene_graph, int count, int seed = 42){
        
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
        

    return generate_random_fusions(gene_graph, fusion_count_per_chrX2 ,seed);
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

auto flatten_gtf( vector<exon> &exons, double min_overlap = 0.90){
//Ensure sorted
//
    std::sort(exons.begin(),exons.end());
    vector<exon> processed;
    map<string, string> mergeindex;

    if(exons.empty()){
        throw std::runtime_error("Empty GTF file!!\n");
    }
    
    exon merger(exons.front());
    int cnt = 0;
    int dupcount = 0;
    for(auto iter = std::next(std::begin(exons)); iter != std::end(exons); ++iter){
        if( merger.reciprocal(*iter) > min_overlap &&
                merger.gene_ref == iter->gene_ref){ //Also check strands and genes

            if(mergeindex.find(iter->exon_id) != mergeindex.end()){
                ++dupcount;
            }
            else{
                merger = exon(merger, *iter);//, string{"E"} + std::to_string(cnt) ); //Write smarter iding
                mergeindex[iter->exon_id] = merger.exon_id;
            }
        }
        else{
            processed.push_back(merger);
            merger = exon(*iter);
            cnt += 1;
        }
    }

    return make_pair(processed, mergeindex);
}

int qmain(){
    tree<int, double> tre(0);


    tre.add_child(1,0.5);
    tre.add_child(2,0.5);
    tree<int,double> &t2 = tre[2];
    t2.add_child(3,1);
    auto &t3 = tre[1];

    t3.add_child(4,0.3);
    t3.add_child(5,0.7);


    vector<int> vals;
    auto lt3 = [&vals] (int d, int val, double e) -> void{
        if(e >0.6){
            vals.push_back(val);
        }
    };
    tre.df_execute(lt3);
    for(int i: vals){
        cout << i << "\t";
    }
    cout << "\n";
    tre.df_execute( [](int d, int val, double e) ->void {
        for(int i=0;i<d;++i){
            cout << "\t";
        }
        cout << val << "\n";
    });
    return 0;
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

class isoform{
    public:
    vector<exon> segments;
    isoform(vector<exon> &segs) : segments(segs){}

};

enum class bpstrategy{
    uniform,
    exon_biased,
};

vector<std::pair<int,int>> generate_fusion_breakpoints( const vector<std::pair<gene,gene>> &gene_pairs, bpstrategy strat, int seed = 1453){

    vector<std::pair<int,int>> breakpoints;
    std::mt19937 rand_gen{std::random_device{}()};
    rand_gen.seed(seed);

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

map<gene, int>  generate_transcript_expression( const tree<exon, int> &exon_tree){

    map<gene,int> expressions;
    for( auto p: exon_tree.children){
        expressions[*p.first.gene_ref]+= p.second.first;
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

void  print_normal_transcript_fasta(const tree<exon, int> &exon_tree,
        const string &out_cdna_path,
        const map<string, string> &transcript2sequence,
        const map<gene, int> &transcriptexpression){

    std::ofstream ost(out_cdna_path, std::ios_base::out);

    for(auto [id, seq] : transcript2sequence){
        ost << ">" << id <<  " depth=0"  << "\n";//Add depth
        ost << seq << "\n";
    }
}

/*
            auto find_sum = [&sum](int depth, const tree<exon, int> *nod) -> void{
                if(depth == 1){
                    sum+= nod->parent->children[nod->data].first;
                }
            };
            child.df_execute2(find_sum);
         */
map<string, vector<isoform>> generate_fusion_isoforms(
        const vector<std::pair<gene, gene>> &fusions, 
        const vector<std::pair<int,int>> &fusion_breakpoints, 
        const vector<int> &fusion_expression,
        const tree<exon, int> &exon_tree,
        int seed = 42){

    std::mt19937 rand_gen{std::random_device{}()};
    rand_gen.seed(seed);
    map<string, vector<isoform>> isoforms;

    map<gene, vector<exon>> gene2firstexons; //reverse_index
    map<gene, int>          gene2count;
    for( const auto &pair : exon_tree.children){
        const exon &e = pair.first;
//        int count = pair.second.first;
//        const tree<exon, int> &child = pair.second.second;
        gene2count[*e.gene_ref]+=  (exon_tree.value(e));
        gene2firstexons[*e.gene_ref].push_back(e);
    }


    cerr << "There are " << gene2firstexons.size()  << " genes!\n";

    auto find_isoform = [&gene2count, &gene2firstexons, &exon_tree, &rand_gen ] (const gene &g, int limit){
        vector<exon> isoform;
    
        const tree<exon, int> *current = nullptr;
        int count = gene2count.at(g);
        std::uniform_int_distribution<> dist(0, count);
        int sum = 0;
        int rand_val = dist(rand_gen);
        for(const exon &e : gene2firstexons[g]){
            const auto &p = exon_tree.children.at(e);
            int count = p.first;
            current = &p.second;
            sum += count;
            if(sum > rand_val){
                break;
            }
        }
/*
        for(const exon &e : gene2firstexons[g]){
            const auto &p = exon_tree.children.at(e);
            int count = p.first;
            //const tree<exon, int> *current = &p.second;
            std::uniform_int_distribution<> dist(0, count);
            int rand_val = dist(rand_gen);
            int sum = 0;
*/


        while(current->data.end < limit){
            if( current->children.size() == 0){
                break;
            }
            int sum = 0;
            int rand_val = dist(rand_gen);
            for(const auto &pc : current->children){
                int cc = pc.second.first;
                sum+=cc;
                if(sum > rand_val){
                    isoform.push_back(current->data);
                    current = &pc.second.second;
                    break;
                }   //Decide
            }
 //           }
        }
        isoform.push_back(current->data);
        return isoform;
    };
    auto bpiter = fusion_breakpoints.begin();
    auto exiter = fusion_expression.begin();
    for( auto fusion : fusions){
        const gene &g1 = fusion.first;
        const gene &g2 = fusion.second;


        auto iso1 = find_isoform(g1, bpiter->first);
        auto iso2 = find_isoform(g2, bpiter->second);
        iso1.reserve(iso1.size()+iso2.size());
        iso1.insert(iso1.end(),iso2.begin(),iso2.end());
        isoforms[g1.gene_id + "_" + g2.gene_id].push_back(isoform{iso1});
        
        ++bpiter;
        ++exiter;
    }

    return isoforms;
}

void reverse_complement(string &seq){

}

void generate_and_print_fusion_fasta(
        const map<string, vector<isoform>> &isoforms,
        const map<string,string> &sequences,
        const string &out_cdna_path,
        int seed = 42){

    std::mt19937 rand_gen{std::random_device{}()};
    rand_gen.seed(seed);

    std::ofstream ost(out_cdna_path);
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
                string seq = sequences.at(e.chr).substr( e.start, e.end - e.start);
                if(! e.plus_strand){
                    reverse_complement(seq);
                }

                ost << e << "\t"  << seq << "\n";
            }
            ++cnt;
        }

    }
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

int main(int argc, char **argv){


    //parameters will be moved to argument parser when done
    string path_to_aligs = argv[1];
    string path_to_gtf   = argv[2];
    string path_to_cdna  = argv[3];
    string path_to_dna {argv[4]};
    string out_cdna_path = argv[5];
    
    int min_dist = 1000;
    int max_dist = 500000;
    
    //

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



/*
    map<string, string> chr2sequence;
    do{
        int r;
        gzFile gzf = gzopen(path_to_dna.c_str(), "r");
        kekseq::kseq<gzFile,gzread,3000000> fr(gzf);//load wholegenome 
        while((r = fr.read()) >= 0){
            chr2sequence[fr.name.s] = string{fr.seq.s};
        }
        gzclose(gzf);

    }while(0);
*/
    map<string, string> chr2contig = read_fasta_fast(path_to_dna);



    vector<mapping> reads = convert_aligs_to_segments(path_to_aligs);   
    auto [gtf_exons, gene_ptrs] = read_gtf_exons(path_to_gtf);

    auto [merged_exons, merge_index] = flatten_gtf(gtf_exons,0.75); // vector<exon>, map<string, string>



    //map<string, IITree<int, size_t>> exon_tree = make_exon_tree(merged_exons);
    map<string, IITree<int, size_t>> exon_tree = make_exon_tree(gtf_exons);


    cerr << "Counting reads on tree\n";
    tree<exon, int> ec = count_reads_on_tree(gtf_exons, exon_tree, reads, min_dist,  max_dist);

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

    map<gene, vector<exon>> gene2firstexons; //reverse_index
    map<gene, int>          gene2count;
    for( const auto &pair : ec.children){
        const exon &e = pair.first;
//        int count = pair.second.first;
//        const tree<exon, int> &child = pair.second.second;
        gene2count[*e.gene_ref]+=  (ec.value(e));
        gene2firstexons[*e.gene_ref].push_back(e);
    }
    cerr << "Building graphs\n";
    graph<gene, double> gene_graph = build_gene_graph(path_to_gtf, gene2count); // I should replace this
    
    cerr << "Building graphs\n";
    vector< std::pair<gene,gene>> fusions = generate_random_fusions(gene_graph,200);
        

    cerr << "Generating breakpoints\n";
    vector<std::pair<int,int>> fusion_breakpoints = generate_fusion_breakpoints( fusions, bpstrategy::uniform);

    map<gene,int> transcript_expression = generate_transcript_expression(  ec);//Convert this to transcript 
    vector<int> fusion_expression = generate_fusion_expression( fusions, ec);


//    cerr << "Printing Fasta\n";
//    print_normal_transcript_fasta(ec, out_cdna_path, transcript2sequence, transcript_expression);
    
    cerr << "Simulating fusionsa\n";
    map<string, vector<isoform>> fusion_isoforms = generate_fusion_isoforms(fusions, fusion_breakpoints, fusion_expression, ec);
    cerr << "Printing Fasta\n";
    generate_and_print_fusion_fasta(fusion_isoforms, chr2contig, out_cdna_path);

    cerr << "Cleaning Up\n";
    //Cleanup
    for( auto &pp : gene_ptrs){
        delete pp;
    }
    return 0;
/** tree exec example   
    ec.df_execute2([](int depth, tree<exon, int> *nod){
        if( nod->parent != nullptr && nod->parent->data != exon{} &&  nod->parent->data.gene_ref != nod->data.gene_ref){
            print_tsv(nod->parent->data, " -- ", nod->parent->value(nod->data), " --> ", nod->data);
        }
    });
    return 0;
    ec.df_execute([](int depth, exon data, int arcval){
        for(int i=0;i<depth;++i){
            cout << "\t";
        }
        cout << data << "\t"  << arcval <<"\n";
    });
**/
    return 0;
}
