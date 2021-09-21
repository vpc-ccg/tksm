

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <fstream>
#include <string_view>
#include <utility>
#include <sstream>
#include <map>

#include <cctype>
#include <cassert>

#include "IITree.h"

using std::string;
using std::vector;
using std::cout;
using std::string_view;
using std::map;

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
};

class ginterval: public interval{
    public:
        string chr;
        ginterval(): chr(""),interval(0,0){}
        ginterval(string chr, int start, int end): 
            chr(chr), interval(start, end){}
};

class annotation: public ginterval{
    string gene;
    string transcript;
    string exon;
    string strand;
    public:
    annotation(string chr, int start, int end, string strand, string gene,
            string transcript, string exon) : ginterval(chr, start, end),
    strand(strand),
    gene(gene), transcript(transcript), exon(exon){}
};

struct segment{
    ginterval tmplt;
    interval query;
    bool reverse_complemented;

    segment(string chr, int start, int end, int s2, int e2, bool reverse_complemented): 
        tmplt(chr, s2, e2), query(start,end), reverse_complemented(reverse_complemented)  {}
};

class read{
    vector<segment> segments;
    public:
    read(const string &paf, int max_skip){

        std::istringstream ps(paf);

        std::string id;
        ps >> id;

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
        while(!ps.eof() && field.substr(0,2) != "cg"){ ps >> field;}

        if( field.substr(0,2) != "cg"){
            std::cerr << "Paf line doesn't include alignment cigar! Exiting!." << std::endl;
            exit(-1);
        }
        cigar cig(field.substr(5));

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
                    aligs.emplace_back(chr,eq,sq-1,st,et-1,strand=="-");
                }
                else{
                    eq = sq + length;
                    aligs.emplace_back(chr,sq,eq-1,st,et-1,strand=="-");
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
        std::vector<segment> merged;
        if( max_skip > 0){

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
                        merged.emplace_back(chr,sq,eq,st,et,strand=="-");
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
                        merged.emplace_back(chr,sq,eq,st,et,strand=="-");
                        auto nxt = std::next(iter);
                        st = nxt->tmplt.start;
                        et = nxt->tmplt.end;
                        sq = nxt->query.start;
                        eq = nxt->query.end;
                    }
                }
            }
            merged.emplace_back(chr,sq,eq,st,et,strand=="-");

        }
        else{
            merged = aligs;
        }
    }
};


template <class N, class A>
class graph{
    vector<N> nodes;
    vector< map<size_t, A> > arcs;

    public:
    void add( N node){
        nodes.push_back(node);
        arcs.emplace_back();
    }

    A &arc(size_t i, size_t j){
        assert(i < nodes.size());
        assert(j < nodes.size());
        return arcs[i][j];
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
            std::pair<N, A> operator *(){
                return std::make_pair(owner.nodes[index->first], index->second);
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
};


string strip_str(const std::string &inpt)
{
    auto start_it = inpt.begin();
    auto end_it = inpt.rbegin();
    while (std::isspace(*start_it))
        ++start_it;
    while (std::isspace(*end_it))
        ++end_it;
    return std::string(start_it, end_it.base());
}

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
            auto it = std::search(source.begin()+known_delimiters.back(),source.end(), std::boyer_moore_searcher(delim.begin(),delim.end()));

            if( it == source.end()){
                known_delimiters.push_back(source.size() + delim.size() - 1);
                assert( index + 1< known_delimiters.size());
                break;
            }
            known_delimiters.push_back(it - source.begin() + delim.size());  
        }


        string_view sv(source);

        sv.remove_prefix(known_delimiters[index]);// + delim.size()-1
        sv.remove_suffix(source.size() - known_delimiters[index+1] + delim.size() );
        if(strip.size()>0){
            sv.remove_prefix(std::min(sv.find_first_not_of(strip),sv.size()));
            sv.remove_suffix( sv.size() - std::min(1+sv.find_last_not_of(strip),sv.size()));
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

vector<annotation> read_gtf_file( string path_to_gtf){

    vector<annotation> annots;

    std::ifstream file(path_to_gtf);

    string str;
    while( std::getline(file, str)){
        if(str[0] == '#'){
            continue;
        }
        lazy_split fields(str, "\t");

        if(fields[2] != "exon"){
            continue;
        }
        string chr{ fields[0]};
        int s = stoi(string{fields[3]});
        int e =  stoi(string{fields[4]});

        string st{ fields[6]};

        string info_str{fields[8]};
        lazy_split info(info_str, ";", " \n\t");


        string gene_id = ""; 
        string transcript_id = ""; 
        string exon_id = ""; 
        int count = 0;
        for( auto iter = info.begin(); iter != info.end(); ++iter){
            if((*iter).size() <= 1){ continue;}
            string f{*iter};
            lazy_split fs(f , " ", "\"");

            if(fs[0] == "gene_id"){
                gene_id = f[1];
                ++ count;
            }
            else if(fs[0] == "transcript_id"){
                transcript_id = f[1];
                ++ count;
            }
            else if(fs[0] == "exon_id"){
                exon_id = f[1];
                ++ count;
            }
            if(count == 3){
                break;
            }
        }

        annots.emplace_back(chr, s, e, st, gene_id, transcript_id, exon_id);
    }
    return annots;
}

IITree<int, size_t> make_exon_tree( vector<annotation> &annots){

    IITree<int, size_t> itree;

    size_t index = 0;
    for(annotation &annot: annots){
        itree.add(annot.start, annot.end, index);
        ++index;
    }
    itree.index();

    return itree;
}

vector<read> convert_aligs_to_segments(const string &path_to_aligs, int max_skip = 25 ){
    vector<read> reads;

    std::ifstream file( path_to_aligs);
    string str;

    //0.180525.0      1636    1192    1592    +       18      80373285        3254010 3256236 56      401     60      tp:A:P  mm:i:2  gn:i:343        go:i:3  cg:Z:17M1D3M341I4M2167N20M1I14M
    while( std::getline(file, str)){
        reads.emplace_back(str, max_skip);
    }
    return reads;
}


int main(int argc, char **argv){

    string path_to_aligs = argv[1];
    string path_to_gtf   = argv[2];

    // 0. Read alignments
    // 1. Convert to segments
/*
    vector<read> reads = convert_aligs_to_segments(path_to_aligs);   

    // 2. Read GTF file

    vector<annotation> gtf_annot = read_gtf_file(path_to_gtf);

    IITree<int, size_t> exon_tree = make_exon_tree(gtf_annot);

    // 3. Annotate segments
*/

    // 4. Create Graph

    // 5. Create Nodes from GTF 

    // 6. Create and count edges using the segments

    // 7. Convert counts to probabilities




    //SPLIT TEST
    //
    //
    string stest = "Hi,  Hello, Howdy!";

    lazy_split ltest(stest, "  ", "!,");



    for( auto q : ltest){
        cout << q << "@\n";   
    }
    for( auto q : ltest.known_delimiters){
        cout << q << "\n";   
    }
    auto it = ltest.begin();
    for(; it != ltest.end(); ++it){

        cout << it.index << " - \n";
    }
    cout << it.index << "\t" << it.end << "\t" << "\n";

    cout << ltest[0] << ".\n";
    cout << ltest[1] << ".\n";
    cout << ltest[2] << ".\n";

    cout << "Q " << (ltest.begin() == ltest.end()) << "\n";


    /*
       for( auto ite = ltest.begin(); ite == ltest.end(); ++ite){
       cout << *ite << "__\n";   
       }


       for( auto q : ltest){
       cout << q << "!\n";   
       }
       for( auto q : ltest.known_delimiters){
       cout << q << "\n";   
       }

*/
    //Graph Test
    graph<string, double> test;

    test.add("Test1");
    test.add("Test2");
    test.add("Test3");


    ++test.arc(0,1);
    ++test.arc(0,1);
    ++test.arc(0,1);
    ++test.arc(0,2);

    ++test.arc(1,2);
    ++test.arc(1,2);
    ++test.arc(1,2);

    auto nei = test.neighbours(1);

    for(auto i : nei){
        cout << i.first << "\t" << i.second << "\n";
    }

    return 0;
}
