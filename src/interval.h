#ifndef INTERVAL_H
#define INTERVAL_H

#include <algorithm>
#include <numeric>
#include <string>
#include <sstream>
#include <iostream>
#include <map>

#include "util.h"

#include "cigar.h"

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
        bool contains(int pos) const {
            return pos > start && pos < end;
        }
};


class ginterval: public interval{
    public:
        std::string chr;
        bool plus_strand;
        ginterval(): interval(0,0), chr(""), plus_strand(true){}
        ginterval(std::string chr, int start, int end, const std::string &strand): 
            interval(start, end), chr(chr), plus_strand(strand=="+"){}
        ginterval( const ginterval &g1, const ginterval &g2) : interval(g1,g2), chr(g1.chr), plus_strand(g1.plus_strand) {} // Merge constructor
        ginterval( const ginterval &g1) = default;
        virtual ~ginterval() {}
        int overlap( const ginterval &other) const{
            if( chr != other.chr){
                return 0;
            }
            return interval::overlap(other);
        } 
        double reciprocal(const ginterval &other) const {
            return static_cast<double>(overlap(other)) / larger_interval(*this, other);
        }
        bool operator<( const ginterval &other) const{
            if( other.chr != chr){
                return chr < other.chr;
            }
            if( start == other.start){
                return end < other.end;
            }
            return start < other.start;
            //}

        }
        bool operator==(const ginterval &other) const {
            return chr == other.chr && start == other.start && end == other.end && plus_strand == other.plus_strand;
        }
};

struct gtf: public ginterval{
    enum class entry_type{
        gene, transcript, exon,
        five_prime_utr, three_prime_utr,
        start_codon, stop_codon,
        CDS, Selenocysteine,
        other
    };
    static entry_type type_from_string( const std::string& type_str){

        if(type_str == "gene"){
            return entry_type::gene;
        }
        else if(type_str == "transcript"){
            return entry_type::transcript;
        }
        else if(type_str == "exon"){
            return entry_type::exon;
        }
        else if(type_str == "five_prime_utr"){
            return entry_type::five_prime_utr;
        }

        else if(type_str == "three_prime_utr"){
            return entry_type::three_prime_utr;
        }
        else if(type_str == "start_codon"){
            return entry_type::start_codon;
        }
        else if(type_str == "stop_codon"){
            return entry_type::stop_codon;
        }
        else if(type_str == "CDS"){
            return entry_type::CDS;
        }
        else if(type_str == "Selenocysteine"){
            return entry_type::Selenocysteine;
        }
        else{
            return entry_type::other;
        }
    }

    entry_type type;
    std::map<std::string, std::string> info;
    gtf( const std::string &gtf_line){
        std::vector<std::string> fields = rsplit(gtf_line, "\t");
        type = type_from_string(fields[2]);
        chr = fields[0];
        start = stoi(fields[3]);
        end   = stoi(fields[4]);
        plus_strand = (fields[6] == "+");
        
        std::string info_str{fields[8]};
        std::vector<std::string> info_vec = rsplit(info_str, ";");
        strip_for_each(info_vec , " ");

        for( auto iter = info_vec.begin(); iter != info_vec.end(); ++iter){
            if((*iter).size() <= 1){ continue;}
            std::string f{*iter};

            std::vector<std::string> fs = rsplit(f , " ");
            strip_for_each(fs, "\"");
            info[fs[0]] = fs[1];
        }
    }
};


struct gene: public ginterval{
    std::string gene_id;
    std::string gene_name;
    gene() : ginterval(),gene_id("NAN"), gene_name("NAN") {}
    gene(const gtf &entry) : ginterval(entry), gene_id(entry.info.at("gene_id")), gene_name(entry.info.at("gene_name")){}
    gene(const std::string &id) : gene_id(id) {} //Mock constructor for map access
    gene(std::string chr, int start, int end, const std::string &strand, const std::string &gene_id,
            const std::string &gene_name) : ginterval(chr, start, end, strand),
    gene_id(gene_id), gene_name(gene_name){}



    friend auto operator ==(const gene &a, const gene &b){
        return a.gene_id == b.gene_id;
    }
    friend auto operator !=(const gene &a, const gene &b){
        return ! (a==b);
    }

    friend auto operator > (const gene &a, const gene &b){
        return a.gene_id > b.gene_id;
    }
    friend auto operator >= (const gene &a, const gene &b){
        return a.gene_id >= b.gene_id;
    }
    friend auto operator < (const gene &a, const gene &b){
        return a.gene_id < b.gene_id;
    }
    friend auto operator <= (const gene &a, const gene &b){
        return a.gene_id <= b.gene_id;
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
    std::string transcript_id;
    gene* gene_ref;
    transcript() :transcript_id{"NULL"}, gene_ref{nullptr}{}
    transcript(const gtf &entry, gene *gref) : ginterval(entry), transcript_id(entry.info.at("transcript_id")), gene_ref{gref}{}
    transcript(std::string chr, int start, int end, const std::string &strand, const std::string &transcript_id, gene *gref):
        ginterval(chr, start, end, strand), transcript_id(transcript_id), gene_ref(gref) {}
};

struct exon: public ginterval{
    std::string exon_id;
    std::string strand;

    virtual ~exon() {}
    gene * gene_ref;

    exon(std::string chr, int start, int end, const std::string &strand, const std::string &exon_id, gene* gref): 
        ginterval(chr, start, end, strand), exon_id(exon_id), gene_ref(gref){}
    exon(const gtf &entry, gene *gref) : ginterval(entry), exon_id(entry.info.at("exon_id")), gene_ref{gref}{}
    exon() : exon_id("NULL") {}

    exon(const exon &g1, const exon &g2, const std::string &id) : ginterval(g1,g2), exon_id(id), strand(g1.strand), gene_ref(g1.gene_ref)  {
    }
    exon(const exon &g1, const exon &g2) : ginterval(g1,g2), exon_id(g1.exon_id), strand(g1.strand), gene_ref(g1.gene_ref){
    } // Merge constructor

    bool operator<( const exon &other) const{
        if(gene_ref != nullptr){
            if( gene_ref-> chr != other.gene_ref->chr){
                return gene_ref->chr < other.gene_ref->chr;
            }
            if( gene_ref != other.gene_ref){
                return *gene_ref < *other.gene_ref;
            }
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



inline std::ostream &operator << ( std::ostream &ost, const interval &ex){
    ost << ex.start << "-" << ex.end;
    return ost;
}

inline std::ostream &operator << ( std::ostream &ost, const ginterval &ex){
    ost << ex.chr << ":" << ex.start << "-" << ex.end;
    return ost;
}


inline std::ostream &operator << ( std::ostream &ost, const exon &ex){
    ost << dynamic_cast<const ginterval&>(ex) << " " << ex.strand << " " << ex.exon_id;
    if(ex.gene_ref != nullptr){
        ost << " " << ex.gene_ref->gene_name << " " << (ex.gene_ref->plus_strand?"+":"-");
    }
    return ost;
}
inline std::ostream &operator << ( std::ostream &ost, const gene &ex){
    ost << (ginterval) ex <<  " " << ex.gene_name ;
    return ost;
}


struct segment{
    ginterval tmplt;
    interval query;
    segment(std::string chr, int start, int end, int s2, int e2, const std::string &strand): 
        tmplt(chr, s2, e2, strand), query(start,end){}
};

struct mapping{
    std::string rid;
    std::vector<segment> segments;
    bool primary;
    mapping(): rid("-1"){}
    mapping(const std::string &paf, int max_skip, int min_segment){

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
        std::string cigar_str = "-1";

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
            std::cerr << "Paf line doesn't include alignment cigar! Exiting!.\n";
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
            cigar::type type = cigar::what_is_this(c);

            if( type == cigar::type::matched){
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
            else if( type == cigar::type::ontemplate){
                st = st + length;
            }
            else if( type == cigar::type::onquery){
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
            std::cerr << "NOALIG " << rid << "\n";
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


class isoform{
    public:
        std::vector<exon> segments;
        int depth;
        std::string gene;

    isoform(std::vector<exon> &segs) : segments(segs), depth(1), gene("NULL"){}
    isoform(std::vector<exon> &segs, int depth) : segments(segs), depth(depth), gene("NULL"){}
    isoform(std::vector<exon> &segs, int depth, const std::string &gene) : segments(segs), depth(depth), gene(gene){}

};

//pcr copy structure that tracks pcr errors introduced
struct pcr_copy{
    std::string id;
    std::vector<ginterval> segments;
    std::vector< std::pair< int, char>> errors_so_far;
    bool reversed;
    int depth;
    std::string comment;

    pcr_copy(){}
    pcr_copy( const std::string &id) : id(id), depth(1){}
    pcr_copy( const std::string &id, const isoform &iso) : id(id), depth(1){
        for(const exon &e :iso.segments){
            segments.push_back(e);
        }
    }
    pcr_copy( const std::string &id, const std::vector<ginterval> &segments, const std::vector< std::pair< int, char>> &errors_so_far) : id(id), segments(segments), errors_so_far(errors_so_far), depth(1) {}
    pcr_copy( const std::string &id, const std::vector<ginterval> &segments, const std::vector< std::pair< int, char>> &errors_so_far, int depth) : id(id), segments(segments), errors_so_far(errors_so_far), depth(depth) {}
    pcr_copy( const std::vector<ginterval> &segments, const std::vector< std::pair< int, char>> &errors_so_far) : id("copy"), segments(segments), errors_so_far(errors_so_far), depth(1) {}
    pcr_copy( const std::vector<ginterval> &segments) : id("copy"), segments(segments), depth(1) {}

    void prepend(const ginterval &g){
        segments.insert(segments.begin(),g);
        for(std::pair<int, char> &errs : errors_so_far){
            errs = std::make_pair(errs.first+g.end-g.start,errs.second);
        }
    }
    size_t size() const{
        return std::accumulate(segments.begin(), segments.end(), 0L, [] (size_t sum_so_far, const ginterval &g) -> size_t { return sum_so_far + g.end - g.start;});
    }
};

//pcr molecule structure that can model paired molecules
struct pcr_molecule{
    std::vector<pcr_copy> paired;

    pcr_molecule() {}

    pcr_molecule(const pcr_molecule &other) : paired(other.paired) {}

    pcr_molecule(const pcr_copy &other) {
        paired.push_back(other);
    }
    pcr_molecule(const std::string &id, const isoform &other) {
        paired.push_back(pcr_copy{id, other});

    }
    pcr_molecule(const pcr_molecule &first, const pcr_molecule &second) : paired(first.paired) {
        paired.reserve(first.paired.size() + second.paired.size());
        paired.insert(paired.end(), second.paired.begin(), second.paired.end());
    }
};

#endif
