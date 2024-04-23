
#pragma once
#ifndef INTERVAL_H
#define INTERVAL_H

#include <fmt/ostream.h>



#include <string>
#include <vector>
#include <iosfwd>
#include <map>


using std::string;
using std::vector;
using std::ostream;

class interval {
public:
    int start;
    int end;

    interval();
    interval(int start, int end);
    interval(const interval &i1, const interval &i2);
    int distance(const interval &other) const;
    int overlap(const interval &other) const;
    friend int larger_interval(const interval &i1, const interval &i2);
    double reciprocal(const interval &other) const;
    constexpr bool contains(int pos) const;
    constexpr int size() const { return end - start; };
};

class contig_str : public string {
    bool number;

public:

    static bool is_number(const string &str);
    contig_str(const string &str);
    contig_str();
    bool operator<(const contig_str &other) const;
    bool operator==(const contig_str &other) const = default;
    bool operator==(const string &other) const;
};

class intra_event : public interval {
public:
    contig_str chr;
    string layout;
    intra_event() : interval(0, 0), chr(""), layout("+") {}
    intra_event(string chr, int start, int end, const string &strand)
        : interval(start, end), chr(chr), layout(strand) {}
    intra_event(const intra_event &g1, const intra_event &g2)
        : interval(g1, g2), chr(g1.chr), layout(g1.layout + g2.layout) {}  // Merge constructor
    intra_event(const intra_event &g1) = default;
    virtual ~intra_event() {}
    int overlap(const intra_event &other) const {
        if (chr != other.chr) {
            return 0;
        }
        return interval::overlap(other);
    }
    double reciprocal(const intra_event &other) const {
        return static_cast<double>(overlap(other)) / larger_interval(*this, other);
    }
    bool operator<(const intra_event &other) const {
        if (other.chr != chr) {
            return chr < other.chr;
        }
        if (start == other.start) {
            return end < other.end;
        }
        return start < other.start;
        //}
    }
    bool operator==(const intra_event &other) const {
        return chr == other.chr && start == other.start && end == other.end && layout == other.layout;
    }
};

class ginterval : public interval {
public:
    contig_str chr;
    bool plus_strand;
    ginterval();
    ginterval(string chr, int start, int end, const string &strand);
    ginterval(string chr, int start, int end, bool plus_strand);
    ginterval(const ginterval &g1, const ginterval &g2);
    ginterval(const ginterval &g1) = default;
    virtual ~ginterval();
    int overlap(const ginterval &other) const;
    double reciprocal(const ginterval &other) const;
    bool operator<(const ginterval &other) const;
    bool operator==(const ginterval &other) const;
};

struct gtf : public ginterval {
    enum class entry_type {
        gene,
        transcript,
        exon,
        five_prime_utr,
        three_prime_utr,
        start_codon,
        stop_codon,
        CDS,
        Selenocysteine,
        other
    };
    static string type_to_string(entry_type type);
    static entry_type type_from_string(const string &type_str);
    entry_type type;
    std::map<string, string> info;
    string source;
    gtf(const string &gtf_line);
    gtf(const gtf &other);

    gtf(const ginterval &pos, entry_type type, const std::map<string, string> &info);
    gtf(const ginterval &pos, entry_type type);
    gtf();
    friend ostream &operator<<(ostream &os, const gtf &g);
    string to_string() const;
};

class transcript : public gtf {
    double abundance;
    vector<gtf> exons;
    string comment;

public:
    transcript(const transcript &other);
    transcript(const gtf &entry);
    transcript(const gtf &entry, double abundance);
    transcript(const gtf &entry, double abundance, const string &comment);

    friend ostream &operator<<(ostream &os, const transcript &t);
    string to_abundance_str() const;
    string to_string() const;
    void add_exon(const gtf &exon);

    void prepend_exon(const gtf &exon);

    void add_exons(const vector<gtf> &exons);
    void sort_exons();
    void set_abundance(double abundance);
    double get_abundance() const;
    vector<gtf> &get_exons();
    const vector<gtf> &cget_exons() const;
    auto get_exons(int start, int end);
    string get_comment() const;
    auto operator==(const transcript &other) const -> bool;
};


inline std::ostream &
operator<<(std::ostream &ost, const interval &ex) {
    ost << ex.start << "-" << ex.end;
    return ost;
}

inline std::ostream &
operator<<(std::ostream &ost, const intra_event &ex) {
    ost << ex.chr << ":" << ex.start << "-" << ex.end << "|" << ex.layout;
    return ost;
}
inline std::ostream &
operator<<(std::ostream &ost, const ginterval &ex) {
    ost << ex.chr << ":" << ex.start << "-" << ex.end;
    return ost;
}

template <>
struct fmt::formatter<interval> : ostream_formatter {};

template <>
struct fmt::formatter<ginterval> : ostream_formatter {};
template <>
struct fmt::formatter<transcript> : ostream_formatter {};

struct segment {
    ginterval tmplt;
    interval query;
    segment(string chr, int start, int end, int s2, int e2, const string &strand)
        : tmplt(chr, s2, e2, strand), query(start, end) {}
};

struct mapping {
    string rid;
    vector<segment> segments;
    bool primary;
    mapping();
    mapping(const string &paf, int max_skip, int min_segment) ;
};


struct base_mod {
    int position;
    char base;
    base_mod(int position, char base) : position(position), base(base) {}
    base_mod(const std::pair<int, char> &error) : position(error.first), base(error.second) {}
    base_mod() : position(-1), base('N') {}
    operator string() const { return std::to_string(position) + base; }
};

class einterval : public ginterval {
    vector<base_mod> errors;

public:
    einterval();
    einterval(const ginterval &gi);

    einterval(string chr, int start, int end, const bool plus_strand);
    einterval(string chr, int start, int end, const string &strand);

    einterval(const einterval &ei, int start, int end);
    void add_error(int position, char base) ;

    void add_error(const std::pair<int, char> &error);
    void add_error(const base_mod &error) ;
    void add_errors(const vector<base_mod> &err);
    void truncate(int start, int end) ;
    void parse_and_add_errors(const string &error_string) ;
    void sort_errors() ;
    auto error_str() const -> string ;
};

struct molecule_descriptor {
    string _id;
    bool _reversed;
    int _depth;
    vector<einterval> _segments;

    std::map<string, vector<string>> meta;
    //    vector<std::pair<int, char>> errors_so_far;

public:
    molecule_descriptor() ;

    molecule_descriptor(const string &id, bool reversed);
    molecule_descriptor(const molecule_descriptor &other);
    molecule_descriptor(const transcript &trans);

    auto get_id() const -> string;
    auto get_depth() const -> int;
    auto cget_segments() const -> const vector<einterval> &;
    auto get_segments() -> vector<einterval> &;
    molecule_descriptor *assign_segments(const vector<einterval> &segments) ; 
    molecule_descriptor *id(const string &id) ;

    // Assumes that special characters are not used in the meta data (',', ';', '=')
    molecule_descriptor *comment(const string &comment) ;
    molecule_descriptor *add_comment(const string &key, const string &value);
    molecule_descriptor *drop_comment(const string &key) ;
    // TODO add key verification
    const vector<string> &get_comment(const string &key) const ;

    molecule_descriptor *depth(int depth) ;
    molecule_descriptor *prepend_segment(const einterval &i) ;
    molecule_descriptor *append_segment(const einterval &i) ;
    molecule_descriptor *add_error(base_mod error) ;
    size_t size() const ;
    string dump_comment() const ;
    molecule_descriptor *concat(const molecule_descriptor &other) ;
    friend std::ostream &operator<<(std::ostream &ost, const molecule_descriptor &md) ;
};

molecule_descriptor
flip_molecule(const molecule_descriptor &md) ;
#endif
