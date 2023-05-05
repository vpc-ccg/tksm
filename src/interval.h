
#pragma once
#include <memory>
#ifndef INTERVAL_H
#define INTERVAL_H

#include <fmt/ostream.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <numeric>
#include <ranges>
#include <sstream>
#include <string>

#include "cigar.h"
#include "util.h"

using std::string;
using std::vector;
class interval {
public:
    int start;
    int end;

    interval() : start(0), end(0) {}
    interval(int start, int end) : start(start), end(end) {}
    interval(const interval &i1, const interval &i2)
        : start(std::min(i1.start, i2.start)), end(std::max(i1.end, i2.end)) {}
    int distance(const interval &other) const {
        if (other.start < this->start) {
            return other.distance(*this);
        }
        return other.start - this->end;
    }
    int overlap(const interval &other) const {
        if (other.end <= start) {  // BEFORE
            return 0;
        }
        else if (other.start >= end) {  // AFTER
            return 0;
        }
        else if (other.start >= start && other.end <= end) {  // IN
            return other.end - other.start;
        }
        else if (other.start < start && other.end > end) {  // AROUND
            return end - start;
        }
        else if (other.start < start && other.end < end && other.end > start) {  // LEFT OVERLAP
            return other.end - start;
        }
        else if (other.start > start && other.start < end && other.end > end) {  // RIGHT OVERLAP
            return end - other.start;
        }
        return 0;
    }
    friend int larger_interval(const interval &i1, const interval &i2) {
        return std::max(i1.end - i1.start, i2.end - i2.start);
    }
    double reciprocal(const interval &other) const {
        return static_cast<double>(overlap(other)) / larger_interval(*this, other);
    }
    bool contains(int pos) const { return pos > start && pos < end; }

    int size() const { return end - start; }
};

class contig_str : public string {
    bool number;

public:
    static bool is_number(const string &str) {
        string _str = str;
        if (_str.empty()) {
            return false;
        }
        if (_str.find("chr") == 0) {
            _str = str.substr(3);
        }
        return std::find_if(_str.begin(), _str.end(), [](unsigned char c) { return !std::isdigit(c); }) == _str.end();
    }

    contig_str(const string &str) : string{str}, number{is_number(str)} {}
    contig_str() : string{} {}
    bool operator<(const contig_str &other) const {
        if (number && other.number) {
            return stoi(*this) < stoi(other);
        }
        else if (number) {
            return true;
        }
        else if (other.number) {
            return false;
        }
        else {
            return static_cast<string>(*this) < static_cast<string>(other);
        }
    }

    bool operator==(const contig_str &other) const = default;
    bool operator==(const string &other) const { return *this == static_cast<contig_str>(other); }
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
    ginterval() : interval(0, 0), chr(""), plus_strand(true) {}

    ginterval(string chr, int start, int end, const string &strand)
        : interval(start, end), chr(chr), plus_strand(strand == "+") {}
    ginterval(string chr, int start, int end, bool plus_strand)
        : interval(start, end), chr(chr), plus_strand(plus_strand) {}
    ginterval(const ginterval &g1, const ginterval &g2)
        : interval(g1, g2), chr(g1.chr), plus_strand(g1.plus_strand) {}  // Merge constructor
    ginterval(const ginterval &g1) = default;
    virtual ~ginterval() {}
    int overlap(const ginterval &other) const {
        if (chr != other.chr) {
            return 0;
        }
        return interval::overlap(other);
    }
    double reciprocal(const ginterval &other) const {
        return static_cast<double>(overlap(other)) / larger_interval(*this, other);
    }
    bool operator<(const ginterval &other) const {
        if (other.chr != chr) {
            return chr < other.chr;
        }
        if (start == other.start) {
            return end < other.end;
        }
        return start < other.start;
        //}
    }
    bool operator==(const ginterval &other) const {
        return chr == other.chr && start == other.start && end == other.end && plus_strand == other.plus_strand;
    }
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
    static entry_type type_from_string(const string &type_str) {
        if (type_str == "gene") {
            return entry_type::gene;
        }
        else if (type_str == "transcript") {
            return entry_type::transcript;
        }
        else if (type_str == "exon") {
            return entry_type::exon;
        }
        else if (type_str == "five_prime_utr") {
            return entry_type::five_prime_utr;
        }

        else if (type_str == "three_prime_utr") {
            return entry_type::three_prime_utr;
        }
        else if (type_str == "start_codon") {
            return entry_type::start_codon;
        }
        else if (type_str == "stop_codon") {
            return entry_type::stop_codon;
        }
        else if (type_str == "CDS") {
            return entry_type::CDS;
        }
        else if (type_str == "Selenocysteine") {
            return entry_type::Selenocysteine;
        }
        else {
            return entry_type::other;
        }
    }
    static string type_to_string(entry_type type) {
        switch (type) {
            case entry_type::gene:
                return "gene";
            case entry_type::transcript:
                return "transcript";
            case entry_type::exon:
                return "exon";
            case entry_type::five_prime_utr:
                return "five_prime_utr";
            case entry_type::three_prime_utr:
                return "three_prime_utr";
            case entry_type::start_codon:
                return "start_codon";
            case entry_type::stop_codon:
                return "stop_codon";
            case entry_type::CDS:
                return "CDS";
            case entry_type::Selenocysteine:
                return "Selenocysteine";
            case entry_type::other:
                return "other";
            default:
                return "other";
        }
    }
    entry_type type;
    std::map<string, string> info;
    string source;
    gtf(const string &gtf_line) {
        vector<string> fields = rsplit(gtf_line, "\t");
        type                  = type_from_string(fields[2]);
        chr                   = fields[0];
        source                = fields[1];
        start                 = stoi(fields[3]) - 1;  // GTF is 1-based
        end                   = stoi(fields[4]);
        plus_strand           = (fields[6] == "+");

        string info_str{fields[8]};
        vector<string> info_vec = rsplit(info_str, ";");
        strip_for_each(info_vec, " ");

        for (auto iter = info_vec.begin(); iter != info_vec.end(); ++iter) {
            if ((*iter).size() <= 1) {
                continue;
            }
            string f{*iter};

            vector<string> fs = rsplit(f, " ");
            strip_for_each(fs, "\"");
            info[fs[0]] = fs[1];
        }
    }

    gtf(const gtf &other) = default;

    gtf(const ginterval &pos, entry_type type, const std::map<string, string> &info)
        : ginterval(pos), type(type), info(info), source{"TKSM"} {}
    gtf(const ginterval &pos, entry_type type) : ginterval(pos), type(type), source{"TKSM"} {}

    gtf() = default;
    friend ostream &operator<<(ostream &os, const gtf &g) {
        os << g.chr << "\t" << g.source << "\t" << type_to_string(g.type) << "\t" << g.start << "\t" << g.end << "\t"
           << (g.plus_strand ? '+' : '-') << "\t"
           << "."
           << "\t";
        for (auto iter = g.info.begin(); iter != g.info.end(); ++iter) {
            os << iter->first << " \"" << iter->second << "\";";
        }
        os << "\n";
        return os;
    }
    string to_string() const {
        std::ostringstream os;
        os << *this;
        return os.str();
    }
};

class transcript : public gtf {
    double abundance;
    vector<gtf> exons;
    string comment;

public:
    transcript(const transcript &other) : gtf{other}, abundance{other.abundance}, exons{other.exons}, comment{other.comment}  {}
    transcript(const gtf &entry) : gtf(entry), abundance(0.0) {}
    transcript(const gtf &entry, double abundance) : gtf(entry), abundance(abundance) {}
    transcript(const gtf &entry, double abundance, const string &comment)
        : gtf(entry), abundance(abundance), comment(comment) {}

    friend ostream &operator<<(ostream &os, const transcript &t) {
        os << static_cast<const gtf &>(t);
        for (auto &exon : t.exons) {
            os << exon;
        }
        return os;
    }
    string to_abundance_str() const {
        std::ostringstream os;
        os << info.at("transcript_id") << "\t" << abundance << "\t" << comment;
        return os.str();
    }
    string to_string() const {
        std::ostringstream os;
        os << *this;
        return os.str();
    }
    void add_exon(const gtf &exon) { exons.push_back(exon); }

    void prepend_exon(const gtf &exon) { exons.insert(exons.begin(), exon); }

    void add_exons(const vector<gtf> &exons) {
        for (auto &exon : exons) {
            add_exon(exon);
        }
    }
    void sort_exons() {
        std::sort(exons.begin(), exons.end(), [](const gtf &a, const gtf &b) { return a.start < b.start; });
    }
    void set_abundance(double abundance) { this->abundance = abundance; }
    double get_abundance() const { return abundance; }
    vector<gtf> &get_exons() { return exons; }
    const vector<gtf> &cget_exons() const { return exons; }
    auto get_exons(int start, int end) {
        return exons | std::ranges::views::filter([start, end](const gtf &exon) -> bool {
                   return exon.start >= start && exon.end <= end;
               });  // clangd ignore
    }

    string get_comment() const { return comment; }

    auto operator==(const transcript &other) const -> bool {
        return exons == other.exons && info.at("transcript_id") == other.info.at("transcript_id") &&
               info.at("gene_id") == other.info.at("gene_id") && info.at("gene_name") == other.info.at("gene_name") &&
               start == other.start && end == other.end && chr == other.chr && plus_strand == other.plus_strand;
    }
};

/*
struct gene : public ginterval {
    string gene_id;
    string gene_name;
    gene() : ginterval(), gene_id("NAN"), gene_name("NAN") {}
    gene(const gtf &entry)
        : ginterval(entry), gene_id(entry.info.at("gene_id")), gene_name(entry.info.at("gene_name")) {}
    gene(const string &id) : gene_id(id) {}  // Mock constructor for map access
    gene(string chr, int start, int end, const string &strand, const string &gene_id, const string &gene_name)
        : ginterval(chr, start, end, strand), gene_id(gene_id), gene_name(gene_name) {}

    friend auto operator==(const gene &a, const gene &b) { return a.gene_id == b.gene_id; }
    friend auto operator!=(const gene &a, const gene &b) { return !(a == b); }

    friend auto operator>(const gene &a, const gene &b) { return a.gene_id > b.gene_id; }
    friend auto operator>=(const gene &a, const gene &b) { return a.gene_id >= b.gene_id; }
    friend auto operator<(const gene &a, const gene &b) { return a.gene_id < b.gene_id; }
    friend auto operator<=(const gene &a, const gene &b) { return a.gene_id <= b.gene_id; }
};

namespace std {
template <>
struct hash<gene> {
    std::size_t operator()(gene const &s) const noexcept { return std::hash<string>{}(s.gene_id); }
};
}  // namespace std

struct transcript : public ginterval {
    string transcript_id;
    gene *gene_ref;
    transcript() : transcript_id{"NULL"}, gene_ref{nullptr} {}
    transcript(const gtf &entry, gene *gref)
        : ginterval(entry), transcript_id(entry.info.at("transcript_id")), gene_ref{gref} {}
    transcript(string chr, int start, int end, const string &strand, const string &transcript_id, gene *gref)
        : ginterval(chr, start, end, strand), transcript_id(transcript_id), gene_ref(gref) {}
};

struct exon : public ginterval {
    string exon_id;
    string strand;
    string transcript_id;

    virtual ~exon() {}
    gene gene_ref;

    exon(string chr, int start, int end, const string &strand, const string &exon_id, const string &transcript_id,
         const gene &g)
        : ginterval(chr, start, end, strand), exon_id(exon_id), transcript_id(transcript_id), gene_ref(g) {}
    exon(string chr, int start, int end, const string &strand, const string &exon_id, const string &transcript_id,
         gene *gref)
        : ginterval(chr, start, end, strand), exon_id(exon_id), transcript_id(transcript_id), gene_ref(*gref) {}
    exon(const gtf &entry, gene *gref)
        : ginterval(entry),
          exon_id(entry.info.at("exon_id")),
          transcript_id(entry.info.at("transcript_id")),
          gene_ref{*gref} {}
    exon() : exon_id("NULL") {}

    exon(const exon &g1, const exon &g2, const string &id)
        : ginterval(g1, g2), exon_id(id), strand(g1.strand), transcript_id(g1.transcript_id), gene_ref(g1.gene_ref) {}
    exon(const exon &g1, const exon &g2)
        : ginterval(g1, g2),
          exon_id(g1.exon_id),
          strand(g1.strand),
          transcript_id(g1.transcript_id),
          gene_ref(g1.gene_ref) {}  // Merge constructor

    bool operator<(const exon &other) const {
        if (gene_ref.chr != other.gene_ref.chr) {
            return gene_ref.chr < other.gene_ref.chr;
        }
        if (gene_ref != other.gene_ref) {
            return gene_ref < other.gene_ref;
        }

        if (start == other.start) {
            return end < other.end;
        }
        return start < other.start;
        //}
    }
    bool operator==(const exon &other) const { return exon_id == other.exon_id; }
};
*/
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

/*
inline std::ostream &
operator<<(std::ostream &ost, const exon &ex) {
    ost << dynamic_cast<const ginterval &>(ex) << " " << ex.strand << " " << ex.exon_id;

    ost << " " << ex.gene_ref.gene_name << " " << (ex.gene_ref.plus_strand ? "+" : "-");

    return ost;
}
inline std::ostream &
operator<<(std::ostream &ost, const gene &ex) {
    ost << (ginterval)ex << " " << ex.gene_name;
    return ost;
}
*/
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
    mapping() : rid("-1") {}
    mapping(const string &paf, int max_skip, int min_segment) {
        std::istringstream ps(paf);

        ps >> rid;

        int rlen;
        int rstart;
        int rend;
        ps >> rlen;
        ps >> rstart;
        ps >> rend;

        string strand;
        ps >> strand;
        bool complemented = strand == "-";
        string chr;
        ps >> chr;

        if (chr.find("chr") != string::npos) {
            chr = chr.substr(3);
        }

        int template_len;
        int template_start;
        int template_end;
        int num_matches;
        int alig_block_len;
        int mapq;
        string field;
        ps >> template_len >> template_start >> template_end >> num_matches >> alig_block_len >> mapq;

        ps >> field;
        string cigar_str = "-1";

        while (!ps.eof()) {
            ps >> field;
            if (field.substr(0, 2) == "tp") {
                primary = (field == "tp:A:P");
            }
            if (field.substr(0, 2) == "cg") {
                cigar_str = field.substr(5);
            }
        }

        if (cigar_str == "-1") {
            std::cerr << "Paf line doesn't include alignment cigar! Exiting!.\n";
            exit(-1);
        }
        cigar cig(cigar_str);

        vector<segment> aligs;

        int st = template_start;
        int sq = rstart;
        int et, eq;
        if (complemented) {
            sq = rend - 1;
        }
        for (auto pair : cig) {
            int length       = pair.first;
            char c           = pair.second;
            cigar::type type = cigar::what_is_this(c);

            if (type == cigar::type::matched) {
                et = st + length;

                if (complemented) {
                    eq = sq - length;
                    aligs.emplace_back(chr, eq, sq - 1, st, et - 1, strand);
                }
                else {
                    eq = sq + length;
                    aligs.emplace_back(chr, sq, eq - 1, st, et - 1, strand);
                }

                sq = eq;
                st = et;
            }
            else if (type == cigar::type::ontemplate) {
                st = st + length;
            }
            else if (type == cigar::type::onquery) {
                if (complemented) {
                    sq = sq - length;
                }
                else {
                    sq = sq + length;
                }
            }
            else {
            }
        }
        if (aligs.size() == 0) {
            std::cerr << "NOALIG " << rid << "\n";
        }
        if (max_skip > 0 && aligs.size() > 0) {
            st = aligs[0].tmplt.start;
            et = aligs[0].tmplt.end;
            sq = aligs[0].query.start;
            eq = aligs[0].query.end;

            for (auto iter = std::begin(aligs); std::next(iter) != std::end(aligs); ++iter) {
                if (complemented) {
                    if (std::next(iter)->tmplt.start - et < max_skip) {
                        sq = std::next(iter)->query.start;
                        et = std::next(iter)->tmplt.end;
                    }
                    else {
                        if (eq - sq >= min_segment && et - st >= min_segment) {
                            segments.emplace_back(chr, sq, eq, st, et, strand);
                        }
                        auto nxt = std::next(iter);
                        st       = nxt->tmplt.start;
                        et       = nxt->tmplt.end;
                        sq       = nxt->query.start;
                        eq       = nxt->query.end;
                    }
                }
                else {
                    if (std::next(iter)->tmplt.start - et < max_skip) {
                        et = std::next(iter)->tmplt.end;
                        eq = std::next(iter)->query.end;
                    }
                    else {
                        if (eq - sq >= min_segment && et - st >= min_segment) {
                            segments.emplace_back(chr, sq, eq, st, et, strand);
                        }
                        auto nxt = std::next(iter);
                        st       = nxt->tmplt.start;
                        et       = nxt->tmplt.end;
                        sq       = nxt->query.start;
                        eq       = nxt->query.end;
                    }
                }
            }
            if (eq - sq >= min_segment && et - st >= min_segment) {
                segments.emplace_back(chr, sq, eq, st, et, strand);
            }
        }
        else {
            segments = aligs;
        }
    }
};
/*
class isoform {
public:
    vector<exon> segments;
    int depth;
    string gene;
    string transcript_id;

    isoform() : segments({}), depth(1), gene("NULL"), transcript_id("NULL") {}
    isoform(const isoform &other)
        : segments(other.segments), depth(other.depth), gene(other.gene), transcript_id(other.transcript_id) {}
    isoform(const vector<exon> &segs) : segments(segs), depth(1), gene("NULL"), transcript_id("NULL") {}
    isoform(const vector<exon> &segs, int depth) : segments(segs), depth(depth), gene("NULL"), transcript_id("NULL") {}
    isoform(const vector<exon> &segs, int depth, const string &gene)
        : segments(segs), depth(depth), gene(gene), transcript_id("NULL") {}
    isoform(const vector<exon> &segs, int depth, const string &gene, const string &tid)
        : segments(segs), depth(depth), gene(gene), transcript_id(tid) {}

    bool operator<(const isoform &other) const {
        if (other.gene != gene) {
            return gene < other.gene;
        }

        return segments < other.segments;
        //}
    }
    bool operator==(const isoform &other) const { return gene == other.gene && segments == other.segments; }
    isoform &operator=(const vector<exon> &segs) {
        this->segments = segs;
        return *this;
    }
};
*/
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
    einterval() : ginterval() {}
    einterval(const ginterval &gi) : ginterval(gi) {}
    einterval(string chr, int start, int end, const string &strand) : ginterval(chr, start, end, strand) {}

    void add_error(int position, char base) { errors.emplace_back(position, base); }

    void add_error(const std::pair<int, char> &error) { errors.emplace_back(error.first, error.second); }
    void add_error(const base_mod &error) { errors.emplace_back(error); }

    void add_errors(const vector<base_mod> &err) { errors.insert(errors.end(), err.begin(), err.end()); }

    void truncate(int start, int end) {
        sort_errors();
        assert(start >= 0);
        assert(end <= this->end - this->start);
        assert(start < end);

        this->start += start;
        this->end = this->start + (end - start);

        if (errors.size() == 0) {
            return;
        }

        errors.erase(std::remove_if(errors.begin(), errors.end(),
                                    [start, end](const base_mod &error) {
                                        return error.position < start || error.position >= end;
                                    }),
                     errors.end());
    }
    void parse_and_add_errors(const string &error_string) {
        auto mutations = rsplit(error_string, ",");
        for (string mutation : mutations) {
            if (mutation == "") {
                continue;
            }
            char target      = mutation[mutation.size() - 1];
            int mutation_pos = stoi(mutation.substr(0, mutation.size() - 1));
            add_error(mutation_pos, target);
        }
    }

    void sort_errors() {
        std::sort(errors.begin(), errors.end(),
                  [](const base_mod &a, const base_mod &b) { return a.position < b.position; });
    }

    auto error_str() const -> string { return join_str(errors.begin(), errors.end(), ","); }
};

struct molecule_descriptor {
    string _id;
    bool _reversed;
    int _depth;
    vector<einterval> _segments;

    std::map<string, vector<string>> meta;
    //    vector<std::pair<int, char>> errors_so_far;

public:
    molecule_descriptor() : _id(""), _reversed{false}, _depth{0}{}

    molecule_descriptor(const string &id, bool reversed) : _id(id), _reversed(reversed), _depth(1) {}

    molecule_descriptor(const molecule_descriptor &other)
        : _id(other._id),
          _reversed(other._reversed),
          _depth(other._depth),
          _segments(other._segments.begin(), other._segments.end()),
          meta{other.meta} {}

    molecule_descriptor(const transcript &trans) : 
        _id(trans.info.at("transcript_id")),
        _reversed(!trans.plus_strand),
        _depth(trans.get_abundance())
    {
        for(const gtf &gi : trans.cget_exons()){
            append_segment(gi);
        }
    }

    /*
    molecule_descriptor(const isoform &iso)
        : _id(iso.transcript_id), _reversed(!iso.segments.back().plus_strand), _depth(iso.depth) {
        for (const exon &e : iso.segments) {
            _segments.push_back(e);
        }
    }
*/
    auto get_id() const { return _id; }

    auto get_depth() const { return _depth; }

    const auto &cget_segments() const { return _segments; }
    auto &get_segments() { return _segments; }
    molecule_descriptor *assign_segments(const vector<einterval> &segments) {
        _segments = segments;
        return this;
    }
    molecule_descriptor *id(const string &id) {
        _id = id;
        return this;
    }

    // Assumes that special characters are not used in the meta data (',', ';', '=')
    molecule_descriptor *comment(const string &comment) {
        vector<std::string_view> fields = splitSV(comment, ";");
        for (const auto &f : fields) {
            if (f.find_first_of("=") == std::string_view::npos) {
                add_comment(string{f}, ".");
            }
            else {
                vector<std::string_view> kv     = splitSV(f, "=");
                const auto &key                 = kv[0];
                const auto &values_string       = kv[1];
                vector<std::string_view> values = splitSV(values_string, ",");
                for (const auto &v : values) {
                    add_comment(string{key}, string{v});
                }
            }
        }
        // TODO
        return this;
    }

    molecule_descriptor *add_comment(const string &key, const string &value) {
        meta[key].push_back(value);
        return this;
    }

    molecule_descriptor *drop_comment(const string &key) {
        meta.erase(key);
        return this;
    }

    // TODO add key verification
    const vector<string> &get_comment(const string &key) const { return meta.at(key); }

    molecule_descriptor *depth(int depth) {
        _depth = depth;
        return this;
    }

    molecule_descriptor *prepend_segment(const einterval &i) {
        _segments.insert(_segments.begin(), i);
        // for (std::pair<int, char> &errs : errors_so_far) {
        //    errs = std::make_pair(errs.first + i.end - i.start, errs.second);
        // }
        return this;
    }

    molecule_descriptor *append_segment(const ginterval &i) {
        _segments.push_back(i);
        return this;
    }
    /*
        molecule_descriptor *update_errors(const vector<std::pair<int, char>> &errors_so_far) {
            this->errors_so_far = errors_so_far;
            return this;
        }
    */

    molecule_descriptor *add_error(base_mod error) {
        auto iter = _segments.begin();
        while (iter->size() <= error.position) {
            error.position -= iter->size();
            iter++;
        }
        iter->add_error(error);
        return this;
    }

    size_t size() const {
        return std::accumulate(
            _segments.begin(), _segments.end(), 0L,
            [](size_t sum_so_far, const ginterval &g) -> size_t { return sum_so_far + g.end - g.start; });
    }

    string dump_comment() const {
        std::ostringstream buffer;
        for (const auto &kv : meta) {
            buffer << kv.first;
            if (kv.second[0] != ".") {
                buffer << "=" << (join_str(kv.second.cbegin(), kv.second.cend(), ","));
            }
            buffer << ";";
        }
        return buffer.str();
    }
    
    molecule_descriptor *concat(const molecule_descriptor &other) {
        _segments.insert(_segments.end(), other._segments.begin(), other._segments.end());
        return this;
    }

    friend std::ostream &operator<<(std::ostream &ost, const molecule_descriptor &md) {
        print_tsv(ost, "+" + md._id, md._depth, md.dump_comment());

        for (const einterval &ival : md._segments) {
            print_tsv(ost, ival.chr, ival.start, ival.end, (ival.plus_strand ? "+" : "-"), ival.error_str());
        }
        return ost;
    }
};

inline
molecule_descriptor flip_molecule(const molecule_descriptor &md){
    molecule_descriptor flipped_md{md.get_id(), md._reversed};
    for(const auto &segment : md.cget_segments() | std::views::reverse){
        auto flipped_segment = segment;
        flipped_segment.plus_strand = !flipped_segment.plus_strand;
        flipped_md.append_segment(flipped_segment);
    }
    flipped_md.meta = md.meta;
    flipped_md.depth(md.get_depth());

    return flipped_md;
}
/*
// pcr copy structure that tracks pcr errors introduced
struct pcr_copy {
    string id;
    vector<ginterval> segments;
    vector<std::pair<int, char>> errors_so_far;
    bool reversed;
    int depth;
    string comment;

    pcr_copy() {}
    pcr_copy(const string &id) : id(id), depth(1) {}

    pcr_copy(const string &id, const isoform &iso) : id(id), depth(iso.depth) {
        for (const exon &e : iso.segments) {
            segments.push_back(e);
        }
    }
    pcr_copy(const isoform &iso) : id(iso.transcript_id), depth(iso.depth) {
        for (const exon &e : iso.segments) {
            segments.push_back(e);
        }
    }
    pcr_copy(const string &id, const vector<ginterval> &segments, const vector<std::pair<int, char>> &errors_so_far)
        : id(id), segments(segments), errors_so_far(errors_so_far), depth(1) {}
    pcr_copy(const string &id, const vector<ginterval> &segments, const vector<std::pair<int, char>> &errors_so_far,
             int depth)
        : id(id), segments(segments), errors_so_far(errors_so_far), depth(depth) {}
    pcr_copy(const vector<ginterval> &segments, const vector<std::pair<int, char>> &errors_so_far)
        : id("copy"), segments(segments), errors_so_far(errors_so_far), depth(1) {}
    pcr_copy(const vector<ginterval> &segments) : id("copy"), segments(segments), depth(1) {}

    void append(const ginterval &g) { segments.push_back(g); }
    void prepend(const ginterval &g) {
        segments.insert(segments.begin(), g);
        for (std::pair<int, char> &errs : errors_so_far) {
            errs = std::make_pair(errs.first + g.end - g.start, errs.second);
        }
    }
    size_t size() const {
        return std::accumulate(
            segments.begin(), segments.end(), 0L,
            [](size_t sum_so_far, const ginterval &g) -> size_t { return sum_so_far + g.end - g.start; });
    }
};

// pcr molecule structure that can model paired molecules
struct pcr_molecule {
    vector<molecule_descriptor> paired;

    pcr_molecule() {}

    pcr_molecule(const pcr_molecule &other) : paired(other.paired) {}

    pcr_molecule(const molecule_descriptor &other) { paired.push_back(other); }
    pcr_molecule(const string &id, const isoform &other) { paired.push_back(*molecule_descriptor{other}.id(id)); }
    pcr_molecule(const pcr_molecule &first, const pcr_molecule &second) : paired(first.paired) {
        paired.reserve(first.paired.size() + second.paired.size());
        paired.insert(paired.end(), second.paired.begin(), second.paired.end());
    }
};
*/
#endif
