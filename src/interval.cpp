#include "interval.h"

#include <fmt/ostream.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <ranges>
#include <sstream>
#include <string>

#include "cigar.h"
#include "util.h"

using std::string;
using std::vector;

template <class B>
inline void
print_tsv(std::ostream &ost, B b) {
    ost << b << "\n";
}

template <class B, class... A>
inline void
print_tsv(std::ostream &ost, B b, A... a) {
    ost << b << "\t";
    print_tsv(ost, a...);
}
interval::interval() : start(0), end(0) {}
interval::interval(int start, int end) : start(start), end(end) {}
interval::interval(const interval &i1, const interval &i2)
    : start(std::min(i1.start, i2.start)), end(std::max(i1.end, i2.end)) {}
int
interval::distance(const interval &other) const {
    if (other.start < this->start) {
        return other.distance(*this);
    }
    return other.start - this->end;
}
int
interval::overlap(const interval &other) const {
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
int
larger_interval(const interval &i1, const interval &i2) {
    return std::max(i1.size(), i2.size());
}
double
interval::reciprocal(const interval &other) const {
    return static_cast<double>(overlap(other)) / larger_interval(*this, other);
}
constexpr bool
interval::contains(int pos) const {
    return pos > start && pos < end;
}

contig_str::contig_str(const string &str) : string{str}, number{is_number(str)} {}
contig_str::contig_str() : string{} {}
bool
contig_str::operator<(const contig_str &other) const {
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

bool
contig_str::operator==(const string &other) const {
    return *this == static_cast<contig_str>(other);
}

bool
contig_str::is_number(const string &str) {
    string _str = str;
    if (_str.empty()) {
        return false;
    }
    if (_str.find("chr") == 0) {
        _str = str.substr(3);
    }
    return std::find_if(_str.begin(), _str.end(), [](unsigned char c) { return !std::isdigit(c); }) == _str.end();
}
ginterval::ginterval() : interval(0, 0), chr(""), plus_strand(true) {}

ginterval::ginterval(string chr, int start, int end, const string &strand)
    : interval(start, end), chr(chr), plus_strand(strand == "+") {}
ginterval::ginterval(string chr, int start, int end, bool plus_strand)
    : interval(start, end), chr(chr), plus_strand(plus_strand) {}
ginterval::ginterval(const ginterval &g1, const ginterval &g2)
    : interval(g1, g2), chr(g1.chr), plus_strand(g1.plus_strand) {}  // Merge constructor
ginterval::~ginterval() {}
int
ginterval::overlap(const ginterval &other) const {
    if (chr != other.chr) {
        return 0;
    }
    return interval::overlap(other);
}
double
ginterval::reciprocal(const ginterval &other) const {
    return static_cast<double>(overlap(other)) / larger_interval(*this, other);
}
bool
ginterval::operator<(const ginterval &other) const {
    if (other.chr != chr) {
        return chr < other.chr;
    }
    if (start == other.start) {
        return end < other.end;
    }
    return start < other.start;
    //}
}
bool
ginterval::operator==(const ginterval &other) const {
    return chr == other.chr && start == other.start && end == other.end && plus_strand == other.plus_strand;
}

gtf::entry_type
gtf::type_from_string(const string &type_str) {
    if (type_str == "gene") {
        return gtf::entry_type::gene;
    }
    else if (type_str == "transcript") {
        return gtf::entry_type::transcript;
    }
    else if (type_str == "exon") {
        return gtf::entry_type::exon;
    }
    else if (type_str == "five_prime_utr") {
        return gtf::entry_type::five_prime_utr;
    }

    else if (type_str == "three_prime_utr") {
        return gtf::entry_type::three_prime_utr;
    }
    else if (type_str == "start_codon") {
        return gtf::entry_type::start_codon;
    }
    else if (type_str == "stop_codon") {
        return gtf::entry_type::stop_codon;
    }
    else if (type_str == "CDS") {
        return gtf::entry_type::CDS;
    }
    else if (type_str == "Selenocysteine") {
        return gtf::entry_type::Selenocysteine;
    }
    else {
        return gtf::entry_type::other;
    }
}
string
gtf::type_to_string(gtf::entry_type type) {
    switch (type) {
        case gtf::entry_type::gene:
            return "gene";
        case gtf::entry_type::transcript:
            return "transcript";
        case gtf::entry_type::exon:
            return "exon";
        case gtf::entry_type::five_prime_utr:
            return "five_prime_utr";
        case gtf::entry_type::three_prime_utr:
            return "three_prime_utr";
        case gtf::entry_type::start_codon:
            return "start_codon";
        case gtf::entry_type::stop_codon:
            return "stop_codon";
        case gtf::entry_type::CDS:
            return "CDS";
        case gtf::entry_type::Selenocysteine:
            return "Selenocysteine";
        case gtf::entry_type::other:
            return "other";
        default:
            return "other";
    }
}

gtf::gtf(const string &gtf_line) {
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

gtf::gtf(const gtf &other) = default;

gtf::gtf(const ginterval &pos, gtf::entry_type type, const std::map<string, string> &info)
    : ginterval(pos), type(type), info(info), source{"TKSM"} {}
gtf::gtf(const ginterval &pos, gtf::entry_type type) : ginterval(pos), type(type), source{"TKSM"} {}

gtf::gtf() = default;
ostream &
operator<<(ostream &os, const gtf &g) {
    os << g.chr << "\t" << g.source << "\t" << gtf::type_to_string(g.type) << "\t" << g.start << "\t" << g.end << "\t"
       << (g.plus_strand ? '+' : '-') << "\t"
       << "."
       << "\t";
    for (auto iter = g.info.begin(); iter != g.info.end(); ++iter) {
        os << iter->first << " \"" << iter->second << "\";";
    }
    os << "\n";
    return os;
}
string
gtf::to_string() const {
    std::ostringstream os;
    os << *this;
    return os.str();
}

transcript::transcript(const transcript &other)
    : gtf{other}, abundance{other.abundance}, exons{other.exons}, comment{other.comment} {}
transcript::transcript(const gtf &entry) : gtf(entry), abundance(0.0) {}
transcript::transcript(const gtf &entry, double abundance) : gtf(entry), abundance(abundance) {}
transcript::transcript(const gtf &entry, double abundance, const string &comment)
    : gtf(entry), abundance(abundance), comment(comment) {}

ostream &
operator<<(ostream &os, const transcript &t) {
    os << static_cast<const gtf &>(t);
    for (auto &exon : t.exons) {
        os << exon;
    }
    return os;
}
string
transcript::to_abundance_str() const {
    std::ostringstream os;
    os << info.at("transcript_id") << "\t" << abundance << "\t" << comment;
    return os.str();
}
string
transcript::to_string() const {
    std::ostringstream os;
    os << *this;
    return os.str();
}
void
transcript::add_exon(const gtf &exon) {
    exons.push_back(exon);
}

void
transcript::prepend_exon(const gtf &exon) {
    exons.insert(exons.begin(), exon);
}

void
transcript::add_exons(const vector<gtf> &exons) {
    for (auto &exon : exons) {
        add_exon(exon);
    }
}
void
transcript::sort_exons() {
    std::sort(exons.begin(), exons.end(), [](const gtf &a, const gtf &b) { return a.start < b.start; });
}
void
transcript::set_abundance(double abundance) {
    this->abundance = abundance;
}
double
transcript::get_abundance() const {
    return abundance;
}
vector<gtf> &
transcript::get_exons() {
    return exons;
}
const vector<gtf> &
transcript::cget_exons() const {
    return exons;
}
auto
transcript::get_exons(int start, int end) {
    return exons | std::ranges::views::filter([start, end](const gtf &exon) -> bool {
               return exon.start >= start && exon.end <= end;
           });  // clangd ignore
}

string
transcript::get_comment() const {
    return comment;
}



auto
transcript::operator==(const transcript &other) const -> bool {
    return exons == other.exons && info.at("transcript_id") == other.info.at("transcript_id") &&
           info.at("gene_id") == other.info.at("gene_id") && info.at("gene_name") == other.info.at("gene_name") &&
           start == other.start && end == other.end && chr == other.chr && plus_strand == other.plus_strand;
}

mapping::mapping() : rid("-1") {}
mapping::mapping(const string &paf, int max_skip, int min_segment) {
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

einterval::einterval() : ginterval() {}
einterval::einterval(const ginterval &gi) : ginterval(gi) {}

einterval::einterval(string chr, int start, int end, const bool plus_strand)
    : ginterval(chr, start, end, plus_strand) {}
einterval::einterval(string chr, int start, int end, const string &strand) : ginterval(chr, start, end, strand) {}

einterval::einterval(const einterval &ei, int start, int end)
    : ginterval(ei.chr, ei.start + start, ei.start + end, ei.plus_strand) {
    for (const base_mod &e : ei.errors) {
        if (e.position < start || e.position > end) {
            continue;
        }
        errors.emplace_back(e.position - start, e.base);
    }
}

void
einterval::add_error(int position, char base) {
    errors.emplace_back(position, base);
}

void
einterval::add_error(const std::pair<int, char> &error) {
    errors.emplace_back(error.first, error.second);
}
void
einterval::add_error(const base_mod &error) {
    errors.emplace_back(error);
}

void
einterval::add_errors(const vector<base_mod> &err) {
    errors.insert(errors.end(), err.begin(), err.end());
}

void
einterval::truncate(int start, int end) {
    sort_errors();
    assert(start >= 0);
    assert(end <= this->end - this->start);
    assert(start < end);

    this->start += start;
    this->end = this->start + (end - start);

    if (errors.size() == 0) {
        return;
    }

    if (start > 0) {
        for (auto &e : errors) {
            e.position -= start;
        }
    }

    errors.erase(std::remove_if(errors.begin(), errors.end(),
                                [start, end](const base_mod &error) {
                                    return error.position < 0 || error.position >= end - start;
                                }),
                 errors.end());
}
void
einterval::parse_and_add_errors(const string &error_string) {
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

void
einterval::sort_errors() {
    std::sort(errors.begin(), errors.end(),
              [](const base_mod &a, const base_mod &b) { return a.position < b.position; });
}

auto
einterval::error_str() const -> string {
    return join_str(errors.begin(), errors.end(), ",");
}

molecule_descriptor::molecule_descriptor() : _id(""), _reversed{false}, _depth{0} {}

molecule_descriptor::molecule_descriptor(const string &id, bool reversed) : _id(id), _reversed(reversed), _depth(1) {}

molecule_descriptor::molecule_descriptor(const molecule_descriptor &other)
    : _id(other._id),
      _reversed(other._reversed),
      _depth(other._depth),
      _segments(other._segments.begin(), other._segments.end()),
      meta{other.meta} {}

molecule_descriptor::molecule_descriptor(const transcript &trans)
    : _id(trans.info.at("transcript_id")), _reversed(!trans.plus_strand), _depth(trans.get_abundance()) {
    for (const gtf &gi : trans.cget_exons()) {
        append_segment(gi);
    }
}

string
molecule_descriptor::get_id() const {
    return _id;
}

int
molecule_descriptor::get_depth() const {
    return _depth;
}

auto
molecule_descriptor::cget_segments() const -> const vector<einterval> & {
    return _segments;
}
auto
molecule_descriptor::get_segments() -> vector<einterval> & {
    return _segments;
}
molecule_descriptor *
molecule_descriptor::assign_segments(const vector<einterval> &segments) {
    _segments = segments;
    return this;
}
molecule_descriptor *
molecule_descriptor::id(const string &id) {
    _id = id;
    return this;
}

// Assumes that special characters are not used in the meta data (',', ';', '=')
molecule_descriptor *
molecule_descriptor::comment(const string &comment) {
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

molecule_descriptor *
molecule_descriptor::add_comment(const string &key, const string &value) {
    meta[key].push_back(value);
    return this;
}

molecule_descriptor *
molecule_descriptor::drop_comment(const string &key) {
    meta.erase(key);
    return this;
}

// TODO add key verification
const vector<string> &
molecule_descriptor::get_comment(const string &key) const {
    return meta.at(key);
}

bool
molecule_descriptor::has_comment(const string &key) const {
    return meta.contains(key);
}
molecule_descriptor *
molecule_descriptor::depth(int depth) {
    _depth = depth;
    return this;
}

molecule_descriptor *
molecule_descriptor::prepend_segment(const einterval &i) {
    _segments.insert(_segments.begin(), i);
    // for (std::pair<int, char> &errs : errors_so_far) {
    //    errs = std::make_pair(errs.first + i.end - i.start, errs.second);
    // }
    return this;
}

molecule_descriptor *
molecule_descriptor::append_segment(const einterval &i) {
    _segments.push_back(i);
    return this;
}

molecule_descriptor *
molecule_descriptor::add_error(base_mod error) {
    auto iter = _segments.begin();
    while (iter->size() <= error.position) {
        error.position -= iter->size();
        iter++;
    }
    iter->add_error(error);
    return this;
}

size_t
molecule_descriptor::size() const {
    return std::accumulate(_segments.begin(), _segments.end(), 0L,
                           [](size_t sum_so_far, const ginterval &g) -> size_t { return sum_so_far + g.size(); });
}

string
molecule_descriptor::dump_comment() const {
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

molecule_descriptor *
molecule_descriptor::concat(const molecule_descriptor &other) {
    _segments.insert(_segments.end(), other._segments.begin(), other._segments.end());
    return this;
}

std::ostream &
operator<<(std::ostream &ost, const molecule_descriptor &md) {
    print_tsv(ost, "+" + md._id, md._depth, md.dump_comment());

    for (const einterval &ival : md._segments) {
        print_tsv(ost, ival.chr, ival.start, ival.end, (ival.plus_strand ? "+" : "-"), ival.error_str());
    }
    return ost;
}

molecule_descriptor
flip_molecule(const molecule_descriptor &md) {
    molecule_descriptor flipped_md{md.get_id(), md._reversed};
    for (const auto &segment : md.cget_segments() | std::views::reverse) {
        auto flipped_segment        = segment;
        flipped_segment.plus_strand = !flipped_segment.plus_strand;
        flipped_md.append_segment(flipped_segment);
    }
    flipped_md.meta = md.meta;
    flipped_md.depth(md.get_depth());

    return flipped_md;
}

