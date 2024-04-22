#ifndef UTIL_H
#define UTIL_H

#include <fmt/core.h>

#include <algorithm>
#include <array>
#include <map>
#include <ostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

#define FMTLOG_HEADER_ONLY
#include <fmtlog/fmtlog.h>

class runtime_error : public std::runtime_error {
public:
    runtime_error(const std::string &msg) : std::runtime_error(msg) { fmtlog::poll(true); }
};

using std::map;
using std::ostream;
using std::string;
using std::vector;

static const string ASCII_ART = R"(
 _______  ___   _  _______  __   __    
|       ||   | | ||       ||  |_|  |   
|_     _||   |_| ||  _____||       |   
  |   |  |      _|| |_____ |       |   
  |   |  |     |_ |_____  ||       |   
  |   |  |    _  | _____| || ||_|| |   
  |___|  |___| |_||_______||_|   |_|
)";

class num2seq {
public:
    std::array<char, 4> table{{'A', 'T', 'C', 'G'}};
    uint64_t max;
    num2seq(uint64_t max) : max(max) {}
    string operator[](uint64_t num) const {
        std::stringstream st;
        int back = 0;
        for (uint64_t mask = 3; mask < max; mask <<= 2) {
            st << table[(num & mask) >> back];
            back += 2;
        }
        return st.str();
    }
};

class fmt2seq {
public:
    std::array<char, 4> table{{'A', 'T', 'C', 'G'}};
    std::array<string, 128> lookup;
    string fmt;
    fmt2seq(const string &fmt) : fmt(fmt) {
#define set_l(A, B) lookup[A[0]] = B
#define set_sl(A, B) set_l(#A, #B)
        set_sl(A, A);
        set_sl(G, G);
        set_sl(T, T);
        set_sl(C, C);
        set_sl(U, U);
        set_sl(R, GA);
        set_sl(Y, TC);
        set_sl(K, GT);
        set_sl(M, AC);
        set_sl(S, GC);
        set_sl(W, AT);
        set_sl(B, GTC);
        set_sl(D, GAT);
        set_sl(H, ACT);
        set_sl(V, GCA);
        set_sl(N, AGCT);
#undef set_l
#undef set_sl
    }
    template <class RANDGEN>
    string operator[](RANDGEN &gen) const {
        std::stringstream st;

        for (char c : fmt) {
            string buffer;

            std::sample(lookup[c].begin(), lookup[c].end(), std::back_inserter(buffer), 1, gen);
            st << buffer;
        }
        return st.str();
    }
};

class LogLevels {
public:
    static std::string log_choices() {
        static map<string, fmtlog::fmtlogT::LogLevel> log_level_map{{"DEBUG", fmtlog::DBG},
                                                                    {"INFO", fmtlog::INF},
                                                                    {"WARN", fmtlog::WRN},
                                                                    {"ERROR", fmtlog::ERR},
                                                                    {"OFF", fmtlog::OFF}};

        std::string str;
        for (auto &p : log_level_map) {
            if (str.size()) {
                str += ", ";
            }
            str += p.first;
        }
        return str;
    }
    static fmtlog::fmtlogT::LogLevel parse_loglevel(const std::string &str) {
        static map<string, fmtlog::fmtlogT::LogLevel> log_level_map{{"DEBUG", fmtlog::DBG},
                                                                    {"INFO", fmtlog::INF},
                                                                    {"WARN", fmtlog::WRN},
                                                                    {"ERROR", fmtlog::ERR},
                                                                    {"OFF", fmtlog::OFF}};
        return log_level_map.at(str);
    }
};

inline void
report_missing_parameter(const std::string &param_name) {
    loge("Missing parameter: {}", param_name);
}

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

template <class ITER>
inline std::string
join_str(const ITER &begin, const ITER &end, const std::string &del) {
    std::string str;
    for (ITER i = begin; i != end; ++i) {
        if (i != begin) {
            str += del;
        }
        str += *i;
    }
    return str;
}

// https://github.com/fenbf/StringViewTests
inline std::vector<std::string_view>
splitSV(std::string_view strv, std::string_view delims = " ") {
    std::vector<std::string_view> output;
    // output.reserve(strv.length() / 4);
    size_t first = 0;

    while (first < strv.size()) {
        const auto second = strv.find_first_of(delims, first);
        if (first != second) {
            output.emplace_back(strv.substr(first, second - first));
        }

        if (second == std::string_view::npos) break;

        first = second + 1;
    }

    return output;
}

inline std::vector<std::string>
rsplit(std::string str, std::string delim) {
    std::vector<std::string> splits;
    size_t p1 = 0;
    size_t p2 = 0;
    while ((p2 = str.find(delim, p1)) != std::string::npos) {
        splits.push_back(str.substr(p1, p2 - p1));
        p1 = p2 + delim.size();
    }
    splits.push_back(str.substr(p1));
    return splits;
}

inline std::string
strip_str(const std::string &inpt, const std::string &chrs) {
    auto frst = inpt.find_first_not_of(chrs);
    auto last = inpt.find_last_not_of(chrs);
    if (frst == std::string::npos || last == std::string::npos) {
        return "";
    }
    return inpt.substr(frst, last - frst + 1);
}
inline void
strip_for_each(std::vector<std::string> &vec, const std::string &chrs = " ") {
    for (std::string &st : vec) {
        st = strip_str(st, chrs);
    }
}

inline void
format_annot_id(std::string &id, bool remove_version = true) {
    if (remove_version) {
        if (id.find_last_of(".") != std::string::npos) {
            id = rsplit(id, ".")[0];
        }
    }
}

#endif
