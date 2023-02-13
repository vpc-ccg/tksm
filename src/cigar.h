#ifndef CIGAR_H
#define CIGAR_H

#include <string>
#include <utility>
#include <vector>

#define CIGAR_CHARACTERS "DHIMNPSX="

class cigar {
    std::vector<std::pair<int, char>> cigarray;

public:
    cigar(const std::string& cgr) {
        size_t prev  = -1;
        size_t index = cgr.find_first_of(CIGAR_CHARACTERS);
        while (index != std::string::npos) {
            int len = std::stoi(cgr.substr(1 + prev, index));
            cigarray.push_back(std::make_pair(len, cgr[index]));
            prev  = index;
            index = cgr.find_first_of(CIGAR_CHARACTERS, index + 1);
        }
    }
    decltype(cigar::cigarray.begin()) begin() { return cigarray.begin(); }

    decltype(cigar::cigarray.end()) end() { return cigarray.end(); }

    auto operator[](size_t index) const { return cigarray[index]; }

    enum class type {
        matched,
        onquery,
        ontemplate,
        hardclip,
        softclip,
        notcigar,
    };

    static type what_is_this(char c) {
        switch (c) {
            case 'M':
            case '=':
            case 'X':
                return type::matched;
            case 'D':
            case 'N':
                return type::ontemplate;
            case 'I':
            case 'P':
                return type::onquery;
            case 'H':
                return type::hardclip;
            case 'S':
                return type::softclip;
            default:
                return type::notcigar;
        }
    }
};

#endif
