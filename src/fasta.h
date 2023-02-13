#ifndef FASTA_H
#define FASTA_H
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using std::cerr;
using std::ifstream;
using std::map;
using std::string;

template <typename K>
inline constexpr int
roundup32(K x) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    ++x;
    return x;
}

inline map<string, string>
read_fasta_fast(const string& fasta_path) {
    string fai_path = fasta_path + ".fai";
    map<string, string> contig2seq;

    string buf;
    if (std::filesystem::exists(fai_path)) {  // This helps saves little time, but keeping it for reference
        std::cerr << "Index exists... Allocating memory ahead of time!\n";
        std::ifstream f(fai_path);

        while (std::getline(f, buf)) {
            std::istringstream str(buf);
            string contig;
            int base_count;
            str >> contig >> base_count;
            contig2seq[contig].reserve(roundup32(base_count));
        }
    }

    std::cerr << "Reading fasta file!\n";
    string contig;
    std::ifstream f(fasta_path);
    while (std::getline(f, buf)) {
        if (buf[0] == '>') {  // new contig
            std::istringstream str(buf.c_str() + 1);
            str >> contig;
        }
        else {
            contig2seq[contig] += buf;
        }
    }
    return contig2seq;
}
#endif

