#include "mdf.h"



generator<molecule_descriptor>
stream_mdf(std::istream &ist, bool unroll) {
    string buffer;
    buffer.reserve(1000);
    std::getline(ist, buffer);
    while (ist) {
        auto fields = rsplit(buffer, "\t");
        string id{fields[0].substr(1)};
        int depth{stoi(fields[1])};
        //        int exon_count{stoi(fields[2])};
        string comment{fields[2]};
        vector<einterval> segments;
        vector<std::pair<int, char>> errors_so_far;

        std::getline(ist, buffer);
        while (ist && buffer[0] != '+') {
            //        for(int i = 0; i< exon_count; ++i){

            auto fields = rsplit(buffer, "\t");
            string chr{fields[0]};
            int start{stoi(fields[1])};
            int end{stoi(fields[2])};
            string strand{fields[3]};
            string error_str;
            if (fields.size() > 4) {
                error_str = fields[4];
            }

            segments.emplace_back(chr, start, end, strand).parse_and_add_errors(error_str);
            std::getline(ist, buffer);
        }
        molecule_descriptor md{id, !segments[0].plus_strand};
        md.depth(depth)->assign_segments(segments)->comment(comment);
        if (unroll && md.get_depth() > 1) {
            molecule_descriptor mdc = md;
            mdc.depth(1);
            for (int i = 0; i < md.get_depth(); ++i) {
                molecule_descriptor mdcx = mdc;
                mdcx.id(md.get_id() + "_" + std::to_string(i));
                co_yield mdcx;
            }
        }
        else {
            co_yield md;
        }
    }
}

generator<molecule_descriptor>
stream_mdf(const string &filename, bool unroll) {
    std::basic_ifstream<char> ist{filename};
    if (!ist) {
        throw std::runtime_error("Could not open file " + filename);
    }

    string buffer;
    buffer.reserve(1000);
    std::getline(ist, buffer);
    while (ist) {
        auto fields = rsplit(buffer, "\t");
        string id{fields[0].substr(1)};
        int depth{stoi(fields[1])};
        //        int exon_count{stoi(fields[2])};
        string comment{fields[2]};
        vector<einterval> segments;
        vector<std::pair<int, char>> errors_so_far;

        std::getline(ist, buffer);
        while (ist && buffer[0] != '+') {
            //        for(int i = 0; i< exon_count; ++i){

            auto fields = rsplit(buffer, "\t");
            string chr{fields[0]};
            int start{stoi(fields[1])};
            int end{stoi(fields[2])};
            string strand{fields[3]};
            string error_str;
            if (fields.size() > 4) {
                error_str = fields[4];
            }

            segments.emplace_back(chr, start, end, strand).parse_and_add_errors(error_str);
            std::getline(ist, buffer);
        }
        molecule_descriptor md{id, !segments[0].plus_strand};
        md.depth(depth)->assign_segments(segments)->comment(comment);
        if (unroll && md.get_depth() > 1) {
            molecule_descriptor mdc = md;
            mdc.depth(1);
            for (int i = 0; i < md.get_depth(); ++i) {
                molecule_descriptor mdcx = mdc;
                mdcx.id(md.get_id() + "_" + std::to_string(i));
                co_yield mdcx;
            }
        }
        else {
            co_yield md;
        }
    }
}
