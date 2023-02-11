#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <ostream>
#include <map>

#define FMT_HEADER_ONLY
#include <fmt/core.h>
#define FMTLOG_HEADER_ONLY
#include <fmtlog/fmtlog.h>
using std::map;
using std::string;
using std::vector;
using std::ostream;

class LogLevels{
    public:

        static std::string log_choices(){
            static map<string, fmtlog::fmtlogT::LogLevel> log_level_map {
                {"DEBUG", fmtlog::DBG},
                {"INFO", fmtlog::INF},
                {"WARN", fmtlog::WRN},
                {"ERROR", fmtlog::ERR},
                {"OFF", fmtlog::OFF}
            };
        
            std::string str;
            for(auto &p: log_level_map){
                if(str.size()){
                    str+=", ";
                }
                str+=p.first;
            }
            return str;
        }
        static fmtlog::fmtlogT::LogLevel parse_loglevel(const std::string &str){   
            static map<string, fmtlog::fmtlogT::LogLevel> log_level_map {
                {"DEBUG", fmtlog::DBG},
                {"INFO", fmtlog::INF},
                {"WARN", fmtlog::WRN},
                {"ERROR", fmtlog::ERR},
                {"OFF", fmtlog::OFF}
            };
            return log_level_map.at(str);
        }
};

template< class B>
inline void print_tsv(std::ostream &ost, B b){
    ost << b << "\n";
}

template <class B, class... A>
inline void print_tsv(std::ostream &ost, B b, A... a){
    ost << b << "\t";
    print_tsv(ost, a...);
}


template<class ITER>
inline std::string join_str(const ITER &begin, const ITER &end, const std::string &del){
    std::string str;
    for(ITER i = begin; i != end; ++i){
        if(i!=begin){
            str+=del;
        }
        str+=*i;
    }
    return str;
}


//https://github.com/fenbf/StringViewTests
inline std::vector<std::string_view>
splitSV(std::string_view strv, std::string_view delims = " ")
{
	std::vector<std::string_view> output;
	//output.reserve(strv.length() / 4);
	size_t first = 0;

	while (first < strv.size())
	{
		const auto second = strv.find_first_of(delims, first);
		//std::cout << first << ", " << second << '\n';
		if (first != second)
		{
			output.emplace_back(strv.substr(first, second-first));
		}

		if (second == std::string_view::npos)
			break;

		first = second + 1;
	}

	return output;
}

inline std::vector<std::string> rsplit(std::string str, std::string delim){
    std::vector<std::string> splits;
    size_t p1 = 0;
    size_t p2 = 0;
    while((p2= str.find(delim,p1)) != std::string::npos){
        splits.push_back(str.substr(p1,p2-p1));
        p1 = p2+delim.size();
    }
    splits.push_back(str.substr(p1));
    return splits;
}

inline std::string strip_str(const std::string &inpt, const std::string &chrs)
{

    auto frst = inpt.find_first_not_of(chrs);
    auto last = inpt.find_last_not_of(chrs);
    if(frst == std::string::npos || last == std::string::npos){
        return "";
    }
    return inpt.substr(frst, last-frst+1);
}
inline void strip_for_each(std::vector<std::string> &vec, const std::string &chrs = " "){

    for( std::string &st : vec){
        st = strip_str(st, chrs);
    }
}


inline void format_annot_id(std::string& id, bool remove_version = true){
    if(remove_version){
        if(id.find_last_of(".") != std::string::npos){
            id = rsplit(id, ".")[0];
        }
    }
}

#endif
