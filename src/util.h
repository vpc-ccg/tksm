#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>


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
#endif
