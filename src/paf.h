#pragma once
#ifndef _PAF_H_
#define _PAF_H_


#include <exception>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>

#include "util.h"
#include "interval.h"

//0046a427-30d0-48b7-a0db-0643ed8154b1    748     78      691     -       16      90338345        85801201        85807002        540     656     60      tp:A:P  mm:i:46 gn:i:70 go:i:41 cg:Z:16M1D38M1D8M3671N9M3I4M2I9M2I8M1D31M1D6M4I50M1D15M1I8M2D23M628N42M1I7M2I19M3D37M1I3M2I14M1D4M2D873N16M1I7M1D4M1I6M2D11M1D5M3D11M1D13M1I7M4D2M2D6M5D29M1D6M2D4M1D6M2I17M1D6M1D9M1D3M1I32M3I4M1D2M2D9M1D20M
struct paf{ 
    std::string qname;
    int qlen;
    int qstart;
    int qend;
    std::string strand;
    std::string tname;
    int tlen;
    int tstart;
    int tend;
    int matchc;
    int alig_len;
    int map_qual;
    std::map<std::string, std::string> info;

    paf(const std::string &line){
//        std::vector<std::string> fields = rsplit(line, "\t");
        std::stringstream ss {line};
        ss >> qname >> qlen >> qstart >> qend >> strand
            >> tname >> tlen >> tstart >> tend >> matchc
            >> alig_len >> map_qual;
       
        std::string buffer;
        while(ss >> buffer){
            std::vector<std::string> fields = rsplit(buffer, ":");
            info[fields[0]] = fields[2];
        }
    }
    bool primary() const{
        auto iter_info = info.find("tp");
        if(iter_info == info.end()){
            std::cerr << "Paf field doesn't have primary/secondary field!\n";
            throw std::exception();
        }
        return iter_info->second == "P";
    }
};

#endif
