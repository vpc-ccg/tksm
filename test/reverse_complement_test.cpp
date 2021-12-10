
#include "ut.hpp"
#include "../src/reverse_complement.h"

#include <iostream>
#include <string>

using namespace std;
using namespace boost::ut;

using namespace reverse_complement;

int main(int argc, char **argv){

    "empty"_test = []{
        string empty = "";
        expect(complement_string(empty) == empty);
    };
    "twice_complement"_test = []{
        string s1 = "AGTCATCGATCGACGACTACG";
        expect(complement_seq(complement_seq(s1)) == s1);
    };
    "singles"_test = [] {
        expect(complement_seq("A")=="T");
        expect(complement_seq("g")=="c");
        expect(complement_seq("U")=="A");
        expect(complement_seq("t")=="a");
        expect(complement_seq("C")=="G");
    };

    "inplace"_test = [] {
        string seq = "AGTCATGC";
        string cop = seq;
        complement_inplace(cop);
        expect(cop == "GCATGACT");
        complement_inplace(cop);
        expect(cop==seq);
    };

    return 0;
}
