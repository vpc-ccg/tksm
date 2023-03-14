
#include "ut.hpp"
#include "../src/truncate.cpp"
//#include "../src/interval.h"

#include <iostream>
#include <string>

using namespace std;
using namespace boost::ut;

int main(int argc, char **argv){
    {
        molecule_descriptor md;

        md   .append_segment({"1",   0, 100, true})
            ->append_segment({"1", 100, 200, true})
            ->append_segment({"1", 200, 300, true})
            ->append_segment({"1", 300, 400, true});

        truncate(md, 200, 100);
        "truncated_size"_test = [&]{
            expect(eq(md.size(), 200UL));
        };
        "truncate_segments"_test = [&]{
            expect(eq(md.cget_segments().size(), 2));
            expect(eq(md.cget_segments()[0].start, 0));
            expect(eq(md.cget_segments()[0].end, 100));
            expect(eq(md.cget_segments()[1].start, 100));
            expect(eq(md.cget_segments()[1].end, 200));
        };
        
        truncate (md, 150, 100);
        "truncated_between_size"_test = [&]{
            expect(eq(md.size(), 150UL));
        };
        "truncated_between_segments"_test = [&]{
            expect(eq(md.cget_segments().size(), 2));
            expect(eq(md.cget_segments()[0].start, 0));
            expect(eq(md.cget_segments()[0].end, 100));
            expect(eq(md.cget_segments()[1].start, 100));
            expect(eq(md.cget_segments()[1].end, 150));
        };
        
        truncate(md, 50, 100);
        "truncated_min_size"_test = [&]{
            expect(eq(md.size(), 100UL));
        };
        "truncated_min_segments"_test = [&]{
            expect(eq(md.cget_segments().size(), 1));
            expect(eq(md.cget_segments()[0].start, 0));
            expect(eq(md.cget_segments()[0].end, 100));
        };
    }

    return 0;
}
