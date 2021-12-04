/* The MIT License of the original algorithm

   Copyright (c) 2008, 2009, 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, andgcor sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Last Modified: 2021-11-10 by f0t1h */

#ifndef KSEQ_H
#define KSEQ_H

#include <ctype.h>
#include <string.h>

#include <algorithm>

#include <type_traits>


namespace kekseq{
    namespace{ // internal
        enum class delimiters{
            space = 0, // isspace(): \t, \n, \v, \f, \r
            tab   = 1, // isspace() && !' '
            line  = 2, // line separator: "\n" (Unix) or "\r\n" (Windows)
        };
        template <typename Enumeration>
        constexpr auto as_integer(Enumeration const value)
            -> typename std::underlying_type<Enumeration>::type
        {
            return static_cast<typename std::underlying_type<Enumeration>::type>(value);
        }

        template<typename K>
        constexpr int roundup32(K x){
            --x;
            x|=x>>1;
            x|=x>>2;
            x|=x>>4;
            x|=x>>8;
            x|=x>>16;
            ++x;
            return x;
        }
    }
    struct kstring {
        size_t l, m;
        char *s;
        kstring (): l(0), m(0), s(0) {}
        kstring (int m): l(0), m(m), s(new char[m]) {}
        ~kstring (){
            delete[] s;
        }
        void expand(int increase){
            m = l + increase; 
            roundup32(m);
            char *tmp = new char[m];
            std::copy_n(s, l, tmp);
            delete[] s;
            s = tmp;
        }
        size_t size() const {
            return l;
        }
        void size(size_t i){
            l = i;
        }
        size_t limit() const {
            return m;
        }
        void limit(size_t i){
            m = i;
        }
        char *seq() const {
            return s;
        }
    };

    template<class FILETYPE, int (READ)(FILETYPE , void *, unsigned), size_t BUFSIZE = 16384 >
    class kstream { 
        using type_t = FILETYPE;

        static const size_t bufsize = BUFSIZE;

        public:
        unsigned char buf[BUFSIZE]; 

        int begin, end, is_eof; 
        type_t f;
        kstream(type_t f): begin(0), end(0), is_eof(0), f(f){
        }
        ~kstream(){

        }
        int err(){
            return end < 0;
        }
        int eof(){
            return is_eof && begin >= end;
        }
        int getc() { 
            if (err()) return -3; 
            if (eof()) return -1; 
            if (begin >= end) { 
                begin = 0; 
                end = READ(f, buf, bufsize); 
                if (end == 0) { is_eof = 1; return -1; } 
                else if (end < 0) {is_eof = 1; return -3; } 
            } 
            return (int)buf[begin++]; 
        }
        int getuntil2(int delimiter, kstring &str, int *dret, int append) 
        { 
            int gotany = 0; 
            if (dret) *dret = 0; 
            str.l = append? str.l : 0; 
            for (;;) { 
                int i; 
                if (err()) return -3; 
                if (begin >= end) { 
                    if (!is_eof) { 
                        begin = 0; 
                        end = READ(f, buf, bufsize); 
                        if (end == 0) { is_eof = 1; break; } 
                        if (end == -1) {is_eof = 1; return -3; } 
                    } else break; 
                } 
                if (delimiter == as_integer(delimiters::line)) { 

//                    unsigned char *sep = (unsigned char*)memchr(buf + begin, '\n', end - begin); 
                    unsigned char *sep = std::find(buf + begin, buf + end, '\n'); 
                    i = sep != NULL ? sep - buf : end; 
                } else if (delimiter > as_integer(delimiters::line)) { 
                    for (i = begin; i < end; ++i) 
                        if (buf[i] == delimiter) break; 
                } else if (delimiter == as_integer(delimiters::space)) { 
                    for (i = begin; i < end; ++i) 
                        if (isspace(buf[i])) break; 
                } else if (delimiter == as_integer(delimiters::tab)) { 
                    for (i = begin; i < end; ++i) 
                        if (isspace(buf[i]) && buf[i] != ' ') break; 
                } else i = 0; /* never come to here! */ 
                if (str.m - str.l < (size_t)(i - begin + 1)) {
                    str.expand( i - begin + 1);
                } 
                gotany = 1; 
//                memcpy(str.s + str.l, buf + begin, i - begin);
                std::copy(buf+begin, buf+i, str.s+str.l);
                str.l = str.l + (i - begin); 
                begin = i + 1; 
                if (i < end) { 
                    if (dret) *dret = buf[i]; 
                    break; 
                } 
            } 
            if (!gotany && eof()) return -1; 
            if (str.s == 0) { 
                str.m = 1;
                str.s = new char[1]();

            } else if (delimiter == as_integer(delimiters::line) && str.l > 1 && str.s[str.l-1] == '\r') --str.l; 
            str.s[str.l] = '\0'; 
            return str.l; 
        } 
        int getuntil(int delimiter, kstring &str, int *dret) 
        {
            return getuntil2(delimiter, str, dret, 0); 
        }

    };


    template<class FILETYPE, int (READ)(FILETYPE , void *, unsigned),size_t BUFSIZE = 16384>
    struct kseq{ 
        kstring name, comment, seq, qual; 
        int last_char, is_fastq; 
        kstream<FILETYPE, READ, BUFSIZE> ks;
       
        kseq(FILETYPE fd): name(128), comment(128), seq(256), qual(256), last_char(0), is_fastq(0), ks(fd){
        }

        ~kseq(){
        }

        void rewind(){
            last_char = ks.is_eof = ks.begin = ks.end;
        }
        
    /* Return value:
       >=0  length of the sequence (normal)
       -1   end-of-file
       -2   truncated quality string
       -3   error reading stream
     */
        int read(){

            int c,r; 

            if (last_char == 0) { /* then jump to the next header line */ 
                while ((c = ks.getc()) >= 0 && c != '>' && c != '@'); 
                if (c < 0) return c; /* end of file or error*/ 
                last_char = c; 
            } /* else: the first header char has been read in the previous call */ 
            comment.l = seq.l = qual.l = 0; /* reset all members */ 
            if ((r=ks.getuntil( 0, name, &c)) < 0) return r;  /* normal exit: EOF or error */ 
            if (c != '\n') ks.getuntil(as_integer(delimiters::line), comment, 0); /* read FASTA/Q comment */ 

            while ((c = ks.getc()) >= 0 && c != '>' && c != '+' && c != '@') { 
                if (c == '\n') continue; /* skip empty lines */ 
                seq.s[seq.l++] = c; /* this is safe: we always have enough space for 1 char */ 
                ks.getuntil2( as_integer(delimiters::line), seq, 0, 1); /* read the rest of the line */ 
            } 
            if (c == '>' || c == '@') last_char = c; /* the first header char has been read */ 
            if (seq.l + 1 >= seq.m) { /* seq.s[seq.l] below may be out of boundary */ 
                seq.expand(2);

            } 
            seq.s[seq.l] = 0;	/* null terminated string */ 
            is_fastq = (c == '+'); 
            if (!is_fastq) return seq.l; /* FASTA */ 
            if (qual.m < seq.m) {	/* allocate memory for qual in case insufficient */ 
                qual.expand(- qual.m + seq.m);
            } 
            while ((c = ks.getc()) >= 0 && c != '\n'); /* skip the rest of '+' line */ 
            if (c == -1) return -2; /* error: no quality string */ 
            while ((c = ks.getuntil2(as_integer(delimiters::line), qual, 0, 1) >= 0 && qual.l < seq.l)); 
            if (c == -3) return -3; /* stream error */ 
            last_char = 0;	/* we have not come to the next header line */ 
            if (seq.l != qual.l) return -2; /* error: qual string is of a different length */ 
            return seq.l; 
        }
    };
}
#endif

