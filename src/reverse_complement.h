#ifndef REVERSE_COMPLEMENT_H
#define REVERSE_COMPLEMENT_H

#include <array>
#include <string>

namespace reverse_complement{
    template<char C>
    struct complement{
        static constexpr char value = C;
    };
    //Specializations for the bases
    #define complement_case(X,Y)            \
        template<>                          \
        struct complement<X>{               \
            static constexpr char value = Y;\
        }
    complement_case('A','T');
    complement_case('T','A');
    complement_case('G','C');
    complement_case('C','G');
    complement_case('a','t');
    complement_case('t','a');
    complement_case('g','c');
    complement_case('c','g');
    complement_case('U','A');
    complement_case('u','a');
    #undef complement_case

    //Lookup table class, inherits array, specialized constructor to make lookuptable compile time
    template< size_t N = 127>
    class complement_lookup : public std::array<char, N>{
        public:
        template<char... args>
        constexpr complement_lookup (std::integer_sequence<char, args...> const&) : 
        std::array<char,N> {complement<args>::value...}
        {}
    };
    //Builds lookup table
    template<size_t N = 127>
    constexpr auto make_complement_table(){
        return complement_lookup<N>{std::make_integer_sequence<char, N> {}};
    }
#ifdef DEBUG
    //This wouldn't compile if make_complement_table is not computed during compilation
    static_assert(make_complement_table() == make_complement_table());
#endif
    inline void complement_inplace(std::string &seq){
        static constexpr complement_lookup table = make_complement_table();    
        std::reverse(seq.begin(), seq.end());
        for( char &c : seq){                 
            c = table[c];          
        }                                    
    }
    //Inheriting string to add specialized complement constructor
    class complement_string: public std::string{
        static constexpr complement_lookup table = make_complement_table();    
        public:
            template<class InputIt>
            constexpr complement_string(InputIt begin, InputIt end){
                reserve(end-begin);
                for(InputIt it = begin; it != end; ++it){
                    push_back(table[*it]);
                }
            }
            complement_string(const complement_string &seq){
                reserve(seq.crend()-seq.crbegin());
                for(auto it = seq.crbegin(); it != seq.crend(); ++it){
                    push_back(table[*it]);
                }
            }
            complement_string(const std::string &seq){
                reserve(seq.crend()-seq.crbegin());
                for(auto it = seq.crbegin(); it != seq.crend(); ++it){
                    push_back(table[*it]);
                }
            }
    };
    inline std::string complement_seq(const std::string &seq){
        return complement_string{seq.rbegin(), seq.rend()};                                    
    }
}

#endif
