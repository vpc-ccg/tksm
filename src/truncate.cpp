#include <cstddef>
#include <exception>
#include <numeric>
#include <stdexcept>
#include <string>
#include <random>
#include <tuple>
#include <variant>
#include <vector>
#include <map>
#include <fstream>
#include <set>
#include <sstream>
#include <climits>

#include <numbers>
#include <functional>
#include <cmath>

#define EPSILON 2.2204460492503131e-16
#define LOG_DBL_EPSILON   (-3.6043653389117154e+01)
#define LOG_DBL_MIN   (-7.0839641853226408e+02)
#define LOG_DBL_MAX    7.0978271289338397e+02

#include <npy/npy.hpp>
#include <cgranges/IITree.h>
#include <cxxopts/cxxopts.hpp>
#include "interval.h"
#include "tree.h"
#include "gtf.h"
#include "mdf.h"
#include "paf.h"

#define FMT_HEADER_ONLY 1
#include <fmt/format.h>


using std::ofstream;
using std::vector;
using std::string;
using std::map;

std::mt19937 rand_gen{std::random_device{}()};

constexpr double approx_digamma(double n){
    return std::log(n) - 1.0 / (2 * n);
}

double approx_trigamma(double n){
    using namespace std::numbers;
    int max_iter = 10000;
    double sum = 0;
    double sump = sum;
    double h = 0.0000001;
    for( int i = 0; i < max_iter; ++i){
        sum += 1.0 / std::pow((n + i),2);
        if( std::abs(sump - sum) < h){
            break;
        }
        sump = sum;
    }
    return sum;
}

/* coefficients for Maclaurin summation in hzeta() Yanked from gsl
 * B_{2j}/(2j)!
 */
static double hzeta_c[15] = {
  1.00000000000000000000000000000,
  0.083333333333333333333333333333,
 -0.00138888888888888888888888888889,
  0.000033068783068783068783068783069,
 -8.2671957671957671957671957672e-07,
  2.0876756987868098979210090321e-08,
 -5.2841901386874931848476822022e-10,
  1.3382536530684678832826980975e-11,
 -3.3896802963225828668301953912e-13,
  8.5860620562778445641359054504e-15,
 -2.1748686985580618730415164239e-16,
  5.5090028283602295152026526089e-18,
 -1.3954464685812523340707686264e-19,
  3.5347070396294674716932299778e-21,
 -8.9535174270375468504026113181e-23
};



double hzeta( double s, double q){
    if( s <= 1 || q <= 0){
        throw std::exception();
    }
    double max_bits = 54;
    double ln_term0 = -s * std::log(q);
   
    if(ln_term0 < LOG_DBL_MIN + 1.0) {
        throw std::exception();
    }
    else if(ln_term0 > LOG_DBL_MAX - 1.0) {
        throw std::exception();
    }

    if((s > max_bits && q < 1) || (s > 0.5 * max_bits && q < 0.25)){
        return std::pow(q, -s );
    }
    else if( s > 0.5 * max_bits && q < 1.0){
        double p1 = std::pow(q, -s);
        double p2 = std::pow(q/(1.0+q), s);
        double p3 = std::pow(q/(2+q), s);
        return p1 * (1 + p2 + p3);
    }
    else{
        int jmax = 12;
        int kmax = 10;
        double pmax = std::pow(kmax + q, -s );
        double scp = s;
        double pcp = pmax / (kmax + q);
        double ans = pmax*((kmax+q)/(s-1) + 0.5);
        for(int k = 0; k < kmax; ++k){
            ans += std::pow(k + q, -s);
        }
        for(int j = 0; j < jmax; ++j){
            double delta = hzeta_c[j+1] * scp *pcp;
            ans += delta;
            if(std::abs(delta/ans) < 0.5 * EPSILON){
                break;
            }
            scp *= (s+2*j+1)*(s+2*j+2);
            pcp /= (kmax+q) *(kmax +q);
        }
        return ans;
    }   
}

double psi_n(int n, double x){
    if( n == 0){
        return approx_digamma(x);
    }
    double z = hzeta(n+1, x);
    double lf = std::lgamma(n+1);
    double e = std::exp(lf) * z;
    if( n%2==0){
        e = -e;
    }
    return e;
}
double psip(double x){
    if( x==0 || x == -1 || x == -2){
        throw std::exception();
    }
    else if( x > 0){
        return psi_n(0, x);
    }
    else if( x > -5){
        int M = - std::floor(x);
        double fx = x + M;
        double sum = 0;


        if(fx == 0){
            throw std::exception();
        }
        
        for(int m = 0; m < M; ++m){
            sum += 1.0 / std::pow(x+m,2);
        }
        
        int stat_psi = psi_n(1, fx);
        stat_psi += sum;
        return stat_psi;
    }
    else{
        using namespace std::numbers;
        double sin_px = std::sin(pi * x);
        double d = pi * pi / (sin_px * sin_px);
        int psi = psi_n(1, 1-x);
        return d - psi;
    }
}
//Yanked from gsl
double polygamma(int n, double x){
    if(n==0){
        return psip(x);//Update
    }
    if(n==1){
        return psi_n(1,x);
    }
    else{ // I dont need atm
        return 0;
    }
}



template <class Distribution, class MeanFoo>
void multi_truncate(molecule_descriptor &md, MeanFoo mu_func, double sigma, Distribution &dist){
    int rand_val = dist(rand_gen);
    
    size_t i = 0;
    int len_so_far = 0;
    auto &segments = md.get_segments();
    for( const ginterval &g : segments){
        len_so_far += (g.end - g.start);
        if(len_so_far >= rand_val){
            break;
        }
        ++i;
    }

    if( i != segments.size()){
        std::stringstream ss;
        ss << segments[i].chr << ':' << segments[i].end - (len_so_far - rand_val) << '-' << segments[i].end;
        md.add_comment("truncated", ss.str());
        segments[i].end -= (len_so_far - rand_val); //Update last segment;
        for(size_t j = i + 1;j < segments.size();++j){
            std::stringstream ss;
            ss << segments[j];
            md.add_comment("truncated", ss.str()); 
        }
        segments.resize(i+1);
    }

    i = 0;
    for( const auto &p : md.errors_so_far){
        if(p.first > rand_val){
            break;
        }
        ++i;
    }
    md.errors_so_far.resize(i);
}



void truncate(molecule_descriptor &md, int rand_val, int min_val = 100){
    if(min_val > rand_val){
        rand_val = min_val;
    }
    size_t i = 0;
    int len_so_far = 0;
    auto &segments = md.get_segments();
    for( const ginterval &g : segments){
        len_so_far += (g.end - g.start);
        if(len_so_far >= rand_val){
            break;
        }
        ++i;
    }
    if( i != segments.size()){

        std::stringstream ss;
        ss << segments[i].chr << ':' << segments[i].end - (len_so_far - rand_val) << '-' << segments[i].end;
        md.add_comment("truncated", ss.str());

        segments[i].end -= (len_so_far - rand_val); //Update last segment;
        for(size_t j = i + 1;j < segments.size();++j){
            std::stringstream ss;
            ss << segments[j];
            md.add_comment("truncated", ss.str());
        }
        segments.resize(i+1);
    }

    i = 0;
    for( const auto &p : md.errors_so_far){
        if(p.first > rand_val){
            break;
        }
        ++i;
    }
    md.errors_so_far.resize(i);
}

template <class Distribution>
void truncate(molecule_descriptor &md, Distribution &dist, int min_val = 100){
    return truncate(md, dist(rand_gen), min_val);
}

struct identity{
    auto operator () (auto val) const {
        return val;
    }
};

template<class Accessor>
auto mean_std(const vector<int> &lengths, Accessor accessor){
    double mean = 0;
    int i = 0;

    for(int val : lengths){
        mean = mean + (accessor(val) - mean) / (i+1);
        ++i;
    }
    double varsum = 0;
    for(int val : lengths){
        varsum = varsum + std::pow( accessor(val)- mean,2) ;
    }
    
    double stdev = std::sqrt(varsum / lengths.size());
    
    return std::make_tuple(mean, stdev);
}

auto mean_std(const vector<int> &lengths){
    return mean_std(lengths, identity{});
}

auto fit(const vector<int> &lengths){

    using namespace std::numbers;

    map<string, std::function<std::tuple<double, double, double> (const vector<int> &)>> ll_computers = {
        { 
            "normal" , [] (const vector<int> &ll){
                auto [_mean, std] = mean_std(ll);
                double mean = _mean; //This is needed because clangd complains 
                size_t n = ll.size();
                double var_cal = std::accumulate(ll.begin(), ll.end(), 0.0L, [mean] (double sum, int val) -> double {
                    return sum + (val - mean) * (val - mean); 
                        });
                double LL = -(n/2.0)*std::log(2 * pi) -(n/2.0) *std::log(std*std) - (1.0/(2*std*std)) * var_cal;
                return std::make_tuple(LL, mean, std);
            }
        },
        { 
            "lognormal" , [] (const vector<int> &ll){
                auto [_mean, std] = mean_std(ll, [](int val){ return std::log(val);});
                double mean = _mean; //This is needed because clangd complains 
                size_t n = ll.size();
                double logsum = std::accumulate(ll.begin(), ll.end(), 0.0L, [] (double sum, int val) -> double {
                    return sum + std::log(val); 
                        });
                double var_cal = std::accumulate(ll.begin(), ll.end(), 0.0L, [mean] (double sum, int val) -> double {
                    return sum + (std::log(val) - mean) * (std::log(val) - mean); 
                        });

                double LL =  - logsum -(n/2.0)*std::log(2 * pi) -(n/2.0) *std::log(std*std) - (1.0/(2*std*std)) * var_cal;
                return std::make_tuple(LL, mean, std);
            }
        },
        {//https://tminka.github.io/papers/minka-gamma.pdf
            "gamma", [] (const vector<int> &ll){
                auto [_mean, std] = mean_std(ll);
                double mean = _mean; //This is needed because clangd complains 
                size_t n = ll.size();
                double logsum = std::accumulate(ll.begin(), ll.end(), 0.0L, [] (double sum, int val) -> double {
                    return sum + std::log(val); 
                        });
                double var_cal = std::accumulate(ll.begin(), ll.end(), 0.0L, [mean] (double sum, int val) -> double {
                    return sum + (val - mean) * (val - mean); 
                        });


                auto l_p = [n,mean,logsum] (double a) ->double {
                    //assert(!std::isnan(std::log(a/mean)));
                    //assert(!std::isnan(approx_digamma(a)));
                    std::cerr << "#\t" << n << " * " << std::log(a/mean) << " - " << n << " * " << polygamma(0,a) << " + " << logsum << "\n";
                    return n * std::log(a/mean) - n * polygamma(0,a) + logsum;
                };
                auto l_pp = [n] (double a) ->double {
                    //assert(!std::isnan(approx_trigamma(a)));
                    std::cerr << "!\t" << n << " / " << a << " - " << n << " * " << polygamma(0,a) << "\n";
                    return n / a - n * approx_trigamma(a);
                };
                
                int max_iter =  100000;
                double a_p = (n * mean * mean) / var_cal; 
                std::cerr << "iter: 0\ta: " << a_p << "\n";
                for(int i =0; i < max_iter; ++i){
                    double der1 = l_p(a_p);
                    double der2 = l_pp(a_p);

                    double a_n = a_p -  der1/der2;
                    std::cerr << "iter: " << i << "\tder1: " << der1 << "\tder2: " << der2 << "\ta: " << a_n << "\t";
                    std::cerr << polygamma(0,a_p) << "\t" << polygamma(0,a_p) << "\t" << logsum << "\n";
                    double convergence = std::abs(a_n - a_p);
                    if(convergence < 0.00000001){
                        break;
                    }
                    if(std::isnan(a_n)){
                        std::cerr << "a is nan at " << i << "'th iteration\n";
                        break;
                    }
                    a_p = a_n;
                }
                double b_p = mean / a_p;
                double LL = (a_p - 1) * logsum - n * std::lgamma(a_p) - n * a_p * std::log(b_p) - n * mean / b_p;
                return std::make_tuple(LL, a_p,b_p);
            }
        }
    };

    string max_dist = "None";
    double max_ll = std::numeric_limits<double>::lowest();
    double chosen_mean = -1;
    double chosen_std  = -1;
    for( auto &p : ll_computers){
        auto [LL, mean, std] = p.second(lengths);
        std::cerr << p.first << ": ";
        std::cerr << "ll: " << LL << " "; 
        std::cerr << "mean: " << mean << " "; 
        std::cerr << "std: " << std << "\n";
        if(LL > max_ll){
            max_dist = p.first;
            max_ll   = LL;
            chosen_mean = mean;
            chosen_std = std;
        }
    }
    return make_tuple(max_dist, max_ll, chosen_mean, chosen_std);
}


auto multi_variate(const vector<int> &rlens, const vector<int> &tlens){
    auto [dist, ll, mean, std] = fit(rlens);
    double mu1 = mean;
    double sig1 = std;

    auto [tmean, tstd] = mean_std(tlens);
    double mu2 = tmean;
    double sig2 = tstd;

    auto tlen_iter = tlens.begin();


    double xy_sum = 0;
    double x_sum = 0;
    double y_sum = 0;
    for( auto rlen_iter = rlens.begin(); rlen_iter != rlens.end(); ++rlen_iter, ++tlen_iter){
        double xdif = (*rlen_iter - mean);
        double ydif = (*tlen_iter - tmean);
        xy_sum += xdif * ydif;
        x_sum  += xdif * xdif;
        y_sum  += ydif * ydif;
    }

    double ro = xy_sum / (std::sqrt(x_sum) * std::sqrt(y_sum));

    auto mean_adjuster = [ro, mu1, mu2, sig1, sig2] (double alpha){
        return mu1 + (sig1/sig2) * ro * (alpha - mu2);
    };
    double adjusted_sig = std::sqrt(( 1 - (ro*ro)) * sig1*sig1);
    return std::make_tuple(dist, mean_adjuster, adjusted_sig);
}



template<class RealType = double, class IndexType = long>
class custom_distribution{
    using result_type = RealType;

    std::uniform_real_distribution<> uniform_dist;

    vector<RealType> pdfv;
    vector<RealType> cdfv;
    vector<IndexType> values;
    public:


    template< class IterType, class IterTypeV>
    custom_distribution (IterType pd_beg, IterType pd_end,
                         IterTypeV val_beg, IterTypeV val_end) :
        uniform_dist{0.0,1},
        pdfv{pd_beg, pd_end},
        cdfv{0},
        values{val_beg, val_end}{
        double sum_pdf = std::accumulate(pd_beg, pd_end, 0.0L); 
        for(double d: pdfv){
            cdfv.push_back(d/sum_pdf+cdfv.back());
        }
    }
    template< class IterType>
    custom_distribution (IterType pd_beg, IterType pd_end,
                         const vector<IndexType> &values) :
        uniform_dist{0.0,1},
        pdfv{pd_beg, pd_end},
        cdfv{0},
        values{values}{
        double sum_pdf = std::accumulate(pd_beg, pd_end, 0.0L); 
        for(double d: pdfv){
            cdfv.push_back(d/sum_pdf+cdfv.back());
        }
    }
    custom_distribution () {}
    custom_distribution (const vector<RealType> &pdf,
                             const vector<RealType> &values) :
            uniform_dist{0.0,1},
            pdfv{pdf},
            cdfv{0},
            values{values}{
            double sum_pdf = std::accumulate(pdf.begin(), pdf.end(), 0.0L); 
            for(double d: pdf){
                cdfv.push_back(d/sum_pdf+cdfv.back());
            }
        }

        template<class Generator>
        result_type operator() ( Generator &g){
            double val = uniform_dist(g);
            auto iter = std::lower_bound(cdfv.cbegin(), cdfv.cend(), val);
            auto val_index = std::distance(cdfv.cbegin(), iter);
            --val_index; // Because we add 0 to cdf

            int range_before = 0;
            int range_after  = 0;
            //Uniformize
            bool not_first = val_index != 0;
            bool not_last =  val_index != ((long int)values.size() - 1);
            
            if(not_first){
                bool odd_before = (values[val_index] - values[val_index-1]) % 2 == 1;

                if( odd_before){
                    range_before = (values[val_index] - values[val_index-1]) / 2; //Rounds down
                }
                else{ // Same if added for clarity
                    range_before = (values[val_index] - values[val_index-1]) / 2;
                }
            }
            if(not_last){
                bool odd_after = (values[val_index+1] - values[val_index]) % 2 == 1;
                if( odd_after){
                    range_after = (values[val_index+1] - values[val_index]) / 2; //Rounds down
                }
                else{
                    range_after = (values[val_index+1] - values[val_index]) / 2 - 1;
                }

            }
            std::uniform_int_distribution <int> smoother_dist(-range_before, range_after);
            double ret = values[val_index] + smoother_dist(g);
//            std::cout << ret << "\t" << values[val_index] << "\t" << range_before << "\t" << range_after << "\n";

            return ret;
        }
        
        const vector<result_type> &cdf() const {
            return cdfv;
        }
        const vector<result_type> &pdf() const {
            return pdfv;
        }
};

template<class RealType = double, class IndexType = long>
class custom_distribution2D{
    using result_type = RealType;
    vector<IndexType> x_axis_index;
    vector<IndexType> y_axis_index;
    vector<custom_distribution<>> distillery;
    public:
        custom_distribution2D (const string &pfile, const string &x_axis_file, const string &y_axis_file){
            vector<unsigned long> shape;

            bool fortran_order;


            npy::LoadArrayFromNumpy(x_axis_file, shape, fortran_order, x_axis_index);
            npy::LoadArrayFromNumpy(y_axis_file, shape, fortran_order, y_axis_index);

            vector<double> pdata;
            npy::LoadArrayFromNumpy(pfile, shape, fortran_order, pdata);

            IndexType width = shape[0];

            if(fortran_order){
                for(size_t i = 0; i < y_axis_index.size(); ++i){
                    vector<double> ps;
                    for(size_t j = 0; j < std::min(shape[1], i+1); ++j){
                        ps.push_back(pdata[j * shape[0] + i]);
                    }
                    distillery.emplace_back(ps.cbegin(), ps.cend() , x_axis_index.begin(), x_axis_index.end()); 
                }
            }
            else{
                for(size_t i = 0; i < y_axis_index.size();++i){
                    distillery.emplace_back(pdata.cbegin() + (i * width), pdata.cbegin() + ((i) * width + i + 1) , x_axis_index.begin(), x_axis_index.end()); 
                }
            }

        }

        template<class Generator>
        result_type operator()( Generator &g, IndexType slice){
            auto iter = std::lower_bound(y_axis_index.begin(), y_axis_index.end(), slice);
            if( iter != y_axis_index.begin() && std::abs(*iter-slice) > std::abs(*(iter-1)-slice)){
                --iter;
            }

            return distillery.at(std::distance(y_axis_index.begin(),iter))(g);
        }

        vector<result_type> cdf( IndexType slice){
            auto iter = std::lower_bound(y_axis_index.begin(), y_axis_index.end(), slice);
            if( iter != y_axis_index.begin() && std::abs(*iter-slice) > std::abs(*(iter-1)-slice)){
                --iter;
            }
            return distillery.at(std::distance(y_axis_index.begin(), iter)).cdf();
        }
        vector<result_type> pdf( IndexType slice){
            auto iter = std::lower_bound(y_axis_index.begin(), y_axis_index.end(), slice);
            if( iter != y_axis_index.begin() && std::abs(*iter-slice) > std::abs(*(iter-1)-slice)){
                --iter;
            }
            return distillery.at(std::distance(y_axis_index.begin(), iter)).pdf();
        }

        friend ostream &operator << ( ostream &ost, const custom_distribution2D<RealType, IndexType> &dist){
            for(IndexType i : dist.x_axis_index){
                ost << i << "\t";
            }
            ost << "\n";

            for(int j = 0; j < dist.y_axis_index.size(); ++j){
                ost << dist.y_axis_index[j] << "\t";
                const custom_distribution<RealType, IndexType> &d1 =dist.distillery[j];
                for(double v : d1.cdf()){
                    ost << v << "\t";
                }
                ost << "\n";
            }
            return ost;
        }
};




inline auto cdf(const vector<double> &pdf) -> vector<double> {
    vector<double> cdf_vec{0};
    double sum_pdf = std::accumulate(pdf.begin(), pdf.end(), 0.0L);
    for(double d: pdf){
        cdf_vec.push_back(d/sum_pdf+cdf_vec.back());
    }

    return cdf_vec;
}


inline auto sample_from_cdf(const vector<double> &cdf, auto &rand_gen)  {
    std::uniform_real_distribution<> dist(cdf.front(), cdf.back());
    double val = dist(rand_gen);
    auto iter = std::lower_bound(cdf.begin(), cdf.end(), val);
    return std::distance(cdf.cbegin(), iter);
}

int main_print_mtx(int argc, char **argv){

    custom_distribution2D<> disko {argv[1], argv[2], argv[3]};
    std::default_random_engine rg{};
    std::cout << disko<< "\n";

    return 0;
}

int main(int argc, char **argv){
    cxxopts::Options options("RNAInfuser Truncate", "Truncate module of RNAInfuser");
    std::cout << approx_trigamma(2) << "\n";
    options.add_options()
        ("i,input", "Molecule description file", cxxopts::value<string>())
        ("g,gtf", "Path to GTF annotation", cxxopts::value<string>())
        ("o,output", "Output path", cxxopts::value<string>())
        ("fit", "Use estimated distribution from the given fast{a,q}", cxxopts::value<string>())
        ("only-single-isoform", "Use only single isoform genes", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("multivariate", "Use estimated distribution from the given paf", cxxopts::value<string>())
        ("kde", "Use Kernel Density Estimation given outputs of the kde.py comma separated", cxxopts::value<vector<string>>())
        ("normal", "Use Normal distribution [μ,σ]", cxxopts::value<vector<double>>())
        ("lognormal", "Use Log-Normal distribution [μ,σ]", cxxopts::value<vector<double>>())
        ("seed", "Random seed", cxxopts::value<int>()->default_value("42"))
        ("h,help", "Help screen")
        ;
    auto args = options.parse(argc, argv);

    if(args.count("help") > 0){
        std::cout << options.help() << std::endl;
        return 0;
    }
    vector<string> mandatory = {"input", "output", "gtf"};

    int missing_parameters = 0;
    for( string &param : mandatory){
        if(args.count(param) == 0){
            std::cerr << param << " is required!\n";
            ++missing_parameters;
        }
    }
    if(missing_parameters  > 0){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    int seed = args["seed"].as<int>();;

    rand_gen.seed(seed);
    vector<string> distributions = {"lognormal", "normal", "fit", "multivariate", "kde"};


    string chosen_dist;
    int arg_dist_count = 0;
    for( const string &dist_name : distributions){
        if( args.count(dist_name) > 0 ){
            ++arg_dist_count;
            chosen_dist = dist_name;
        }
    }

    if(arg_dist_count != 1){
        std::cerr << "Please provide one distribution type!\n";
        return 1;
    }


    string path_to_gtf   {args["gtf"].as<string>()};
    ifstream gtfile( path_to_gtf);

    string buffer;
    map<string, vector<string>> gene2transcripts;
    while(std::getline(gtfile, buffer)){
        if(buffer[0] == '#'){
            continue;
        }
        gtf entry{buffer};

        if(entry.type == gtf::entry_type::transcript && entry.info["transcript_biotype"] == "protein_coding"){
            string gid = entry.info["gene_id"];
            string tid = entry.info["transcript_id"];
            gene2transcripts[gid].push_back(tid);
        }
    }

    std::set<string> single_isoform_transcripts;

    for( const auto &kv : gene2transcripts){
        if(kv.second.size() == 1){
            single_isoform_transcripts.insert(kv.second.back());
        }

    }
    vector<double> dist_params;


    if( chosen_dist == "kde"){

        const vector<string> &kdefiles = args["kde"].as<vector<string>>();
        custom_distribution2D<> disko {kdefiles[0], kdefiles[1], kdefiles[2]};
        string mdf_file_path {args["input"].as<string>()};
        std::ifstream mdf_file {mdf_file_path};
        auto streamer = stream_mdf(mdf_file, true);
        string outfile_name = args["output"].as<string>();
        std::ofstream outfile{outfile_name};
        
        while(streamer){
            auto md = streamer();
            truncate(md, disko(rand_gen, md.size()));
            outfile << md;
        }
        return 0;
    }
    else if( chosen_dist == "multivariate"){
        string dist_fit_file_path {args["multivariate"].as<string>()};
        std::ifstream dist_fit_file(dist_fit_file_path);

        string buffer;
        vector<int> lens;
        vector<int> tlens;

        
        while(std::getline(dist_fit_file, buffer)){ 
            paf pf{buffer};
            if(!pf.primary()){
                continue;
            }
            if(args["only-single-isoform"].as<bool>() && single_isoform_transcripts.find(pf.tname) == single_isoform_transcripts.end()){
                continue;
            }
            lens.push_back( pf.qend - pf.qstart);
            tlens.push_back( pf.tlen);

        }

        auto [distname, mu_func, sigma] = multi_variate(lens, tlens);
        
        std::cerr << distname << " is chosen!\n";
        string mdf_file_path {args["input"].as<string>()};
        

        string outfile_name = args["output"].as<string>();
        std::ofstream outfile{outfile_name};

        std::ifstream mdf_file {mdf_file_path};
        auto streamer = stream_mdf(mdf_file, true);

        while(streamer){
            auto md = streamer();           
            std::variant<
                std::normal_distribution<>,
                std::lognormal_distribution<>,
                std::gamma_distribution<>> dist;
            std::normal_distribution<> normal_dist(mu_func(md.size()), sigma);
            std::lognormal_distribution<> lognormal_dist(mu_func(md.size()), sigma);
            if( distname == "normal"){
                dist = std::normal_distribution<>(mu_func(md.size()), sigma);
            }
            else if( distname == "lognormal"){
                dist = std::lognormal_distribution<>(mu_func(md.size()), sigma);
            }
            else if( distname == "gamma"){
                dist = std::gamma_distribution<>(mu_func(md.size()), sigma);
            }
            else{
                std::cerr << distname << " not implemented!\n";
                return -1;
            }
            std::visit([&md](auto dist){
                    truncate(md, dist);
            }, dist);
            
            outfile << md;
        }

        return 0;   
    }
    else if( chosen_dist == "fit"){
        string dist_fit_file_path {args["fit"].as<string>()};
        std::ifstream dist_fit_file(dist_fit_file_path);

        string buffer;
        vector<int> lens;
        

        
        while(std::getline(dist_fit_file, buffer)){ 
            paf pf{buffer};
            if(!pf.primary()){
                continue;
            }
            if(args["only-single-isoform"].as<bool>() && single_isoform_transcripts.find(pf.tname.substr(0,15)) == single_isoform_transcripts.end()){
                continue;
            }
            lens.push_back( pf.qend - pf.qstart);


        }

        auto [name, ll, mean, std] = fit(lens);
        std::cerr << fmt::format("{} is chosen with {} log-likelihood\n", name, ll);
        chosen_dist = name;
        dist_params = {mean, std};

    }
    else{
        dist_params = args[chosen_dist].as<vector<double>>();
        if( dist_params.size() != 2){
            std::cerr << fmt::format("{} requires 2 parameters! Make sure there is no space between them.\n", chosen_dist);
            return 1;
        }

    }


    std::normal_distribution<> normal_dist(dist_params[0], dist_params[1]);

    std::lognormal_distribution<> lognormal_dist(dist_params[0], dist_params[1]);
    std::gamma_distribution<> gamma_dist(dist_params[0], dist_params[1]);

    string outfile_name = args["output"].as<string>();
    std::ofstream outfile{outfile_name};

    string mdf_file_path {args["input"].as<string>()};
    std::ifstream mdf_file {mdf_file_path};
    auto streamer = stream_mdf(mdf_file, true);
    while(streamer){
        auto md = streamer();
        if(chosen_dist == "normal"){
            truncate(md, normal_dist);
        }
        else if(chosen_dist == "lognormal"){
            truncate(md, lognormal_dist);
        }
        else if(chosen_dist == "gamma"){
            truncate(md, gamma_dist);
        }
        outfile << md;
    }
    return 0;
}
