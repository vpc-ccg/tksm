#ifndef _TRUNCATE_H_
#define _TRUNCATE_H_

#include "cxxopts/cxxopts.hpp"
#include "interval.h"
#include <npy/npy.hpp>
#include <random>




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
/*
template <class Distribution>
void truncate(molecule_descriptor &md, Distribution &dist, int min_val = 100){
    return truncate(md, dist(rand_gen), min_val);
}
*/

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

            for(size_t j = 0; j < dist.y_axis_index.size(); ++j){
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


class Truncate_module : public tksm_module{

    cxxopts::ParseResult parse(int argc, char **argv){
        options.add_options("main")
            ("i,input", "input mdf file", cxxopts::value<string>())
            ("o,output", "output mdf file", cxxopts::value<string>())
            ("kde", "Comma separated files generated by kde module (grid,labelx,labely)", cxxopts::value<vector<string>>())    
            ("normal", "Use Normal distribution [μ,σ]", cxxopts::value<vector<double>>())
            ("lognormal", "Use Log-Normal distribution [μ,σ]", cxxopts::value<vector<double>>())
            ;
        return  options.parse(argc, argv);
    }
    
    cxxopts::ParseResult args;
    std::mt19937 rand_gen;
    public:
    Truncate_module( int argc, char **argv) : tksm_module{"truncate", "Truncate module"}, args(parse(argc, argv)){
    }

    int  validate_arguments(){
        std::vector<string> mandatory = {"input","output"};
        int missing_parameters = 0;
        for( string &param : mandatory){
            if(args.count(param) == 0){
                std::cerr << param << " is required!\n";
                ++missing_parameters;
            }
        }

        if( args.count("kde") == 0 && args.count("normal") == 0 && args.count("lognormal") == 0){
            std::cerr << "One of kde, normal or lognormal is required!\n";
            ++missing_parameters;
        }
        if( args.count("kde") + args.count("normal") + args.count("lognormal") > 1){
            std::cerr << "Only one of kde, normal or lognormal is allowed!\n";
            ++missing_parameters;
        }
        if(missing_parameters  > 0){
            std::cerr << options.help() << std::endl;
            return 1;
        }
        return 0;
    }
    int run(){
        fmtlog::setLogLevel(LogLevels::parse_loglevel(args["verbosity"].as<string>()));
        fmtlog::flushOn(fmtlog::DBG);

        if(help_or_version_is_used(args)){
            return 0;
        }

        if(validate_arguments()){
            return 1;
        }
        describe_program();

        int seed = args["seed"].as<int>();;
        rand_gen.seed(seed);

        string mdf_file_path {args["input"].as<string>()};
        std::ifstream mdf_file {mdf_file_path};
        auto streamer = stream_mdf(mdf_file);

        string output_file_path {args["output"].as<string>()};
        std::ofstream output_file {output_file_path};

        if ( args["kde"].count() > 0){
            std::vector<std::string> kdefiles = args["kde"].as<vector<string>>();

            custom_distribution2D<> disko {kdefiles[0], kdefiles[1], kdefiles[2]};
            while(streamer){
                molecule_descriptor md = streamer();
                truncate(md, disko(rand_gen, md.size()));
                output_file << md;
            }
        }
        else if( args["normal"].count() > 0){
            std::vector<double> normal_params = args["normal"].as<vector<double>>();
            std::normal_distribution<> disko {normal_params[0], normal_params[1]};
            while(streamer){
                molecule_descriptor md = streamer();
                truncate(md, disko(rand_gen));
                output_file << md;
            }
        }
        else if( args["lognormal"].count() > 0){
            std::vector<double> lognormal_params = args["lognormal"].as<vector<double>>();
            std::lognormal_distribution<> disko {lognormal_params[0], lognormal_params[1]};
            while(streamer){
                molecule_descriptor md = streamer();
                truncate(md, disko(rand_gen));
                output_file << md;
            }
        }

        return 0;   
    }

    void describe_program(){  
        logi("Running Single-cell barcoding module");
        logi("Input MDF: {}", args["input"].as<string>());
        logi("Output MDF: {}", args["output"].as<string>());
        if(args.count("normal")){
            logi("Normal distribution: μ={}, σ={}", args["normal"].as<vector<double>>()[0], args["normal"].as<vector<double>>()[1]);
        }
        if(args.count("lognormal")){
            logi("Log-Normal distribution: μ={}, σ={}", args["lognormal"].as<vector<double>>()[0], args["lognormal"].as<vector<double>>()[1]);
        }
        if(args.count("kde")){
            logi("KDE files: {}", fmt::join(args["kde"].as<vector<string>>(), ", "));
        }

        logi("Seed: {}", args["seed"].as<int>());
        logi("Verbosity: {}", args["verbosity"].as<string>());
        
        fmtlog::poll(true);
    }
};

#endif
