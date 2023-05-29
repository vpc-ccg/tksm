#include "truncate.h"

#include <cxxopts.hpp>
#include <npy/npy.hpp>
#include <random>
#include <variant>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "pimpl.h"

inline void
truncate(molecule_descriptor &md, int truncated_length, int min_val = 100) {
    if (min_val > truncated_length) {
        truncated_length = min_val;
    }
    size_t i       = 0;
    int len_so_far = 0;
    auto &segments = md.get_segments();
    for (const ginterval &g : segments) {
        len_so_far += (g.end - g.start);
        if (len_so_far >= truncated_length) {
            break;
        }
        ++i;
    }
    if (i != segments.size()) {
        std::stringstream ss;
        ss << segments[i].chr << ':' << segments[i].end - (len_so_far - truncated_length) << '-' << segments[i].end;
        md.add_comment("truncated", ss.str());
        segments[i].truncate(0, segments[i].end - segments[i].start - (len_so_far - truncated_length));
        for (size_t j = i + 1; j < segments.size(); ++j) {
            std::stringstream ss;
            ss << segments[j];

            md.add_comment("truncated", ss.str());
        }
        logd("Before resize: {}", md.cget_segments().size());
        segments.resize(i + 1);
        logd("After resize: {}", md.cget_segments().size());
    }
    int sum = 0;
    for (const auto &g : md.cget_segments()) {
        sum += g.end - g.start;
    }
    logd("sum: {} truncated_len: {}", sum, truncated_length);
}

template <class RealType = double, class IndexType = long>
class custom_distribution {
    using result_type = RealType;

    std::uniform_real_distribution<> uniform_dist;

    vector<RealType> pdfv;
    vector<RealType> cdfv;
    vector<IndexType> bins;
    vector<std::uniform_int_distribution<int>> smoother_distros;

public:
    // clang-format off
    template <class IterType, class IterTypeV>
    custom_distribution(IterType pdf_beg, IterType pdf_end, IterTypeV bins_beg, IterTypeV bins_end) :
            uniform_dist{0.0, 1},
            pdfv{pdf_beg, pdf_end},
            cdfv{0},
            bins{bins_beg, bins_end} {
        double sum_pdf = std::accumulate(pdf_beg, pdf_end, 0.0L);
        for (double d : pdfv) {
            cdfv.push_back(d / sum_pdf + cdfv.back());
        }
        smoother_distros.reserve(bins.size());
        smoother_distros.emplace_back(0, bins.front());
        for (auto it = bins.begin(); std::next(it) != bins.end(); ++it){
            int current = *it;
            int next = *(it + 1);
            smoother_distros.emplace_back(current, next);
        }
    }
    // clang-format on
    template <class Generator>
    result_type operator()(Generator &g) {
        auto cdf_index_iter = std::lower_bound(cdfv.cbegin(), cdfv.cend(), uniform_dist(g));
        auto cdf_index      = std::distance(cdfv.cbegin(), cdf_index_iter);
        auto bin_index      = cdf_index - 1;  // Because we add 0 to cdf
        return smoother_distros[bin_index](g);
    }

    const vector<result_type> &cdf() const { return cdfv; }
    const vector<result_type> &pdf() const { return pdfv; }
};

template <class RealType = double, class IndexType = long>
class custom_distribution2D {
    using result_type = RealType;
    vector<IndexType> x_axis_index;
    vector<IndexType> y_axis_index;
    vector<custom_distribution<>> distillery;

public:
    custom_distribution2D(const string &pfile, const string &x_axis_file, const string &y_axis_file) {
        vector<unsigned long> shape;

        bool fortran_order;

        npy::LoadArrayFromNumpy(x_axis_file, shape, fortran_order, x_axis_index);
        npy::LoadArrayFromNumpy(y_axis_file, shape, fortran_order, y_axis_index);

        vector<double> pdata;
        npy::LoadArrayFromNumpy(pfile, shape, fortran_order, pdata);

        IndexType width = shape[0];

        if (fortran_order) {
            for (size_t i = 0; i < y_axis_index.size(); ++i) {
                vector<double> ps;
                for (size_t j = 0; j < std::min(shape[1], i + 1); ++j) {
                    ps.push_back(pdata[j * shape[0] + i]);
                }
                distillery.emplace_back(ps.cbegin(), ps.cend(), x_axis_index.begin(), x_axis_index.end());
            }
        }
        else {
            for (size_t i = 0; i < y_axis_index.size(); ++i) {
                distillery.emplace_back(pdata.cbegin() + (i * width), pdata.cbegin() + ((i)*width + i + 1),
                                        x_axis_index.begin(), x_axis_index.end());
            }
        }
    }

    template <class Generator>
    result_type operator()(Generator &g, IndexType slice) {
        auto iter = std::lower_bound(y_axis_index.begin(), y_axis_index.end(), slice);
        if (iter != y_axis_index.begin() && std::abs(*iter - slice) > std::abs(*(iter - 1) - slice)) {
            --iter;
        }
        logd("Slice: {}, index: {}", slice, std::distance(y_axis_index.begin(), iter));
        return distillery.at(std::distance(y_axis_index.begin(), iter))(g);
    }

    vector<result_type> cdf(IndexType slice) {
        auto iter = std::lower_bound(y_axis_index.begin(), y_axis_index.end(), slice);
        if (iter != y_axis_index.begin() && std::abs(*iter - slice) > std::abs(*(iter - 1) - slice)) {
            --iter;
        }
        return distillery.at(std::distance(y_axis_index.begin(), iter)).cdf();
    }
    vector<result_type> pdf(IndexType slice) {
        auto iter = std::lower_bound(y_axis_index.begin(), y_axis_index.end(), slice);
        if (iter != y_axis_index.begin() && std::abs(*iter - slice) > std::abs(*(iter - 1) - slice)) {
            --iter;
        }
        return distillery.at(std::distance(y_axis_index.begin(), iter)).pdf();
    }

    friend ostream &operator<<(ostream &ost, const custom_distribution2D<RealType, IndexType> &dist) {
        for (IndexType i : dist.x_axis_index) {
            ost << i << "\t";
        }
        ost << "\n";

        for (size_t j = 0; j < dist.y_axis_index.size(); ++j) {
            ost << dist.y_axis_index[j] << "\t";
            const custom_distribution<RealType, IndexType> &d1 = dist.distillery[j];
            for (double v : d1.cdf()) {
                ost << v << "\t";
            }
            ost << "\n";
        }
        return ost;
    }
};

// Visit helpers
template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;
// End visit helpers
//
class Truncate_module::impl : public tksm_module {
    cxxopts::ParseResult parse(int argc, char **argv) {
        // clang-format off
        options.add_options("main")
            (
                "i,input",
                "input mdf file",
                cxxopts::value<string>()
            )(
                "o,output",
                "output mdf file",
                cxxopts::value<string>()
            )(
                "kde-model",
                "Comma separated files generated by model-truncate utility (grid,labelx,labely)",
                cxxopts::value<vector<string>>())   
                (
                "normal",
                "Use Normal distribution [μ,σ]",
                cxxopts::value<vector<double>>()
            )(
                "lognormal",
                "Use Log-Normal distribution [μ,σ]",
                cxxopts::value<vector<double>>()
            );
        // clang-format on
        return options.parse(argc, argv);
    }

    cxxopts::ParseResult args;

public:
    impl(int argc, char **argv) : tksm_module{"truncate", "Truncate module"}, args(parse(argc, argv)) {}
    ~impl() = default;

    int validate_arguments() {
        std::vector<string> mandatory = {"input", "output"};
        int missing_parameters        = 0;
        for (string &param : mandatory) {
            if (args.count(param) == 0) {
                report_missing_parameter(param);
                ++missing_parameters;
            }
        }

        if (args.count("kde-model") == 0 && args.count("normal") == 0 && args.count("lognormal") == 0) {
            std::cerr << "One of kde-model, normal or lognormal is required!\n";
            ++missing_parameters;
        }
        if (args.count("kde-model") + args.count("normal") + args.count("lognormal") > 1) {
            std::cerr << "Only one of kde-model, normal or lognormal is allowed!\n";
            ++missing_parameters;
        }
        if (missing_parameters > 0) {
            std::cerr << options.help() << std::endl;
            return 1;
        }
        return 0;
    }

    auto get_dist()
        -> std::variant<custom_distribution2D<>, std::normal_distribution<>, std::lognormal_distribution<>> {
        if (args["kde-model"].count() > 0) {
            std::vector<std::string> kdefiles = args["kde-model"].as<vector<string>>();
            return custom_distribution2D<>{kdefiles[0], kdefiles[1], kdefiles[2]};
        }
        else if (args["normal"].count() > 0) {
            std::vector<double> normal_params = args["normal"].as<vector<double>>();
            return std::normal_distribution<>{normal_params[0], normal_params[1]};
        }
        else {  //(args["lognormal"].count() > 0){ // skipped because we already check this and lambda needs to always
                // return
            std::vector<double> lognormal_params = args["lognormal"].as<vector<double>>();
            return std::lognormal_distribution<>{lognormal_params[0], lognormal_params[1]};
        }
    }

    static auto truncate_transformer(auto &disko, auto &rand_gen) {
        return std::ranges::views::transform([&](auto &md) {
            auto truncate_length =
                std::visit(overloaded{[&](auto &arg) { return arg(rand_gen); },
                                      [&](custom_distribution2D<> &arg) { return arg(rand_gen, md.size()); }},
                           disko);
            truncate(md, truncate_length);
            return md;
        });
    }

    auto operator()() {
        auto dist = get_dist();
        return truncate_transformer(dist, rand_gen);
    }

    int run() {
        if (process_utility_arguments(args)) {
            return 0;
        }

        if (validate_arguments()) {
            return 1;
        }
        describe_program();

        string output_file_path{args["output"].as<string>()};
        std::ofstream output_file{output_file_path};
        auto disko = get_dist();

        //        std::ranges::copy( stream_mdf(args["input"].as<string>(), true) | truncate_transformer(disko),
        //                std::ostream_iterator<molecule_descriptor>{output_file});
        for (const auto &md : stream_mdf(args["input"].as<string>(), true) | truncate_transformer(disko, rand_gen)) {
            output_file << md;
        }

        return 0;
    }

    void describe_program() {
        logi("Running Truncate module");
        logi("Input MDF: {}", args["input"].as<string>());
        logi("Output MDF: {}", args["output"].as<string>());
        if (args.count("normal")) {
            logi("Normal distribution: μ={}, σ={}", args["normal"].as<vector<double>>()[0],
                 args["normal"].as<vector<double>>()[1]);
        }
        if (args.count("lognormal")) {
            logi("Log-Normal distribution: μ={}, σ={}", args["lognormal"].as<vector<double>>()[0],
                 args["lognormal"].as<vector<double>>()[1]);
        }
        if (args.count("kde-model")) {
            logi("KDE files: {}", fmt::join(args["kde-model"].as<vector<string>>(), ", "));
        }

        logi("Seed: {}", args["seed"].as<int>());
        logi("Verbosity: {}", args["verbosity"].as<string>());

        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Truncate_module);
