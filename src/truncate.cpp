#include "truncate.h"

#include <fmt/compile.h>
#include <fmt/core.h>

#include <concepts>
#include <cxxopts.hpp>
#include <npy/npy.hpp>
#include <random>
#include <variant>

#include "interval.h"
#include "mdf.h"
#include "module.h"
#include "pimpl.h"

inline void
truncate(molecule_descriptor &md, int post_truncation_length, int min_val = 100) {
    if (post_truncation_length == (int)md.size()) {
        return;
    }
    if (min_val > post_truncation_length) {
        post_truncation_length = min_val;
    }
    size_t i        = 0;
    int kept_so_far = 0;
    auto &segments  = md.get_segments();
    for (const ginterval &g : segments) {
        if (kept_so_far + g.size() >= post_truncation_length) {
            break;
        }
        kept_so_far += g.size();
        ++i;
    }
    if (i != segments.size()) {
        int segment_kept_len = (post_truncation_length - kept_so_far);

        int trunc_start, trunc_end;
        if (segments[i].plus_strand) {
            trunc_start = segments[i].start + segment_kept_len;
            trunc_end   = segments[i].end;
            segments[i].truncate(0, segment_kept_len);
        }
        else {
            trunc_start = segments[i].start;
            trunc_end   = segments[i].end - segment_kept_len;
            segments[i].truncate(segments[i].size() - segment_kept_len, segments[i].size());
        }

        md.add_comment("truncated", fmt::format(FMT_COMPILE("{}:{}-{}"), segments[i].chr, trunc_start, trunc_end));
        for (size_t j = i + 1; j < segments.size(); ++j) {
            md.add_comment("truncated",
                           fmt::format(FMT_COMPILE("{}:{}-{}"), segments[j].chr, segments[j].start, segments[j].end));
        }
        logd("Before resize: {}", md.cget_segments().size());
        segments.resize(i + 1);
        logd("After resize: {}", md.cget_segments().size());
    }
}

template <class N>  //Continiuous
struct dist_picker {
    using type = std::uniform_real_distribution<N>;
};

template <std::integral N> // Discrete
struct dist_picker<N> {
    using type = std::uniform_int_distribution<N>;
};

template <class RealType = double, class IndexType = long>
class custom_distribution {
    using result_type     = RealType;
    using index_type      = IndexType;
    using smoot_dist_type = dist_picker<index_type>::type;
    std::uniform_real_distribution<> uniform_dist;

    vector<RealType> pdfv;
    vector<RealType> cdfv;
    vector<IndexType> bins;
    vector<smoot_dist_type> smoother_distros;

public:
    // clang-format off
    template <class IterType, class IterTypeV>
    custom_distribution(IterType pdf_beg, IterType pdf_end, IterTypeV bins_beg, IterTypeV bins_end) :
            uniform_dist{0.0, 1},
            pdfv{pdf_beg, pdf_end},
            cdfv{0},
            bins{bins_beg, bins_end} {
        double sum_pdf = std::accumulate(pdf_beg, pdf_end, RealType{0});
        for (double d : pdfv) {
            cdfv.push_back(d / sum_pdf + cdfv.back());
        }
        smoother_distros.reserve(bins.size());
        smoother_distros.emplace_back(0, bins.front());
        for (auto it = bins.begin(); std::next(it) != bins.end(); ++it){
            index_type current = *it;
            index_type next = *(it + 1);

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

    template <class Generator>
    result_type operator()(Generator &g, double uniform_random_value) {
        auto cdf_index_iter = std::lower_bound(cdfv.cbegin(), cdfv.cend(), uniform_random_value);
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
    using index_type  = IndexType;
    vector<IndexType> x_axis_index;
    vector<IndexType> y_axis_index;
    vector<custom_distribution<>> distillery;
    std::uniform_real_distribution<> uniform_dist;

public:
    custom_distribution2D(const string &pfile, const string &x_axis_file, const string &y_axis_file) :
        uniform_dist{0.0, 1}{
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
    result_type operator()(Generator &g, IndexType slice, bool smoothed= true) {
        auto iter = std::lower_bound(y_axis_index.begin(), y_axis_index.end(), slice);
        if (iter != y_axis_index.begin() && std::abs(*iter - slice) > std::abs(*(iter - 1) - slice)) {
            --iter;
        }
        logd("Slice: {}, index: {}", slice, std::distance(y_axis_index.begin(), iter));
        size_t didx = std::distance(y_axis_index.begin(), iter);
        double uidx = uniform_dist(g);
        result_type val = distillery.at(didx)(g, uidx);
        if(smoothed && didx+1 < distillery.size()){
            val = (val + distillery.at(didx+1)(g,uidx))/2;
        }
        return val;
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
                "Comma separated files generated by model-truncate utility (grid,labelx,labely,sider)",
                cxxopts::value<vector<string>>()
            )(
                "always-end",
                "Ignore sider.tsv and always 3' end truncate",
                cxxopts::value<bool>()->default_value("false")->implicit_value("true")
            )(
                "kde-models-length",
                "KDE models read length instead of truncation length",
                cxxopts::value<bool>()->default_value("false")->implicit_value("true")
            )(
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
    static auto truncate_transformer_kde(auto &disko, auto &sider_decider, auto &rand_gen, bool always_end,
                                         bool models_length) {
        return std::ranges::views::transform([&, always_end, models_length](auto &md) {
            auto truncate_length = [models_length, &disko, &rand_gen](auto &md) {
                if (models_length) {
                    return md.size() - disko(rand_gen, md.size());
                }
                else {
                    return disko(rand_gen, md.size());
                }
            }(md);

            auto side_ratio = [always_end, &sider_decider, &rand_gen]() {
                if (always_end) {
                    sider_decider(rand_gen);  // Move random seed forward;
                    return 1.0;
                }
                else {
                    return sider_decider(rand_gen);
                }
            }();
            truncate(md, md.size() - truncate_length * side_ratio);
            auto md_reversed = flip_molecule(md);
            truncate(md_reversed, md_reversed.size() - truncate_length * (1 - side_ratio));
            md = flip_molecule(md_reversed);
            md.add_comment("TR", fmt::format("{},{:.2f}", truncate_length, side_ratio));

            return md;
        });
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
        if (args["kde-model"].count() > 0) {
            std::ifstream side_pdf_tsv(args["kde-model"].as<vector<string>>().back());
            vector<double> pdfs;
            vector<double> bins;
            string buffer;
            while (std::getline(side_pdf_tsv, buffer)) {
                std::istringstream bufst{buffer};
                double pdf_val, bin_val;
                bufst >> pdf_val >> bin_val;
                pdfs.push_back(pdf_val);
                bins.push_back(bin_val * 1.0);
            }

            custom_distribution<double, double> sider_decider{pdfs.begin(), pdfs.end(), bins.begin(), bins.end()};

            bool always_end      = args["always-end"].as<bool>();
            bool model_is_length = args["kde-models-length"].as<bool>();
            for (const auto &md : stream_mdf(args["input"].as<string>(), true) |
                                      truncate_transformer_kde(std::get<custom_distribution2D<>>(disko), sider_decider,
                                                               rand_gen, always_end, model_is_length)) {
                output_file << md;
            }
        }
        else {
            for (const auto &md :
                 stream_mdf(args["input"].as<string>(), true) | truncate_transformer(disko, rand_gen)) {
                output_file << md;
            }
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
        logi("Always end truncate: {}", args["always-end"].as<bool>());
        logi("Seed: {}", args["seed"].as<int>());
        logi("Verbosity: {}", args["verbosity"].as<string>());

        fmtlog::poll(true);
    }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Truncate_module);
