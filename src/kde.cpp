#include "kde.h"
#include "pimpl_impl.h"
#include <fmt/format.h>
#include <truncate_kde.h>

#include <cstdlib>
#include <string>

#include "module.h"
#include "python_runner.h"
#include "util.h"

using std::string;


class KDE_module::impl : public tksm_module {
    MAKE_PYTHON_RUNNER(static, kde, py_truncate_kde_py, NO_MODULES);

    int argc;
    char **argv;

public:
    impl(int argc, char **argv) : tksm_module("kde", "Model truncation using KDE"), argc(argc), argv(argv) {}
    int run() {
        if (validate_arguments()) {
            return 1;
        }
        describe_program();


        return run_kde(argc, argv);
    }
    void describe_program() {}
    int validate_arguments() { return 0; }
};

MODULE_IMPLEMENT_PIMPLE_CLASS(KDE_module);
