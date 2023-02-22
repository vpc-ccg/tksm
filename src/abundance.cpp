#include "abundance.h"
#include "pimpl_impl.h"
#include <fmt/format.h>

#include <transcript_abundance.h>
#include <cstdlib>
#include <string>

#include "module.h"
#include "python_runner.h"
#include "util.h"

using std::string;

#ifndef TKSM_MODELS_PATH
#define TKSM_MODELS_PATH ""
#endif

class Abundance_module::impl : public tksm_module {
    MAKE_PYTHON_RUNNER(static, run_abundance, py_transcript_abundance_py, NO_MODULES);

    int argc;
    char **argv;

public:
    impl(int argc, char **argv) : tksm_module("Abundance", "Abundance Module"), argc(argc), argv(argv) {}
    int run() {
        if (validate_arguments()) {
            return 1;
        }
        describe_program();
        return run_abundance(argc, argv);
    }
    void describe_program() {}
    int validate_arguments() { return 0; }
};

MODULE_IMPLEMENT_PIMPLE_CLASS(Abundance_module);
