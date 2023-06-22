#include "sequence.h"

#include <fmt/format.h>
#include <sequence.h>
#include <tksm_badread.h>

#include <cstdlib>
#include <string>

#include "module.h"
#include "pimpl.h"
#include "python_runner.h"
#include "util.h"

using std::string;

#ifndef TKSM_MODELS_PATH
#define TKSM_MODELS_PATH ""
#endif

class Sequencer_module::impl : public tksm_module {
    MAKE_PYTHON_RUNNER(static, run_sequencer, py_sequence_py, (ModuleDict{{"tksm_badread", py_tksm_badread_py}}));

    int argc;
    char **argv;

public:
    impl(int argc, char **argv) : tksm_module("sequence", "Sequencer Module"), argc(argc), argv(argv) {}
    ~impl() = default;
    int run() {
        if (validate_arguments()) {
            return 1;
        }
        describe_program();

        // Update envrionment for the pytthon modules directory using the macro definition defined in the makefile

        char *env = getenv("TKSM_MODELS");
        if (env == NULL) {
            if (string{""} == TKSM_MODELS_PATH) {
                loge("TKSM_MODELS not set and TKSM_MODELS_PATH not defined in makefile");
                // Let the python runner handle the error
            }
            logd("TKSM_MODELS not set, setting to {}", std::string(TKSM_MODELS_PATH));
            setenv("TKSM_MODELS", TKSM_MODELS_PATH, 1);
        }
        else {
            logd("TKSM_MODELS was set to {}", std::string(env));
            logd("Appending {}", std::string(TKSM_MODELS_PATH));
            string new_env = fmt::format("{}:{}", TKSM_MODELS_PATH, env);
            setenv("TKSM_MODELS", new_env.c_str(), 1);
        }
        return run_sequencer(argc, argv);
    }
    void describe_program() {}
    int validate_arguments() { return 0; }
};

MODULE_IMPLEMENT_PIMPL_CLASS(Sequencer_module);
