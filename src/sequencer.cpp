#include "sequencer.h"
#include "pimpl_impl.h"
#include <fmt/format.h>
#include <sequence.h>

#include <cstdlib>
#include <string>

#include "module.h"
#include "python_runner.h"
#include "util.h"

using std::string;

#ifndef TKSM_MODELS_PATH
#define TKSM_MODELS_PATH ""
#endif

class Sequencer_module::impl : public tksm_module {
    MAKE_PYTHON_RUNNER(static, sequencer, py_sequence_py, NO_MODULES);

    int argc;
    char **argv;

public:
    impl(int argc, char **argv) : tksm_module("sequencer", "Sequencer Module"), argc(argc), argv(argv) {}
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
            setenv("TKSM_MODELS", TKSM_MODELS_PATH, 0);
        }
        else {
            logd("TKSM_MODELS was set to {}", std::string(env));
            logd("Appending {}", std::string(TKSM_MODELS_PATH));
            string new_env = fmt::format("{}:{}", TKSM_MODELS_PATH, env);
            setenv("TKSM_MODELS", new_env.c_str(), 0);
        }
        return run_sequencer(argc, argv);
    }
    void describe_program() {}
    int validate_arguments() { return 0; }
};

MODULE_IMPLEMENT_PIMPLE_CLASS(Sequencer_module);
