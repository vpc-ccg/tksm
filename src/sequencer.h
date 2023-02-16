#ifndef _SEQUENCER_H_
#define _SEQUENCER_H_

#include <sequence.h>

#include "python_runner.h"
#include "module.h"
#include <cstdlib>
#include <string>

using std::string;


#ifndef TKSM_MODELS_PATH
#define TKSM_MODELS_PATH ""
#endif

class Sequencer_module : public tksm_module{
    MAKE_PYTHON_RUNNER(static, sequencer, py_sequence_py);

    int argc;
    char **argv;

    public:
        Sequencer_module(int argc, char **argv) : tksm_module("sequencer", "Sequencer Module"),
            argc(argc), argv(argv){}
        int run(){
            if(validate_arguments()){
                return 1;
            }
            describe_program();
            

            //Update envrionment for the pytthon modules directory using the macro definition defined in the makefile
            
            char *env = getenv("TKSM_MODELS");
            if(env == NULL){
                if(string{""} == TKSM_MODELS_PATH){
                    loge("TKSM_MODELS not set and TKSM_MODELS_PATH not defined in makefile");
                    return 1;
                }
                logd("TKSM_MODELS not set, setting to {}", std::string(TKSM_MODELS_PATH));
                setenv("TKSM_MODELS", TKSM_MODELS_PATH, 0);
            }else{
                logd("TKSM_MODELS was set to {}", std::string(env));
                logd("Appending {}", std::string(TKSM_MODELS_PATH));
                string new_env = string{env} + ":" + TKSM_MODELS_PATH;
                setenv("TKSM_MODELS", new_env.c_str(), 0);
            }
            run_sequencer(argc, argv);
            return 0;
        }
        void describe_program() {
            //TODO we need to have a python call for this
        }
        int validate_arguments(){
            //TODO we need to have a python call for this
            return 0;
        }

};



#endif
