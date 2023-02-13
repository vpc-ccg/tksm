
#include "sequencer.h"

#ifndef INSTALL_PATH
#define INSTALL_PATH "."
#endif


int main(int argc, char **argv){
    const char *path = "--badread-model-path=" INSTALL_PATH "/badread_models";
    char **argv_cpy = (char **)malloc((argc+2)*sizeof(char *));
    argv_cpy[0] = argv[0];
    for(int i=1; i<argc; i++){
        argv_cpy[i + 1] = argv[i];
    }
    argv_cpy[1] = (char *)path;
    return run_sequencer(argc+1,  argv_cpy);
}
