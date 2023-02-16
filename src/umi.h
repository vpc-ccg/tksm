#ifndef _UMI_H
#define _UMI_H
#include <memory>

using std::unique_ptr;

class UMI_module{
    class impl;
    unique_ptr<impl> pimpl;
    public:
    UMI_module(int argc, char **argv);
    ~UMI_module();
    int run();
};
#endif
