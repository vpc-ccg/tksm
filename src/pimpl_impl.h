#ifndef _PIMPLE_IMP_H
#define _PIMPLE_IMP_H



#define MODULE_IMPLEMENT_PIMPLE_CLASS(CLASS_NAME) \
    CLASS_NAME::CLASS_NAME(int argc, char **argv) : pimpl{std::make_unique<impl>(argc, argv)} {} \
    CLASS_NAME::~CLASS_NAME() = default; \
    int CLASS_NAME::run() { return pimpl->run(); }

#endif
