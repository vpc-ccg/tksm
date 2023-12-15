#ifndef _PIMPL_H
#define _PIMPL_H
#include <memory>

#define MODULE_IMP_PIMPL(CLASS_NAME)                                                             \
    CLASS_NAME::CLASS_NAME(int argc, char **argv) : pimpl{std::make_unique<impl>(argc, argv)} {} \
    CLASS_NAME::~CLASS_NAME() = default;                                                         \
    int CLASS_NAME::run() { return pimpl->run(); }

#ifndef MULTI_BINARY
#define MODULE_IMPLEMENT_PIMPL_CLASS(CLASS_NAME) MODULE_IMP_PIMPL(CLASS_NAME)
#else
#define MODULE_IMPLEMENT_PIMPL_CLASS(CLASS_NAME) \
    MODULE_IMP_PIMPL(CLASS_NAME)                 \
    int main(int argc, char **argv) { return CLASS_NAME{argc, argv}.run(); }
#endif

#define MODULE_DECLARE_PIMPL_CLASS(CLASS_NAME) \
    class CLASS_NAME {                         \
        class impl;                            \
        std::unique_ptr<impl> pimpl;           \
                                               \
    public:                                    \
        CLASS_NAME(int argc, char **argv);     \
        ~CLASS_NAME();                         \
        int run();                             \
    };                                         \
    static_assert(true)

#endif
