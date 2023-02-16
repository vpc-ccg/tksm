#ifndef _PIMPLE_H
#define _PIMPLE_H
#include <memory>

#define MODULE_DECLARE_PIMPLE_CLASS(CLASS_NAME)     \
    class CLASS_NAME {                     \
        class impl;                                 \
        std::unique_ptr<impl> pimpl;                \
                                                    \
    public:                                         \
        CLASS_NAME(int argc, char **argv); \
        ~CLASS_NAME();                     \
                                                    \
        int run();                                  \
    };                                              \
    static_assert(true)

#endif
