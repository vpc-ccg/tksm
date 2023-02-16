#ifndef __PYTHON_RUNNER_H__
#define __PYTHON_RUNNER_H__
#define PY_SSIZE_T_CLEAN

extern "C" {
#include <Python.h>
}

inline const wchar_t *
GetWC(const char *c) {
    const size_t cSize = strlen(c) + 1;
    wchar_t *wc        = new wchar_t[cSize];
    mbstowcs(wc, c, cSize);

    return wc;
}
#define MAKE_PYTHON_RUNNER(FUNC_NAME, ARRAY_NAME)                                  \
    inline int run_##FUNC_NAME(int argc, char **argv_s) {                          \
        wchar_t **argv = new wchar_t *[argc];                                      \
        for (int i = 0; i < argc; i++) argv[i] = (wchar_t *)GetWC(argv_s[i]);      \
        wchar_t *program = Py_DecodeLocale("", NULL);                              \
        if (program == NULL) {                                                     \
            fprintf(stderr, "Fatal error: cannot decode argv[0]");                 \
            exit(1);                                                               \
        }                                                                          \
        Py_SetProgramName(program);                                                \
        Py_Initialize();                                                           \
        PySys_SetArgvEx(argc, argv, 0);                                            \
        if (PyRun_SimpleString(reinterpret_cast<const char *>(ARRAY_NAME)) != 0) { \
            PyErr_Print();                                                         \
            return 1;                                                              \
        }                                                                          \
        if (Py_FinalizeEx() < 0) {                                                 \
            return 120;                                                            \
        }                                                                          \
        PyMem_RawFree(program);                                                    \
        return 0;                                                                  \
    }                                                                              \
    static_assert(true, "")

#endif
