#ifndef __PYTHON_RUNNER_H__
#define __PYTHON_RUNNER_H__
#define PY_SSIZE_T_CLEAN
#include <string.h>

#include <map>
#include <string>
using std::map;
using std::string;
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

class ModuleDict {
    map<string, const char *> modules;

public:
    ModuleDict(const std::initializer_list<std::pair<const char *, const unsigned char *>> &modules) {
        for (const auto &module : modules) {
            this->modules[module.first] = reinterpret_cast<const char *>(module.second);
        }
    }
    ModuleDict() {}

    operator map<string, const char *>() const { return modules; }

    auto begin() const { return modules.begin(); }
    auto end() const { return modules.end(); }
    auto cbegin() const { return modules.cbegin(); }
    auto cend() const { return modules.cend(); }
};

#define ARRAY_LENGTH(X) X##_len

#define NO_MODULES \
    ModuleDict {}
#define MAKE_PYTHON_RUNNER(SCOPE, FUNC_NAME, ARRAY_NAME, MODULE_MAP)                                   \
    SCOPE int FUNC_NAME(int argc, char **argv_s) {                                                     \
        extern unsigned int ARRAY_LENGTH(ARRAY_NAME);                                                  \
        char program_string[ARRAY_LENGTH(ARRAY_NAME) + 1];                                             \
        strncpy(program_string, reinterpret_cast<const char *>(ARRAY_NAME), ARRAY_LENGTH(ARRAY_NAME)); \
        program_string[ARRAY_LENGTH(ARRAY_NAME)] = 0;                                                  \
        wchar_t **argv                           = new wchar_t *[argc];                                \
        for (int i = 0; i < argc; i++) argv[i] = (wchar_t *)GetWC(argv_s[i]);                          \
        wchar_t *program = Py_DecodeLocale("", NULL);                                                  \
        if (program == NULL) {                                                                         \
            fprintf(stderr, "Fatal error: cannot decode argv[0]");                                     \
            exit(1);                                                                                   \
        }                                                                                              \
        Py_SetProgramName(program);                                                                    \
        Py_Initialize();                                                                               \
        for (const auto module_pair : MODULE_MAP) {                                                    \
            PyObject *module     = PyImport_AddModule(module_pair.first.c_str());                      \
            PyObject *moduleDict = PyModule_GetDict(module);                                           \
            PyRun_String(module_pair.second, Py_file_input, moduleDict, moduleDict);                   \
        }                                                                                              \
        PySys_SetArgvEx(argc, argv, 0);                                                                \
        if (PyRun_SimpleString(program_string) != 0) {                                                 \
            PyErr_Print();                                                                             \
            return 1;                                                                                  \
        }                                                                                              \
        if (Py_FinalizeEx() < 0) {                                                                     \
            return 120;                                                                                \
        }                                                                                              \
        PyMem_RawFree(program);                                                                        \
        return 0;                                                                                      \
    }                                                                                                  \
    static_assert(true, "")

#define MAKE_PYTHON_RUNNER_NO_MODULES(SCOPE, FUNC_NAME, ARRAY_NAME) \
    MAKE_PYTHON_RUNNER(SCOPE, FUNC_NAME, ARRAY_NAME, NO_MODULES)

#endif
