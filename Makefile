


SRC_PATH = src
BUILD_PATH = build
PY_HEADER_PATH = py_header
EXTERN_HEADER_PATH = extern/include
BIN_PATH = ${BUILD_PATH}/bin

INSTALL_PREFIX ?= ${BUILD_PATH}
CXXFLAGS += -DINSTALL_PATH=\"${INSTALL_PREFIX}/bin\"

GIT_VERSION:=$(shell git describe --dirty --always --tags)
ifneq ($(GIT_VERSION),"")
	CXXFLAGS += -DVERSION=\"${GIT_VERSION}\"
endif
CXX?=g++


CXX_OPT ?= -O2
CXX_DBG ?=
CXX_STD ?=c++20
CXXFLAGS += -std=$(CXX_STD) -Wall -Werror $(CXX_OPT) $(CXX_DBG) 
LDFLAGS += -lz -lpthread -lstdc++fs  -lfmt


PY_CXXFLAGS = $(shell python3-config --cflags --embed)
PY_LDFLAGS = $(shell python3-config --ldflags --embed)

MAIN = $(SRC_PATH)/tksm.cpp
EXEC = $(BIN_PATH)/tksm

PY_FILES = $(wildcard py/*.py)
PY_HEADERS = $(PY_FILES:py/%.py=py_header/%.h)

$(EXEC): $(MAIN) $(PY_HEADERS) install.sh
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(PY_CXXFLAGS) -I$(PY_HEADER_PATH) -I$(EXTERN_HEADER_PATH) -o $@ $< $(PY_LDFLAGS) $(LDFLAGS)

$(BIN_PATH)/%: $(SRC_PATH)/%.cpp $(PY_HEADERS) 
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(PY_CXXFLAGS) -I$(PY_HEADER_PATH) -I$(EXTERN_HEADER_PATH) -o $@ $< $(PY_LDFLAGS) $(LDFLAGS)

py_header/%.h: py/%.py
	@mkdir -p py_header
	@xxd -i $< > $@

install.sh: ${EXEC_FILES} Makefile
	@echo mkdir -p ${INSTALL_PREFIX}/bin > $@
	@echo cp ${EXEC_FILES} ${INSTALL_PREFIX}/bin >> $@
	@echo cp py/badread_models ${INSTALL_PREFIX}/bin -r >> $@
	chmod +x $@

.PHONY: clean
clean:
	rm -rf ${OBJD} ${BUILD_PATH} ${TSTB} ${PY_HEADER_PATH} ${BIN_PATH} install.sh

compile_flags.txt: compile_flags.txt.pre
	@cat $< | grep -v -- "-flto-partition=none" | sort | uniq > $@
	rm $<

compile_flags.txt.pre: Makefile
	@echo -stdlib=libc++ > $@
	@echo -xc++ >> $@
	@echo -std=$(CXX_STD) >> $@
	@echo -Wall >> $@
	@echo -Werror >> $@
	@echo -I$(PY_HEADER_PATH) >> $@
	@echo -I$(EXTERN_HEADER_PATH) >> $@
	@echo $(CXXFLAGS) | sed 's/isystem /I/g' | sed -e "s/\s/\n/g" >> $@
	@echo $(LDFLAGS)  | sed 's/isystem /I/g' | sed -e "s/\s/\n/g" >> $@
	@echo $(PY_CXXFLAGS) | sed 's/isystem /I/g'  | sed -e "s/\s/\n/g" >> $@
	@echo $(PY_LDFLAGS)  | sed 's/isystem /I/g' | sed -e "s/\s/\n/g" >> $@
	@echo -fexperimental-library >> $@


# Testing 
TSTD = test
TSTB = ${TSTD}/build
TEST_SOURCES = $(wildcard ${TSTD}/*.cpp)
TEST_BINARIES = $(TEST_SOURCES:${TSTD}/%.cpp=${TSTB}/%)

$(TSTB)/%: $(TSTD)/%.cpp
	@mkdir -p ${TSTB}
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

check: ${TEST_BINARIES}
	for t in $^; do echo $$t && $$t; done