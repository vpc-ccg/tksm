

INSTALL_PREFIX ?= /usr
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
LDFLAGS += -lz -lpthread -lstdc++fs 


PY_CXXFLAGS = $(shell python3-config --cflags --embed)
PY_LDFLAGS = $(shell python3-config --ldflags --embed)


SRC_PATH = src
BUILD_PATH = build
PY_HEADER_PATH = py_header
EXTERN_HEADER_PATH = extern/include
BIN_PATH = bin

MAIN = $(SRC_PATH)/tksm.cpp
EXEC = $(BIN_PATH)/tksm
PY_FILES = $(wildcard py/*.py)
PY_HEADERS = $(PY_FILES:py/%.py=py_header/%.h)

$(EXEC): $(MAIN) $(PY_HEADERS) install.sh
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(PY_CXXFLAGS) -I$(PY_HEADER_PATH) -I$(EXTERN_HEADER_PATH) -o $@ $< $(PY_LDFLAGS)

py_header/%.h: py/%.py
	@mkdir -p py_header
	@xxd -i $< > $@

install.sh: ${EXEC_FILES} Makefile
	@echo mkdir -p ${INSTALL_PREFIX}/bin > $@
	@echo cp ${BIN_PATH}/* ${INSTALL_PREFIX}/bin >> $@
	@echo cp py/badread_models ${INSTALL_PREFIX}/bin -r >> $@
	chmod +x $@

.PHONY: clean
clean:
	rm -rf ${OBJD} ${BUILDD} ${TSTB} ${PY_HEADER_PATH} ${BIN_PATH} install.sh

# Testing 
#

TSTD = test
TSTB = ${TSTD}/build
${TSTB}/%:${TSTD}/%.cpp
	@mkdir -p ${TSTB}
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

reverse_complement_test: ${TSTB}/reverse_complement_test
	./${TSTB}/reverse_complement_test

check: reverse_complement_test
