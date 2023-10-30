# Check lld exists

LLD_EXISTS := $(shell command -v ld.lld 2> /dev/null)

ifdef LLD_EXISTS
	CXXFLAGS += -fuse-ld=lld
endif

SRC_PATH = src
BUILD_PATH = build
OBJ_PATH = ${BUILD_PATH}/obj
PY_HEADER_PATH = py_header
EXTERN_HEADER_PATH = extern/include
BIN_PATH = ${BUILD_PATH}/bin

INSTALL_PREFIX ?= ${BUILD_PATH}
CXXFLAGS += -DINSTALL_PATH=\"${INSTALL_PREFIX}/bin\"

GIT_VERSION:=$(shell git describe --dirty --always --tags)
ifneq ($(GIT_VERSION),"")
	CXXFLAGS += -DVERSION=\"${GIT_VERSION}\"
endif
TKSM_MODELS_PATH:=${INSTALL_PREFIX}/bin/tksm_models
ifneq ($(TKSM_MODELS_PATH),"")
	CXXFLAGS += -DTKSM_MODELS_PATH=\"${TKSM_MODELS_PATH}\"
endif
CXX?=g++

ifneq ($(DEBUG),1)
	CXX_OPT ?= -O3
else
	CXX_OPT ?= -O0 -g
endif
CXX_STD ?=c++20
CXXFLAGS += -std=$(CXX_STD) -Wall  $(CXX_OPT) $(CXX_DBG) 
LDFLAGS += -lz -lpthread -lfmt $(LD_DBG) -flto=auto


ifneq ($(DEBUG),1)
	PY_CXXFLAGS = $(shell python3-config --cflags --embed | sed 's|-flto-partition=none||g')
	PY_LDFLAGS = $(shell python3-config --ldflags --embed | sed 's|-flto-partition=none||g')
else
	PY_CXXFLAGS = $(shell python3-config --cflags --embed | sed 's|-flto-partition=none||g' | sed 's|-O3|-O0|g')
	PY_LDFLAGS = $(shell python3-config --ldflags --embed | sed 's|-flto-partition=none||g' | sed 's|-O3|-O0|g')
endif



MAIN = $(SRC_PATH)/tksm.cpp
EXEC = $(BIN_PATH)/tksm


SRC_FILES = tksm.cpp tag.cpp truncate.cpp transcribe.cpp scb.cpp sequence.cpp polyA.cpp pcr.cpp model_truncation.cpp abundance.cpp strand_man.cpp filter.cpp random_wgs.cpp shuffle.cpp unsegment.cpp mutate.cpp

#Append SRC_PATH to SRC_FILES
SRC_FILES := $(addprefix $(SRC_PATH)/,$(SRC_FILES))

OBJECTS = $(SRC_FILES:$(SRC_PATH)/%.cpp=$(OBJ_PATH)/%.o)

PCH_HEADER = $(SRC_PATH)/headers.h
PCH_OBJECT = $(PCH_HEADER:$(SRC_PATH)/%.h=$(OBJ_PATH)/%.gch)

PY_FILES = $(wildcard py/*.py)
PY_HEADERS = $(PY_FILES:py/%.py=py_header/%.h)

$(EXEC): $(OBJECTS) $(PY_HEADERS) install.sh $(PCH_OBJECT)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(PY_CXXFLAGS) -I$(PY_HEADER_PATH) -I$(EXTERN_HEADER_PATH) -o $@ $(OBJECTS) $(PY_LDFLAGS) $(LDFLAGS) -include $(PCH_OBJECT)

$(BIN_PATH)/%: $(SRC_PATH)/%.cpp $(PY_HEADERS) 
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(PY_CXXFLAGS) -I$(PY_HEADER_PATH) -I$(EXTERN_HEADER_PATH) -o $@ $< $(PY_LDFLAGS) $(LDFLAGS)

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.cpp $(PY_HEADERS) 
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(PY_CXXFLAGS) -I$(EXTERN_HEADER_PATH) -I$(PY_HEADER_PATH) -c -o $@ $< 

# Compile PCH headers.h
$(OBJ_PATH)/%.gch: $(SRC_PATH)/%.h
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(PY_CXXFLAGS) -I$(EXTERN_HEADER_PATH) -I$(PY_HEADER_PATH) -c -o $@ $<

py_header/%.h: py/%.py
	@mkdir -p py_header
	@echo "#pragma once" > $@
	@xxd -i $< | sed 's|unsigned|inline\ unsigned|g' >> $@

install.sh:  Makefile
	@echo mkdir -p ${INSTALL_PREFIX}/bin > $@
	@echo cp ${EXEC} ${INSTALL_PREFIX}/bin >> $@
	@echo cp py/tksm_models ${INSTALL_PREFIX}/bin -r >> $@
	chmod +x $@

.PHONY: clean
clean:
	rm -rf ${OBJ_PATH} ${BUILD_PATH} ${TSTB} ${PY_HEADER_PATH} ${BIN_PATH} install.sh

compile_flags.txt: compile_flags.txt.pre
	@cat $< | grep -v -- "-flto-partition=none" | sort | uniq > $@
	rm $<

compile_flags.txt.pre: Makefile
	@echo -stdlib=libstdc++ > $@
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
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I$(EXTERN_HEADER_PATH) $(PY_LDLFAGS) $(PY_CXXFLAGS) -o $@ $^

check: ${TEST_BINARIES}
	for t in $^; do echo $$t && $$t; done
