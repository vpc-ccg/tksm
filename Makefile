

INSTALL_PREFIX ?= /usr
CXXFLAGS += -DINSTALL_PATH=\"${INSTALL_PREFIX}/bin\"

CXX ?= g++
OPT ?= -O2
CXXFLAGS += ${OPT} -Wall -std=c++20 -Iextern/include -Ipy_header
LDFLAGS += -lz -lpthread 
LDFLAGS += $(shell python3-config --ldflags --embed)
CXXFLAGS += $(shell python3-config --cflags --embed)

#GIT_VERSION:=$(shell git describe --dirty --always --tags)
GIT_VERSION:=""
ifneq ($(GIT_VERSION),"")
	CXXFLAGS += -DVERSION=\"${GIT_VERSION}\"
endif

SRCD=src
OBJD=obj
TSTD=test
TSTB=test/binaries
BUILDD=build
PCH_LIB = libs.h
PCH_HEADERS = tree.h graph.h interval.h reverse_complement.h cigar.h extern/IITree.h extern/cxxopts.h kde.h
SOURCE_FILES =  pcr.cpp  sequencer.cpp splicer.cpp polyA.cpp truncate.cpp umi.cpp single-cell-barcoder.cpp kde.cpp

SOURCE_PATH = $(SOURCE_FILES:%.cpp=${SRCD}/%.cpp)
EXEC_FILES = $(SOURCE_FILES:%.cpp=${BUILDD}/%)
OBJECTS = $(SOURCE_FILES:%.cpp=${OBJD}/%.o)
PCH_OUT = $(PCH_LIB:%.h=${OBJD}/%.pch)
HEADER_PATH = $(PCH_HEADERS:%.h=${SRCD}/%.h)
PY_FILES = py/truncate_kde.py py/transcript_abundance.py py/sequence.py
PY_HEADERS = $(PY_FILES:py/%.py=py_header/%.h)


all: $(EXEC_FILES) install.sh

$(BUILDD)/%:$(SRCD)/%.cpp $(PY_HEADERS) 
	@mkdir -p ${BUILDD}
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(PCH_OUT): $(PCH_LIB) $(PCH_HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<

${OBJD}/%.o:$(SRCD)/%.cpp 
	@mkdir -p ${OBJD}
	$(CXX) $(CXXFLAGS) -o $@ $<

${TSTB}/%:${TSTD}/%.cpp
	@mkdir -p ${TSTB}
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

py_header/%.h: py/%.py
	@mkdir -p py_header
	@xxd -i $< > $@

reverse_complement_test: ${TSTB}/reverse_complement_test
	./${TSTB}/reverse_complement_test

check: reverse_complement_test
	
all:
	@echo ${SOURCE_PATH}
	@echo ${OBJECTS}
	@echo ${PCH_LIB}
	@echo ${PCH_OUT}
	@echo ${HEADER_PATH}
	@echo ${EXEC_FILES}

.PHONY: clean

install.sh:
	@echo mkdir -p ${INSTALL_PREFIX}/bin > $@
	@echo cp ${BUILDD}/* ${INSTALL_PREFIX}/bin >> $@
	@echo cp py/badread_models ${INSTALL_PREFIX}/bin -r >> $@
	chmod +x $@
clean:
	@rm -f ${BUILDD}/*
	@rm -f ${OBJD}/*.o
	@rm -f ${PY_HEADERS}
