
CXX ?= g++
OPT ?= -O2
CXXFLAGS += ${OPT} -Wall -std=c++20 -Iextern/include
LDFLAGS += -lz -lpthread 

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
PCH_HEADERS = tree.h graph.h interval.h reverse_complement.h cigar.h extern/IITree.h extern/cxxopts.h
SOURCE_FILES = fusion.cpp pcr.cpp  sequencer.cpp splicer.cpp polyA.cpp truncate.cpp umi.cpp single-cell-barcoder.cpp


SOURCE_PATH = $(SOURCE_FILES:%.cpp=${SRCD}/%.cpp)
EXEC_FILES = $(SOURCE_FILES:%.cpp=${BUILDD}/%)
OBJECTS = $(SOURCE_FILES:%.cpp=${OBJD}/%.o)
PCH_OUT = $(PCH_LIB:%.h=${OBJD}/%.pch)
HEADER_PATH = $(PCH_HEADERS:%.h=${SRCD}/%.h)

all: $(EXEC_FILES)

$(BUILDD)/%:$(SRCD)/%.cpp
	mkdir -p ${BUILDD}
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(PCH_OUT): $(PCH_LIB) $(PCH_HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<

${OBJD}/%.o:$(SRCD)/%.cpp 
	mkdir -p ${OBJD}
	$(CXX) $(CXXFLAGS) -o $@ $<

${TSTB}/%:${TSTD}/%.cpp
	mkdir -p ${TSTB}
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

py_header/%.h: py/%.py
	mkdir -p py_header
	xxd -i $< > $@

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

clean:
	rm -f ${BUILDD}/*
	rm -f ${OBJD}/*.o
