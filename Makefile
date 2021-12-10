
CXX ?= g++
CXXFLAGS += -O2  -Wall -std=c++20 
LDFLAGS += -lz -lpthread

SRCD=src
OBJD=obj
BUILDD=build
PCH_LIB = libs.h
PCH_HEADERS = tree.h graph.h interval.h reverse_complement.h cigar.h extern/IITree.h extern/cxxopts.h
SOURCE_FILES = fusion.cpp pcr.cpp truncation_analysis.cpp sequencer.cpp
SOURCE_PATH = $(SOURCE_FILES:%.cpp=${SRCD}/%.cpp)
EXEC_FILES = $(SOURCE_FILES:%.cpp=${BUILDD}/%)
OBJECTS = $(SOURCE_FILES:%.cpp=${OBJD}/%.o)
PCH_OUT = $(PCH_LIB:%.h=${OBJD}/%.pch)
HEADER_PATH = $(PCH_HEADERS:%.h=${SRCD}/%.h)

all: $(EXEC_FILES)

$(BUILDD)/%:$(SRCD)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(PCH_OUT): $(PCH_LIB) $(PCH_HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<

${OBJD}/%.o:$(SRCD)/%.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $<


all:
	@echo ${SOURCE_PATH}
	@echo ${OBJECTS}
	@echo ${PCH_LIB}
	@echo ${PCH_OUT}
	@echo ${HEADER_PATH}
	@echo ${EXEC_FILES}