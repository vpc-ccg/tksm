CXX ?= g++
CXXFLAGS += -O2  -Wall -std=c++20 
LDFLAGS += -lz
OBJ = main.o #more objects here eventually I would think!

PCH_LIB = libs.h
PCH_HEADERS = tree.h graph.h interval.h reverse_complement.h IITree.h cigar.h
PCH_OUT = $(PCH_LIB:.h=.pch)
#PCH_INC = $(addprefix -include , ${PCH_HEADERS})

main: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ ${LDFLAGS}

# Compiles your PCH
$(PCH_OUT):  $(PCH_LIB) $(PCH_HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<

# the -include flag instructs the compiler to act as if lib.hpp
# were the first header in every source file
%.o: %.cpp $(PCH_OUT) $(OBJ_DIR)
	$(CXX) $(CXXFLAGS)  -include $(PCH_LIB) -c -o $@ $<

.PHONY:  clean
clean:
	${RM} *.o
	${RM} main
	${RM} *.pch
