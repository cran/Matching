ifeq ($(CXX),g++)
#CXXFLAGS = -O3 -funroll-loops -fstrict-aliasing -fsched-interblock -falign-loops=16 -falign-jumps=16 -falign-functions=16 -falign-jumps-max-skip=15 -falign-loops-max-skip=15 -malign-natural -ffast-math -mpowerpc-gpopt -force_cpusubtype_ALL -fstrict-aliasing -mtune=G5 -mcpu=G5 -mpowerpc64
#CXXFLAGS = -O3 -ffast-math -funroll-loops -fexpensive-optimizations -ftree-vectorize -fdump-tree-vect-stats -msse -msse2 -march=pentium4 -mtune=pentium4
#CXXFLAGS = -O3 -ansi -pedantic -Wall
CXXFLAGS = -O3 -ffast-math -funroll-loops -fexpensive-optimizations 
endif

SOURCES = matching.cc scythematrix.cc 

HEADERS = matching.h scythematrix.h 

OBJS = $(SOURCES:.cc=.o)

PROGRAM = Matching${SHLIB_EXT}

#$(PROGRAM):	$(OBJS) $(LIBS)
#		$(SHLIB_CXXLD) $(DYLIB_LDFLAGS) $(OBJS) -o $(PROGRAM)
#		@echo library $(PROGRAM) ready.

ALL_CXXFLAGS = $(R_XTRA_CXXFLAGS) $(PKG_CXXFLAGS) $(CXXPICFLAGS) $(SHLIB_CXXFLAGS) $(CXXFLAGS)

$(PROGRAM):	$(OBJS) 
		$(SHLIB_CXXLD) $(SHLIB_CXXLDFLAGS) $(ALL_CXXFLAGS) $(OBJS) -o $(PROGRAM)
		@echo library $(PROGRAM) ready.

clean:
	rm -f core*
	rm -f *.o





