#
# Makefile for RSchr-5.x program for GCC compiler
#

RSCHRBIN = rschr.x

.PRECIOUS: $(RSCHRBIN)
# .SILENT:
#


PETSCLIB  = -L../lib/petsc-3.3-p5/aa-aa/lib -lpetsc
PETSCINC  = -I../lib/petsc-3.3-p5/aa-aa/include -I../lib/petsc-3.3-p5/include

SLEPCLIB  = -L../lib/slepc-3.3-p3/aa-aa/lib -lslepc
SLEPCINC  = -I../lib/slepc-3.3-p3/aa-aa/include -I../lib/slepc-3.3-p3/include


PARMETISLIB  = -L../lib/parmetis-4.0.2/build/aa-aa/libparmetis -lparmetis -L../lib/parmetis-4.0.2/build/aa-aa/libmetis -lmetis
PARMETISINC  = -I../lib/parmetis-4.0.2/include -I../lib/parmetis-4.0.2/metis/include


#PARMETISLIB = -L../lib/ParMetis-3.2.0/ -lparmetis -lmetis
#PARMETISINC = -I../lib/ParMetis-3.2.0/ -I../lib/ParMetis-3.2.0/ParMETISLib/

BLASLIB    = -L../lib/lapack-3.4.1 -lblas
LAPACKLIB  = -L../lib/lapack-3.4.1 -llapack



CXX       = mpicxx
CXXFLAGS  = -DNDEBUG -O3 -Wall
# CXXFLAGS  = -DNDEBUG -O3 -Wall -march=native
#CXXFLAGS  = -DNDEBUG -O2 -funroll-all-loops -ffast-math -Wall
#CXXFLAGS  =  -g -O -Wall
CXXINCL   = $(SLEPCINC) $(PETSCINC) $(PARMETISINC)
CXXLIB    = $(SLEPCLIB) $(PETSCLIB) $(PARMETISLIB) $(LAPACKLIB) $(BLASLIB) -lgfortran -lm -ldl




SOURCE = bspline.cpp   \
	cubat.cpp          \
	eigsol.cpp         \
	graph.cpp          \
	graphbspl.cpp      \
	main.cpp           \
	ovrkinelt.cpp      \
	paramdb.cpp        \
	potcachehash.cpp   \
	potcachetree.cpp   \
	potelt.cpp         \
	poten.cpp          \
	recmtx.cpp         \
	schr.cpp           \
	util.cpp

OBJECT = $(SOURCE:.cpp=.o)


all : $(RSCHRBIN)


#
#.depend: depend
#
#depend: $(SOURCE)
#	rm -f ./.depend
#	$(CXX) $(CXXFLAGS) -MM $^  >> ./.depend;
#
#include .depend
#

$(RSCHRBIN): $(OBJECT)
	$(CXX) $^ $(CXXFLAGS) $(CXXLIB) -o $@
	strip ./$(RSCHRBIN)
	mv ./$(RSCHRBIN) ../bin/$(RSCHRBIN)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(CXXINCL) -c $*.cpp


clean:
	rm -f *.o $(RSCHRBIN) ../bin/$(RSCHRBIN)

