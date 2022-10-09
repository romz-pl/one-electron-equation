#
# Makefile for RSchr-5.x program for XL compiler on Bleu Gene/P NOTOS
#

RSCHRBIN = rschr.x

.PRECIOUS: $(RSCHRBIN)
# .SILENT:
#



PETSCLIB  = -L../lib/petsc-3.3-p5/aa-aa/lib/  -lpetsc 
PETSCINC  = -I../lib/petsc-3.3-p5/aa-aa/include   -I../lib/petsc-3.3-p5/include

SLEPCLIB  = -L../lib/slepc-3.3-p3/aa-aa/lib/  -lslepc
SLEPCINC  = -I../lib/slepc-3.3-p3/aa-aa/include   -I../lib/slepc-3.3-p3/include


#PETSCLIB  = -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc
#PETSCINC  = -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include

#SLEPCLIB  = -L${SLEPC_DIR}/${PETSC_ARCH}/lib -lslepc
#SLEPCINC  = -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${SLEPC_DIR}/include

PARMETISLIB  = -L../lib/parmetis-4.0.2/build/aa-aa/libparmetis -lparmetis -L../lib/parmetis-4.0.2/build/aa-aa/libmetis -lmetis
PARMETISINC  = -I../lib/parmetis-4.0.2/include -I../lib/parmetis-4.0.2/metis/include

#BLASLIB    = -L../libs/lapack-3.4.1/  -lblas
BLASLIB    = -lessl
LAPACKLIB  = -L../lib/lapack-3.4.1 -llapack

CXX       = mpCC_r
CXXFLAGS  = -q64 -DNDEBUG -O3 -qarch=pwr7 -qtune=pwr7 -qmaxmem=-1 -qsimd=auto -qhot=level=2 -qcache=auto -qsuppress=1500-036
#CXXFLAGS  =  -g -O -Wall
CXXINCL   = $(SLEPCINC) $(PETSCINC) $(PARMETISINC)
CXXLIB    = $(SLEPCLIB) $(PETSCLIB) $(PARMETISLIB) $(BLASLIB) $(LAPACKLIB) $(BLASLIB) -lm -lxlf90 -lX11 -lxlsmp




SOURCE = bspline.cpp   \
	cubat.cpp      \
	eigsol.cpp     \
	graph.cpp      \
	graphbspl.cpp  \
	main.cpp       \
	ovrkinelt.cpp  \
	paramdb.cpp    \
	potcachehash.cpp   \
	potcachetree.cpp   \
	potelt.cpp     \
	poten.cpp      \
	recmtx.cpp     \
	schr.cpp       \
	util.cpp

OBJECT = $(SOURCE:.cpp=.o)


all : $(RSCHRBIN)


# $(OBJECT) : Makefile


$(RSCHRBIN): $(OBJECT)
	$(CXX) -q64 $^ $(CXXLIB) $(SYSLIB) -o $@
	mv ./$(RSCHRBIN) ../bin/$(RSCHRBIN)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(CXXINCL) -c $*.cpp


clean:
	rm -f *.o $(RSCHRBIN) ../bin/$(RSCHRBIN)

