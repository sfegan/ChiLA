include ../Makefile.common
CHILADIR = ..

LIBOBJECTS = VSAMath.o VSAAlgebra.o VSAAtmosphere.o                   \
	VSAReconstruction.o                                           \
	VSAReconstructionMethod1.o VSAReconstructionMethod2.o         \
	VSAReconstructionMethod3.o VSALinearLeastSquares.o            \
	VSASVD.o VSAFunction.o VSAQuadrature.o VSANonlinearFitting.o  \
	VSAStatistics.o

OBJECTS = $(LIBOBJECTS)

TARGETS = libVSAnalysis.a

LIBS =  -lPhysics -lVSUtility -lVSCommon -lmysqlclient -lz

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

libVSAnalysis.a: $(LIBOBJECTS)
	$(AR) r $@ $^

doc: Reconstuction.pdf

Reconstuction.pdf: Reconstruction.tex Geometry.pdf Coherent.pdf
	pdflatex $<
	pdflatex $<

.PHONY: clean

clean:
	$(RM) \
	$(TARGETS) $(OBJECTS) *~ test \
	Geometry.pdf Coherent.pdf Reconstruction.pdf \
	Reconstruction.aux Reconstruction.log Reconstruction.out \
	x.log

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

%.pdf: %.fig
	fig2dev -Lpdf $< $@
