include ../Makefile.common
CHILADIR = ..

# -----------------------------------------------------------------------------
# Root
# -----------------------------------------------------------------------------

ifdef ROOTSYS
ROOTBIN := $(ROOTSYS)/bin/
else
ROOTBIN :=
endif

ROOTCONFIG = $(ROOTBIN)root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := 
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs) -lMathCore
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs) -lMathCore

DEFFLAGS += -DNOVDB

CXXFLAGS += $(VBFCFLAGS) $(ROOTCFLAGS)
LDFLAGS += $(VBFLDFLAGS)

LIBS =  -lVSRootAnalysis -lVSStage3 -lVSDataReduction -lVSNSpace -lVSAnalysis \
	-lPhysics -lVSUtility -lSEphem -lVSCommon -lVBF -lm -lz -lbz2 -ldl \
	-lpthread $(ROOTGLIBS) $(MYSQLLIBS)

LIBOBJECTS	= VSRHistogram.o VSRHistogram1D.o VSRHistogram2D.o \
	VSRGraph1D.o VSRHistogramCamera.o VSRH5Loader.o VSRText.o \
	VSRCanvas.o VSRCanvasElement.o VSRTraceVisitor.o VSRPalette.o \
	VSREventVisitor.o VSRVBFVisitor.o

TARGETS         = h5plot traceplot eventplot h5root

all: $(TARGETS)

libVSRootAnalysis.a: $(LIBOBJECTS)
	$(AR) r $@ $^

testMinimization: testMinimization.o libVSRootAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

plotpsf: plotpsf.o libVSRootAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS) -lcfitsio

h5plot: h5plot.o libVSRootAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

h5root: h5root.o libVSRootAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

eventplot: eventplot.o libVSRootAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

pedplot: pedplot.o libVSRootAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

nonlinear_fit: nonlinear_fit.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

nonlinear_fit2: nonlinear_fit2.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

nonlinear_fit3: nonlinear_fit3.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

nonlinear_fit4: nonlinear_fit4.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

lrtest: lrtest.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

extract_channel_diagnostics: extract_channel_diagnostics.o libVSRootAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

traceplot: traceplot.o libVSRootAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

avg_trace_plot: avg_trace_plot.o libVSRootAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)


.PHONY: clean

clean:
	$(RM) *.o *.a $(TARGETS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
