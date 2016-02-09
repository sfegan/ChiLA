include ../Makefile.common
CHILADIR = ..

DEFFLAGS += -DNOVDB

CXXFLAGS += $(DEFFLAGS) $(VBFCFLAGS)
LDFLAGS += $(VBFLDFLAGS)

LIBS =  -lVSMiscAnalysis -lVSStage3 -lVSDataReduction -lVSNSpace -lVSAnalysis \
	-lPhysics -lSEphem -lVBF -lVSUtility -lVSCommon -lm -lz -lbz2 -ldl \
	-lpthread

LIBOBJECTS = VSCutOptimizer.o VSSimpleCutOptimizer.o VSNSpaceCutOptimizer.o \
	VSSpectrumBiasCalc.o

TARGETOBJECTS = optimize_msc.o test_msc.o polyfit.o stat.o vbf_compress.o

OBJECTS = $(LIBOBJECTS) $(TARGETOBJECTS)

TARGETS = $(TARGETOBJECTS:.o=)

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

libVSMiscAnalysis.a: $(LIBOBJECTS)
	$(AR) r $@ $^

vbf_compress: vbf_compress.o VSPixelCleaningVisitor.o
	$(CXX) $(LDFLAGS) -o $@ $< VSPixelCleaningVisitor.o $(LIBS) \
	 $(MYSQLLIBS) $(H5LIBS)

optimize_msc: optimize_msc.o libVSMiscAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(MYSQLLIBS) $(H5LIBS)

extract_cuts: extract_cuts.o libVSMiscAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

test_msc: test_msc.o libVSMiscAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

polyfit: polyfit.o
	$(CXX) $(LDFLAGS) -o $@ $< -lVSAnalysis -lVSUtility -lVSCommon

stat: stat.o
	$(CXX) $(LDFLAGS) -o $@ $< -lVSAnalysis -lVSUtility -lVSCommon

spbias: spbias.o libVSMiscAnalysis.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

testreg: testreg.o libVSMiscAnalysis.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

.PHONY: clean

clean:
	$(RM) \
	$(TARGETS) $(OBJECTS) *~ test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
