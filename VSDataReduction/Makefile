include ../Makefile.common
CHILADIR = ..

DEFFLAGS += -DNOVDB

CXXFLAGS += $(DEFFLAGS) $(VBFCFLAGS)
LDFLAGS += $(VBFLDFLAGS)

LIBS =  -lVSDataReduction -lVSNSpace -lVSAnalysis -lPhysics -lVSUtility \
        -lSEphem -lVSCommon -lVBF -lm -lz -lbz2 -ldl -lpthread 

LIBOBJECTS = VBFChargeIntegration.o VBFSimplePeds.o VSTimingCalc.o           \
	VSSimpleVBF.o VBFTOffsetCalc.o VSSimplePedData.o VSLaserData.o       \
	VSPointing.o VSCleaner.o VBFRunInfo.o VBFTimeDepPeds.o               \
	VBFTimeDepNSB.o VBFDumper.o VBFPrintFrequency.o VSAnalysisStage1.o   \
	VSTimeDepSupData.o VSTimeDepPedData.o VSAnalysisStage2.o             \
	VBFAnalysisStage2.o VBFAnalysisStage2Diagnostics.o                   \
	VBFAnalysisStage2IED.o VBFAnalysisStage2MergedCal.o                  \
	VSAnalysisData.o VBFSimpleRunInfoExtractor.o VSCentralizedDBAccess.o \
	VSMiscellaneousDBData.o VSRingFitter.o VSMuonAnalysis.o              \
	VSMuonAnalysisData.o VSEventData.o VSDiagnosticsData.o               \
	VSMergedCalibrationData.o VBFLaserCalc.o VSCutsData.o                \
	VSNSpaceCutsCalc.o VSSimulationData.o VSSimpleCutsCalc.o             \
	VSCutsCalc.o  VSScaledParameterLibrary.o VSChannelMap.o              \
	VSStarCatalog.o VSCutsEvaluator.o VSHiLoData.o                       \
	VSEnergyCalcLT.o VSScaledParameterCalc.o VSNSpaceEventData.o         \
	VSSimCoordTransform.o VSSimEnergyWeightCalc.o VBFHiLoCalc.o          \
	VSLTLibraryCalc.o 

LIBOBJECTS_STAGE3 = VSAnalysisStage3.o VSResultsSimData.o                    \
	VSAnalysisStage3Visitor.o VSEventDataVisitor.o VSSourceModel.o       \
	VSBkgndModel.o VSDataModel.o VSSourceInjector.o                      \
	VSAnalysisStage3Data.o VSStage3SimCalc.o VSLTLibraryVisitor.o        \
	VSEffectiveAreaCalc.o VSPSFCalc.o VSInstrumentResponseCalc.o         \
	VSSpectrumCalc.o                                                     \
	VSRingBackgroundAnalysis.o VSIntegralAnalysis.o VSApertureAnalysis.o \
	VSIntegralAnalysisData.o VSMLMAnalysis.o VSReflectedRegionAnalysis.o 

TARGETOBJECTS = eventprint.o collect_muon.o dump_muon.o \
		nspace.o stage1.o stage2.o stage3.o laser.o logain.o \
		make_sp_tables.o make_egy_tables.o \
		extract_diagnostics.o combine_sp_tables.o s2opt.o \
		make_effarea_library.o

OBJECTS = $(LIBOBJECTS) $(LIBOBJECTS_STAGE3) $(TARGETOBJECTS)

TARGETS = libVSDataReduction.a libVSStage3.a $(TARGETOBJECTS:.o=)

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

libVSDataReduction.a: $(LIBOBJECTS)
	$(AR) r $@ $^

libVSStage3.a: $(LIBOBJECTS_STAGE3)
	$(AR) r $@ $^

$(TARGETOBJECTS): libVSDataReduction.a

VSStarCatalog.o: star_catalog_radec_j2000_Vmag_sorted.h

eventprint: eventprint.o libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

toffset: toffset.o libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(MYSQLLIBS) $(H5LIBS)

laser: laser.o libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(MYSQLLIBS) $(H5LIBS)

logain: logain.o libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(MYSQLLIBS) $(H5LIBS)

flatfield: flatfield.o libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(MYSQLLIBS) $(H5LIBS)

stage1: stage1.o libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(MYSQLLIBS) \
	$(H5LIBS) -ldl

stage2: stage2.o libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< \
	$(LIBS) $(MYSQLLIBS) $(H5LIBS) -ldl

stage3: stage3.o libVSStage3.a libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< \
	-lVSStage3 $(LIBS) $(MYSQLLIBS) $(H5LIBS) -ldl

nspace: nspace.o ../VSUtility/libVSUtility.a \
		 ../VSNSpace/libVSNSpace.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS) 

collect_muon: collect_muon.o libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

dump_muon: dump_muon.o libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(H5LIBS)

make_sp_tables: make_sp_tables.o libVSStage3.a libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< -lVSStage3 \
	$(LIBS) $(MYSQLLIBS) $(H5LIBS) -ldl

make_egy_tables: make_egy_tables.o libVSStage3.a libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< -lVSStage3 \
	$(LIBS) $(MYSQLLIBS) $(H5LIBS) -ldl

make_effarea_library: make_effarea_library.o libVSStage3.a \
	libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< -lVSStage3 \
	$(LIBS) $(MYSQLLIBS) $(H5LIBS) -ldl

extract_diagnostics: extract_diagnostics.o
	$(CXX) $(LDFLAGS) -o $@ $< \
	$(LIBS) $(MYSQLLIBS) $(H5LIBS) -ldl

combine_sp_tables: combine_sp_tables.o libVSDataReduction.a 
	$(CXX) $(LDFLAGS) -o $@ $< \
	$(LIBS) $(MYSQLLIBS) $(H5LIBS) -ldl

s2opt: s2opt.o
	$(CXX) $(LDFLAGS) -o $@ $< \
	$(LIBS) $(MYSQLLIBS) $(H5LIBS) -ldl

.PHONY: clean

clean:
	$(RM) \
	$(TARGETS) $(OBJECTS) *~ test

clean_stage3:
	$(RM) $(LIBOBJECTS_STAGE3) *~ test 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
