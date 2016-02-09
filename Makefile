ANLTARGETS = Physics SEphem VSUtility VSOptics VSSimDB VSAnalysis \
	VSNSpace VSDataReduction

SIMTARGETS = VSShower VSBin VSElectronics

TARGETS = $(ANLTARGETS) $(SIMTARGETS) VSRootAnalysis VSMiscAnalysis \
	VSMatlabAnalysis VQMultiCamera

CLEANTARGETS = $(TARGETS) VSCommon

analysis: $(ANLTARGETS)

sim: analysis $(SIMTARGETS)

root: analysis VSRootAnalysis

matlab: analysis VSMatlabAnalysis

misc: analysis VSMiscAnalysis

camera: VQMultiCamera

all: analysis sim root matlab misc camera

static: STATICFLAG=static
static: analysis sim

analysis-static: STATICFLAG=static
analysis-static: analysis

LIBDEPS = Physics/libPhysics.a VSCommon/libVSCommon.a VSOptics/libVSOptics.a \
	 VSShower/libVSShower.a VSSimDB/libVSSimDB.a VSUtility/libVSUtility.a \
	 VSAnalysis/libVSAnalysis.a VSNSpace/libVSNSpace.a

VSBin: $(LIBDEPS)

Physics/libPhysics.a:		Physics
VSCommon/libVSCommon.a:		VSCommon
VSOptics/libVSOptics.a:		VSOptics
VSShower/libVSShower.a: 	VSShower
VSSimDB/libVSSimDB.a: 		VSSimDB
VSUtility/libVSUtility.a: 	VSUtility
VSAnalysis/libVSAnalysis.a: 	VSAnalysis
VSNSpace/libVSNSpace.a: 	VSNSpace
VSDataReduction/libSimpleVBF.a: VSDataReduction
VQMultiCamera/libVQMultiCamera.a:	VQMultiCamera

$(TARGETS): 		VSCommon

VSOptics:		VSUtility Physics
VSShower:		VSOptics
VSSimDB:		VSOptics
VSElectronics:          VSUtility VSOptics VSSimDB Physics
VSDataReduction:	VSUtility VSNSpace VSAnalysis Physics SEphem
VSMiscAnalysis:		VSDataReduction
VSRootAnalysis:		VSDataReduction
VQMultiCamera:		VSUtility Physics

clean: $(addsuffix -clean,$(CLEANTARGETS))
	$(RM) *~ test

docs:
	doxygen Docs/doxy.conf

.PHONY: $(TARGETS) $(addsuffix -clean,$(CLEANTARGETS)) VSCommon	

$(TARGETS) VSCommon: 
	$(MAKE) -C $@ $(STATICFLAG)

$(addsuffix -clean,$(CLEANTARGETS)):
	$(MAKE) -C $(@:-clean=) clean

