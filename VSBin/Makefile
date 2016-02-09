include ../Makefile.common
CHILADIR = ..

LIBOBJECTS = 

BINOBJECTS = dump_corsika.o dump_optics.o manage_simdb.o make_optics.o \
	make_simdb.o make_steering.o raytrace_corsika_to_simdb.o \
	raytrace_psf.o make_density_hist.o

OBJECTS = $(LIBOBJECTS) $(BINOBJECTS)

TARGETS = $(BINOBJECTS:.o=)

CORSIKA_DATA_FILES = cors_nucnuccs.h cors_venusdat.h \
        cors_egsdat5_005.h cors_egsdat5_015.h cors_egsdat5_025.h \
        cors_egsdat5_100.h cors_egsdat5_300.h cors_egsdat5_040.h \
        cors_atmabs_dat.h cors_quanteff_dat.h cors_mirreff_corsika_dat.h \
	cors_mirreff_falcone_dat.h cors_mirreff_falcone_squared_dat.h \
        cors_glaubtar_dat.h cors_nuclear_bin.h cors_sectnu.h cors_qgsdat.h

FLUKA_DATA_FILES = fluka_nuclear_bin.h fluka_sigmapi_bin.h fluka_elasct_bin.h

BERNLOHR_DATA_FILES = bern_atmprof1_dat.h bern_atmprof2_dat.h \
	bern_atmprof3_dat.h bern_atmprof4_dat.h bern_atmprof5_dat.h \
	bern_atmprof6_dat.h bern_atmprof9_dat.h

DATA_FILES = $(CORSIKA_DATA_FILES) $(FLUKA_DATA_FILES) $(BERNLOHR_DATA_FILES)

LIBDEPS = ../Physics/libPhysics.a ../VSOptics/libVSOptics.a \
	../VSShower/libVSShower.a ../VSSimDB/libVSSimDB.a \
	../VSUtility/libVSUtility.a ../VSCommon/libVSCommon.a

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

dump_corsika: dump_corsika.o $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $< $(VSLIBS)

dump_optics: dump_optics.o $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $< $(VSLIBS) $(MYSQLLIBS)

manage_simdb: manage_simdb.o ../Physics/libPhysics.a ../VSSimDB/libVSSimDB.a \
        ../VSUtility/libVSUtility.a ../VSCommon/libVSCommon.a
	$(CXX) $(LDFLAGS) -o $@ $^ $(MYSQLLIBS) $(H5LIBS) -lpthread

make_optics: make_optics.o $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $< $(VSLIBS) $(MYSQLLIBS)

make_simdb: make_simdb.o $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(VSLIBS) $(MYSQLLIBS)

make_steering: make_steering.o $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(VSLIBS) $(MYSQLLIBS)

raytrace_corsika_to_simdb: raytrace_corsika_to_simdb.o $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(VSLIBS) $(MYSQLLIBS) $(H5LIBS)

raytrace_psf: raytrace_psf.o $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(VSLIBS) $(MYSQLLIBS)

make_simdb.o: db_config_template.h

make_steering.o: $(DATA_FILES)

make_density_hist: make_density_hist.o $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $^ -lVSNSpace $(VSLIBS) $(MYSQLLIBS) $(H5LIBS)

.PHONY: clean

clean:
	$(RM) \
	$(TARGETS) $(OBJECTS) $(DATA_FILES) \
	db_config_template.h \
	*~ test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

# Data files included in binary files

cors_nucnuccs.h: $(CORSIKA_DATA_DIR)/NUCNUCCS
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_NUCNUCCS = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_venusdat.h: $(CORSIKA_DATA_DIR)/VENUSDAT
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_VENUSDAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_egsdat5_005.h: $(CORSIKA_DATA_DIR)/EGSDAT5_.05
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_EGSDAT5_005 = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_egsdat5_015.h: $(CORSIKA_DATA_DIR)/EGSDAT5_.15
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_EGSDAT5_015 = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_egsdat5_025.h: $(CORSIKA_DATA_DIR)/EGSDAT5_.25
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_EGSDAT5_025 = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_egsdat5_040.h: $(CORSIKA_DATA_DIR)/EGSDAT5_.4
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_EGSDAT5_040 = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_egsdat5_100.h: $(CORSIKA_DATA_DIR)/EGSDAT5_1.
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_EGSDAT5_100 = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_egsdat5_300.h: $(CORSIKA_DATA_DIR)/EGSDAT5_3.
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_EGSDAT5_300 = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_atmabs_dat.h: $(CORSIKA_DATA_DIR)/atmabs.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_ATMABS_DAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_quanteff_dat.h: $(CORSIKA_DATA_DIR)/quanteff.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_QUANTEFF_DAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_mirreff_corsika_dat.h: $(CORSIKA_DATA_DIR)/mirreff.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_MIRREFF_CORSIKA_DAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_mirreff_falcone_dat.h: ../Data/falcone_mirror.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_MIRREFF_FALCONE_DAT = { \"mirreff.dat\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_mirreff_falcone_squared_dat.h: ../Data/falcone_double_mirror.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_MIRREFF_FALCONE_SQUARED_DAT = { \"mirreff.dat\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_glaubtar_dat.h: $(CORSIKA_DATA_DIR)/GLAUBTAR.DAT
	perl -MPOSIX -ne 'BEGIN{print "DataFile CORS_GLAUBTAR_DAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

cors_nuclear_bin.h: $(CORSIKA_DATA_DIR)/NUCLEAR.BIN
	perl -MPOSIX -MFileHandle -e 'print "DataFile CORS_NUCLEAR_BIN = { \"$(notdir $<)\", \""; $$fp=new FileHandle("<".$$ARGV[0]); while($$fp->read($$_,100)) { $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); }; print "\", ",$$ll, " };\n";' $< > $@

#cors_sectnu.h: $(CORSIKA_DATA_DIR)/sectnu-II-03
#	perl -MPOSIX -MFileHandle -e 'print "DataFile CORS_SECTNU = { \"$(notdir $<)\", \""; $$fp=new FileHandle("<".$$ARGV[0]); while($$fp->read($$_,100)) { $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); }; print "\", ",$$ll, " };\n";' $< > $@
#
#cors_qgsdat.h: $(CORSIKA_DATA_DIR)/qgsdat-II-03
#	perl -MPOSIX -MFileHandle -e 'print "DataFile CORS_QGSDAT = { \"$(notdir $<)\", \""; $$fp=new FileHandle("<".$$ARGV[0]); while($$fp->read($$_,100)) { $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); }; print "\", ",$$ll, " };\n";' $< > $@
#
cors_sectnu.h: $(CORSIKA_DATA_DIR)/SECTNU
	perl -MPOSIX -MFileHandle -e 'print "DataFile CORS_SECTNU = { \"$(notdir $<)\", \""; $$fp=new FileHandle("<".$$ARGV[0]); while($$fp->read($$_,100)) { $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); }; print "\", ",$$ll, " };\n";' $< > $@

cors_qgsdat.h: $(CORSIKA_DATA_DIR)/QGSDAT01
	perl -MPOSIX -MFileHandle -e 'print "DataFile CORS_QGSDAT = { \"$(notdir $<)\", \""; $$fp=new FileHandle("<".$$ARGV[0]); while($$fp->read($$_,100)) { $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); }; print "\", ",$$ll, " };\n";' $< > $@

fluka_nuclear_bin.h: $(FLUKA_DIR)/nuclear.bin
	perl -MPOSIX -MFileHandle -e 'print "DataFile FLUKA_NUCLEAR_BIN = { \"$(notdir $<)\", \""; $$fp=new FileHandle("<".$$ARGV[0]); while($$fp->read($$_,100)) { $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); }; print "\", ",$$ll, " };\n";' $< > $@

fluka_sigmapi_bin.h: $(FLUKA_DIR)/sigmapi.bin
	perl -MPOSIX -MFileHandle -e 'print "DataFile FLUKA_SIGMAPI_BIN = { \"$(notdir $<)\", \""; $$fp=new FileHandle("<".$$ARGV[0]); while($$fp->read($$_,100)) { $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); }; print "\", ",$$ll, " };\n";' $< > $@

fluka_elasct_bin.h: $(FLUKA_DIR)/elasct.bin
	perl -MPOSIX -MFileHandle -e 'print "DataFile FLUKA_ELASCT_BIN = { \"$(notdir $<)\", \""; $$fp=new FileHandle("<".$$ARGV[0]); while($$fp->read($$_,100)) { $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); }; print "\", ",$$ll, " };\n";' $< > $@

bern_atmprof1_dat.h: $(BERNLOHR_DIR)/atmprof1.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile BERN_ATMPROF1_DAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

bern_atmprof2_dat.h: $(BERNLOHR_DIR)/atmprof2.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile BERN_ATMPROF2_DAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

bern_atmprof3_dat.h: $(BERNLOHR_DIR)/atmprof3.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile BERN_ATMPROF3_DAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

bern_atmprof4_dat.h: $(BERNLOHR_DIR)/atmprof4.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile BERN_ATMPROF4_DAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

bern_atmprof5_dat.h: $(BERNLOHR_DIR)/atmprof5.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile BERN_ATMPROF5_DAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

bern_atmprof6_dat.h: $(BERNLOHR_DIR)/atmprof6.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile BERN_ATMPROF6_DAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

bern_atmprof9_dat.h: $(BERNLOHR_DIR)/atmprof9.dat
	perl -MPOSIX -ne 'BEGIN{print "DataFile BERN_ATMPROF9_DAT = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@

db_config_template.h: DB_CONFIG_TEMPLATE
	perl -MPOSIX -ne 'BEGIN{print "DataFile DB_CONFIG_TEMPLATE = { \"$(notdir $<)\", \"";}; $$ll+=length($$_); print map { ((isgraph(chr)||($$_==32))&&($$_!=92)&&($$_!=34)&&($$_!=63))?chr:sprintf("\\%03o", $$_) } unpack("C*",$$_); END{print "\", ",$$ll, " };\n"}' $< > $@
