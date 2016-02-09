include ../Makefile.common
CHILADIR = ..

BERNLOHR_OBJECTS = $(BERNLOHR_DIR)/$(BERNLOHR_IACT_PREFIX)io_simtel.o \
        $(BERNLOHR_DIR)/$(BERNLOHR_BERN_PREFIX)fileopen.o \
        $(BERNLOHR_DIR)/$(BERNLOHR_BERN_PREFIX)straux.o \
        $(BERNLOHR_DIR)/$(BERNLOHR_IACT_PREFIX)eventio.o \
        $(BERNLOHR_DIR)/$(BERNLOHR_IACT_PREFIX)warning.o

LIBOBJECTS = VSCORSIKAEvent.o VSReconstruction.o \
	VSTargeting.o VSRayTracingVisitor.o VSSimDBVisitor.o 

BINOBJECTS = 

OBJECTS = $(LIBOBJECTS) $(BINOBJECTS)

TARGETS = libVSShower.a $(BINOBJECTS:.o=)

CXXFLAGS += -I$(BERNLOHR_DIR) -O3

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

libVSShower.a: $(LIBOBJECTS) $(BERNLOHR_OBJECTS)
	$(AR) r $@ $^

VSReconstruction_fast.o: VSReconstruction_fast.cpp
#	icc -I/usr/local/intel/cc/include/c++  -Wall -g -I../Physics -vec_report2 -xN -O3 -c $<
#	icc -Wall -g -I../Physics -vec_report2 -xP -O3 -c $<
	g++ -m32 -g -I../Physics -O3 -c $<

.PHONY: clean

clean:
	$(RM) \
	$(TARGETS) $(OBJECTS) $(CORSIKA_DATA_FILES) $(FLUKA_DATA_FILES) *~ test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
