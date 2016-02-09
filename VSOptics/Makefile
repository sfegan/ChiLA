include ../Makefile.common
CHILADIR = ..

LIBOBJECTS = VSOArrayParameters.o \
	     VSOMirror.o VSOPixel.o VSOTelescope.o VSOTelescopeArray.o \
             VSORayTracer.o

BINOBJECTS = test_raytrace.o

OBJECTS = $(LIBOBJECTS) $(BINOBJECTS)

TARGETS = libVSOptics.a $(BINOBJECTS:.o=)

LIBS = -lVSOptics -lVSUtility -lPhysics -lVSCommon $(MYSQLLIBS) -lpthread 

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

libVSOptics.a: $(LIBOBJECTS)
	$(AR) r $@ $^

test_raytrace: test_raytrace.o libVSOptics.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

VSOArrayParameters.o: VSOArrayParameters.cpp VSOArrayParameters.hpp

VSOTelescopeArray.o: VSOArrayParameters.hpp

VSOArrayParameters.cpp: ARRAY_INI_TEMPLATE code_gen_array_ini_cpp.pl
	perl code_gen_array_ini_cpp.pl $< > $@

VSOArrayParameters.hpp: ARRAY_INI_TEMPLATE code_gen_array_ini_hpp.pl
	perl code_gen_array_ini_hpp.pl $< > $@

.PHONY: clean

clean:
	$(RM) \
	VSOArrayParameters.cpp VSOArrayParameters.hpp \
	$(TARGETS) $(OBJECTS) *~ test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
