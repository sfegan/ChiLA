include ../Makefile.common
CHILADIR = ..

OBJECTS = Vec3D.o Vec4D.o Particle.o Constants.o \
	  RandomNumbers_TOS.o RandomNumbers_TNG.o xytohex.o xytosquare.o

#CXXFLAGS = -xN -O3 -Wall
#LDFLAGS = -xN -O3 -Wall
#CXX = /opt/intel_cc_80/bin/icc
#CXXFLAGS = -g -O3 -Wall
#LDFLAGS = -g -O3 -Wall 
#CXX = g++

#%CXXFLAGS = -O3
#%LDFLAGS = -O3
#%CXX = icc

all: libPhysics.a

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

libPhysics.a: $(OBJECTS)
	$(AR) r $@ $^

.PHONY: clean

clean:
	$(RM) libPhysics.a $(OBJECTS) *~ test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
