include ../Makefile.common

LIBOBJECTS = 

BINOBJECTS = simple_reconstruct.o simple_reconstruct_new.o

OBJECTS = $(LIBOBJECTS) $(BINOBJECTS)

TARGETS = $(BINOBJECTS:.o=)

LIBDEPS = ../Physics/libPhysics.a ../VSOptics/libVSOptics.a \
	../VSShower/libVSShower.a ../VSSimDB/libVSSimDB.a \
	../VSUtility/libVSUtility.a ../VSAnalysis/libVSAnalysis.a

all: $(TARGETS)

#CXXFLAGS += -DNO_OUTPUT

simple_reconstruct: simple_reconstruct.o $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $< $(VSLIBS) $(MYSQLLIBS)

simple_reconstruct_new: simple_reconstruct_new.o $(LIBDEPS)
	$(CXX) $(LDFLAGS) -o $@ $< $(VSLIBS) $(MYSQLLIBS)

.PHONY: clean

clean:
	$(RM) \
	$(TARGETS) $(OBJECTS) $(DATA_FILES) \
	*~ test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

