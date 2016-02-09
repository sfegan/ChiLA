include ../Makefile.common
CHILADIR = ..

LIBOBJECTS = VSNSpace.o VSNSpaceOctaveH5IO.o VSNSpaceUtility.o \
		VSNSpaceDataSource.o VSNSpaceOptimizer.o VSGenSpace.o

OBJECTS = $(LIBOBJECTS)

TARGETS = libVSNSpace.a

LIBS =  -lPhysics -lVSUtility -lVSCommon $(MYSQLLIBS)

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

libVSNSpace.a: $(LIBOBJECTS)
	$(AR) r $@ $^

.PHONY: clean

clean:
	$(RM) \
	$(TARGETS) $(OBJECTS) *~ test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

%.pdf: %.fig
	fig2dev -Lpdf $< $@
