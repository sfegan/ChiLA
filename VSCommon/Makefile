include ../Makefile.common
CHILADIR = ..

LIBOBJECTS = VSAssert.o

BINOBJECTS = 

OBJECTS = $(LIBOBJECTS) $(BINOBJECTS)

TARGETS = libVSCommon.a $(BINOBJECTS:.o=)

LIBS = 

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

libVSCommon.a: $(LIBOBJECTS)
	$(AR) r $@ $^

.PHONY: clean

clean:
	$(RM) \
	$(TARGETS) $(OBJECTS) *~ test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
