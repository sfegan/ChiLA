include ../Makefile.common
CHILADIR = ..

LIBOBJECTS = VSLineTokenizer.o VSOptions.o \
             VSDatabase.o VSDBMySQLBase.o  VSDBMySQL3X.o  VSDBMySQL41.o \
	     VSDBFactory.o VSDBParameterTable.o VSTestTimer.o \
	     VSOctaveIO.o VSOctaveH5Reader.o VSOctaveH5Writer.o VATime.o \
	     VSFileUtility.o VSFileLock.o VSH5Library.o

BINOBJECTS = octaveio.o

OBJECTS = $(LIBOBJECTS) $(BINOBJECTS)

TARGETS = libVSUtility.a $(BINOBJECTS:.o=)

LIBS = -lVSUtility -lVSCommon $(MYSQLLIBS)

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

libVSUtility.a: $(LIBOBJECTS)
	$(AR) r $@ $^

mysql2latex: mysql2latex.o libVSUtility.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

octaveio: octaveio.o libVSUtility.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) -lhdf5 -lpthread -lz

test_ps: VSDBMySQL41.cpp libVSUtility.a $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -DTEST_MAIN -o $@ VSDBMySQL41.cpp \
		VSDBMySQL3X.o VSDBMySQLBase.o VSOptions.o VSDatabase.o $(LIBS)

vatime: vatime.o libVSUtility.a
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) -lVSUtility -lVSCommon

.PHONY: clean

clean:
	$(RM) \
	$(TARGETS) $(OBJECTS) *~ test test_ps mysql2latex.o mysql2latex \
				vatime.o vatime

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
