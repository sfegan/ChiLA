include ../Makefile.common
CHILADIR = ..

LIBOBJECTS = VSSimDB.o VSSimDBCORSIKADatasets.o VSEventSupplier.o            \
	     VSEventSupplierDBTable.o VSEventSupplierTextFile.o              \
             VSEventSupplierLaser.o VSEventSupplierPedestal.o                \
	     VSHDFEventStore.o VSEventSupplierChainedTables.o                \
	     VSEventSupplierOneTableCommon.o VSEventSupplierOneTableDB.o     \
	     VSEventSupplierOneTableHDF.o VSEventSupplierFactory.o           \
	     VSEventSupplierManyTables.o VSChannelReMapping.o                \
	     VSChannelReMappingDBTable.o VSHDFEventStoreManyTables.o

OBJECTS = $(LIBOBJECTS)

TARGETS = libVSSimDB.a

LIBS =  -lPhysics -lVSUtility -lVSCommon $(MYSQLLIBS)

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

libVSSimDB.a: $(LIBOBJECTS)
	$(AR) r $@ $^

createarray: createarray.o libVSOptics.a 
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)


.PHONY: clean

clean:
	$(RM) \
	$(TARGETS) $(OBJECTS) *~ test

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
