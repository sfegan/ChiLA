include ../Makefile.common
CHILADIR = ..
ifndef MATLABDIR
MATLABDIR = /veritas/matlab/extern
endif

MEXCFLAGS = -I${MATLABDIR}/include -DMATLAB_MEX_FILE -DMX_COMPAT_32 -fPIC -fno-omit-frame-pointer -ansi -D_GNU_SOURCE -pthread -O -DNDEBUG -fexceptions

MEXLDFLAGS = -pthread -shared -Wl,--version-script,${MATLABDIR}/lib/${MATLABARCH}/mexFunction.map -Wl,-rpath-link,${MATLABDIR}/../bin/${MATLABARCH}

OPTFLAGS = 
CFLAGS = $(CHILADEFFLAGS) -I../VSUtility -I../VSCommon $(MEXCFLAGS) 
CXXFLAGS = $(CHILADEFFLAGS) -I../VSUtility -I../VSCommon $(MEXCFLAGS)
LDFLAGS = $(MEXLDFLAGS) -L$(HDF5DIR)/lib

ifdef MATLABHDF5DIR
CFLAGS += -I$(MATLABHDF5DIR)/include
CXXFLAGS += -I$(MATLABHDF5DIR)/include
else 
CFLAGS += -I$(HDF5DIR)/include
CXXFLAGS += -I$(HDF5DIR)/include
endif

LIBS = -lhdf5

all: a64

glx: loadh5.mexglx 
glx: MATLABARCH = glnx86

a64: loadh5.mexa64
a64: MATLABARCH = glnxa64

loadh5.mexglx: loadh5.o mexversion.o VSOctaveH5Reader.o VSOctaveH5Writer.o \
		VSOctaveIO.o VSFileUtility.o VSDataConverter.o
	${CXX} ${LDFLAGS} -o $@  $^ ${LIBS}
#	chcon -t texrel_shlib_t $@

loadh5.mexa64: loadh5.o mexversion.o VSOctaveH5Reader.o VSOctaveH5Writer.o \
		VSOctaveIO.o VSFileUtility.o VSDataConverter.o
	${CXX} ${LDFLAGS} -o $@  $^ ${LIBS}
#	chcon -t texrel_shlib_t $@

VSOctaveH5Reader.o: ../VSUtility/VSOctaveH5Reader.cpp
	${CXX} $(CXXFLAGS) -c $<

VSOctaveH5Writer.o: ../VSUtility/VSOctaveH5Writer.cpp
	${CXX} $(CXXFLAGS) -c $<

VSOctaveIO.o: ../VSUtility/VSOctaveIO.cpp
	${CXX} $(CXXFLAGS) -c $<

VSFileUtility.o: ../VSUtility/VSFileUtility.cpp
	${CXX} $(CXXFLAGS) -c $<

VSDataConverter.o: ../VSUtility/VSDataConverter.cpp
	${CXX} $(CXXFLAGS) -c $<

loadh5.o: loadh5.cpp
	${CXX} $(CXXFLAGS) -c $<

mexversion.o: ${MATLABDIR}/src/mexversion.c
	${CC} $(CFLAGS) -c $<

clean:
	${RM} loadh5.mexglx loadh5.mexa64 loadh5.o mexversion.o *~ \
		VSOctaveH5Reader.o VSOctaveH5Writer.o \
		VSOctaveIO.o VSFileUtility.o VSDataConverter.o
