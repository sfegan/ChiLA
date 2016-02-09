CHILADIR = ..
include ../Makefile.common

ifndef QTDIR
QTDIR = /usr/local/qt4
endif

QLIBOBJECTS = VQMultiCamera.o

LIBOBJECTS = VSimpleCameraLayout.o VQMCChannelRenderer.o VQColormap.o\
	$(QLIBOBJECTS) $(addprefix moc_,$(QLIBOBJECTS))

BINOBJECTS = drawcamera.o

QOBJECTS = $(QLIBOBJECTS)

OBJECTS = $(LIBOBJECTS) $(BINOBJECTS)

TARGETS = libVQMultiCamera.a $(BINOBJECTS:.o=)

LIBS = -L$(QTDIR)/lib -lVQMultiCamera -lPhysics -lVSUtility -lVSCommon \
	$(MYSQLLIBS) -lQtGui

INCLUDE_DIRS += -I$(QTDIR)/include/Qt -I$(QTDIR)/include

ifndef MOC
MOC = $(QTDIR)/bin/moc
endif	

all: $(TARGETS)

static: LDFLAGS += -static 
static: MYSQLLIBS = $(MYSQLLIBS_STATIC)
static: all

drawcamera: drawcamera.o libVQMultiCamera.a	
	$(CXX) $(subst -static,,$(LDFLAGS)) -o $@ $< $(LIBS)

libVQMultiCamera.a: $(LIBOBJECTS)
	$(AR) r $@ $^

.PHONY: clean

$(addprefix moc_,$(QOBJECTS:.o=.cpp)): moc_%.cpp: %.hpp
	$(MOC) -o $@ $<

clean:
	$(RM) \
	$(TARGETS) $(OBJECTS) $(addprefix moc_,$(QOBJECTS:.o=.cpp)) \
	*~ test test_ps

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<
