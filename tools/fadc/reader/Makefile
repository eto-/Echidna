# Makefile for BOREXINO FADC ROOT-files production tool

PROG	= fadc_reader
SOURCES = *.cc
OBJECTS = $(patsubst %.cc, %.o, $(wildcard *.cc))
HEADERS = *.h

ECHIDNA_CYCLE	= 16
ECHIDNA_PATH	= $(HOME)/Echidna_c$(ECHIDNA_CYCLE)
LIBECHIDNA		= $(ECHIDNA_PATH)/rootechidna.so

ROOT_FILES_PATH	= /home/lukyanch/fadc/rootfiles/

ROOT_LIBS		= `root-config --glibs` -lSpectrum
ROOT_FLAGS  	= `root-config --cflags`

INCLUDES		= -I$(ECHIDNA_PATH) $(ROOT_FLAGS)
LIBS			= $(ROOT_LIBS) $(LIBECHIDNA) 

# ugly hack for bxmaster additional libs
ifeq ($(shell hostname),bxmaster)
	LIBS	+= ./libz.a -lssp
else
	LIBS	+= -lz
endif

$(PROG): $(OBJECTS) $(LIBECHIDNA)
	$(CXX) $(OBJECTS) -o $(PROG) $(LIBS) -Wl,-rpath=$(ECHIDNA_PATH)

%.o: %.cc $(HEADERS)
	$(CXX) -O $(INCLUDES) $(patsubst %.o, %.cc, $@) -o $@ -c \
	-DCYCLE=$(ECHIDNA_CYCLE) -DROOT_FILES_PATH=\"$(ROOT_FILES_PATH)\"

clean:
	rm -f *.o $(PROG)
