# Directories where code is (add subfolders here)
DIRS := event modules framework interface particleid mach4

# Include path defined here
INCLUDEFLAGS := $(patsubst %, -I%, $(DIRS)) -I. $(LOCAL_BOREX_INCLUDE)
#CXXFLAGS := $(INCLUDEFLAGS) -Wall -Wno-non-virtual-dtor -m32
CXXFLAGS := $(INCLUDEFLAGS) -Wall -Wno-deprecated

# Add first libz
#LDFLAGS := -lz -lfftw3 -lm $(LOCAL_BOREX_LIBS) -m32
LDFLAGS := -lz -lfftw3 -lm $(LOCAL_BOREX_LIBS) 

OPTIMIZE_ECHIDNA := 1
#GCC_MARCH := "i686"
GCC_MARCH := "native"
NOXML := 1	

# Postgres libs
ifdef PGLIB
LDFLAGS += -L$(PGLIB) -lpq
else
LDFLAGS += -lpq
endif
ifdef PGINCLUDE
CXXFLAGS += -I$(PGINCLUDE)
else
CXXFLAGS += -I/usr/include/postgresql -I/usr/include/pgsql
endif

ifdef NOXML
CXXFLAGS += -DNOXML
else
CXXFLAGS += $(shell xmlrpc-c-config c++ libwww-client --cflags)
LDFLAGS += $(shell xmlrpc-c-config c++ libwww-client --libs)
endif

ifdef CMAP_HARD_DEBUG
CXXFLAGS += -DCMAP_ABORT
endif

# Root libs
ifeq ($(shell which root-config),)
$(error root is not configured properly. Define ROOTSYS variable!)
endif
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS += $(shell root-config --glibs) -lMinuit
ifneq ($(strip $(wildcard $(shell root-config --libdir)/libSpectrum.so)),)
LDFLAGS += -lSpectrum
endif

# Optimize
ifdef OPTIMIZE_ECHIDNA
ifndef GCC_MARCH
GCC_MARCH = $(shell uname -m)
endif
LDFLAGS += -O3 -Wall -march=$(GCC_MARCH) #-mfpmath=sse -msse2
CXXFLAGS += -O3 -Wall -fomit-frame-pointer -march=$(GCC_MARCH) #-mfpmath=sse -msse2
else
LDFLAGS += -g -Wall
CXXFLAGS += -g
endif

# Define CC files
LINKDEF := $(shell find $(DIRS) -maxdepth 1 -name '*LinkDef.hh')
LIBCC := $(shell find $(DIRS) -maxdepth 1 -name '*.cc') $(patsubst %LinkDef.hh,%_dict.cc,$(LINKDEF)) constants.cc event/Mach4Event.cc
ROOTCC := $(patsubst %LinkDef.hh,%.cc,$(LINKDEF)) $(patsubst %LinkDef.hh,%_dict.cc,$(LINKDEF))

# Define OBJS
LIBOBJ := $(patsubst %.cc, %.lo, $(LIBCC))
ROOTOBJ := $(patsubst %.cc, %_root.lo, $(ROOTCC))

# Get the Echidna version
LIBVERSION := $(shell grep "define ECHIDNA_VERSION" bx_rec_general.hh | cut -f 2 -d \")

# Default targets
all: status echidna rootechidna tools/fix_db
rootechidna: rootechidna.so.$(LIBVERSION)
rootmap: rootechidna.rootmap

# Real targets
status:
	@if [ -z "$(OPTIMIZE_ECHIDNA)" ]; then echo "Building debugging version of Echidna"; else echo "Building optimized version of Echidna for $(GCC_MARCH) architecture"; fi 

libechidna.so.$(LIBVERSION): $(LIBOBJ)
	@echo "  [LD] $@"
	@$(CCACHE) $(CXX) -shared -o $@ $^  -Wl,-soname -Wl,$@ 
	@ln -sf $@ `basename $@ .$(LIBVERSION)`

echidna: main.o libechidna.so.$(LIBVERSION)
	@echo "  [LD] $@"
	@$(CCACHE) $(CXX) -o $@ $< -L. -lechidna $(LDFLAGS) -Wl,--rpath -Wl,./

rootechidna.so.$(LIBVERSION): $(ROOTOBJ)
	@echo "  [LD] $@"
	@$(CCACHE) $(CXX) -shared -o $@ $^ -Wl,-soname -Wl,`basename $@` 
	@cd `dirname $@` && ln -sf `basename $@` `basename $@ .$(LIBVERSION)` 

rootechidna.rootmap: rootechidna.so.$(LIBVERSION)
	@echo "  [RLIBMAP] $@"
	@rm -f $@; for i in $(LINKDEF); do rlibmap -f -l rootechidna.so -c $$i >> $@ ; done
	@sed -i -e /vector.float/d -e /vector.int/d $@

online: echidna tools/omon_read
	@tools/install.sh

tools/omon_read: tools/omon_read.c
	@echo "  [CC] $@"
	@$(CC) $< -o $@ -Wall

event/Mach4Event.cc event/BxEvent.cc: event/Mach4Event.hh

event/Mach4Event.hh: mach4/Mach4Event.hh-scheme mach4/Mach4Event.cc-scheme mach4/M4_variables.txt mach4/build_m4event.sh
	@echo "  [SH] $@"
	@mach4/build_m4event.sh

# Internal build pattern rules
%.o: %.cc
	@echo "  [CXX] $@"
	@$(CCACHE) $(CXX) -c $< -o $@ $(CPPFLAGS) $(CXXFLAGS)

%_dict.lo: %_dict.cc
	@echo "  [CXX] $@"
	@$(CCACHE) $(CXX) -c $< -o $@ $(CPPFLAGS) $(CXXFLAGS) -fPIC

%.lo: %.cc
	@echo "  [CXX] $@"
	@$(CCACHE) $(CXX) -c $< -o $@ $(CPPFLAGS) $(CXXFLAGS) -fPIC

%_dict_root.lo: %_dict.cc
	@echo "  [CXX] $@"
	@$(CCACHE) $(CXX) -c $< -o $@ $(CPPFLAGS) $(CXXFLAGS) -fPIC -D_ECHIDNA_ROOTLIB_

%_root.lo: %.cc
	@echo "  [CXX] $@"
	@$(CCACHE) $(CXX) -c $< -o $@ $(CPPFLAGS) $(CXXFLAGS) -fPIC -D_ECHIDNA_ROOTLIB_

%_dict.cc: %.hh %LinkDef.hh
	@echo "  [ROOTCINT] $@"
	@rootcint -f $@ -c $(INCLUDEFLAGS) $^

# Ask to make to preserv dictionary intermediate files
#.PRECIOUS: %_dict.cc

# Include file specific dependencies
-include .deps

tools/fix_db: tools/fix_db.cc
	$(MAKE) -C tools fix_db

# Make dependencies table
deps: .deps
.deps: $(filter _dict.cc, $(LIBCC)) $(filter _dict.cc, $(ROOTCC)) $(shell find . $(DIRS) -maxdepth 1 -name '*.hh')
	@touch .deps
	@makedepend $(filter-out -pthread, $(CXXFLAGS)) -o.lo -f .deps $(LIBCC) 2>/dev/null
	@makedepend $(filter-out -pthread, $(CXXFLAGS)) -o_root.lo -f .deps -a $(ROOTCC) 2>/dev/null
	@makedepend $(filter-out -pthread, $(CXXFLAGS)) -f .deps -a main.cc 2>/dev/null
	@for i in $(LINKDEF); do\
	  j=`grep extra_include $$i| grep ^#pragma | cut -f2 -d\"`;\
	  if [ "$$j" ]; then\
	    echo >> .deps;\
	    echo "$$i: $$j" | sed s,\.\/,, | xargs >> .deps;\
	  fi;\
	done

# Some clean targets
clean:
	find . -name \*.lo -exec rm {} \;
	find . -name \*_dict\* -exec rm {} \;
	rm -f *.log main.o 

distclean: clean
	rm -f echidna libechidna* rootechidna* event/Mach4Event.hh event/Mach4Event.cc .deps .deps.bak

ctags: clean
	ctags -R --c++-kinds=+p --fields=+iaS --extra=+q
