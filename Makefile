# fortran 90/95 compiler
FC=gfortran

# awk and sed program
AWK=gawk
SED=/usr/bin/sed

# archiver and flags
AR=ar
AFLAGS=rcv

# command for preprocessing (should be empty);
# used for parallel profiling with scalasca (see below)
PREP=

# some files
EXAMPLES=$(basename $(wildcard examples/*.f90))
ifndef TARGET
	TARGET=$(EXAMPLES)
endif
SOURCES=fosite.f90
OBJECTS=$(SOURCES:.f90=.o)

# some directories
BASEDIR=$(shell pwd)
SUBDIRS=numtools common mesh physics boundary fluxes sources io timedisc
INCDIRS=$(foreach dir,$(SUBDIRS),-I$(BASEDIR)/$(dir))
LIBDIRS=$(foreach dir,$(SUBDIRS),-L$(BASEDIR)/$(dir))
ALLSRCS=$(SOURCES) $(foreach dir,$(SUBDIRS), $(wildcard $(dir)/*.f90))
ALLOBJS=$(ALLSRCS:.f90=.o)

# compiler dependent variables
FCFLAGS_ALL= -cpp -fdefault-real-8     -DFORTRAN_STREAMS    $(INCDIRS)
FCFLAGS_OPT= -O3 -finline-functions
FCFLAGS_DBG= -g -O2 -DDEBUG -fcheck=all
FCFLAGS_PROF=-pg
FCFLAGS_MPI= 
LDFLAGS_ALL=       
LDFLAGS_OPT= -O3
LDFLAGS_DBG= -g -fcheck=all
LDFLAGS_PROF=-pg
LDFLAGS_MPI= 

# default compiler flags for target "all"
FCFLAGS=$(FCFLAGS_ALL) $(FCFLAGS_OPT)
LDFLAGS=$(LDFLAGS_ALL) $(LDFLAGS_OPT)

# target dependent compiler flags 
debug : FCFLAGS=$(FCFLAGS_ALL) $(FCFLAGS_DBG)
debug : LDFLAGS=$(LDFLAGS_ALL) $(LDFLAGS_DBG)
parallel : FCFLAGS=$(FCFLAGS_ALL) $(FCFLAGS_OPT) $(FCFLAGS_MPI)
parallel : LDFLAGS=$(LDFLAGS_ALL) $(LDFLAGS_OPT) $(LDFLAGS_MPI)
prof : FCFLAGS=$(FCFLAGS_ALL) $(FCFLAGS_OPT) $(FCFLAGS_PROF)
prof : LDFLAGS=$(LDFLAGS_ALL) $(LDFLAGS_OPT) $(LDFLAGS_PROF)
parprof : FCFLAGS=$(FCFLAGS_ALL) $(FCFLAGS_OPT) $(FCFLAGS_MPI)
parprof : LDFLAGS=$(LDFLAGS_ALL) $(LDFLAGS_OPT) $(LDFLAGS_MPI)
parprof : PREP=
baldr : FCFLAGS=$(FCFLAGS_ALL) $(FCFLAGS_OPT) -fPIC
baldr : LDFLAGS=$(LDFLAGS_ALL) $(LDFLAGS_OPT)

export FC FCFLAGS LDFLAGS AR AFLAGS AWK SED PREP
# variable definitions end here

# compile rules
all : $(addsuffix .out,$(TARGET))

debug parallel prof parprof: all

%.o : %.f90
	$(PREP) $(FC) $(FCFLAGS) -o $@ -c $<

subdirs : $(SUBDIRS)

fosite.a : subdirs $(OBJECTS)
	$(AR) $(AFLAGS) $@ $(ALLOBJS)

%.out : fosite.a
	$(PREP) $(FC) $(FCFLAGS) $(@:.out=.f90) -o $@ fosite.a $(LDFLAGS)

$(SUBDIRS) :
	$(MAKE) -C $@

baldr : fosite.a
	$(MAKE) -C $@

common : numtools
mesh : common
physics : common mesh
boundary : common physics
fluxes : common physics mesh
sources : common physics boundary fluxes mesh
timedisc : common physics boundary fluxes mesh sources io
io: common physics fluxes mesh sources

fosite.o : subdirs 

clean :
	for dir in $(SUBDIRS); do \
	  $(MAKE) clean -C $$dir; \
	done
	$(MAKE) clean -C baldr
	rm -f $(OBJECTS) $(TARGET) fosite.a fosite.so *.mod *.o $(addsuffix .out,$(EXAMPLES)) $(addsuffix .o,$(EXAMPLES))

distclean : clean
	for dir in $(SUBDIRS); do \
	  $(MAKE) distclean -C $$dir; \
	done
	$(MAKE) distclean -C baldr
	rm -f $(OBJECTS) $(TARGET) *.mod *.bak *.dat *.vts *.pvd *.bin *.nc *.log *~
	rm -rf autom4te.cache epik*

.SUFFIXES:

.PHONY: all debug parallel prof parprof subdirs $(SUBDIRS) clean distclean

