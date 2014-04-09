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
TARGET=fosite
SOURCES=$(wildcard *.f90)
OBJECTS=$(SOURCES:.f90=.o)

# some directories
BASEDIR=$(shell pwd)
SUBDIRS=numtools common physics boundary fluxes mesh sources timedisc io
INCDIRS=$(foreach dir,$(SUBDIRS),-I$(BASEDIR)/$(dir))
LIBDIRS=$(foreach dir,$(SUBDIRS),-L$(BASEDIR)/$(dir))
ALLSRCS=$(SOURCES) $(foreach dir,$(SUBDIRS), $(wildcard $(dir)/*.f90))
ALLOBJS=$(ALLSRCS:.f90=.o)

# compiler dependent variables
FCFLAGS_ALL= -cpp    -DFORTRAN_STREAMS   $(INCDIRS)
FCFLAGS_OPT= -O2 
FCFLAGS_DBG=-fcheck=all
FCFLAGS_PROF=-pg
FCFLAGS_MPI= 
LDFLAGS_ALL=     
LDFLAGS_OPT= -O2
LDFLAGS_DBG=-fcheck=all
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

export FC FCFLAGS LDFLAGS AR AFLAGS AWK SED PREP
# variable definitions end here


# compile rules
all : $(TARGET)

debug parallel prof parprof : all

%.o : %.f90
	$(PREP) $(FC) $(FCFLAGS) -c $<

subdirs : $(SUBDIRS)

$(TARGET) : subdirs $(OBJECTS)
	$(PREP) $(FC) $(ALLOBJS) -o $(TARGET) $(LDFLAGS)

$(SUBDIRS) :
	$(MAKE) -C $@

common : numtools
physics : common
boundary : common physics
fluxes : common physics
mesh : common fluxes
sources : common physics boundary fluxes mesh
timedisc : common physics boundary fluxes mesh sources
io: common physics fluxes mesh

init.o : subdirs
main.o : subdirs init.o

clean :
	for dir in $(SUBDIRS); do \
	  $(MAKE) clean -C $$dir; \
	done
	rm -f $(OBJECTS) $(TARGET) *.mod

distclean :
	for dir in $(SUBDIRS); do \
	  $(MAKE) distclean -C $$dir; \
	done
	rm -f $(OBJECTS) $(TARGET) *.mod *.bak *.dat *.vts *.bin *.nc *.log *~
	rm -rf autom4te.cache epik*

.SUFFIXES:

.PHONY: all debug parallel prof parprof subdirs $(SUBDIRS) clean distclean

