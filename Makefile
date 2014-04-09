# fortran 90/95 compiler
FC=gfortran

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

# archiver and flags
AR=ar
AFLAGS=rcv

# command for preprocessing (should be empty);
# used for parallel profiling with scalasca (see below)
PREP=

# compiler dependent variables
FCFLAGS_ALL= -x f95-cpp-input -fdefault-real-8  -DFORTRAN_STREAMS $(INCDIRS)
FCFLAGS_OPT= -O3 -finline-functions
FCFLAGS_DBG= -g -O2 -eC
FCFLAGS_PROF=-pg
FCFLAGS_MPI= -DPARALLEL
LDFLAGS_ALL= 
LDFLAGS_OPT=
LDFLAGS_DBG= -g -eC
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

export FC FCFLAGS LDFLAGS AR AFLAGS PREP
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
	rm -f $(OBJECTS) $(TARGET) *.mod *.bak *.dat *.bin *.nc *.log *~
	rm -rf epik*
	rm -rf autom4te.cache

.SUFFIXES:

.PHONY: all debug parallel prof parprof subdirs $(SUBDIRS) clean distclean

