# fortran 90/95 compiler
FC=mpif90

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
FCFLAGS_ALL= -x f95-cpp-input -I/astro/home/tillense/lib/netcdf/include  -DHAVE_HDF5 -DHAVE_NETCDF $(INCDIRS)
FCFLAGS_OPT= -O3 
FCFLAGS_DBG= -g -O2
FCFLAGS_PROF=
FCFLAGS_MPI= -DPARALLEL
LDFLAGS_ALL= -L/astro/home/tillense/lib/netcdf/lib -L/astro/home/tillense/lib/hdf//lib -lnetcdf -lhdf5_hl -lz  -lhdf5
LDFLAGS_OPT=
LDFLAGS_DBG= -g
LDFLAGS_PROF=
LDFLAGS_MPI=   -lpthread -lrt

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
parprof : PREP=scalasca -instrument

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
sources : common physics fluxes mesh
timedisc : common physics boundary fluxes mesh sources
io: common physics

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
	rm -f $(OBJECTS) $(TARGET) *.mod *.bak *.dat *.bin *.nc *~
	rm -rf epik*

.SUFFIXES:

.PHONY: all debug parallel prof parprof subdirs $(SUBDIRS) clean distclean

