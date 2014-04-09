# fortran 90/95 compiler
FC=ifort

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

# compiler dependent variables
FCFLAGS_ALL= -cpp -r8 $(INCDIRS) 
FCFLAGS_OPT=-O3 -ip -mcpu=pentium4
FCFLAGS_DBG=-g -check all
FCFLAGS_PROF=-p
FCFLAGS_MPI= -I/astro/home/tillense/lib/mpich2/include -DPARALLEL
LDFLAGS_ALL=
LDFLAGS_OPT=
LDFLAGS_DBG=-g
LDFLAGS_PROF=-p
LDFLAGS_MPI= -L/astro/home/tillense/lib/mpich2/lib -L/astro/home/tillense/lib/pvfs2/lib -lmpich   -lpthread -lrt -lcrypto -lpvfs2

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

export FC FCFLAGS LDFLAGS AR AFLAGS
# variable definitions end here


# compile rules
all : $(TARGET)

debug parallel prof : all

%.o : %.f90
	$(FC) $(FCFLAGS) -c $<

subdirs : $(SUBDIRS)

$(TARGET) : subdirs $(OBJECTS)
	$(FC) $(ALLOBJS) -o $(TARGET) $(LDFLAGS)

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
	rm -f $(OBJECTS) $(TARGET) *.mod *.bak *.dat *.dx *~

.SUFFIXES:

.PHONY: all debug parallel subdirs $(SUBDIRS) clean distclean

