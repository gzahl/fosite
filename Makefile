# some files
TARGET=fosite
SOURCES=$(wildcard *.f90)
OBJECTS=$(SOURCES:.f90=.o)

#some directories
BASEDIR=$(shell pwd)
SUBDIRS=numtools common physics boundary fluxes mesh sources timedisc io
INCDIRS=$(foreach dir,$(SUBDIRS),-I$(BASEDIR)/$(dir))
#LIBDIRS=$(foreach dir,$(SUBDIRS),-L$(BASEDIR)/$(dir))
ALLSRCS=$(SOURCES) $(foreach dir,$(SUBDIRS), $(wildcard $(dir)/*.f90))
ALLOBJS=$(ALLSRCS:.f90=.o)

# compiler and flags
FC = ifort
FFLAGS = $(INCDIRS) -W1 -O
#LDFLAGS= $(LIBDIRS) -lphysics -lmesh -lfluxes -lboundary -lsources \
#	-lio -lcommon -lnt
LDFLAGS= $(LIBDIRS) -lsvml -i-static -static

# archiver and flags
AR = ar
AFLAGS = rcv

# target dependent variables
gfortran : FC = gfortran
gfortran : FFLAGS = $(INCDIRS) -Wall -O3 -pg
gfortran : FFLAGS_NOVEC = $(FFLAGS)
gfortran : LDFLAGS= $(LIBDIRS) -pg

g95 : FC = g95
g95 : FFLAGS = $(INCDIRS) -Wall -g -ftrace=full
g95 : FFLAGS_NOVEC = $(FFLAGS)
g95 : LDFLAGS= $(LIBDIRS)

gprof : FC = ifort
gprof : FFLAGS = $(INCDIRS) -W1 -O3 -xP -ipo -qp 
gprof : FFLAGS_NOVEC = $(FFLAGS)
gprof : LDFLAGS= $(LIBDIRS) -qp -lsvml

debug : FFLAGS = $(INCDIRS) -g -warn all -CB
	FFLAGS_NOVEC = $(FFLAGS)

fast : FFLAGS = $(INCDIRS) -W1 -ipo -O3 -xP
#fast :	FFLAGS = $(INCDIRS) -W1 -ipo -O3 -xB
fast : FFLAGS_NOVEC = $(INCDIRS) -W1 -O3

double : FFLAGS = $(INCDIRS) -W1 -O3 -r8
double : FFLAGS_NOVEC = $(FFLAGS)

doublefast : FFLAGS = $(INCDIRS) -W1 -ipo -xP -O3 -r8
#doublefast : FFLAGS = $(INCDIRS) -W1 -ipo -xB -O3 -r8 
doublefast : FFLAGS_NOVEC = $(INCDIRS) -W1 -O3 -r8

parallel : FFLAGS = $(INCDIRS) -O -W0 -r8 -msse3
parallel : FFLAGS_NOVEC = $(INCDIRS) -W1 -O3 -r8

profgen : FFLAGS = $(INCDIRS) -O -prof-gen -prof-dir=$(BASEDIR)/profdata
profuse : FFLAGS = $(INCDIRS) -prof-use -ipo -prof-dir=$(BASEDIR)/profdata

export FC FFLAGS FFLAGS_NOVEC AR AFLAGS


.SUFFIXES:

.PHONY: all debug parallel subdirs $(SUBDIRS) clean distclean

all : $(TARGET)

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

# compile without vectorization
init.o : init.f90
	$(FC) $(FFLAGS_NOVEC) -c $<

main.o : main.f90
	$(FC) $(FFLAGS_NOVEC) -c $<

debug double doublefast fast gfortran g95 gprof parallel profgen profuse : all

subdirs : $(SUBDIRS)

$(TARGET) : subdirs $(OBJECTS)
	$(FC) $(ALLOBJS) -o $(TARGET) $(LDFLAGS)

$(SUBDIRS) :
	$(MAKE) -C $@

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

