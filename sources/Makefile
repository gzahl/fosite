# some files
TARGET=libsources.a
SOURCES=$(wildcard *.f90)
OBJECTS=$(SOURCES:.f90=.o)

.SUFFIXES: .f90

.PHONY: clean distclean

%.o : %.f90
	$(FC) $(FCFLAGS) -c $<

$(TARGET) : $(OBJECTS)
	$(AR) $(AFLAGS) $@ $(OBJECTS)

clean :
	rm -f $(OBJECTS) $(TARGET) *.mod

distclean :
	rm -f $(OBJECTS) $(TARGET) *.mod *.bak *~

sources_pointmass.o : ../common/sources_common.o
sources_diskthomson.o : sources_pointmass.o
sources_viscosity.o : sources_pointmass.o
sources_cooling.o : sources_pointmass.o
sources_boundary.o :
sources_selfgravitation.o : sources_pointmass.o sources_boundary.o
sources_c_accel.o : sources_pointmass.o
sources_generic.o : sources_pointmass.o sources_diskthomson.o \
	sources_c_accel.o sources_cooling.o sources_viscosity.o \
        sources_selfgravitation.o
