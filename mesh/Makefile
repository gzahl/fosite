# some files
TARGET=libmesh.a
SOURCES=$(wildcard *.f90)
OBJECTS=$(SOURCES:.f90=.o)

.SUFFIXES: .f90

.PHONY: clean distclean

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

$(TARGET) : $(OBJECTS)
	$(AR) $(AFLAGS) $@ $(OBJECTS)

clean :
	rm -f $(OBJECTS) $(TARGET) *.mod

distclean :
	rm -f $(OBJECTS) $(TARGET) *.mod *.bak *~

geometry_generic.o : geometry_cartesian.o geometry_cylindrical.o \
	geometry_polar.o geometry_spherical.o \
	geometry_oblatespheroidal.o
mesh_generic.o : geometry_generic.o mesh_midpoint.o mesh_trapezoidal.o
