# some files
TARGET=libmesh.a
SOURCES=$(wildcard *.f90)
OBJECTS=$(SOURCES:.f90=.o)

.SUFFIXES: .f90

.PHONY: clean distclean

%.o : %.f90
	$(PREP) $(FC) $(FCFLAGS) -c $<

$(TARGET) : $(OBJECTS)
	$(AR) $(AFLAGS) $@ $(OBJECTS)

clean :
	rm -f $(OBJECTS) $(TARGET) *.mod

distclean :
	rm -f $(OBJECTS) $(TARGET) *.mod *.bak *~

geometry_cartesian.o : ../common/geometry_common.o
geometry_cylindrical.o : geometry_cartesian.o
geometry_polar.o : geometry_cartesian.o
geometry_logpolar.o : geometry_oblatespheroidal.o
geometry_spherical.o : geometry_cartesian.o
geometry_oblatespheroidal.o : geometry_cartesian.o
geometry_generic.o : geometry_cartesian.o geometry_cylindrical.o \
	geometry_polar.o geometry_logpolar.o geometry_spherical.o \
	geometry_oblatespheroidal.o
mesh_midpoint.o : ../common/mesh_common.o geometry_generic.o
mesh_trapezoidal.o : mesh_midpoint.o geometry_generic.o
mesh_generic.o : mesh_midpoint.o mesh_trapezoidal.o
