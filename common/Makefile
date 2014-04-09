# some files
TARGET=libcommon.a
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

boundary_common.o : common_types.o
constants_common.o : common_types.o
geometry_common.o : common_types.o
logio_common.o : common_types.o
output_common.o : common_types.o
reconstruction_common.o : common_types.o
sources_common.o : common_types.o
timedisc_common.o : common_types.o
physics_common.o : common_types.o constants_common.o sources_common.o
fluxes_common.o : common_types.o reconstruction_common.o
mesh_common.o : boundary_common.o geometry_common.o
