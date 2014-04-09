# some files
TARGET=libfluxes.a
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

reconstruction_generic.o : reconstruction_constant.o reconstruction_linear.o
fluxes_midpoint.o : reconstruction_generic.o
fluxes_trapezoidal.o : fluxes_midpoint.o reconstruction_generic.o
fluxes_generic.o : fluxes_midpoint.o fluxes_trapezoidal.o