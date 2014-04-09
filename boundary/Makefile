# some files
TARGET=libboundary.a
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

boundary_nogradients.o : ../common/boundary_common.o
boundary_periodic.o : boundary_nogradients.o
boundary_noslip.o : boundary_nogradients.o
boundary_fixed.o : boundary_nogradients.o
boundary_noh.o : boundary_nogradients.o boundary_fixed.o
boundary_custom.o : boundary_nogradients.o boundary_fixed.o
boundary_farfield.o : boundary_nogradients.o boundary_fixed.o
boundary_extrapolation.o : boundary_nogradients.o
boundary_reflecting.o : boundary_nogradients.o
boundary_axis.o : boundary_reflecting.o
boundary_folded.o : boundary_reflecting.o
boundary_generic.o : boundary_nogradients.o boundary_periodic.o \
	boundary_reflecting.o boundary_folded.o boundary_axis.o \
	boundary_fixed.o boundary_extrapolation.o boundary_noh.o\
        boundary_noslip.o boundary_custom.o boundary_farfield.o

