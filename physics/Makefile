# some files
TARGET=libphysics.a
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

physics_euler3Drotsym.o : physics_euler2D.o
physics_euler3Drotamt.o : physics_euler2D.o physics_euler3Drotsym.o
physics_generic.o : physics_euler2D.o physics_euler3Drotsym.o \
	physics_euler3Drotamt.o constants_generic.o
constants_generic.o : constants_SI.o constants_geometrical.o
