# some files
TARGET=libio.a
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

fileio_binary.o : fileio_gnuplot.o 
fileio_hdf5.o : fileio_gnuplot.o 
fileio_generic.o : fileio_binary.o fileio_hdf5.o fileio_gnuplot.o 
