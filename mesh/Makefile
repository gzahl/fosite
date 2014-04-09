# some files
WORKDIR=$(shell pwd | $(SED) -e 's/.*\///')
TARGET=lib$(WORKDIR).a
DEPFILE=Makefile.dep
SOURCES=$(wildcard *.f90)
OBJECTS=$(SOURCES:.f90=.o)

.SUFFIXES: .f90

.PHONY: clean distclean

%.o : %.f90
	$(PREP) $(FC) $(FCFLAGS) -c $<

$(TARGET) : $(OBJECTS)
	$(AR) $(AFLAGS) $@ $(OBJECTS)

$(DEPFILE) : $(SOURCES)
	$(AWK) -f ../makedeps -v objfiles="$(OBJECTS)" $(SOURCES) > $@

-include $(DEPFILE)

clean :
	rm -f $(OBJECTS) $(TARGET) *.mod

distclean :
	rm -f $(OBJECTS) $(TARGET) *.mod *.bak *~


