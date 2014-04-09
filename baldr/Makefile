# some files
BASEDIR=$(shell pwd)/..
WORKDIR=$(shell pwd | $(SED) -e 's/.*\///')
TARGET=baldr_impl.so

.SUFFIXES: .f90

.PHONY: clean distclean

all : $(TARGET)

baldr_impl.f90 : baldr_impl.fpp
	$(shell gfortran -E -cpp $< | $(SED) 's/@@/\n/g;s/##//g' > $@)

#On Pluto set "--fcompiler=gnu95", because otherwise it defaults to g77 and fails to compile therefore
$(TARGET) : baldr_impl.f90
	LDFLAGS="$(LDFLAGS) -shared" f2py -DF2PY_REPORT_ON_ARRAY_COPY=1 --f90exec=$(shell which $(FC)) --f90flags="$(FCFLAGS) $(LIBDIRS) $(INCDIRS) -I./.. $(LDFLAGS) -shared" -c baldr_impl.f90 $(BASEDIR)/fosite.a -m baldr_impl $(LIBDIRS) -L$(BASEDIR) $(INCDIRS)

clean :
	rm -f baldr_impl.f90 baldr_impl.s baldr_impl.so

distclean : clean
	rm -f *.bak *~ *.pyc *.vts *.pvd


