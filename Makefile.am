#                                               -*- Automake -*-
# Process this file with automake --add-missing to produce a Makefile.
#############################################################################
#                                                                           #
# fosite - 2D hydrodynamical simulation program                             #
# configuration file: Makefile.am                                           #
#                                                                           #
# Copyright (C) 2012-2013                                                   #
# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
#                                                                           #
# This program is free software; you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 2 of the License, or (at     #
# your option) any later version.                                           #
#                                                                           #
# This program is distributed in the hope that it will be useful, but       #
# WITHOUT ANY WARRANTY; without even the implied warranty of                #
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, GOOD TITLE or        #
# NON INFRINGEMENT.  See the GNU General Public License for more            #
# details.                                                                  #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with this program; if not, write to the Free Software               #
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
#                                                                           #
#############################################################################

include common.am
SUBDIRS=numtools common mesh physics boundary fluxes sources io timedisc . tests doc

SOURCES=fosite.f90
OBJECTS=$(SOURCES:.f90=.o)
SUBOBJS=$(foreach dir,$(SRCDIRS), $(wildcard $(top_builddir)/$(dir)/*.o))

noinst_LIBRARIES = libfosite.a
libfosite_a_SOURCES = $(SOURCES)
libfosite_a_LIBADD = $(SUBOBJS)

BINDIRS = @BINDIRS@ @EXAMPLES@ $EXAMPLES
include binaries.am

EXTRA_DIST = $(BINARIES_SRC) bootstrap.sh
