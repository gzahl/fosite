#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#############################################################################
#                                                                           #
# fosite - 2D hydrodynamical simulation program                             #
# configuration file: configure.ac                                          #
#                                                                           #
# Copyright (C) 2008-2012                                                   #
# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
# Marc Junker <maj@astrophysik.uni-kiel.de>                                 #
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

AC_PREREQ(2.62)
AC_INIT(fosite, 0.5.1, tillense@astrophysik.uni-kiel.de)

# ===========================================================================
# define command line arguments 
# ===========================================================================

# add some libraries
linkClibs=no
AC_ARG_WITH([mpi],
   AS_HELP_STRING([--with-mpi=[DIR]],[
      Enable MPI support and check for installation in directory [DIR]; 
      required to compile parallel version of fosite @<:@DEFAULT=no@:>@]),[dnl
   AS_CASE([$withval],[yes],[dnl
      linkClibs=yes],[no],[],[dnl
      linkClibs=yes
      MPIDIR=$withval
      with_mpi=yes
      # check if the given directory exists
      AS_IF([test -d "$MPIDIR"],[dnl
         mpi_bindirs="m4_combine([ ], [$MPIDIR], [/], [bin], [sbin])"
         mpi_libdirs="m4_combine([ ], [$MPIDIR], [/], [lib], [lib64])"
         mpi_moddirs="m4_combine([ ], [$MPIDIR], [/], [include],
                                      [mod], [modules], [lib], [lib64])"
         # set PATH variable for MPI programs
         psep=""
         for dir in $mpi_bindirs; do
            AS_IF([test -d "$dir"],[dnl
               AS_VAR_APPEND([MPIPATH],[$psep$dir])
               psep=$PATH_SEPARATOR
            ])
         done
      ],[dnl
         AC_MSG_WARN([MPI installation directory $MPIDIR doesn't exist])
      ])
   ])
],[dnl
  with_mpi=no
]) # with_mpi

AC_ARG_WITH([pvfs2],
   AS_HELP_STRING([--with-pvfs2=[DIR]],[
      Check for PVFS2 installation in directory [DIR];
      recommended for parallel I/O @<:@DEFAULT=no@:>@]),[dnl
   AS_CASE([$withval],[yes],[dnl
      linkClibs=yes],[no],[],[dnl
      linkClibs=yes
      PVFSDIR=$withval
      with_pvfs2=yes
      # check if the given directory exists
      AS_IF([test -d "$PVFSDIR"],[dnl
         pvfs2_libdirs="m4_combine([ ], [$PVFSDIR], [/], [lib], [lib64])"
      ],[dnl
         AC_MSG_WARN([PVFS2 installation directory $PVFSDIR doesn't exist])
      ])
   ])
],[dnl
  with_pvfs2=no
]) # with_pvfs2

AC_ARG_WITH([netcdf],
   AS_HELP_STRING([--with-netcdf=[DIR]],[
      Check for NetCDF installation in directory [DIR];
      required for NetCDF output @<:@DEFAULT=no@:>@]),[dnl
   AS_CASE([$withval],[yes],[dnl
      linkClibs=yes],[no],[],[dnl
      linkClibs=yes
      NETCDFDIR=$withval
      with_netcdf=yes
      # check if the given directory exists
      AS_IF([test -d "$NETCDFDIR"],[dnl
         netcdf_libdirs="m4_combine([ ], [$NETCDFDIR], [/], [lib], [lib64])"
         netcdf_moddirs="m4_combine([ ], [$NETCDFDIR], [/], [include],
                                         [mod], [modules], [lib], [lib64])"
      ],[dnl
         AC_MSG_WARN([NetCDF installation directory $NETCDFDIR doesn't exist])
      ])
   ])
],[dnl
  with_netcdf=no
]) # with_netcdf

AC_ARG_WITH([hdf5],
   AS_HELP_STRING([--with-hdf5=[DIR]],[
      Check for HDF5 installation in directory [DIR];
      required for parallel NetCDF output @<:@DEFAULT=no@:>@]),[dnl
   AS_CASE([$withval],[yes],[dnl
      linkClibs=yes],[no],[],[dnl
      linkClibs=yes
      HDF5DIR=$withval
      with_hdf5=yes
      # check if the given directory exists
      AS_IF([test -d "$HDF5DIR"],[dnl
         hdf5_libdirs="m4_combine([ ], [$HDF5DIR], [/], [lib], [lib64])"
      ],[dnl
         AC_MSG_WARN([HDF5 installation directory $HDF5DIR doesn't exist])
      ])
   ])
],[dnl
  with_hdf5=no
]) # with_hdf5

AC_ARG_WITH([npy],
	AS_HELP_STRING([--with-npy=[DIR]],[
		Check for libnpy (>=0.5) installation in directory [DIR];
		required for numpy npy binary file output @<:@DEFAULT=no@:>@]),[dnl
	AS_CASE([$withval],[yes],[dnl
		linkClibs=yes],[no],[],[dnl
		linkClibs=yes
		NPYDIR=$withval
		with_npy=yes
		# check if the given directory exists
		AS_IF([test -d "$NPYDIR"],[dnl
			   npy_libdirs="m4_combine([ ], [$NPYDIR], [/], [lib], [lib64])"
			   npy_moddirs="m4_combine([ ], [$NPYDIR], [/], [include])"
		],[dnl
		   AC_MSG_WARN([NPY installation directory $NPYDIR doesn't exist])
		])
	])
],[dnl
   with_npy=no
]) # with_npy

# enable additional features
AC_ARG_ENABLE([static],
	AC_HELP_STRING([--enable-static],
		[Create a statically linked executable file.
                 @<:@DEFAULT=no@:>@ ]),
		[arg_static=$enableval])

AC_ARG_ENABLE([double],
	AC_HELP_STRING([--enable-double],
		[turn on autodouble compiler flag if available. 
                 @<:@DEFAULT=yes@:>@ ]),[dnl
 		 arg_autodouble=$enableval],[dnl
		 arg_autodouble=yes
])
AC_ARG_ENABLE([cbind],
	AC_HELP_STRING([--enable-cbind],
		[turn on ISO_C_BINDING compiler flag if available. 
                 @<:@DEFAULT=no@:>@ ]),[dnl
 		 arg_cbind=$enableval],[dnl
		 arg_cbind=no
])


AC_ARG_ENABLE([fast],
   AC_HELP_STRING([--enable-fast=option],
      [turn on standard compiler optimizations;
      e.g. --enable-fast=O2 activates optimization level 2
      @<:@DEFAULT=yes@:>@ ]),[dnl
   arg_fast=$enableval],[dnl
   arg_fast=yes
])

AC_ARG_ENABLE([inlining],
   AC_HELP_STRING([--enable-inlining],
      [turn on inline expansion of functions;@<:@DEFAULT=no@:>@ ]),[dnl
   arg_inline=$enableval
])

AC_ARG_ENABLE([mpi-module],
   AC_HELP_STRING([--enable-mpi-module],
      [use Fortran 90 module file for MPI;@<:@DEFAULT=no@:>@ ]),[dnl
   arg_mpi_module=$enableval
],[dnl
   arg_mpi_module=no
])

AC_ARG_VAR([PREP],[preprocessor command for parallel profiling])

AC_ARG_VAR([MPIFC],[MPI fortran compiler command])


# ===========================================================================
# main configure script starts here!
# ===========================================================================

# determine the target host architecture we are building for
AC_CANONICAL_HOST

# standard archiver
AR="ar"

# host specific settings
AS_CASE([$host_cpu],
  [sx8],[dnl
     arg_fast_default="Cvopt"
	 AR="sxar"
     FCFLAGS="$FCFLAGS -$host_cpu"
     LDFLAGS="$LDFLAGS -$host_cpu"
     DEFS="$DEFS -DNECSX8"],
  [sx9],[dnl
     arg_fast_default="Cvopt"
	 AR="sxar"
     FCFLAGS="$FCFLAGS -$host_cpu"
     LDFLAGS="$LDFLAGS -$host_cpu"
     DEFS="$DEFS -DNECSX9"],
  [arg_fast_default="O2"]
)

# default programming language
AC_LANG(Fortran)


# ===========================================================================
# Fortran 90/95 compiler checks
# ===========================================================================
AS_VAR_IF([with_mpi],[no],[],[dnl
   ac_prog_fc_mpi_FC="sxmpif90 mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c"
   AS_VAR_SET_IF([FC],[MPIFC=$FC],[dnl
      # if _not_ set check for MPI compiler command in MPIPATH
      AS_VAR_SET_IF([MPIPATH],[dnl
         AC_PATH_PROGS([MPIFC],[$MPIFC $ac_prog_fc_mpi_FC],[],[$MPIPATH])
         AS_VAR_IF([MPIFC],[""],[dnl
            AC_MSG_WARN([no MPI Fortran compiler found in $MPIPATH])
         ],[dnl
            FC=$MPIFC
         ])
      ])
   ])
   AS_VAR_IF([MPIFC],[""],[dnl
      # if _not_ set check these MPI compiler commands (see below)
      MPIFC=$ac_prog_fc_mpi_FC
   ])
])
# check all other compilers if MPI is disabled or
# detection of MPI compiler failed
ac_prog_fc_FCFLAGS_save=$FCFLAGS
AC_PROG_FC([$MPIFC sxf90 f90 f95 gfortran ifort g95])

# debug flags and normal compile flags
AS_VAR_IF([ac_cv_prog_fc_g],[yes],[dnl
   FCFLAGS=$ac_prog_fc_FCFLAGS_save
   FCFLAGS_DBG="$FCFLAGS -g"
   LDFLAGS_DBG="$LDFLAGS -g"
   AS_VAR_IF([ac_cv_fc_compiler_gnu],[yes],[dnl
      FCFLAGS_DBG="$FCFLAGS_DBG -O2"])
])

# default file name extension for Fortran 90
AC_FC_SRCEXT(f90)

# ensure that the Fortran compiler accepts free-form
AC_FC_FREEFORM

# search for C compiler
AS_VAR_IF([with_mpi],[no],[],[dnl
   ac_prog_cc_mpi_CC="sxmpicc mpicc hcc mpxlc_r mpxlc mpcc cmpicc"
   AS_VAR_SET_IF([CC],[MPICC=$CC],[dnl
      # if _not_ set check for MPI compiler command in MPIPATH
      AS_VAR_SET_IF([MPIPATH],[dnl
         AC_PATH_PROGS([MPICC],[$MPICC $ac_prog_cc_mpi_CC],[],[$MPIPATH])
         AS_VAR_IF([MPICC],[""],[dnl
            AC_MSG_WARN([no MPI C compiler found in $MPIPATH])
         ],[dnl
            CC=$MPICC
         ])
      ])
   ])
   AS_VAR_IF([MPICC],[""],[dnl
      # if _not_ set check these MPI compiler commands (see below)
      MPICC=$ac_prog_cc_mpi_CC
   ])
])
AC_PROG_CC([$MPICC sxcc cc xlC xlc pgcc icc pathcc gcc])

# determine name mangling scheme of the Fortran compiler
# this makro needs a working C compiler, thus we check
# for ANSI C89 compliance
#AS_VAR_IF([ac_cv_prog_cc_c89],["no"],[dnl
#   AC_MSG_WARN([C compiler doesn't support ANSI C89])],[dnl
#   AC_FC_WRAPPERS
#])

# check preprocessor flag
AC_FC_CPPFLAG

# check Fortran 2003 streams
AC_FC_STREAMS([DEFS="$DEFS -DFORTRAN_STREAMS"])

# check whether the Fortran compiler supports ISO C bindings
AS_VAR_IF([arg_cbind],[yes],[
   AC_FC_ISO_C_BINDING([
      DEFS="$DEFS -DHAVE_ISO_C_BINDING"
      AS_VAR_IF([ac_cv_fc_cbind],[none],[],[dnl	
           FCFLAGS="$FCFLAGS $ac_cv_fc_cbind"])
   ])
])

# enable advanced debugging features if available
AC_FC_ADVANCED_DEBUG

# check for runtime profiling
AC_FC_PROFILING

# check for autodouble flag if requested
# ATTENTION: autdouble may break MPI and/or VTK , because some
# autodouble compiler flags do not only effect real numbers and
# promote 4 byte integers to 8 byte as well.
AS_VAR_IF([arg_autodouble],[yes],[dnl
   AC_FC_AUTODOUBLE
   AS_VAR_IF([arg_cbind],[yes],[
      AS_VAR_IF([ac_cv_fc_cbind],[unknown],[],[dnl
         AC_MSG_WARN([ISO C binding with autodouble may break MPI and VTK.])
      ])
   ])
])

# create a statically linked binary
AS_VAR_IF([arg_static],[yes],[AC_FC_STATIC])dnl

# check inline expansion flag
AS_VAR_IF([arg_inline],[yes],[AC_FC_INLINE])dnl

# check optimization flag
# ATTENTION: This check cannot test if the flag really turns on optimizations!
# Several compilers only give a warning message and don't abort with an error
# code if the given flag is not recognized. Those compilers would pass the
# test although the given flag is invalid.
AS_VAR_IF([arg_fast],[no],[],[dnl
   ac_cv_fc_FCFLAGS_save="$FCFLAGS"
   AC_CACHE_CHECK([for Fortran compiler optimizations],
      [ac_cv_fc_fast],[dnl
      ac_cv_fc_fast=unknown
      AS_VAR_IF([arg_fast],[yes],[arg_fast="$arg_fast_default"])
      m4_foreach([ac_opt_flag],[[$arg_fast],[$arg_fast_default]],[dnl
         AS_VAR_IF([ac_cv_fc_fast],[unknown],[dnl
            FCFLAGS="$ac_cv_fc_FCFLAGS_save -ac_opt_flag"
            AC_LINK_IFELSE([AC_LANG_PROGRAM([],[])],[dnl
               ac_cv_fc_fast="-ac_opt_flag"
               FCFLAGS_OPT="$FCFLAGS_OPT -ac_opt_flag"
               LDFLAGS_OPT="$LDFLAGS_OPT -ac_opt_flag"
            ])
         ])
      ])
   ])
   FCFLAGS="$ac_cv_fc_FCFLAGS_save"
   AS_VAR_IF([ac_cv_fc_fast],[unknown],[dnl
      AC_MSG_FAILURE([Invalid compiler flag.])
   ])
])

# special Fortran compiler flags to prevent name mangling
AS_VAR_IF([linkClibs],[yes],[dnl
   ac_cv_fc_FCFLAGS_save=$FCFLAGS
   ac_cv_fc_linkClibs=""
   AC_FC_NOUNDERSCORE([dnl
      ac_cv_fc_linkClibs="$ac_cv_fc_linkClibs $ac_cv_fc_nous"
      FCFLAGS="$ac_cv_fc_FCFLAGS_save $ac_cv_fc_nous"],[dnl
      AC_MSG_WARN([Cannot remove trailling underscores from function names])
      ac_cv_fc_nous="unknown"
   ])
   AC_FC_PRESERVECASE([dnl
      ac_cv_fc_linkClibs="$ac_cv_fc_linkClibs $ac_cv_fc_preservecase"],[dnl
      AC_MSG_WARN([Cannot link C programs with mixed upper/lower case functions])
      ac_cv_fc_preservecase="unknown"
   ])
   # search for C compiler
   AS_IF([test "x$ac_cv_fc_nous" = "xunknown" \
            -o "x$ac_cv_fc_preservecase" = "xunknown"],[dnl
      AC_MSG_NOTICE([trying fallback solution for linking against C libraries])
   ])
   FCFLAGS=$ac_cv_fc_FCFLAGS_save
])


# ===========================================================================
# check for some other programs / command line tools
# ===========================================================================
# text processing tool awk
AC_PROG_AWK

# stream editor sed
AC_PROG_SED

# check for scalasca (http://www.scalasca.org)
# for automatic instrumentation and profiling of
# the parallel version
AC_CHECK_PROG([PREP],[scalasca],[scalasca -instrument])


# ===========================================================================
# MPI library and module/header tests
# ===========================================================================
ac_mpi_FCFLAGS=""
ac_mpi_LDFLAGS=""
ac_mpi_LIBS=""
ac_mpi_DEFS=""
AS_VAR_IF([with_mpi],[no],[],[dnl
   ac_mpi_FCFLAGS_save="$FCFLAGS"
   ac_mpi_LDFLAGS_save="$LDFLAGS"
   ac_mpi_LIBS_save="$LIBS"
   ac_mpi_DEFS_save="$DEFS"
   # add MPI library path
   for dir in $mpi_libdirs; do
      AS_IF([test -d "$dir"],[ac_mpi_LDFLAGS="$ac_mpi_LDFLAGS -L$dir"])
   done

   # check some libraries possibly needed for MPI
   # set special compiler flags for linking against C libraries
   LDFLAGS="$ac_mpi_LDFLAGS_save $ac_mpi_LDFLAGS"
   LIBS=""
   AS_VAR_IF([ac_cv_fc_nous],[unknown],[dnl
      AC_LANG_PUSH(C)
      AC_SEARCH_LIBS(pthread_key_create,pthread,[dnl
         ac_mpi_LIBS="$LIBS $ac_mpi_LIBS"],[],[$ac_mpi_LIBS_save $ac_mpi_LIBS])
      LIBS=""
      AC_SEARCH_LIBS(aio_init,rt,[dnl
         ac_mpi_LIBS="$LIBS $ac_mpi_LIBS"],[],[$ac_mpi_LIBS_save $ac_mpi_LIBS])
      AC_LANG_POP(C)
   ],[dnl
      LDFLAGS="$ac_cv_fc_linkClibs $LDFLAGS"
      AC_SEARCH_LIBS(pthread_key_create,pthread,[dnl
         ac_mpi_LIBS="$LIBS $ac_mpi_LIBS"],[],[$ac_mpi_LIBS_save $ac_mpi_LIBS])
      LIBS=""
      AC_SEARCH_LIBS(aio_init,rt,[dnl
         ac_mpi_LIBS="$LIBS $ac_mpi_LIBS"],[],[$ac_mpi_LIBS_save $ac_mpi_LIBS])
   ])

   # check for PVFS2 libraries (possibly needed for MPI-IO)
   AS_VAR_IF([with_pvfs2],[no],[],[dnl
      # add PVFS2 library path
      for dir in $pvfs2_libdirs; do
      	  AS_IF([test -d "$dir"],[ac_mpi_LDFLAGS="$ac_mpi_LDFLAGS -L$dir"]) 
      done
      # set special Fortran compiler flags for linking against C libraries
      # if possible, otherwise use C compiler to check the linking
      LDFLAGS="$ac_mpi_LDFLAGS_save $ac_mpi_LDFLAGS"
      LIBS=""
      AS_VAR_IF([ac_cv_fc_preservecase],[unknown],[dnl
         AC_LANG_PUSH(C)
         AC_SEARCH_LIBS(BIO_new,crypto,[ac_mpi_LIBS="$LIBS $ac_mpi_LIBS"],
            [],[$ac_mpi_LIBS_save $ac_mpi_LIBS])
         AC_SEARCH_LIBS(PVFS_sys_io,pvfs2,[ac_mpi_LIBS="$LIBS $ac_mpi_LIBS"],
            [],[$ac_mpi_LIBS_save $ac_mpi_LIBS])
         AC_LANG_POP(C)
      ],[dnl
         LDFLAGS="$ac_cv_fc_linkClibs $LDFLAGS"
         AC_SEARCH_LIBS(BIO_new,crypto,[ac_mpi_LIBS="$LIBS $ac_mpi_LIBS"],
            [],[$ac_mpi_LIBS_save $ac_mpi_LIBS])
         AC_SEARCH_LIBS(PVFS_sys_io,pvfs2,[ac_mpi_LIBS="$LIBS $ac_mpi_LIBS"],
            [],[$ac_mpi_LIBS_save $ac_mpi_LIBS])
      ])
   ])

   # test linking against MPI libraries
   LDFLAGS="$ac_mpi_LDFLAGS_save $ac_mpi_LDFLAGS"
   LIBS=""
   AC_CACHE_CHECK([for library providing MPI_Init],
      ac_cv_search_MPI_Init,[
      ac_cv_search_MPI_Init=unknown
      # generate a list of various MPI library combinations
m4_define([ac_mpi_list],[m4_flatten([,]dnl
m4_combine([,],m4_combine([,],[[],[-lmpichf90]],[ ],[],[-lfmpich]),[ ],[-lmpich])[,]
m4_combine([,],m4_combine([,],[[],[-lmpi_f90]],[ ],[],[-lmpi_f77]),[ ],[-lmpi])
[,-lfmpi -lmpi])])
     # search for MPI libraries
     m4_foreach([ac_libs],ac_mpi_list,[dnl
         AS_VAR_IF([ac_cv_search_MPI_Init],[unknown],[dnl
            LIBS="ac_libs $ac_mpi_LIBS $ac_mpi_LIBS_save"
            AC_LINK_IFELSE(AC_LANG_CALL([],[MPI_Init]),[dnl
               m4_if(ac_libs,[],[dnl
                  ac_cv_search_MPI_Init="none required"],[dnl
                  ac_cv_search_MPI_Init="ac_libs"
               ])
            ])           
         ])
      ])
      AS_VAR_IF([ac_cv_search_MPI_Init],["none required"],[],[dnl
         ac_mpi_LIBS="$ac_cv_search_MPI_Init $ac_mpi_LIBS"
      ])
   ])
   AS_VAR_IF([ac_cv_search_MPI_Init],[unknown],[dnl
      AC_MSG_FAILURE([cannot find library providing MPI_Init])
   ])

   # check MPI Fortran 90 module file
   AS_VAR_IF([arg_mpi_module],[yes],[dnl
      AC_FC_MODULE([mpi],[],[$mpi_moddirs],[dnl
         AS_VAR_IF([ac_cv_module_mpi_moddir],[""],[],[dnl
            ac_mpi_FCFLAGS="$ac_mpi_FCFLAGS -I$ac_cv_module_mpi_moddir"
         ])
         ac_mpi_DEFS="$ac_mpi_DEFS -DHAVE_MPI_MOD"
      ],[dnl
        ac_cv_module_mpi_moddir=unknown
      ])
   ])

   # look for "mpif.h" header file
   AC_CACHE_CHECK([whether we have to include mpif.h for MPI],
      [ac_cv_header_mpif_h],[dnl
      ac_cv_header_mpif_h=unknown
      FCFLAGS="$ac_mpi_FCFLAGS_save $ac_mpi_FCFLAGS $ac_mpi_DEFS"
      AC_COMPILE_IFELSE([dnl
         AC_LANG_PROGRAM([],[
#ifdef HAVE_MPI_MOD
use mpi
#endif
implicit none
integer, parameter :: comm = MPI_COMM_WORLD
         ])],[ac_cv_header_mpif_h=no
      ],[dnl
         # else check if mpif.h works
         ac_include_dirs="m4_combine([ ],[[$MPIDIR],[/usr],[/usr/local]],[/],[include])"
         for ac_incdir in "" $ac_include_dirs
         do
            AS_VAR_IF([ac_incdir],[""],[],[dnl
               FCFLAGS="$ac_mpi_FCFLAGS_save $ac_mpi_FCFLAGS -I$ac_incdir"])
            AC_COMPILE_IFELSE([dnl
               AC_LANG_PROGRAM([],[
implicit none
include 'mpif.h'
integer, parameter :: comm = MPI_COMM_WORLD
               ])],[dnl
               ac_cv_header_mpif_h=yes
               AS_VAR_IF([ac_incdir],[""],[],[dnl
                  ac_mpi_FCFLAGS="$ac_mpi_FCFLAGS -I$ac_incdir"])
               ac_mpi_DEFS="$ac_mpi_DEFS -DHAVE_MPIF_H"; break
            ])
         done         
      ])
   ])
   rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
   AS_VAR_IF([ac_cv_header_mpif_h],[unknown],[dnl
      AC_MSG_FAILURE([cannot compile MPI programs])
   ])

   # look for "mpiof.h" header file
   AC_CACHE_CHECK([whether we have to include mpiof.h for MPI-IO],[ac_cv_header_mpiof_h],[dnl
      ac_cv_header_mpiof_h=unknown
      FCFLAGS="$ac_mpi_FCFLAGS_save $ac_mpi_FCFLAGS $ac_mpi_DEFS"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
#ifdef HAVE_MPI_MOD
use mpi
#endif
implicit none
#ifdef HAVE_MPIF_H
include 'mpif.h'
#endif
integer, parameter :: mode = MPI_MODE_RDONLY
])],[
      ac_cv_header_mpiof_h=no],[dnl
         AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
#ifdef HAVE_MPI_MOD
use mpi
#endif
implicit none
#ifdef HAVE_MPIF_H
include 'mpif.h'
#endif
include 'mpiof.h'
integer, parameter :: mode = MPI_MODE_RDONLY
])],[
            ac_cv_header_mpiof_h=yes
            ac_mpi_DEFS="$ac_mpi_DEFS -DHAVE_MPIOF_H"
         ])
      ])
   ])
   rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
   AS_VAR_IF([ac_cv_header_mpiof_h],[unknown],[dnl
      AC_MSG_FAILURE([cannot compile MPI programs using MPI-IO])
   ])
   ac_mpi_DEFS="-DPARALLEL $ac_mpi_DEFS"
   DEFS="$ac_mpi_DEFS_save"
   LIBS="$ac_mpi_LIBS_save"
   FCFLAGS="$ac_mpi_FCFLAGS_save"
   LDFLAGS="$ac_mpi_LDFLAGS_save"
]) # end MPI checks


# ==========================================================================
# HDF5 library tests (possibly necessary for NetCDF output)
# ===========================================================================
ac_hdf5_FCFLAGS=""
ac_hdf5_LDFLAGS=""
ac_hdf5_LIBS=""
ac_hdf5_DEFS=""
AS_VAR_IF([with_hdf5],[no],[],[dnl
   ac_hdf5_FCFLAGS_save="$FCFLAGS"
   ac_hdf5_LDFLAGS_save="$LDFLAGS"
   ac_hdf5_LIBS_save="$LIBS"
   ac_hdf5_DEFS_save="$DEFS"
   for dir in $hdf5_libdirs; do
      AS_IF([test -d "$dir"],[ac_hdf5_LDFLAGS="$ac_hdf5_LDFLAGS -L$dir"])
   done

   # check additional libraries required for linking
   # against hdf5 libraries
   # set special Fortran compiler flags for linking against C libraries
   # if possible, otherwise use C compiler to check the linking
   LDFLAGS="$ac_hdf5_LDFLAGS_save $ac_hdf5_LDFLAGS $ac_mpi_LDFLAGS"
   LIBS=""
   AS_VAR_IF([ac_cv_fc_nous],[unknown],[dnl
      AC_LANG_PUSH(C)
      AC_SEARCH_LIBS(ceil,m,[ac_hdf5_LIBS="$LIBS $ac_hdf5_LIBS"],[],
         [$ac_hdf5_LIBS_save $ac_mpi_LIBS $ac_hdf5_LIBS])
      AC_SEARCH_LIBS(inflate,z,[ac_hdf5_LIBS="$LIBS $ac_hdf5_LIBS"],[],
         [$ac_hdf5_LIBS_save $ac_mpi_LIBS $ac_hdf5_LIBS])
      AC_SEARCH_LIBS(sz_inflate,sz,[ac_hdf5_LIBS="$LIBS $ac_hdf5_LIBS"],[],
         [$ac_hdf5_LIBS_save $ac_mpi_LIBS $ac_hdf5_LIBS])
     AC_LANG_POP(C)
   ],[dnl
      LDFLAGS="$ac_cv_fc_linkClibs $LDFLAGS"
      AC_SEARCH_LIBS(ceil,m,[ac_hdf5_LIBS="$LIBS $ac_hdf5_LIBS"],[],
         [$ac_hdf5_LIBS_save $ac_mpi_LIBS $ac_hdf5_LIBS])
      AC_SEARCH_LIBS(inflate,z,[ac_hdf5_LIBS="$LIBS $ac_hdf5_LIBS"],[],
         [$ac_hdf5_LIBS_save $ac_mpi_LIBS $ac_hdf5_LIBS])
      AC_SEARCH_LIBS(sz_inflate,sz,[ac_hdf5_LIBS="$LIBS $ac_hdf5_LIBS"],[],
         [$ac_hdf5_LIBS_save $ac_mpi_LIBS $ac_hdf5_LIBS])
   ])

   # set special Fortran compiler flags for linking against C libraries
   # if possible, otherwise use C compiler to check the linking
   LDFLAGS="$ac_hdf5_LDFLAGS_save $ac_hdf5_LDFLAGS $ac_mpi_LDFLAGS"
   LIBS=""
   AS_VAR_IF([ac_cv_fc_preservecase],[unknown],[dnl
      AC_LANG_PUSH(C)
      AC_SEARCH_LIBS(H5open,hdf5,[dnl
         ac_hdf5_LIBS="$LIBS $ac_hdf5_LIBS"
         ac_hdf5_DEFS="$ac_hdf5_DEFS -DHAVE_HDF5"],[],
         [$ac_hdf5_LIBS_save $ac_mpi_LIBS $ac_hdf5_LIBS])
      LIBS=""
      AC_SEARCH_LIBS(H5DSattach_scale,hdf5_hl,[dnl
         ac_hdf5_LIBS="$LIBS $ac_hdf5_LIBS"],[],
         [$ac_hdf5_LIBS_save $ac_mpi_LIBS $ac_hdf5_LIBS])
      AC_LANG_POP(C)
   ],[dnl
      LDFLAGS="$ac_cv_fc_linkClibs $LDFLAGS"
      AC_SEARCH_LIBS(H5open,hdf5,[dnl
         ac_hdf5_LIBS="$LIBS $ac_hdf5_LIBS"
         ac_hdf5_DEFS="$ac_hdf5_DEFS -DHAVE_HDF5"],[],
         [$ac_hdf5_LIBS_save $ac_mpi_LIBS $ac_hdf5_LIBS])
      LIBS=""
      AC_SEARCH_LIBS(H5DSattach_scale,hdf5_hl,[dnl
         ac_hdf5_LIBS="$LIBS $ac_hdf5_LIBS"],[],
         [$ac_hdf5_LIBS_save $ac_mpi_LIBS $ac_hdf5_LIBS])
   ])

   FCFLAGS="$ac_hdf5_FCFLAGS_save"
   LDFLAGS="$ac_hdf5_LDFLAGS_save"
   LIBS="$ac_hdf5_LIBS_save"
   DEFS="$ac_hdf5_DEFS_save"
]) # end HDF5 checks


# ===========================================================================
# NetCDF library and module/header file tests
# ===========================================================================
ac_netcdf_FCFLAGS=""
ac_netcdf_LDFLAGS=""
ac_netcdf_LIBS=""
ac_netcdf_DEFS=""
AS_VAR_IF([with_netcdf],[no],[],[dnl
   ac_netcdf_FCFLAGS_save="$FCFLAGS"
   ac_netcdf_LDFLAGS_save="$LDFLAGS"
   ac_netcdf_LIBS_save="$LIBS"
   ac_netcdf_DEFS_save="$DEFS"
   for dir in $netcdf_libdirs; do
      AS_IF([test -d "$dir"],[ac_netcdf_LDFLAGS="$ac_netcdf_LDFLAGS -L$dir"])
   done

   # check NetCDF Fortran 90 module
   AC_FC_MODULE([netcdf],[
implicit none
integer, parameter :: nfflag = nf90_noclobber
],[$netcdf_moddirs],[dnl
      AS_VAR_IF([ac_cv_module_netcdf_moddir],[""],[],[dnl
         ac_netcdf_FCFLAGS="$ac_netcdf_FCFLAGS -I$ac_cv_module_netcdf_moddir"
      ])
      # NetCDF library checks
      AC_CACHE_CHECK([for library providing nf90_create],
         ac_cv_search_nf90_create,[
         ac_cv_search_nf90_create=unknown
         # search for NetCDF library
         LDFLAGS="$ac_netcdf_LDFLAGS_save $ac_mpi_LDFLAGS $ac_hdf5_LDFLAGS $ac_netcdf_LDFLAGS"
         for ac_libs in "m4_combine([" "],[[],[-lnetcdff],[-lnetcdff90]],[],[],[ -lnetcdf])"
         do
            AS_VAR_IF([ac_cv_search_nf90_create],[unknown],[dnl
               LIBS="$ac_libs $ac_hdf5_LIBS $ac_mpi_LIBS $ac_netcdf_LIBS_save"
               AC_LINK_IFELSE([
                  AC_LANG_PROGRAM([],[
use netcdf
implicit none
integer :: status,ncid
status = nf90_create("test.nc",nf90_noclobber,ncid)
])
               ],[dnl
                  m4_if(m4_normalize([$ac_libs]),[],[dnl
                     ac_cv_search_nf90_create="none required"],[dnl
                     ac_cv_search_nf90_create=m4_normalize([$ac_libs])
                     ac_netcdf_LIBS="$ac_cv_search_nf90_create $ac_netcdf_LIBS"
                  ])
               ])
            ])
         done
      ])
   ])

   AS_IF([test "x$ac_cv_module_netcdf" = "xno" \
         -o "x$ac_cv_search_nf90_create" = "xunknown"],[dnl
      AC_MSG_WARN([NetCDF support disabled.])
   ],[dnl
      ac_netcdf_DEFS="$ac_netcdf_DEFS -DHAVE_NETCDF"
      AS_VAR_IF([with_hdf5],[yes],[dnl
         ac_cv_netcdf_hdf5="no"
         AC_MSG_CHECKING([whether the NetCDF library was build with HDF5 support])
         AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
use netcdf
implicit none
integer, parameter :: ncfmt = nf90_netcdf4
])],[ac_cv_netcdf_hdf5=yes])
         AC_MSG_RESULT([$ac_cv_netcdf_hdf5])
      ])
      AS_VAR_IF([ac_cv_netcdf_hdf5],[yes],[dnl
         AS_VAR_IF([with_mpi],[yes],[dnl
            # this is a workaround for a bug in the NetCDF 4
            # Fortran 90 interface; both constants are missing
            # in older versions
	    ac_netcdf_mpi_def=yes
            AC_MSG_CHECKING([whether we have to define NF90_INDEPENDENT])
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
use netcdf
implicit none
integer, parameter :: comm = NF90_INDEPENDENT
])],[ac_netcdf_mpi_def=no],[ac_netcdf_DEFS="$ac_netcdf_DEFS -DNF90_INDEPENDENT=0"])
	    AC_MSG_RESULT([$ac_netcdf_mpi_def])
	    ac_netcdf_mpi_def=yes
            AC_MSG_CHECKING([whether we have to define NF90_COLLECTIVE])
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
use netcdf
implicit none
integer, parameter :: comm = NF90_COLLECTIVE
])],[ac_netcdf_mpi_def=no],[ac_netcdf_DEFS="$ac_netcdf_DEFS -DNF90_COLLECTIVE=1"])
	    AC_MSG_RESULT([$ac_netcdf_mpi_def])
         ])
      ],[dnl
         # disable all hdf5 flags
         ac_hdf5_FCFLAGS=""
         ac_hdf5_LDFLAGS=""
         ac_hdf5_LIBS=""
         ac_hdf5_DEFS=""         
         AS_VAR_IF([with_mpi],[yes],[dnl
            AC_MSG_WARN([HDF5 is required for parallel NetCDF I/O.])
         ])
      ])
   ])
   FCFLAGS="$ac_netcdf_FCFLAGS_save"
   LDFLAGS="$ac_netcdf_LDFLAGS_save"
   LIBS="$ac_netcdf_LIBS_save"
   DEFS="$ac_netcdf_DEFS_save"
]) # end NetCDF checks


# ==========================================================================
# NPY library tests
# ===========================================================================
ac_npy_FCFLAGS=""
ac_npy_LDFLAGS=""
ac_npy_LIBS=""
ac_npy_DEFS=""
AS_VAR_IF([with_npy],[no],[],[dnl
   ac_npy_FCFLAGS_save="$FCFLAGS"
   ac_npy_LDFLAGS_save="$LDFLAGS"
   ac_npy_LIBS_save="$LIBS"
   ac_npy_DEFS_save="$DEFS"
   for dir in $npy_libdirs; do
      AS_IF([test -d "$dir"],[ac_npy_LDFLAGS="$ac_npy_LDFLAGS -L$dir"])
   done

   # set special Fortran compiler flags for linking against C libraries
   # if possible, otherwise use C compiler to check the linking
   LDFLAGS="$ac_npy_LDFLAGS_save $ac_npy_LDFLAGS"
   LIBS=""
   AS_VAR_IF([ac_cv_fc_preservecase],[unknown],[dnl
      AC_LANG_PUSH(C)
      AC_SEARCH_LIBS(npy_save_double,npy,[dnl
         ac_npy_LIBS="$LIBS $ac_npy_LIBS"
		 ac_npy_FCFLAGS="$ac_npy_FCFLAGS -I$npy_moddirs"
         ac_npy_DEFS="$ac_npy_DEFS -DHAVE_NPY"],[],
         [$ac_npy_LIBS_save $ac_npy_LIBS])
      AC_LANG_POP(C)
   ],[dnl
      LDFLAGS="$ac_cv_fc_linkClibs $LDFLAGS"
      AC_SEARCH_LIBS(npy_save_double,npy,[dnl
         ac_npy_LIBS="$LIBS $ac_npy_LIBS"
		 ac_npy_FCFLAGS="$ac_npy_FCFLAGS -I$npy_moddirs"
         ac_npy_DEFS="$ac_npy_DEFS -DHAVE_NPY"],[],
         [$ac_npy_LIBS_save $ac_npy_LIBS])
   ])

   FCFLAGS="$ac_npy_FCFLAGS_save"
   LDFLAGS="$ac_npy_LDFLAGS_save"
   LIBS="$ac_npy_LIBS_save"
   DEFS="$ac_npy_DEFS_save"
]) # end NPY checks


# ===========================================================================
# set Makefile variables
# ===========================================================================
AC_SUBST(FCFLAGS_ALL,[m4_normalize("$FCFLAGS $ac_hdf5_FCFLAGS \
                                   $ac_netcdf_FCFLAGS $ac_npy_FCFLAGS")])
AC_SUBST(DEFS_ALL,[m4_normalize("$DEFS $ac_hdf5_DEFS $ac_netcdf_DEFS $ac_npy_DEFS")])
AC_SUBST(FCFLAGS_OPT,[m4_normalize("$FCFLAGS_OPT $ac_cv_fc_inline")])
AC_SUBST(FCFLAGS_DBG,[m4_normalize("$FCFLAGS_DBG -DDEBUG $ac_cv_fc_advanced_debug")])
AC_SUBST(FCFLAGS_PROF,[m4_normalize("$ac_cv_fc_profiling")])
AC_SUBST(LDFLAGS_ALL,[m4_normalize("$LDFLAGS $ac_hdf5_LDFLAGS $ac_netcdf_LDFLAGS \
                                    $ac_npy_LDFLAGS \
                                    $LIBS $ac_netcdf_LIBS $ac_hdf5_LIBS \
		 		    $ac_npy_LIBS")])
AC_SUBST(LDFLAGS_OPT)
AC_SUBST(LDFLAGS_PROF,[m4_normalize("$ac_cv_fc_profiling")])
AC_SUBST(LDFLAGS_DBG,[m4_normalize("$LDFLAGS_DBG $ac_cv_fc_advanced_debug")])

AC_SUBST(FCFLAGS_MPI,[m4_normalize("$ac_mpi_FCFLAGS")])
AC_SUBST(DEFS_PARALLEL,[m4_normalize("$ac_mpi_DEFS")])
AC_SUBST(LDFLAGS_MPI,[m4_normalize("$ac_mpi_LDFLAGS")])
AC_SUBST(MPILIBS,[m4_normalize("$ac_mpi_LIBS")])
AC_SUBST(PREP)
AC_SUBST(AR)
AC_SUBST(AWK)
AC_SUBST(SED)

# generate Makefile
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
