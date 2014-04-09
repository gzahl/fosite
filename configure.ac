#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#############################################################################
#                                                                           #
# fosite - 2D hydrodynamical simulation program                             #
# configuration file: configure.ac                                          #
#                                                                           #
# Copyright (C) 2008-2010                                                   #
# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
# Marc Junker <maj@astrophysik.uni-kiel.de>                                 #
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

AC_PREREQ(2.61)
AC_INIT(fosite, 0.3.1, tillense@astrophysik.uni-kiel.de)

# set up configure switches.
linkClibs="no"
AC_ARG_WITH([mpi],
    AS_HELP_STRING([--with-mpi=[DIR]],[
	Enable MPI support and check for installation in directory [DIR]; 
        required to compile parallel version of fosite @<:@DEFAULT=no@:>@]),[
    AS_IF([test "x$with_mpi" = "xyes"],[MPIDIR=""],[
    	test "x$with_mpi" != "xno"],[MPIDIR=$with_mpi])],[
    with_mpi="no"])

AC_ARG_WITH([pvfs2],
    AS_HELP_STRING([--with-pvfs2=[DIR]],[
	Check for PVFS2 installation in directory [DIR]; 
        recommended for parallel I/O.]),[
    PVFSDIR=$with_pvfs2
    test "x$PVFSDIR" = "xyes" && PVFSDIR=""
    linkClibs="yes"
    with_pvfs2="yes"],[with_pvfs2="no"])

AC_ARG_WITH([netcdf],
    AS_HELP_STRING([--with-netcdf=[DIR]],[
	Check for NetCDF installation in directory [DIR];
        required for NetCDF output @<:@DEFAULT=yes@:>@]),[
    AS_IF([test "x$with_netcdf" = "xyes"],[NETCDFDIR=""],
          [test "x$with_netcdf" = "xno"], [NETCDFDIR=""],[
          NETCDFDIR=$with_netcdf
          with_netcdf="yes"])],[with_netcdf="yes"])

AC_ARG_WITH([hdf5],
    AS_HELP_STRING([--with-hdf5=[DIR]],[
	Check for HDF5 installation in directory [DIR];
        required for parallel NetCDF output]),[
    AS_IF([test "x$with_hdf5" = "xyes"],[HDF5DIR=""],
          [test "x$with_netcdf" = "xno"], [HDF5DIR=""],[
          HDF5DIR=$with_hdf5
          linkClibs="yes"
          with_hdf5="yes"])],[
    with_hdf5="no"
])

AC_ARG_ENABLE([static],
	AC_HELP_STRING([--enable-static],
		[Create a statically linked executable file.
                 @<:@DEFAULT=no@:>@ ]),
		[arg_static=$enableval])

AC_ARG_ENABLE([double],
	AC_HELP_STRING([--enable-double],
		[turn on autodouble compiler flag if available. 
                 @<:@DEFAULT=no@:>@ ]),
 		[arg_autodouble=$enableval])

AC_ARG_ENABLE([fast],
	AC_HELP_STRING([--enable-fast=option],
		[turn on standard compiler optimizations;
		e.g. --enable-fast=O2 activates optimization level 2
                 @<:@DEFAULT=O@:>@ ]),[
        AS_IF([test "x$enableval" = "xyes"],[arg_fast=O],
	      [test "x$enableval" != "xno"],[arg_fast=$enableval])],[
	arg_fast=""])

AC_ARG_ENABLE([inlining],
	AC_HELP_STRING([--enable-inlining],
		[turn on inline expansion of functions;
                 @<:@DEFAULT=no@:>@ ]),
 		[arg_inline=$enableval])

AC_ARG_VAR([PREP],[preprocessor command for parallel profiling])

# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_PRESERVECASE([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# Looks for a compiler flag to prevent the Fortran (FC) compiler from
# transforming subroutine names to upper/lower case and adds it to
# FCFLAGS.  Call ACTION-IF-SUCCESS (defaults to nothing) if successful 
# (i.e. can compile code using the specific compiler flag) and
# ACTION-IF-FAILURE (defaults to failing with an error message) if not.
#
# The known flags are:
#           -names as_is: Intel compiler (ifort, ifc)
# -fsource-case-preserve: gfortran and g95 compiler
# --------------------------------------------------------------------------
AC_DEFUN_ONCE([AC_FC_PRESERVECASE],[
    AC_REQUIRE([AC_FC_NOUNDERSCORE])dnl
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to preserve upper/lower cases in names],
    	ac_cv_fc_preservecase,[
	ac_cv_fc_preservecase=unknown
	ac_fc_preservecase_FCFLAGS_save="$FCFLAGS"
	for ac_flag in none -fsource-case-preserve "-names as_is" 
	do
		test "x$ac_flag" != xnone && \
		     FCFLAGS="$ac_fc_preservecase_FCFLAGS_save $ac_flag"
    		AC_LINK_IFELSE(AC_LANG_PROGRAM,[
			AC_LINK_IFELSE(AC_LANG_CALL([],[maLLoc]),[],[
				ac_cv_fc_preservecase=$ac_flag; break])])
	done
	rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
	FCFLAGS=$ac_fc_preservecase_FCFLAGS_save
    ])
    if test "x$ac_cv_fc_preservecase" = xunknown; then
        ac_cv_fc_preservecase=""
        m4_default([$2],[
            AC_MSG_ERROR([Cannot preserve case of external C function names])
        ])
    else
	test "x$ac_cv_fc_preservecase" != xnone && \
	    FCFLAGS="$FCFLAGS $ac_cv_fc_preservecase"
	$1
    fi
    AC_LANG_POP(Fortran)dnl
])# AC_FC_PRESERVECASE

# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_NOUNDERSCORE([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# Looks for a compiler flag to prevent the Fortran (FC) compiler from
# adding trailing underscores to external subroutine names and adds it to
# FCFLAGS.  Call ACTION-IF-SUCCESS (defaults to nothing) if successful 
# (i.e. can compile code using the specific compiler flag) and
# ACTION-IF-FAILURE (defaults to failing with an error message) if not.
#
# The known flags are:
#               -nus: Intel compiler (ifort, ifc)
#  -fno-underscoring: gfortran and g95 compiler
# --------------------------------------------------------------------------
AC_DEFUN_ONCE([AC_FC_NOUNDERSCORE],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to remove trailing underscores],
    	ac_cv_fc_nous,[
	ac_cv_fc_nous=unknown
	ac_fc_nous_FCFLAGS_save="$FCFLAGS"
	for ac_flag in none -nus -fno-underscoring
	do
	    test "x$ac_flag" != xnone && \
	    	FCFLAGS="$ac_fc_nous_FCFLAGS_save $ac_flag"
            AC_LINK_IFELSE(AC_LANG_CALL([],[malloc]),[
	        ac_cv_fc_nous=$ac_flag; break])
        done
	rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
	FCFLAGS=$ac_fc_nous_FCFLAGS_save
    ])
    if test "x$ac_cv_fc_nous" = xunknown; then
        ac_cv_fc_nous=""
        m4_default([$2],[
            AC_MSG_ERROR([Cannot remove trailling underscores from external C function names])])
    else
        test "x$ac_cv_fc_nous" != xnone && \
            FCFLAGS="$FCFLAGS $ac_cv_fc_nous"
        $1
    fi
    AC_LANG_POP(Fortran)dnl
])# AC_FC_NOUNDERSCORE

# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_ISO_C_BINDING([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# Checks whether the Fortran compiler supports ISO_C_BINDING for
# linking against C functions and sets the variable ac_cv_fc_bind to
# either "yes" or "no". 
#
# REMARK: This is _not_ Fortran 90/95 but Fortran 2003 standard!
#
# --------------------------------------------------------------------------
AC_DEFUN_ONCE([AC_FC_ISO_C_BINDING],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([whether $FC supports ISO C binding],
    	ac_cv_fc_bind,[
	ac_cv_fc_bind=no
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
use iso_c_binding
call foo
contains
    subroutine foo() bind(c)
    implicit none
    print *,'dummy output'
    end subroutine foo
])],[
           ac_cv_fc_bind=yes
        ])
    ])
    rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
    AS_IF([test "x$ac_cv_fc_bind" = xno],[m4_default([$2],[AC_MSG_ERROR(
	[Fortran compiler does not support ISO C bindings])])],[$1])
    AC_LANG_POP(Fortran)dnl
])# AC_FC_ISO_C_BINDING

# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_CPPFLAG([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# Looks for a compiler flag to make the Fortran (FC) compiler process
# preprocessor directives, and adds it to FCFLAGS.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile code using the specific compiler flag) and ACTION-IF-FAILURE
# (defaults to failing with an error message) if not.
#
# The known flags are:
#               -cpp: Intel compiler (ifort, ifc) and g95 compiler
#   -x f95-cpp-input: GNU Fortran compiler (gfortran)
#                -Ep: NEC SX-8 compiler (sxf90)
# --------------------------------------------------------------------------
AC_DEFUN_ONCE([AC_FC_CPPFLAG],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to process preprocessor directives],
    	ac_cv_fc_cppflag,[
	ac_cv_fc_cppflag=unknown
	ac_fc_cppflag_FCFLAGS_save="$FCFLAGS"
	for ac_flag in none -Ep -cpp "-x f95-cpp-input"
	do
	    test "x$ac_flag" != xnone && \
	    	 FCFLAGS="$ac_fc_cppflag_FCFLAGS_save $ac_flag -DFPPTEST"
	    AC_COMPILE_IFELSE([dnl
		AC_LANG_PROGRAM([],[
#ifdef FPPTEST
      print *, 'this code should be compiled'
#else
      this not!!!
#endif
])],[           ac_cv_fc_cppflag=$ac_flag; break])
	done
	rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
	FCFLAGS=$ac_fc_cppflag_FCFLAGS_save
    ])
    if test "x$ac_cv_fc_cppflag" = xunknown; then
       m4_default([$2],[dnl
           AC_MSG_ERROR([Fortran does not process preprocessor directives])])
    else
       test "x$ac_cv_fc_cppflag" != xnone && \
	   FCFLAGS="$FCFLAGS $ac_cv_fc_cppflag"
	   $1
       $1
    fi
    AC_LANG_POP(Fortran)dnl
])# AC_FC_CPPFLAG


# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_AUTODOUBLE([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# Looks for a compiler flag to turn all real variables automatically
# into double precision variables and adds the flag to FCFLAGS.
# On success the variable ac_cv_fc_autodouble holds the compiler flag
# otherwise the string "none". Call ACTION-IF-SUCCESS
# (defaults to nothing) if successful (i.e. can compile code using
# the specific compiler flag) and ACTION-IF-FAILURE (defaults to 
# failing with an error message) if not.
#
# The known flags are:
#                -r8: Intel compiler (ifort, ifc) and g95 compiler
#   -fdefault-real-8: GNU Fortran compiler (gfortran)
#                -dW: NEC SX-8 compiler (sxf90)
#      -Wf"-A idbl4": NEC SX-9 compiler with MPI (sxmpif90)
# --------------------------------------------------------------------------
AC_DEFUN_ONCE([AC_FC_AUTODOUBLE],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to autodouble real numbers],
        ac_cv_fc_autodouble,[
	ac_cv_fc_autodouble=unknown
	ac_fc_autodouble_FCFLAGS_save="$FCFLAGS"
	for ac_flag in -r8 -fdefault-real-8 -dW '-Wf"-A idbl4"'
	do
	    test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_autodouble_FCFLAGS_save $ac_flag"
	    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
real(kind=8) :: y
y=1.0D0
call testdbl(y)
contains
    subroutine testdbl(x)
    real :: x
    print *,x
    end subroutine testdbl
])],[             ac_cv_fc_autodouble=$ac_flag; break
            ])
        done
        rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
        FCFLAGS=$ac_fc_autodouble_FCFLAGS_save
    ])
    if test "x$ac_cv_fc_autodouble" = xunknown; then
        m4_default([$2],[
            AC_MSG_ERROR([Fortran compiler does not support autodouble])])
    else
        test "x$ac_cv_fc_autodouble" != xnone && \
            FCFLAGS="$FCFLAGS $ac_cv_fc_autodouble"
        $1
    fi
    AC_LANG_POP(Fortran)dnl
])# AC_FC_AUTODOUBLE

# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_STATIC([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# Looks for a compiler flag to statically link executables
# and adds the flag to FCFLAGS.
# On success the variable ac_cv_fc_static holds the compiler flag
# otherwise the string "none". Call ACTION-IF-SUCCESS
# (defaults to nothing) if successful (i.e. can compile code using
# the specific compiler flag) and ACTION-IF-FAILURE (defaults to 
# failing with an error message) if not.
#
# The known flags are:
#            -static: most compilers
#        -Wl,-static: g95,gfortran
#   -Xlinker -static: g95,gfortran
# --------------------------------------------------------------------------
AC_DEFUN_ONCE([AC_FC_STATIC],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to create statically linked executables],
        ac_cv_fc_static,[
	ac_cv_fc_static=unknown
	ac_fc_static_FCFLAGS_save="$FCFLAGS"
	for ac_flag in -static "-Wl,-static" "-Xlinker -static"
	do
	    test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_static_FCFLAGS_save $ac_flag"
	    AC_LINK_IFELSE(AC_LANG_PROGRAM,[
	        ac_cv_fc_static=$ac_flag ; break])
	done
        rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
        FCFLAGS=$ac_fc_static_FCFLAGS_save
    ])
    if test "x$ac_cv_fc_static" = xunknown; then
        m4_default([$2],[
            AC_MSG_ERROR([Fortran compiler does not support static linking])])
    else
        test "x$ac_cv_fc_static" != xnone && \
            FCFLAGS="$FCFLAGS $ac_cv_fc_static"
        $1
    fi
    AC_LANG_POP(Fortran)dnl
])# AC_FC_STATIC

# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_INLINE([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# Looks for a compiler flag to activate inline expansion of functions.
# On success the variable ac_cv_fc_inline holds the compiler flag
# otherwise the string "none". Call ACTION-IF-SUCCESS
# (defaults to nothing) if successful (i.e. can compile code using
# the specific compiler flag) and ACTION-IF-FAILURE (defaults to 
# failing with an error message) if not.
#
# The known flags are:
# -finline-functions: gfortran, g95
#                -ip: ifort
#         -pi incdir: sxf90
# --------------------------------------------------------------------------
AC_DEFUN_ONCE([AC_FC_INLINE],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to expand functions inline],
        ac_cv_fc_inline,[
	ac_cv_fc_inline=unknown
	ac_fc_inline_FCFLAGS_save="$FCFLAGS"
	for ac_flag in -finline-functions -ip "-pi incdir"
	do
	    test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_inline_FCFLAGS_save $ac_flag"
	    AC_LINK_IFELSE(AC_LANG_PROGRAM,[
	        ac_cv_fc_inline=$ac_flag ; break])
	done
        rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
        FCFLAGS=$ac_fc_inline_FCFLAGS_save
    ])
    if test "x$ac_cv_fc_inline" = xunknown; then
        m4_default([$2],[
            AC_MSG_ERROR([Fortran compiler does not support function inlining])])
#    else
#        test "x$ac_cv_fc_inline" != xnone && \
#            FCFLAGS="$FCFLAGS $ac_cv_fc_inline"
#        $1
    fi
    AC_LANG_POP(Fortran)dnl
])# AC_FC_INLINE


# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_ADVANCED_DEBUG([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# Looks for a compiler flag for advanced debugging features, e.g. runtime
# checks for out-of-bounds array subscripts. On success the variable 
# ac_cv_fc_advanced_debug holds the compiler flag otherwise the string "none".
# Call ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can 
# compile code using the specific compiler flag) and ACTION-IF-FAILURE
# (defaults to nothing) if not.
#
# The known flags are:
#         -check all: ifort
#       -ftrace=full: g95
#                -eC: sxf90
# --------------------------------------------------------------------------
AC_DEFUN_ONCE([AC_FC_ADVANCED_DEBUG],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([whether $FC supports advanced debugging],
        ac_cv_fc_advanced_debug,[
	ac_cv_fc_advanced_debug=no
	ac_fc_advanced_debug_FCFLAGS_save="$FCFLAGS"
	for ac_flag in "-check all" "-ftrace=full" "-eC"
	do
	    test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_advanced_debug_FCFLAGS_save $ac_flag"
	    AC_LINK_IFELSE(AC_LANG_PROGRAM,[
	        ac_cv_fc_advanced_debug=$ac_flag ; break])
	done
        rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
        FCFLAGS="$ac_fc_advanced_debug_FCFLAGS_save"
    ])
    AS_IF([test "x$ac_cv_fc_advanced_debug" != xno],[$1],[$2])
    AC_LANG_POP(Fortran)dnl
])# AC_FC_ADVANCED_DEBUG

# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_STREAMS([ACTION-IF-SUPPORTED [, ACTION-IF-NOT-SUPPORTED]])
#
# Checks whether the compiler supports Fortran 2003 streams. This is
# required for VDK io. Call ACTION-IF-SUCCESS (defaults to nothing) if
# successful (i.e. can compile code) and ACTION-IF-FAILURE (defaults to 
# nothing) if not.
#
# --------------------------------------------------------------------------
#	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
#open(10,FILE='test',ACCESS='STREAM',ACTION='READWRITE')
AC_DEFUN_ONCE([AC_FC_STREAMS],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([whether $FC supports Fortran 2003 streams],
    	ac_cv_fc_streams,[
	ac_cv_fc_streams=no
        AC_RUN_IFELSE(AC_LANG_PROGRAM([],[
implicit none
character :: c
open(UNIT=10,FILE='conftest.txt',ACCESS='STREAM',ACTION='WRITE')
write (UNIT=10) "A"
close(10)
open(UNIT=10,FILE='conftest.txt',ACCESS='STREAM',ACTION='WRITE',POSITION='APPEND')
write (UNIT=10) "B"
close(10)
open(UNIT=10,FILE='conftest.txt',ACCESS='STREAM',ACTION='READ',POSITION='REWIND')
read (UNIT=10) c
if (c.NE."A") call exit(1)
]),[ac_cv_fc_streams=yes],[],[
		AC_MSG_WARN([cannot run test program while cross compiling])
 	])
    ])
    rm -f conftest.err conftest.$ac_objext conftest.$ac_ext conftest.txt
    AS_IF([test "x$ac_cv_fc_streams" = xyes],[$1],[$2])
    AC_LANG_POP(Fortran)dnl
])# AC_FC_STREAMS


# ===========================================================================
# main configure script starts here!
# ===========================================================================

# determine the target host architecture we are building for
AC_CANONICAL_HOST

# host specific settings
AS_IF([test "x$host_cpu" = "xsx8"],[
    FCFLAGS="$FCFLAGS -$host_cpu"
    LDFLAGS="$LDFLAGS -$host_cpu"
    DEFS_ALL="$DEFS_ALL -DNECSX8"]
)

AS_IF([test "x$host_cpu" = "xsx9"],[
    FCFLAGS="$FCFLAGS -$host_cpu"
    LDFLAGS="$LDFLAGS -$host_cpu"
    DEFS_ALL="$DEFS_ALL -DNECSX9"]
)
# default programming language
AC_LANG(Fortran)

# Fortran 90/95 compiler checks
ac_prog_fc_FCFLAGS_save=$FCFLAGS
MPIFC=""
AS_IF([test "x$with_mpi" != "xno"],[
    AS_IF([test "x$MPIDIR" != "x"],[
        AS_IF([test "x$FC" = "x"],[
	    # check for MPI compiler command in $MPIDIR/bin
            ac_mpifc_PATH_save=$PATH
            PATH=$MPIDIR/bin
	    AC_PATH_PROGS([MPIFC],[mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c])
	    AS_IF([test "x$MPIFC" = "x"],[
	        AC_MSG_WARN([no MPI compiler found in $PATH])],[
		FC=$MPIFC
            ])
            PATH=$ac_mpifc_PATH_save
	])
    ])
    AS_IF([test "x$MPIFC" = "x"],[
        MPIFC="mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c"
    ])
])
# check all other compilers if MPI is disabled or
# detection of MPI compiler failed
AC_PROG_FC([$MPIFC sxf90 f90 f95 gfortran ifort g95])

# default file name extension for Fortran 90
AC_FC_SRCEXT(f90)

# debug flags and normal compile flags
AS_IF([test "x$ac_cv_prog_fc_g" = "xyes"],[
   FCFLAGS=$ac_prog_fc_FCFLAGS_save
   FCFLAGS_DBG="$FCFLAGS -g"
   LDFLAGS_DBG="$LDFLAGS -g"
   AS_IF([test "x$ac_cv_fc_compiler_gnu" = "xyes"],[
      FCFLAGS_DBG="$FCFLAGS_DBG -O2"])
])

# enable advanced debugging features if available
AC_FC_ADVANCED_DEBUG([FCFLAGS_DBG="$FCFLAGS_DBG $ac_cv_fc_advanced_debug"
	LDFLAGS_DBG="$LDFLAGS_DBG $ac_cv_fc_advanced_debug"])

# ensure that the Fortran compiler accepts free-form
AC_FC_FREEFORM

# check whether the Fortran compiler supports ISO C bindings
#AC_FC_ISO_C_BINDING

# check preprocessor flag
AC_FC_CPPFLAG

# check Fortran 2003 streams
AC_FC_STREAMS([DEFS_ALL="$DEFS_ALL -DFORTRAN_STREAMS"])

# check for autodouble flag if requested
AS_IF([test "x$arg_autodouble" = "xyes"],[
   AC_FC_AUTODOUBLE
])

# create a statically linked binary
AS_IF([test "x$arg_static" = "xyes"],[
   AC_FC_STATIC
])

# check inline expansion flag
AS_IF([test "x$arg_inline" = "xyes"],[
   AC_FC_INLINE
])

# check optimization flag
AS_IF([test "x$arg_fast" != "x"],[
   # FIXME add check command here
   ac_cv_fc_fast="-$arg_fast"
])

# Fortran 90 library flags, i. e.  $FCLIBS
#AC_FC_LIBRARY_LDFLAGS

# special Fortran compiler flags to prevent name mangling
AS_IF([test "x$linkClibs" = "xyes"],[
   ac_cv_fc_FCFLAGS_save=$FCFLAGS
   AC_FC_NOUNDERSCORE([],[
      AC_MSG_WARN([Cannot remove trailling underscores from function names])
      ac_cv_fc_nous=unknown
   ])
   AS_IF([test "x$ac_cv_fc_nous" != "xnone" \
	    -a "x$ac_cv_fc_nous" != "xunknown"],[
      ac_cv_fc_FCFLAGS_linkClibs="$ac_cv_fc_FCFLAGS_linkClibs $ac_cv_fc_nous"
   ])
   AC_FC_PRESERVECASE([],[
      AC_MSG_WARN([Cannot link C programs with mixed upper/lower case functions])
      ac_cv_fc_preservecase=unknown
   ])
   AS_IF([test "x$ac_cv_fc_preservecase" != "xnone" \
            -a "x$ac_cv_fc_preservecase" != "xunknown"],[     
      ac_cv_fc_FCFLAGS_linkClibs="$ac_cv_fc_FCFLAGS_linkClibs $ac_cv_fc_preservecase"
   ])
   # search for C compiler
   AS_IF([test "x$ac_cv_fc_nous" = "xunknown" \
            -o "x$ac_cv_fc_preservecase" = "xunknown"],[
      AC_MSG_NOTICE([trying fallback solution for linking against C libraries])
      # search for MPI compiler if requested
      MPICC=""
      AS_IF([test "x$with_mpi" != "xno"],[
         AS_IF([test "x$MPIDIR" != "x"],[
            AS_IF([test "x$CC" = "x"],[
	       # check for MPI compiler command in $MPIDIR/bin
               ac_mpicc_PATH_save=$PATH
               PATH=$MPIDIR/bin
	       AC_PATH_PROGS([MPICC],[mpicc hcc mpxlc_r mpxlc mpcc cmpicc])
	       AS_IF([test "x$MPICC" = "x"],[
	          AC_MSG_WARN([no MPI compiler found in $PATH])],[
		  CC=$MPICC
               ])
               PATH=$ac_mpicc_PATH_save
	    ])
         ])
         AS_IF([test "x$MPICC" = "x"],[
            MPICC="mpicc"
         ])
      ])
      AC_PROG_CC([$MPICC cc xlC xlc pgcc icc pathcc gcc])
   ])
   FCFLAGS=$ac_cv_fc_FCFLAGS_save
])

# MPI library and header checks
AS_IF([test "x$with_mpi" != "xno"],[
   MPILIBS=""
   ac_mpi_LDFLAGS_save=$LDFLAGS
   ac_mpi_LIBS_save=$LIBS
   ac_mpi_FCFLAGS_save=$FCFLAGS
   AS_IF([test "x$MPIDIR" != "x"],[
       	   LDFLAGS_MPI="$LDFLAGS_MPI -L$MPIDIR/lib"
	   LDFLAGS="$LDFLAGS $LDFLAGS_MPI"       
   ])
   # set special compiler flags for linking against C libraries
   FCFLAGS="$FCFLAGS $ac_cv_fc_FCFLAGS_linkClibs"
   # some libraries possibly needed for MPI
   AC_CHECK_LIB(pthread, pthread_key_create,[
           MPILIBS=$MPILIBS" -lpthread"])
   AC_CHECK_LIB(rt, aio_init,[
           MPILIBS=$MPILIBS" -lrt"])
   LIBS="$ac_mpi_LIBS_save $MPILIBS"
   # some libraries possibly needed for MPI-IO
   # (if compiled with pvfs support)
   AS_IF([test "x$with_pvfs2" != "xno"],[
       AS_IF([test "x$PVFSDIR" != "x"],[
       	   LDFLAGS_MPI="$LDFLAGS_MPI -L$PVFSDIR/lib"
	   LDFLAGS="$LDFLAGS $LDFLAGS_MPI"
       ])
       AS_IF([test "x$ac_cv_fc_preservecase" = "xunknown"],[
          AC_LANG_PUSH([C])])
       AC_CHECK_LIB(crypto, BIO_new,[
           MPILIBS=$MPILIBS" -lcrypto"])
       LIBS="$ac_mpi_LIBS_save $MPILIBS"
       AC_CHECK_LIB(pvfs2, PVFS_sys_io,[
           MPILIBS=$MPILIBS" -lpvfs2"])
       LIBS="$ac_mpi_LIBS_save $MPILIBS"
       AS_IF([test "x$ac_cv_fc_preservecase" = "xunknown"],[
          AC_LANG_POP([C])])
   ])
   FCFLAGS=$ac_mpi_FCFLAGS_save
   # test linking against MPI libraries
   AC_SEARCH_LIBS([MPI_Init],[mpi fmpi mpich mpichf90],[MPILIBS=$LIBS],[
       AC_MSG_FAILURE([cannot find library providing MPI_Init])
   ])
   # look for "mpif.h" header file
   AS_IF([test x != x"$MPILIBS"],[
       AC_MSG_CHECKING([for mpif.h])
       ac_mpi_FCFLAGS_mpifh_save=$FCFLAGS
       ac_mpi_INCDIR_mpifh=unknown
       for ac_incdir in none "-I$MPIDIR/include"
       do
           test "x$ac_incdir" != xnone && \
       	       FCFLAGS="$ac_mpi_FCFLAGS_mpifh_save $ac_incdir"
           AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      include 'mpif.h'])],[
	       ac_mpi_INCDIR_mpifh=$ac_incdir
	       AC_MSG_RESULT([yes]); break
           ])
       done
       rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
       FCFLAGS=$ac_mpi_FCFLAGS_mpifh_save
       AS_IF([test "x$ac_mpi_INCDIR_mpifh" = "xunknown"],[
           AC_MSG_RESULT([no])
	   AC_MSG_FAILURE([cannot compile MPI programs])],
       [test "x$ac_mpi_INCDIR_mpifh" != "xnone"],[
	   FCFLAGS_MPI="$FCFLAGS_MPI $ac_mpi_INCDIR_mpifh"
       ])
   ])
   LIBS=$ac_mpi_LIBS_save
   FCFLAGS=$ac_mpi_FCFLAGS_save
   LDFLAGS=$ac_mpi_LDFLAGS_save
   # check for scalasca for parallel profiling
   AC_CHECK_PROG([PREP],[scalasca],[scalasca -instrument])
])
# NetCDF library and module file checks
AS_IF([test "x$with_netcdf" != "xno"],[
   ac_netcdf_FCFLAGS_save=$FCFLAGS
   ac_netcdf_LDFLAGS_save=$LDFLAGS
   ac_netcdf_LIBS_save=$LIBS
   AS_IF([test "x$NETCDFDIR" != "x"],[
      LDFLAGS="$LDFLAGS -L$NETCDFDIR/lib"
   ])
   # HDF5 Fortran library linking   
   AS_IF([test "x$with_hdf5" != "xno"],[
      ac_hdf5_FCFLAGS_save=$FCFLAGS
      ac_hdf5_LDFLAGS_save=$LDFLAGS
      ac_hdf5_LIBS_save=$LIBS
      AS_IF([test "x$HDF5DIR" != "x"],[
      	   LDFLAGS="$LDFLAGS -L$HDF5DIR/lib"
      ])
      ac_hdf5_LDFLAGS_nompi=$LDFLAGS
      AS_IF([test "x$LDFLAGS_MPI" != "x"],[
           LDFLAGS="$LDFLAGS $LDFLAGS_MPI"
      ])
      # set special compiler flags for linking against C libraries
      AS_IF([test "x$ac_cv_fc_nous" = "xunknown" \
               -o "x$ac_cv_fc_preservecase" = "xunknown"],[
         AC_LANG_PUSH(C)
         AC_CHECK_LIB(z,inflate)
         AC_CHECK_LIB(hdf5,H5open,[
	    DEFS_ALL="$DEFS_ALL -DHAVE_HDF5"
            LIBS="$LIBS -lhdf5"],[],[$MPILIBS -lm])
         AC_CHECK_LIB(hdf5_hl,H5DSattach_scale,[],[],[$MPILIBS -lm])
         AC_LANG_POP(C)],[
         FCFLAGS="$FCFLAGS $ac_cv_fc_FCFLAGS_linkClibs"
         AC_CHECK_LIB(z,inflate)
         AC_CHECK_LIB(hdf5,H5open,[
	    DEFS_ALL="$DEFS_ALL -DHAVE_HDF5"
            LIBS="$LIBS -lhdf5"],[],[$MPILIBS])
         AC_CHECK_LIB(hdf5_hl,H5DSattach_scale,[],[],[$MPILIBS])
      ])
      FCFLAGS=$ac_hdf5_FCFLAGS_save
      LDFLAGS=$ac_hdf5_LDFLAGS_nompi
   ])
   # NetCDF Fortran bindings
   # check for netcdf.mod
   AC_CACHE_CHECK([for NetCDF Fortran module],
      ac_cv_netcdf_module,[
      ac_cv_netcdf_module=no
      for ac_incdir in "$NETCDFDIR/include" none /usr/include /usr/include/netcdf \
          /usr/local/include /usr/local/include/netcdf
      do
          AS_IF([test "x$ac_incdir" != "xnone"],[
             FCFLAGS="$ac_netcdf_FCFLAGS_save -I$ac_incdir"])
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[use netcdf])],[
             ac_cv_netcdf_module=yes; break])
      done
      rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
   ])
   AS_IF([test "x$ac_cv_netcdf_module" = "xno"],[
      # no netcdf module -> check for netcdf.inc
      AC_CACHE_CHECK([for netcdf.inc],
         ac_cv_netcdf_include,[
         ac_cv_netcdf_include=no
         for ac_incdir in "$NETCDFDIR/include" none /usr/include /usr/include/netcdf \
             /usr/local/include /usr/local/include/netcdf
         do
         AS_IF([test "x$ac_incdir" != "xnone"],[
            FCFLAGS="$ac_netcdf_FCFLAGS_save -I$ac_incdir"])
      # gfortran version 4.4 doesn't care !!!
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
include 'netcdf.inc'
])],[
                ac_cv_netcdf_include=yes; break])
         done
         rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
      ])
   ])
   # check for netcdf libraries with Fortran 90 bindings
   AC_CACHE_CHECK([for nf90_create in -lnetcdf...],
      ac_cv_lib_netcdf_nf90_create,[
      ac_cv_lib_netcdf_nf90_create=no
      ac_netcdf_LDFLAGS_nompi=$LDFLAGS
      ac_netcdf_LIBS_nompi=$LIBS
      AS_IF([test "x$ac_cv_netcdf_module" = "xyes"],[
         FCFLAGS="$FCFLAGS -DHAVE_NETCDF_MOD"])
      for ac_netcdf_LIBS in none "-lnetcdf" "-lnetcdf -lnetcdff"
      do
         AS_IF([test "x$ac_netcdf_LIBS" != "xnone"],[
            LIBS="$ac_netcdf_LIBS $ac_netcdf_LIBS_nompi"])
         AS_IF([test "x$LDFLAGS_MPI" != "x"],[
            FCFLAGS="$LDFLAGS $LDFLAGS_MPI"
            LIBS="$LIBS $MPILIBS"
         ])
         AC_LINK_IFELSE([AC_LANG_PROGRAM([],[
#ifdef HAVE_NETCDF_MOD
use netcdf
#else
include 'netcdf.inc'
#endif
integer :: status,ncid
status = nf90_create("test.nc",nf90_noclobber,ncid)
])],[
	       ac_cv_lib_netcdf_nf90_create=yes; break])
      done
      rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
   ])
   AS_IF([test "x$ac_cv_lib_netcdf_nf90_create" = "xno"],[
      FCFLAGS=$ac_netcdf_FCFLAGS_save
      LDFLAGS=$ac_netcdf_LDFLAGS_save
      LIBS=$ac_netcdf_LIBS_save
      AC_MSG_WARN([NetCDF support disabled])],[
      LDFLAGS=$ac_netcdf_LDFLAGS_nompi
      LIBS="$ac_netcdf_LIBS $ac_netcdf_LIBS_nompi"
      DEFS_ALL="$DEFS_ALL -DHAVE_NETCDF"
      AS_IF([test "x$ac_cv_module_netcdf" = "xyes"],[
         DEFS_ALL="$DEFS_ALL -DHAVE_NETCDF_MOD"
      ])
   ])
])


# compiler and target specific flags
case $FC in
	sxf90)
		FCFLAGS_PROF="-ftrace"
		LDFLAGS_PROF="-ftrace"
		;;
	sxmpif90)
                FCFLAGS_PROF="-ftrace"
                LDFLAGS_PROF="-ftrace"
                ;;
	ifort)
		FCFLAGS_PROF="-p"
		LDFLAGS_PROF="-p"
		;;
	g95)
		FCFLAGS_PROF="-pg"
		LDFLAGS_PROF="-pg"
		;;
	gfortran)
		FCFLAGS_PROF="-pg"
		LDFLAGS_PROF="-pg"
		;;
	*)
		FCFLAGS_PROF=""
		LDFLAGS_PROF=""
esac

AC_SUBST(FCFLAGS_ALL,[m4_normalize($FCFLAGS)])
AC_SUBST(DEFS_ALL,[m4_normalize($DEFS_ALL)])
AC_SUBST(FCFLAGS_OPT,[m4_normalize("$FCFLAGS_OPT $ac_cv_fc_fast $ac_cv_fc_inline")])
AC_SUBST(FCFLAGS_DBG)
AC_SUBST(FCFLAGS_PROF)
AC_SUBST(LDFLAGS_ALL,[m4_normalize("$LDFLAGS $LIBS")])
AC_SUBST(LDFLAGS_OPT)
AC_SUBST(LDFLAGS_PROF)
AC_SUBST(LDFLAGS_DBG)

AC_SUBST(FCFLAGS_MPI)
AC_SUBST(DEFS_PARALLEL,["-DPARALLEL"])
AC_SUBST(LDFLAGS_MPI)
AC_SUBST(MPILIBS)
AC_SUBST(PREP)

# generate Makefile
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
