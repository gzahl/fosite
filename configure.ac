#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#############################################################################
#                                                                           #
# fosite - 2D hydrodynamical simulation program                             #
# configuration file: configure.ac                                          #
#                                                                           #
# Copyright (C) 2008 Tobias Illenseer <tillense@astrophysik.uni-kiel.de>,   #
#                    Marc Junker <maj@astrophysik.uni-kiel.de>              #
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
AC_INIT(fosite, VERSION, tillense@astrophysik.uni-kiel.de)

# set up configure switches.
AC_ARG_WITH([mpi],
    AS_HELP_STRING([--with-mpi=[DIR]],[
	Enable MPI support and check for installation in directory [DIR]; 
        required to compile parallel version of fosite @<:@DEFAULT=no@:>@]),[
    MPIDIR=$with_mpi
    test "x$MPIDIR" = "xyes" && MPIDIR=""
    with_mpi="yes"],[with_mpi="no"])

AC_ARG_WITH([pvfs2],
    AS_HELP_STRING([--with-pvfs2=[DIR]],[
	Check for PVFS2 installation in directory [DIR]; 
        recommended for parallel I/O.]),[
    PVFSDIR=$with_pvfs2
    test "x$PVFSDIR" = "xyes" && MPIDIR=""
    with_pvfs2="yes"],[with_pvfs2="no"])

#AC_ARG_WITH([intel],
#	AC_HELP_STRING([--with-intel=[DIR]],
#		[Check for Intel Compiler installation in directory [DIR]; 
#                 required if some libraries are linked against the
#		 intel support libraries]),
#		[intel_dir=$withval])

AC_ARG_ENABLE([static],
	AC_HELP_STRING([--enable-static],
		[Create a statically linked executable file.
                 @<:@DEFAULT=no@:>@ ]),
		[link_static=$enableval])

AC_ARG_ENABLE([double],
	AC_HELP_STRING([--enable-double],
		[This option enables the autodouble compiler flag if available. 
                 @<:@DEFAULT=no@:>@ ]),
 		[arg_autodouble=$enableval])

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
    		AC_LINK_IFELSE(AC_LANG_CALL([],[maLLoc]),
			[], [ac_cv_fc_preservecase=$ac_flag; break])
	done
	rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
	FCFLAGS=$ac_fc_preservecase_FCFLAGS_save
    ])
    if test "x$ac_cv_fc_preservecase" = xunknown; then
        m4_default([$2],[
            AC_MSG_ERROR([Cannot link Fortran program against C libraries])
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
        m4_default([$2],[
            AC_MSG_ERROR([Cannot link Fortran program against C libraries])])
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
#                -ew: NEC SX-8 compiler (sxf90)
# --------------------------------------------------------------------------
AC_DEFUN_ONCE([AC_FC_AUTODOUBLE],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to autodouble real numbers],
        ac_cv_fc_autodouble,[
	ac_cv_fc_autodouble=unknown
	ac_fc_autodouble_FCFLAGS_save="$FCFLAGS"
	for ac_flag in -r8 -fdefault-real-8 -ew
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
            AC_MSG_FAILURE([Fortran compiler does not support autodouble])])
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
	for ac_flag in -static "-Wl,-static"
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
            AC_MSG_FAILURE([Fortran compiler does not support static linking])])
    else
        test "x$ac_cv_fc_static" != xnone && \
            FCFLAGS="$FCFLAGS $ac_cv_fc_static"
        $1
    fi
    AC_LANG_POP(Fortran)dnl
])# AC_FC_STATIC


# ===========================================================================
# main configure script starts here!
# ===========================================================================

# determine the target host architecture we are building for
AC_CANONICAL_HOST

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
AC_PROG_FC([$MPIFC sxf90 f90 f95 ifort g95 gfortran])

# default file name extension for Fortran 90
AC_FC_SRCEXT(f90)

# debug flags and normal compile flags
AS_IF([test "x$ac_cv_prog_fc_g" = "xyes"],[
   FCFLAGS_DBG=$FCFLAGS
   FCFLAGS=$ac_prog_fc_FCFLAGS_save
   LDFLAGS_DBG="-g"
])

# ensure that the Fortran compiler accepts free-form
AC_FC_FREEFORM

# check preprocessor flag
AC_FC_CPPFLAG

# check for autodouble flag if requested
AS_IF([test "x$arg_autodouble" = "xyes"],[
   AC_FC_AUTODOUBLE
])

# create a statically linked binary
AS_IF([test "x$link_static" = "xyes"],[
   AC_FC_STATIC
])

# Fortran 90 library flags, i. e.  $FCLIBS
#AC_FC_LIBRARY_LDFLAGS

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
   AC_FC_NOUNDERSCORE
   AC_FC_PRESERVECASE
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
       AC_CHECK_LIB(crypto, BIO_new,[
           MPILIBS=$MPILIBS" -lcrypto"])
       LIBS="$ac_mpi_LIBS_save $MPILIBS"
       AC_CHECK_LIB(pvfs2, PVFS_sys_io,[
           MPILIBS=$MPILIBS" -lpvfs2"])
       LIBS="$ac_mpi_LIBS_save $MPILIBS"
   ])
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
])

# compiler and target specific flags
case $FC in
	sxf90)
		FCFLAGS_OPT="-pi incdir -$host_cpu"
		FCFLAGS_PROF="-ftrace"
		LDFLAGS_PROF="-ftrace"
		;;
	ifort)
		FCFLAGS_OPT="-O3 -ip"
		if test "y$host_cpu" = "yx86_64" ; then
		   FCFLAGS_OPT=$FCFLAGS_OPT" -mcpu=pentium4"
		fi
# FIXME: optimization for intel processors
# doesn't work on AMD x86_64 compatible processors
#		if test "y$host_cpu" = "yx86_64" ; then
#		   FCFLAGS_OPT=$FCFLAGS_OPT" -msse3"
#		fi
		FCFLAGS_DBG=$FCFLAGS_DBG" -check all"
		FCFLAGS_PROF="-p"
		LDFLAGS_PROF="-p"
		;;
	g95)
		FCFLAGS_OPT="-O3"
		FCFLAGS_DBG=$FCFLAGS_DBG" -ftrace=full"
		FCFLAGS_PROF="-pg"
		LDFLAGS_PROF="-pg"
		;;
	gfortran)
		FCFLAGS_OPT="-O3"
		FCFLAGS_PROF="-pg"
		LDFLAGS_PROF="-pg"
		;;
	*)	FCFLAGS_OPT="-O"
		FCFLAGS_PROF=""
		LDFLAGS_PROF=""
esac

AC_SUBST(FCFLAGS_ALL,$FCFLAGS)
AC_SUBST(FCFLAGS_OPT)
AC_SUBST(FCFLAGS_DBG)
AC_SUBST(FCFLAGS_PROF)
AC_SUBST(LDFLAGS_ALL,$ac_cv_fc_libs_all)
AC_SUBST(LDFLAGS_OPT)
AC_SUBST(LDFLAGS_PROF)
AC_SUBST(LDFLAGS_DBG)

# MPI environment checks
#AC_MSG_NOTICE([Performing MPI compiler and library checks...])
#ACX_MPI([],dnl
#[AC_MSG_WARN([no MPI environment found; parallel version won't compile])
#])

AC_SUBST(FCFLAGS_MPI)
AC_SUBST(LDFLAGS_MPI)
AC_SUBST(MPILIBS)

# generate Makefile
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
