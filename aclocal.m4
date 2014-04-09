#                                               -*- Autoconf -*-
# This file is included by configure.ac
#############################################################################
#                                                                           #
# fosite - 2D hydrodynamical simulation program                             #
# additional autoconf makros: aclocal.m4                                    #
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
# -fsource-case-preserve: old gfortran and g95 compiler
# --------------------------------------------------------------------------
AC_DEFUN([AC_FC_PRESERVECASE],[
    AC_REQUIRE([AC_FC_NOUNDERSCORE])dnl
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to preserve upper/lower cases in names],
    	ac_cv_fc_preservecase,[
	ac_cv_fc_preservecase=unknown
	ac_fc_preservecase_FCFLAGS_save="$FCFLAGS"
	for ac_flag in none -fsource-case-preserve "-names as_is" 
	do
           AS_VAR_IF([ac_flag],[none],[],[dnl
              FCFLAGS="$ac_fc_preservecase_FCFLAGS_save $ac_flag"])
	   AC_LINK_IFELSE([AC_LANG_PROGRAM],[dnl
              AC_LINK_IFELSE(AC_LANG_CALL([],[maLLoc]),[],[dnl
                 ac_cv_fc_preservecase=$ac_flag; break])
           ])
	done
	rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
	FCFLAGS=$ac_fc_preservecase_FCFLAGS_save
    ])
    AS_VAR_IF([ac_cv_fc_preservecase],[unknown],[dnl
        ac_cv_fc_preservecase=""
        m4_default([$2],[dnl
            AC_MSG_ERROR([Cannot preserve case of external C function names])])
    ],[dnl
        AS_VAR_IF([ac_cv_fc_preservecase],[none],[dnl
           ac_cv_fc_preservecase=""])
        m4_default([$1],[dnl
           FCFLAGS="$FCFLAGS $ac_cv_fc_preservecase"])
    ])
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
AC_DEFUN([AC_FC_NOUNDERSCORE],[
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
    AS_VAR_IF([ac_cv_fc_nous],[unknown],[dnl
        ac_cv_fc_nous=""
        m4_default([$2],[
            AC_MSG_ERROR([Cannot remove trailling underscores \
                          from external C function names])])
    ],[dnl
        AS_VAR_IF([ac_cv_fc_nous],[none],[dnl
           ac_cv_fc_nous=""])
        m4_default([$1],[
           FCFLAGS="$FCFLAGS $ac_cv_fc_nous"])
    ])
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
AC_DEFUN([AC_FC_ISO_C_BINDING],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([whether the Fortran compiler supports ISO C binding],
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
    AS_VAR_IF([ac_cv_fc_bind],[no],[dnl
        m4_default([$2],[dnl
           AC_MSG_ERROR([Fortran compiler does not support ISO C bindings])])
    ],[$1])
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
#       -Mpreprocess: Portland group compiler (pgf90)
#                -Ep: NEC SX-8 compiler (sxf90)
# --------------------------------------------------------------------------
AC_DEFUN([AC_FC_CPPFLAG],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to process preprocessor directives],
    	ac_cv_fc_cppflag,[
	ac_cv_fc_cppflag=unknown
	ac_fc_cppflag_FCFLAGS_save="$FCFLAGS"
	for ac_flag in none -Ep -cpp "-x f95-cpp-input" -Mpreprocess
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
    AS_VAR_IF([ac_cv_fc_cppflag],[unknown],[dnl
       m4_default([$2],[dnl
           AC_MSG_ERROR([Fortran does not process preprocessor directives])])
    ],[dnl
       AS_VAR_IF([ac_cv_fc_cppflag],[none],[dnl
          ac_cv_fc_cppflag=""])
       m4_default([$1],[dnl
           FCFLAGS="$FCFLAGS $ac_cv_fc_cppflag"])
    ])
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
AC_DEFUN([AC_FC_AUTODOUBLE],[
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
    AS_VAR_IF([ac_cv_fc_autodouble],[unknown],[dnl
        ac_cv_fc_autodouble=""
        m4_default([$2],[
            AC_MSG_ERROR([Fortran compiler does not support autodouble])])
    ],[dnl
        AS_VAR_IF([ac_cv_fc_autodouble],[none],[dnl
           ac_cv_fc_autodouble=""])
        m4_default([$1],[dnl
	    FCFLAGS="$FCFLAGS $ac_cv_fc_autodouble"])
    ])
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
AC_DEFUN([AC_FC_STATIC],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to create statically linked executables],
        ac_cv_fc_static,[
	ac_cv_fc_static=unknown
	ac_fc_static_FCFLAGS_save="$FCFLAGS"
	for ac_flag in "-Wl,-static" "-Xlinker -static" -static
	do
	    test "x$ac_flag" != xnone && FCFLAGS="$ac_fc_static_FCFLAGS_save $ac_flag"
	    AC_LINK_IFELSE(AC_LANG_PROGRAM,[
	        ac_cv_fc_static=$ac_flag ; break])
	done
        rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
        FCFLAGS=$ac_fc_static_FCFLAGS_save
    ])
    AS_VAR_IF([ac_cv_fc_static],[unknown],[dnl
        ac_cv_fc_static=""
        m4_default([$2],[
            AC_MSG_ERROR([Fortran compiler does not support static linking])])
    ],[dnl
        AS_VAR_IF([ac_cv_fc_static],[none],[dnl
           ac_cv_fc_static=""])
        m4_default([$1],[dnl
           FCFLAGS="$FCFLAGS $ac_cv_fc_static"
	   LDFLAGS="$LDFLAGS $ac_cv_fc_static"])
    ])
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
#           -Minline: pgf90
#         -pi incdir: sxf90
# --------------------------------------------------------------------------
AC_DEFUN([AC_FC_INLINE],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to expand functions inline],
        ac_cv_fc_inline,[
	ac_cv_fc_inline=unknown
	ac_fc_inline_FCFLAGS_save="$FCFLAGS"
	for ac_flag in -finline-functions -ip "-pi incdir" -Minline
	do
	    FCFLAGS="$ac_fc_inline_FCFLAGS_save $ac_flag"
	    AC_LINK_IFELSE(AC_LANG_PROGRAM,[
	        ac_cv_fc_inline=$ac_flag ; break])
	done
        rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
        FCFLAGS=$ac_fc_inline_FCFLAGS_save
    ])
    AS_VAR_IF([ac_cv_fc_inline],[unknown],[dnl
        m4_default([$2],[dnl
            AC_MSG_ERROR([Fortran compiler does not support function inlining])])],[$1])
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
#        -fbacktrace: gfortran
#        -fcheck=all: gfortran
#       -ftrace=full: g95
#         -check all: ifort
#                -eC: sxf90
# --------------------------------------------------------------------------
AC_DEFUN([AC_FC_ADVANCED_DEBUG],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([whether the Fortran compiler supports advanced debugging],
        ac_cv_fc_advanced_debug,[
	ac_cv_fc_advanced_debug=no
	ac_fc_advanced_debug_FCFLAGS_save="$FCFLAGS"
	for ac_flag in "-fcheck=all" "-fbacktrace" "-ftrace=full" "-check all" "-eC"
	do
	    FCFLAGS="$ac_fc_advanced_debug_FCFLAGS_save $ac_flag"
	    AC_LINK_IFELSE(AC_LANG_PROGRAM,[
	        ac_cv_fc_advanced_debug=$ac_flag ; break])
	done
        rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
        FCFLAGS="$ac_fc_advanced_debug_FCFLAGS_save"
    ])
    AS_VAR_IF([ac_cv_fc_advanced_debug],[no],[$2],[$1])
    AC_LANG_POP(Fortran)dnl
])# AC_FC_ADVANCED_DEBUG


# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_PROFILING([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# Looks for a compiler flag for profiling features. On success the variable 
# ac_cv_fc_profiling holds the compiler flag otherwise the string "none".
# Call ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can 
# compile code using the specific compiler flag) and ACTION-IF-FAILURE
# (defaults to nothing) if not.
#
# The known flags are:
#                -pg: gfortran, g95
#                 -p: ifort
#            -ftrace: sxf90, sxmpif90
# --------------------------------------------------------------------------
AC_DEFUN([AC_FC_PROFILING],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([whether the Fortran compiler supports runtime profiling],
        ac_cv_fc_profiling,[
	ac_cv_fc_profiling=no
	ac_fc_profiling_FCFLAGS_save="$FCFLAGS"
	for ac_flag in "-pg" "-p" "-ftrace"
	do
            FCFLAGS="$ac_fc_profiling_FCFLAGS_save $ac_flag"
	    AC_LINK_IFELSE(AC_LANG_PROGRAM,[
	        ac_cv_fc_profiling=$ac_flag ; break])
	done
        rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
        FCFLAGS="$ac_fc_profiling_FCFLAGS_save"
    ])
    AS_VAR_IF([ac_cv_fc_profiling],[no],[$2],[$1])
    AC_LANG_POP(Fortran)dnl
])# AC_FC_PROFILING


# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_STREAMS([ACTION-IF-SUPPORTED [, ACTION-IF-NOT-SUPPORTED]])
#
# Checks whether the compiler supports Fortran 2003 streams. This is
# required for VDK io. Call ACTION-IF-SUPPORTED (defaults to nothing) if
# successful (i.e. can compile code) and ACTION-IF-NOT-SUPPORTED (defaults
# to nothing) if not.
#
# --------------------------------------------------------------------------
AC_DEFUN([AC_FC_STREAMS],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([whether the Fortran compiler supports Fortran 2003 streams],
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
    AS_VAR_IF([ac_cv_fc_streams],[yes],[$1],[$2])
    AC_LANG_POP(Fortran)dnl
])# AC_FC_STREAMS


# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AC_FC_MODULE(MODULE [, FORTRAN-CODE [, MODULE-DIR [, ACTION-IF-FOUND 
#                [, ACTION-IF-NOT-FOUND]]]])
#
# Checks for the Fortran 90 module MODULE compiling the program code
# given in FORTRAN-CODE (if specified). It looks for the module file in 
# directory MODULE-DIR (if given) than in some standard directories.
# Call ACTION-IF-FOUND (defaults to nothing) if the module exists and
# ACTION-IF-NOT-FOUND (defaults to nothing) if not.
#
# --------------------------------------------------------------------------
AC_DEFUN([AC_FC_MODULE],[dnl
ac_module_dirs="$3 m4_combine([ ],[[/usr],[/usr/local]],[/],[include],[lib],[lib64])"
AS_VAR_PUSHDEF([ac_Module], [ac_cv_module_$1])dnl
AS_VAR_PUSHDEF([ac_ModDir], [ac_cv_module_$1_moddir])dnl
AC_CACHE_CHECK([for $1 Fortran 90 module],[ac_Module],[dnl
   ac_Module=no
   ac_Module_FCFLAGS_save=$FCFLAGS
   # check for $1 in each directory given by ac_module_dirs
   for ac_dir in "" $ac_module_dirs
   do
      ac_ModDir=$ac_dir
      AS_VAR_IF([ac_dir],[""],[],[dnl
         FCFLAGS="$ac_Module_FCFLAGS_save -I$ac_dir"
      ])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
use $1
$2])],[
         ac_Module=yes; break],[dnl
         FCFLAGS="$ac_Module_FCFLAGS_save"
         ac_ModDir=""
      ])
   done
   rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
])
AS_VAR_IF([ac_Module],[yes],[$4],[dnl
   m4_default([$5],[dnl
      AC_MSG_ERROR([cannot compile programs using Fortran 90 module $1])])
])
AS_VAR_PUSHDEF([ac_ModDir])dnl
AS_VAR_POPDEF([ac_Module])dnl
]) # AC_F90_MODULE
