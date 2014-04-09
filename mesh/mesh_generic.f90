!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: mesh_generic.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2012                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!#                                                                           #
!# This program is free software; you can redistribute it and/or modify      #
!# it under the terms of the GNU General Public License as published by      #
!# the Free Software Foundation; either version 2 of the License, or (at     #
!# your option) any later version.                                           #
!#                                                                           #
!# This program is distributed in the hope that it will be useful, but       #
!# WITHOUT ANY WARRANTY; without even the implied warranty of                #
!# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, GOOD TITLE or        #
!# NON INFRINGEMENT.  See the GNU General Public License for more            #
!# details.                                                                  #
!#                                                                           #
!# You should have received a copy of the GNU General Public License         #
!# along with this program; if not, write to the Free Software               #
!# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
!#                                                                           #
!#############################################################################

!----------------------------------------------------------------------------!
! generic mesh module
!----------------------------------------------------------------------------!
MODULE mesh_generic
  USE mesh_midpoint, InitMesh_common => InitMesh, CloseMesh_common => CloseMesh
  USE mesh_trapezoidal
  USE geometry_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: MIDPOINT     = 1
  INTEGER, PARAMETER :: TRAPEZOIDAL  = 2
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Mesh_TYP, &
       Selection_TYP, &
       ! constants
       PI, &
#ifdef PARALLEL
       DEFAULT_MPI_REAL, &
#endif
       MIDPOINT, TRAPEZOIDAL, &
       CARTESIAN, POLAR, LOGPOLAR, TANPOLAR, SINHPOLAR, &
       CYLINDRICAL, TANCYLINDRICAL, SPHERICAL, SINHSPHERICAL, &
       OBLATE_SPHEROIDAL, &
       ! methods
       InitMesh, &
       Convert2Cartesian, &
       Convert2Curvilinear, &
       Divergence, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error, &
       CloseMesh
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitMesh(this,meshtype,geometry,inum,jnum,xmin,xmax,ymin,ymax, &
                      gparam)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    INTEGER           :: meshtype
    INTEGER           :: geometry
    INTEGER           :: inum,jnum
    REAL              :: xmin,xmax,ymin,ymax
    REAL, OPTIONAL    :: gparam
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: xres,yres
    !------------------------------------------------------------------------!
    INTENT(IN)        :: meshtype,geometry,inum,jnum,xmin,xmax,ymin,ymax, &
                         gparam
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
         
    SELECT CASE(meshtype)
    CASE(MIDPOINT)
       CALL InitMesh_midpoint(this,meshtype,geometry,inum,jnum,xmin,xmax,ymin, &
                              ymax,gparam)
    CASE(TRAPEZOIDAL)
       CALL InitMesh_trapezoidal(this,meshtype,geometry,inum,jnum,xmin,xmax,ymin, &
                                 ymax,gparam)
    CASE DEFAULT
       CALL Error(this,"InitMesh", "Unknown mesh type.")
    END SELECT

    ! compute cartesian coordinates for bary center values
    CALL Convert2Cartesian(this%geometry,this%bcenter,this%bccart)
    ! compute cartesian coordinates for corner positions
    CALL Convert2Cartesian(this%geometry,this%cpos,this%cpcart)

    ! print some information
    CALL Info(this, " MESH-----> quadrature rule:   " // TRIM(GetName(this)))
    WRITE (xres, '(I0)') this%INUM    ! this is just for better looking output
    WRITE (yres, '(I0)') this%JNUM
    CALL Info(this, "            resolution:        " // TRIM(xres) // " x " // TRIM(yres))
    WRITE (xres, '(ES9.2,A,ES9.2)') this%xmin, " ..", this%xmax
    WRITE (yres, '(ES9.2,A,ES9.2)') this%ymin, " ..", this%ymax
    CALL Info(this, "            computat. domain:  x=" // TRIM(xres) // ACHAR(10)  &
                 // "                               y="  // TRIM(yres))
#ifdef PARALLEL
    WRITE (xres, '(I0)') this%dims(1)
    WRITE (yres, '(I0)') this%dims(2)
    CALL Info(this, "            MPI partition:     " // TRIM(xres) // " x " // TRIM(yres))
#endif
    ! print warning message if radial coordinate is negative
    SELECT CASE(geometry)
    CASE(POLAR,LOGPOLAR,TANPOLAR,SINHPOLAR,SPHERICAL,SINHSPHERICAL)
       IF(this%xmin < 0.) THEN
          CALL Warning(this,"InitMesh","xmin < 0 could crash the code and won't" // &
               ACHAR(10) // "yield the expected results.")
       END IF
    CASE(CYLINDRICAL,TANCYLINDRICAL)
       IF(this%ymin < 0.) THEN
          CALL Warning(this,"InitMesh","ymin < 0 could crash the code and won't" // &
               ACHAR(10) // "yield the expected results.")
       END IF
    CASE DEFAULT
       ! do nothing
    END SELECT
  END SUBROUTINE InitMesh


  SUBROUTINE CloseMesh(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(this)) &
        CALL Error(this,"CloseMesh","not initialized")
!CDIR IEXPAND
    SELECT CASE(GetType(this))
    CASE(MIDPOINT)
       CALL CloseMesh_midpoint(this)
    CASE(TRAPEZOIDAL)
       CALL CloseMesh_trapezoidal(this)
    CASE DEFAULT
       CALL Error(this,"CloseMesh", "Unknown mesh type.")
    END SELECT
  END SUBROUTINE CloseMesh


END MODULE mesh_generic
