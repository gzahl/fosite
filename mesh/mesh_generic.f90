!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: mesh_generic.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006-2008                                                   #
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
  USE mesh_midpoint, InitMesh_basic => InitMesh
  USE mesh_trapezoidal
  USE geometry_generic
  USE fluxes_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Mesh_TYP, &
       ! constants
#ifdef PARALLEL
       DEFAULT_MPI_REAL, &
#endif
       CARTESIAN, POLAR, LOGPOLAR, CYLINDRICAL, SPHERICAL, OBLATE_SPHEROIDAL, &
       ! methods
       InitMesh, &
       Convert2Cartesian, &
       Convert2Curvilinear, &
       GetType, &
       GetRank, &
       Info, &
       Warning, &
       Error, &
       CloseMesh
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitMesh(this,Fluxes,geometry,inum,jnum,xmin,xmax,ymin,ymax,gparam)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    TYPE(Fluxes_TYP)  :: Fluxes
    INTEGER           :: geometry
    INTEGER           :: inum,jnum
    REAL              :: xmin,xmax,ymin,ymax
    REAL, OPTIONAL    :: gparam
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: xres,yres
    REAL              :: gparam_default
    !------------------------------------------------------------------------!
    INTENT(IN)        :: geometry,inum,jnum,xmin,xmax,ymin,ymax,gparam
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

    ! default geometry parameter
    IF (PRESENT(gparam)) THEN
       gparam_default = gparam
    ELSE
       gparam_default = 1.0
    END IF

    SELECT CASE(GetType(Fluxes))
    CASE(MIDPOINT)
       CALL InitMesh_midpoint(this,geometry,inum,jnum,xmin,xmax,ymin,ymax, &
            gparam_default)
    CASE(TRAPEZOIDAL)
       CALL InitMesh_trapezoidal(this,geometry,inum,jnum,xmin,xmax,ymin,ymax, &
            gparam_default)
    CASE DEFAULT
       CALL Error(this,"InitMesh", "Unknown flux type.")
    END SELECT

    ! compute cartesian coordinates for bary center values
    CALL Convert2Cartesian(this%geometry,this%bcenter,this%bccart)

    ! print some information
    WRITE (xres, '(I0)') this%INUM    ! this is just for better looking output
    WRITE (yres, '(I0)') this%JNUM
    CALL Info(this, " MESH-----> resolution:        " // TRIM(xres) // " x " // TRIM(yres))
    WRITE (xres, '(ES8.1,A,ES8.1)') this%xmin, " ..", this%xmax
    WRITE (yres, '(ES8.1,A,ES8.1)') this%ymin, " ..", this%ymax    
    CALL Info(this, "            computat. domain:  x=" // TRIM(xres) // ACHAR(10)  &
                 // "                               y="  // TRIM(yres))
  END SUBROUTINE InitMesh


  SUBROUTINE CloseMesh(this,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    TYPE(Fluxes_TYP)  :: Fluxes
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this
    !------------------------------------------------------------------------!

    SELECT CASE(GetType(Fluxes))
    CASE(MIDPOINT)
       CALL CloseMesh_midpoint(this)
    CASE(TRAPEZOIDAL)
       CALL CloseMesh_trapezoidal(this)
    CASE DEFAULT
       CALL Error(this,"InitMesh", "Unknown flux type.")
    END SELECT
  END SUBROUTINE CloseMesh


END MODULE mesh_generic
