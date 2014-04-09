!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: mesh_generic.f90                                                  #
!#                                                                           #
!# Copyright (C) 2006 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
  USE mesh_common, InitMesh_common => InitMesh, CloseMesh_common => CloseMesh
  USE mesh_midpoint
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
       CARTESIAN, POLAR, CYLINDRICAL, SPHERICAL, OBLATE_SPHEROIDAL, &
       ! methods
       InitMesh, &
       Convert2Cartesian, &
       Convert2Curvilinear, &
       GetType, &
       CloseMesh
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitMesh(this,Fluxes,geometry,inum,jnum,xmin,xmax,ymin,ymax,gparam)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    TYPE(Fluxes_TYP)  :: Fluxes
    INTEGER           :: geometry
    INTEGER           :: nvar
    INTEGER           :: inum,jnum
    REAL              :: xmin,xmax,ymin,ymax
    REAL, OPTIONAL    :: gparam
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,dec
    CHARACTER(LEN=32) :: fmt
    !------------------------------------------------------------------------!
    INTENT(IN)        :: geometry,inum,jnum,xmin,xmax,ymin,ymax,gparam
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

    ! basic mesh initialization
    CALL InitMesh_common(this,inum,jnum,xmin,xmax,ymin,ymax)

    ! initialize geometry
    IF (PRESENT(gparam)) THEN
       CALL InitGeometry(this%geometry,geometry,gparam)
    ELSE
       CALL InitGeometry(this%geometry,geometry)
    END IF

    SELECT CASE(GetType(Fluxes))
    CASE(MIDPOINT)
       CALL InitMesh_midpoint(this,inum,jnum,xmin,xmax,ymin,ymax)
    CASE(TRAPEZOIDAL)
       CALL InitMesh_trapezoidal(this,inum,jnum,xmin,xmax,ymin,ymax)
    CASE DEFAULT
       PRINT *, "ERROR in InitMesh: unknown flux type"
       STOP
    END SELECT

    ! count number of decimals for resolution;
    ! this is just for better looking output
    dec = this%INUM
    i=0
    DO WHILE (dec.GT.0)
       dec = dec / 10 
       i=i+1
    END DO
    dec = this%JNUM
    j=0
    DO WHILE (dec.GT.0)
       dec = dec / 10
       j=j+1
    END DO
    ! output format string
    WRITE (fmt, '(A,I1,A,I1,A)') "(A,I", i, ",A,I", j, ")"

    ! print some information
    PRINT fmt, " MESH-----> resolution:        ", &
         this%INUM, " x ", this%JNUM
    PRINT "(A,ES8.1,A,ES8.1)", "            computat. domain:  x=", &
         this%xmin, " ..", this%xmax
    PRINT "(A,ES8.1,A,ES8.1)", "                               y=", &
         this%ymin, " ..", this%ymax
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
       PRINT *, "ERROR in CloseMesh: unknown flux type"
       STOP
    END SELECT
    ! call basic mesh deconstructor
    CALL CloseMesh_common(this)
  END SUBROUTINE CloseMesh


END MODULE mesh_generic
