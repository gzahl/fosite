!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: sources_pointmass.f90                                             #
!#                                                                           #
!# Copyright (C) 2007 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
! source terms module for gravitational acceleration due to
! a point mass at the center of the coordinate system
!----------------------------------------------------------------------------!
MODULE sources_pointmass
  USE sources_common, InitSources_common => InitSources
  USE fluxes_common, ONLY : Fluxes_TYP
  USE physics_generic
  USE mesh_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  CHARACTER(LEN=32), PARAMETER :: source_name = "central point mass"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Sources_TYP, &
       ! methods
       InitSources, &
       InitSources_pointmass, &
       GetType, &
       GetName, &
       ExternalSources_pointmass, &
       CloseSources_pointmass
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitSources(this,stype,sname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    INTEGER           :: stype
    CHARACTER(LEN=32) :: sname
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: newsrc, tmpsrc
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: stype,sname
    !------------------------------------------------------------------------!
    ! allocate memory for new source term
    ALLOCATE(newsrc,STAT=err)
    IF (err.NE.0) THEN
       PRINT *, "ERROR in InitSources: Unable to allocate memory"
       STOP
    END IF
    
    ! basic initialization
    CALL InitSources_common(newsrc,stype,sname)

    ! add new source term to beginning of
    ! list of source terms
    IF (ASSOCIATED(this).EQV..FALSE.) THEN
       this => newsrc
    ELSE
       tmpsrc => this
       this => newsrc
       this%next => tmpsrc
    END IF

  END SUBROUTINE InitSources


  SUBROUTINE InitSources_pointmass(this,Mesh,Physics,stype,mass)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Fluxes_TYP)  :: Fluxes
    TYPE(Physics_TYP) :: Physics
    INTEGER           :: stype
    REAL              :: mass
    !------------------------------------------------------------------------!
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) :: cart, accel
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) :: tanphi
    INTEGER           :: err
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,stype,mass
    !------------------------------------------------------------------------!
    CALL InitSources(this,stype,source_name)

    ! set mass
    this%mass = mass
    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         STAT = err)
    IF (err.NE.0) THEN
       PRINT *, "ERROR in InitSources_pointmass: Unable allocate memory!"
       STOP
    END IF

    ! convert to cartesian coordinates
    CALL Convert2Cartesian(Mesh%geometry,Mesh%bcenter,cart)

    ! initialize gravitational acceleration
    FORALL (i=Mesh%IGMIN:Mesh%IGMAX, j=Mesh%JGMIN:Mesh%JGMAX)
       ! calculate distance ratio
       tanphi(i,j) = cart(i,j,2)/cart(i,j,1)
       ! calculate cartesian components
       accel(i,j,1) = -SIGN(1.0,cart(i,j,1))* Physics%constants%GN &
            * (this%mass / cart(i,j,1)**2 ) &
            / (SQRT(1. + tanphi(i,j)**2)**3 + TINY)
       accel(i,j,2) = accel(i,j,1) * tanphi(i,j)
    END FORALL

    ! convert to curvilinear vector components
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,accel,this%accel)

    ! for advanced time step control: 
    ! minimal free fall time scale on the mesh
    this%dtmin = SQRT(MIN(MINVAL(ABS(Mesh%dlx(:,:)/this%accel(:,:,1))), &
         MINVAL(ABS(Mesh%dly(:,:)/this%accel(:,:,2)))))
  END SUBROUTINE InitSources_pointmass


  PURE SUBROUTINE ExternalSources_pointmass(this,Mesh,Physics,pvar,cvar,sterm)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    TYPE(Fluxes_TYP)  :: Fluxes
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: cvar,pvar,sterm
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,cvar,pvar
    INTENT(INOUT)     :: this
    INTENT(OUT)       :: sterm
    !------------------------------------------------------------------------!
    ! gravitational source terms
    CALL ExternalSources(Physics,Mesh,this%accel,pvar,cvar,sterm)
  END SUBROUTINE ExternalSources_pointmass

 
  SUBROUTINE CloseSources_pointmass(this,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Sources_TYP) :: this
    TYPE(Fluxes_TYP)  :: Fluxes
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Fluxes
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%accel)
  END SUBROUTINE CloseSources_pointmass

END MODULE sources_pointmass
