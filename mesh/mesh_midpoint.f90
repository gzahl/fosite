!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: mesh_midpoint.f90                                                 #
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
! mesh module for midpoint quadrature rule
!----------------------------------------------------------------------------!
MODULE mesh_midpoint
  USE geometry_generic
  USE mesh_common, ONLY : Mesh_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  REAL, PARAMETER :: TINY = 1.0E-30              ! to avoid division by 0    !
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Mesh_TYP, &
       ! methods
       InitMesh_midpoint, &
       CloseMesh_midpoint
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitMesh_midpoint(this,inum,jnum,xmin,xmax,ymin,ymax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    INTEGER           :: inum,jnum
    REAL              :: xmin,xmax,ymin,ymax    
    !------------------------------------------------------------------------!
    INTEGER           :: err,i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: inum,jnum,xmin,xmax,ymin,ymax
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

    ! allocate memory for pointers that are specific for midpoint fluxes
    ALLOCATE(this%dAx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%dAy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%dAxdy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%dAydx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%cyxy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%cxyx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%czxz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%czyz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         STAT=err)
    IF (err.NE.0) THEN
       PRINT *, "ERROR in InitMesh_midpoint: Unable to allocate memory!"
       STOP
    END IF

    ! get geometrical scale factors
    CALL ScaleFactors(this%geometry,this%fpos,this%fhx,this%fhy,this%fhz)   ! faces
    CALL ScaleFactors(this%geometry,this%center,this%bhx,this%bhy,this%bhz) ! centers

    ! surface elements
    this%dAx(:,:,1) = this%fhz(:,:,1)*this%fhy(:,:,1)*this%dy ! perpendicular to x-direction
    this%dAy(:,:,1) = this%fhz(:,:,3)*this%fhx(:,:,3)*this%dx ! perpendicular to y-direction

    ! surface elements divided by dx or dy
    this%dAxdy(:,:,1) = this%fhz(:,:,1)*this%fhy(:,:,1)       ! perpendicular to x-direction
    this%dAydx(:,:,1) = this%fhz(:,:,3)*this%fhx(:,:,3)       ! perpendicular to y-direction

    ! volume elements
    this%volume(:,:) = this%bhx(:,:)*this%bhy(:,:)*this%bhz(:,:)*this%dx*this%dy

    ! inverse volume elements multiplied by dx or dy
    this%dxdV(:,:) = 1./(this%bhx(:,:)*this%bhy(:,:)*this%bhz(:,:)*this%dy + TINY)
    this%dydV(:,:) = 1./(this%bhx(:,:)*this%bhy(:,:)*this%bhz(:,:)*this%dx + TINY)

    ! cell bary centers
    this%bcenter(:,:,:)  = this%center(:,:,:)

    ! commutator coefficients at cell centers
    this%cyxy(:,:,1) = this%bhz(:,:)*(this%fhy(:,:,2)-this%fhy(:,:,1)) * this%dydV(:,:)
    this%cxyx(:,:,1) = this%bhz(:,:)*(this%fhx(:,:,4)-this%fhx(:,:,3)) * this%dxdV(:,:)
    this%czxz(:,:,1) = this%bhy(:,:)*(this%fhz(:,:,2)-this%fhz(:,:,1)) * this%dydV(:,:)
    this%czyz(:,:,1) = this%bhx(:,:)*(this%fhz(:,:,4)-this%fhz(:,:,3)) * this%dxdV(:,:)

    ! center line elements
    this%dlx(:,:) = this%bhx(:,:)*this%dx
    this%dly(:,:) = this%bhy(:,:)*this%dy
  END SUBROUTINE InitMesh_midpoint


  SUBROUTINE CloseMesh_midpoint(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    !------------------------------------------------------------------------!

    DEALLOCATE(this%dAx,this%dAy,this%dAydx,this%dAxdy, &
         this%cyxy,this%cxyx,this%czxz,this%czyz)
  END SUBROUTINE CloseMesh_midpoint


END MODULE mesh_midpoint
