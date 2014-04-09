!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: mesh_common.f90                                                   #
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
! basic mesh module
!----------------------------------------------------------------------------!
MODULE mesh_common
  USE geometry_common, ONLY : Geometry_TYP
  USE boundary_common, ONLY : Boundary_TYP
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! mesh data structure
  TYPE Mesh_TYP
     TYPE(Geometry_TYP)         :: Geometry        ! geometrical properties  !
     TYPE(Boundary_TYP), DIMENSION(4) :: Boundary  ! one for each boundary   !
     INTEGER                    :: GNUM            ! num. ghost cells        !
     INTEGER                    :: INUM,JNUM       ! resolution              !
     INTEGER                    :: IMIN,IMAX       ! min. & max. index in x- !
     INTEGER                    :: JMIN,JMAX       ! and y-direction         !
     INTEGER                    :: IGMIN,IGMAX     ! same with ghost cells   !
     INTEGER                    :: JGMIN,JGMAX     !                         !
     REAL                       :: xmin, xmax      ! comput. domain in x-    !
     REAL                       :: ymin, ymax      ! and y-direction         !
     REAL                       :: dx,dy           ! width of cells in x-&   !
     REAL                       :: invdx, invdy    ! y-direction; inverse    !
     REAL, POINTER              :: center(:,:,:)   ! cell geometr. centers   !
     REAL, POINTER              :: bcenter(:,:,:)  ! cell bary centers       !
     REAL, POINTER              :: fpos(:,:,:,:)   ! face centered positions !
     REAL, POINTER              :: cpos(:,:,:,:)   ! corner positions        !
     REAL, POINTER              :: fhx(:,:,:)      ! face centered scale     !
     REAL, POINTER              :: fhy(:,:,:)      !  factors                !
     REAL, POINTER              :: fhz(:,:,:)      !                         !
     REAL, POINTER              :: chx(:,:,:)      ! corner scale factors    !
     REAL, POINTER              :: chy(:,:,:)      !                         !
     REAL, POINTER              :: chz(:,:,:)      !                         !
     REAL, POINTER              :: bhx(:,:)        ! bary center scale       !
     REAL, POINTER              :: bhy(:,:)        !  factors                !
     REAL, POINTER              :: bhz(:,:)        !                         !
     REAL, POINTER              :: cyxy(:,:,:)     ! commutator coefficients !
     REAL, POINTER              :: cxyx(:,:,:)     ! for geometrical sources !
     REAL, POINTER              :: czxz(:,:,:)
     REAL, POINTER              :: czyz(:,:,:)
     REAL, POINTER              :: sqrtg(:,:)      ! sqrt(det(g_ij))         !
     REAL, POINTER              :: weights(:,:,:,:)! interpolation weights   !
     REAL, POINTER              :: dlx(:,:)        ! line elements           !
     REAL, POINTER              :: dly(:,:)        !   at cell centers       !
     REAL, POINTER              :: dAx(:,:,:)      ! surface elements        !
     REAL, POINTER              :: dAy(:,:,:)      !   on cell faces         !
     REAL, POINTER              :: dAxdy(:,:,:)    ! surface elements divided!
     REAL, POINTER              :: dAydx(:,:,:)    !   by dx or dy           !
     REAL, POINTER              :: volume(:,:)     ! cell volumes            !
     REAL, POINTER              :: dxdV(:,:)       ! inverse volume elements !
     REAL, POINTER              :: dydV(:,:)       ! multiplied with dx or dy!
  END TYPE Mesh_TYP
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Mesh_TYP, &
       ! methods
       InitMesh, &
       CloseMesh
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitMesh(this,inum,jnum,xmin,xmax,ymin,ymax)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    INTEGER           :: inum,jnum
    REAL              :: xmin,xmax,ymin,ymax    
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: inum,jnum,xmin,xmax,ymin,ymax
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

    ! number of ghost rows/columns
    this%GNUM = 2

    ! resolution and index ranges
    this%INUM = inum
    this%JNUM = jnum
    this%IMIN = 1
    this%IMAX = this%INUM
    this%JMIN = 1
    this%JMAX = this%JNUM

    ! index ranges with ghost cells
    this%IGMIN = this%IMIN - this%GNUM
    this%IGMAX = this%IMAX + this%GNUM
    this%JGMIN = this%JMIN - this%GNUM
    this%JGMAX = this%JMAX + this%GNUM

    ! coordinate domain
    this%xmin = xmin
    this%xmax = xmax
    this%ymin = ymin
    this%ymax = ymax

    ! coordinate differences in each direction
    this%dx = (this%xmax - this%xmin) / this%INUM
    this%dy = (this%ymax - this%ymin) / this%JNUM
    
    ! inverse coordinate differences
    this%invdx = 1./this%dx
    this%invdy = 1./this%dy

    ! allocate memory for all pointers that are independent of fluxtype
    ALLOCATE(this%center(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%bcenter(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%fpos(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,2), &
         this%bhx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%bhy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%bhz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%fhx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%fhy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%fhz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%volume(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%dxdV(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%dydV(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%dlx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%dly(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         STAT=err)
    IF (err.NE.0) THEN
       PRINT *, "ERROR in InitMesh: Can't allocate memory!"
       STOP
    END IF

    ! calculate geometrical cell centers
    FORALL (i=this%IGMIN:this%IGMAX,j=this%JGMIN:this%JGMAX)
       this%center(i,j,1) = this%xmin + (i - this%IMIN + 0.5)*this%dx
       this%center(i,j,2) = this%ymin + (j - this%IMIN + 0.5)*this%dy
    END FORALL

    ! calculate face centered positions
    FORALL (i=this%IGMIN:this%IGMAX,j=this%JGMIN:this%JGMAX)
       this%fpos(i,j,1,1) = this%center(i,j,1) - 0.5*this%dx  ! western x coord
       this%fpos(i,j,1,2) = this%center(i,j,2)                ! western y coord
       this%fpos(i,j,2,1) = this%center(i,j,1) + 0.5*this%dx  ! eastern x coord
       this%fpos(i,j,2,2) = this%center(i,j,2)                ! eastern y coord
       this%fpos(i,j,3,1) = this%center(i,j,1)                ! southern x coord
       this%fpos(i,j,3,2) = this%center(i,j,2) - 0.5*this%dy  ! southern y coord
       this%fpos(i,j,4,1) = this%center(i,j,1)                ! northern x coord
       this%fpos(i,j,4,2) = this%center(i,j,2) + 0.5*this%dy  ! northern y coord
    END FORALL
  END SUBROUTINE InitMesh


  SUBROUTINE CloseMesh(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%center,this%bcenter,this%fpos, &
         this%bhx,this%bhy,this%bhz, &
         this%fhx,this%fhy,this%fhz, &
         this%volume,this%dxdV,this%dydV,this%dlx,this%dly)
  END SUBROUTINE CloseMesh

END MODULE mesh_common
