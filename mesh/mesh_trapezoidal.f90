!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: mesh_trapezoidal.f90                                              #
!#                                                                           #
!# Copyright (C) 2007-2012                                                   #
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
! mesh module for trapezoidal quadrature rule
!----------------------------------------------------------------------------!
MODULE mesh_trapezoidal
  USE mesh_midpoint
  USE geometry_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: mesh_name = "trapezoidal"
  ! precision for Newton-Raphson (see CalculateWeights)
  REAL, PARAMETER :: EPS  = 1.0D-04
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       InitMesh_trapezoidal, &
       CloseMesh_trapezoidal
  !--------------------------------------------------------------------------!

CONTAINS


  SUBROUTINE InitMesh_trapezoidal(this,meshtype,geometry,inum,jnum,xmin,xmax, &
                                  ymin,ymax,gparam)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    INTEGER           :: meshtype
    INTEGER           :: geometry
    INTEGER           :: inum,jnum
    REAL              :: xmin,xmax,ymin,ymax    
    REAL, OPTIONAL    :: gparam
    !------------------------------------------------------------------------!
    INTEGER           :: err
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: meshtype,geometry,inum,jnum,xmin,xmax,ymin,ymax, &
                         gparam
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!

    ! basic mesh and geometry initialization
    CALL InitMesh(this,meshtype,mesh_name,geometry,inum,jnum,xmin,xmax,ymin, &
                  ymax,gparam)

    CALL Warning(this,"InitMesh_trapezoidal", &
         "this module is experimental and may produce false results")

    ! allocate memory for pointers that are specific for trapezoidal fluxes
    ALLOCATE(this%dAx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%dAy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%dAxdy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%dAydx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%cyxy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%cxyx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%czxz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%czyz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%chx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%chy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%chz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%sqrtg(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%weights(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2,2), &
         STAT=err)
    IF (err.NE.0) &
         CALL Error(this,"InitMesh_trapezoidal", "Unable to allocate memory.")

    ! corner scale factors
    CALL ScaleFactors(this%geometry,this%cpos,this%chx,this%chy,this%chz)
    
       
    ! scale factors at face centered positions
    CALL ScaleFactors(this%geometry,this%fpos,this%fhx,this%fhy,this%fhz)

!CDIR COLLAPSE
    DO j=this%JGMIN,this%JGMAX
!CDIR NODEP
       DO i=this%IGMIN,this%IGMAX
          ! surface elements
          ! perpendicular to x-direction
          this%dAx(i,j,1) = this%chz(i,j,1)*this%chy(i,j,1)*this%dy   ! south-west
          this%dAx(i,j,2) = this%chz(i,j,3)*this%chy(i,j,3)*this%dy   ! north-west
          ! perpendicular to y-direction
          this%dAy(i,j,1) = this%chz(i,j,1)*this%chx(i,j,1)*this%dx   ! south-west
          this%dAy(i,j,2) = this%chz(i,j,2)*this%chx(i,j,2)*this%dx   ! south-east

          ! surface elements devided by dx or dy
          ! perpendicular to x-direction
          this%dAxdy(i,j,1) = this%chz(i,j,1)*this%chy(i,j,1)         ! south-west
          this%dAxdy(i,j,2) = this%chz(i,j,3)*this%chy(i,j,3)         ! north-west
          ! perpendicular to y-direction
          this%dAydx(i,j,1) = this%chz(i,j,1)*this%chx(i,j,1)         ! south-west
          this%dAydx(i,j,2) = this%chz(i,j,2)*this%chx(i,j,2)         ! south-east

          ! volume elements
          this%volume(i,j) = 0.25*this%dx*this%dy*SUM(this%chx(i,j,:) &
               *this%chy(i,j,:)*this%chz(i,j,:))
          ! inverse volume elements multiplied by dx or dy
          this%dxdV(i,j) = 1./(0.25*this%dy*SUM(this%chx(i,j,:) &
               *this%chy(i,j,:)*this%chz(i,j,:)) + TINY(1.0))
          this%dydV(i,j) = 1./(0.25*this%dx*SUM(this%chx(i,j,:) &
               *this%chy(i,j,:)*this%chz(i,j,:)) + TINY(1.0))

          ! cell bary centers
          this%bcenter(i,j,1)  = 0.25*this%dx * this%dydV(i,j) * (&
               this%cpos(i,j,1,1) * (this%chx(i,j,1)*this%chy(i,j,1)*this%chz(i,j,1) &
               + this%chx(i,j,3)*this%chy(i,j,3)*this%chz(i,j,3)) + &
               this%cpos(i,j,2,1) * (this%chx(i,j,2)*this%chy(i,j,2)*this%chz(i,j,2) &
               + this%chx(i,j,4)*this%chy(i,j,4)*this%chz(i,j,4)))
          
          this%bcenter(i,j,2)  = 0.25*this%dx * this%dydV(i,j) * (&
               this%cpos(i,j,2,2) * (this%chx(i,j,1)*this%chy(i,j,1)*this%chz(i,j,1) &
               + this%chx(i,j,2)*this%chy(i,j,2)*this%chz(i,j,2)) + &
               this%cpos(i,j,3,2) * (this%chx(i,j,3)*this%chy(i,j,3)*this%chz(i,j,3) &
               + this%chx(i,j,4)*this%chy(i,j,4)*this%chz(i,j,4)))

          ! commutator coefficients
          ! south-west
          this%cyxy(i,j,1) = 0.25 * this%dydV(i,j) * &
               this%chz(i,j,1)*(this%chy(i,j,2)-this%chy(i,j,1))
          this%cxyx(i,j,1) = 0.25 * this%dxdV(i,j) * &
               this%chz(i,j,1)*(this%chx(i,j,3)-this%chx(i,j,1))
          this%czxz(i,j,1) = 0.25 * this%dydV(i,j) * &
               this%chy(i,j,1)*(this%chz(i,j,2)-this%chz(i,j,1))
          this%czyz(i,j,1) = 0.25 * this%dxdV(i,j) * &
               this%chx(i,j,1)*(this%chz(i,j,3)-this%chz(i,j,1))
          ! south-east
          this%cyxy(i,j,2) = 0.25 * this%dydV(i,j) * & 
               this%chz(i,j,2)*(this%chy(i,j,2)-this%chy(i,j,1))
          this%cxyx(i,j,2) = 0.25 * this%dxdV(i,j) * &
               this%chz(i,j,2)*(this%chx(i,j,4)-this%chx(i,j,2))
          this%czxz(i,j,2) = 0.25 * this%dydV(i,j) * &
               this%chy(i,j,2)*(this%chz(i,j,2)-this%chz(i,j,1))
          this%czyz(i,j,2) = 0.25 * this%dxdV(i,j) * &
               this%chx(i,j,2)*(this%chz(i,j,4)-this%chz(i,j,2))
          ! north-west
          this%cyxy(i,j,3) = 0.25 * this%dydV(i,j) * &
               this%chz(i,j,3)*(this%chy(i,j,4)-this%chy(i,j,3))
          this%cxyx(i,j,3) = 0.25 * this%dxdV(i,j) * &
               this%chz(i,j,3)*(this%chx(i,j,3)-this%chx(i,j,1))
          this%czxz(i,j,3) = 0.25 * this%dydV(i,j) * &
               this%chy(i,j,3)*(this%chz(i,j,4)-this%chz(i,j,3))
          this%czyz(i,j,3) = 0.25 * this%dxdV(i,j) * &
               this%chx(i,j,3)*(this%chz(i,j,3)-this%chz(i,j,1))
          ! north-east
          this%cyxy(i,j,4) = 0.25 * this%dydV(i,j) * &
               this%chz(i,j,4)*(this%chy(i,j,4)-this%chy(i,j,3))
          this%cxyx(i,j,4) = 0.25 * this%dxdV(i,j) * &
               this%chz(i,j,4)*(this%chx(i,j,4)-this%chx(i,j,2))
          this%czxz(i,j,4) = 0.25 * this%dydV(i,j) * &
               this%chy(i,j,4)*(this%chz(i,j,4)-this%chz(i,j,3))
          this%czyz(i,j,4) = 0.25 * this%dxdV(i,j) * &
               this%chx(i,j,4)*(this%chz(i,j,4)-this%chz(i,j,2))
          
          ! square root of determinant of metric (upper left corner)
          this%sqrtg(i,j) = this%chx(i,j,4)*this%chy(i,j,4)*this%chz(i,j,4)
          
          ! line elements at cell bary centers
          this%dlx(i,j) = 0.5 * this%dx * (this%fhx(i,j,1)+this%fhx(i,j,2)) 
          this%dly(i,j) = 0.5 * this%dy * (this%fhy(i,j,3)+this%fhy(i,j,4))
       END DO
    END DO

    ! scale factors at bary center
    CALL ScaleFactors(this%geometry,this%bcenter,this%bhx,this%bhy,this%bhz)
    
    ! calculate weights for 2D interpolation
    CALL CalculateWeights(this)
  END SUBROUTINE InitMesh_trapezoidal


  SUBROUTINE CalculateWeights(this)
    !------------------------------------------------------------------------!
    ! calculate weights; they are used for interpolation of cell             !
    ! (bary-) center values to corner values                                 !
    ! the algorithm uses a bilinear mapping scheme to map the coordinates of !
    ! the four bary centers surrounding one corner onto a unit square;       !
    ! the Newton-Raphson method solves the non-linear set of 2 equations for !
    ! the coordinates of the mapped corner                                   !
    !------------------------------------------------------------------------!
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)         :: this
    !------------------------------------------------------------------------!
    REAL, DIMENSION(2,2,2) :: x
    REAL, DIMENSION(2,2)   :: df,dfinv
    REAL, DIMENSION(2)     :: y,z,dz,f,a0,a1,a2,a3
    REAL                   :: det,norm_dz
    INTEGER                :: i,j,n
    !------------------------------------------------------------------------!
    INTENT(INOUT)          :: this
    !------------------------------------------------------------------------!

    DO j=this%JGMIN,this%JGMAX-1
       DO i=this%IGMIN,this%IGMAX-1
          ! extracting bary center coordinates of the four neighbouring cells
          x(:,:,:) = this%bcenter(i:i+1,j:j+1,:)
          ! coordinates of the upper right cell corner
          y(:) = this%cpos(i,j,4,:)

          ! Newton-Raphson; solution vector: z(:)
          ! some coefficients
          a0(:) = x(1,1,:)
          a1(:) = x(1,1,:) + x(2,2,:) - x(2,1,:) - x(1,2,:)
          a2(:) = x(2,1,:) - x(1,1,:)
          a3(:) = x(1,2,:) - x(1,1,:)
          z(:) = 0.5         ! initial guess -> center of the unit square
          DO n=1,50
             ! function vector for root finding
             f(:) = a0(:) + a1(:)*z(1)*z(2) + a2(:)*z(1) + a3(:)*z(2) - y(:)
             ! jacobian
             df(:,1) = a1(:)*z(2) + a2(:)
             df(:,2) = a1(:)*z(1) + a3(:)
             ! determinant
             det = df(1,1)*df(2,2) - df(1,2)*df(2,1)
             IF (det.EQ.0) THEN
                CALL Error(this,"CalculateWeights", "det=0 while inverting matrix")
             END IF
             ! invert the jacobian
             dfinv(1,1) = df(2,2) / det
             dfinv(2,2) = df(1,1) / det
             dfinv(1,2) = - df(2,1) / det
             dfinv(2,1) = - df(1,2) / det
             ! update the solution vector
             dz(:) = - MATMUL(dfinv,f)
             z(:) = z(:) + dz(:)
             norm_dz = SQRT(dz(1)*dz(1)+dz(2)*dz(2))
             IF (norm_dz.LE.EPS) EXIT
          END DO

          IF (n.GE.10) THEN
             CALL Error(this,"CalculateWeights", "Newton-Raphson not convergent")
          END IF

          ! set the weights for the corner enclosed by
          ! cells (i,j) (i+1,j) (i,j+1) (i+1,j+1)
          this%weights(i,j,2,2)     = (1.0 - z(1))*(1.0 - z(2)) ! upper right in cell (i,j)
          this%weights(i+1,j,1,2)   = z(1)*(1.0 - z(2))         ! upper left in cell (i+1,j)
          this%weights(i,j+1,2,1)   = (1.0 - z(1))*z(2)         ! lower right in cell (i,j+1)
          this%weights(i+1,j+1,1,1) = z(1)*z(2)                 ! lower left in cell (i+1,j+1)
       END DO
    END DO

  END SUBROUTINE CalculateWeights


  SUBROUTINE CloseMesh_trapezoidal(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%dAx,this%dAy,this%dAydx,this%dAxdy, &
         this%cyxy,this%cxyx,this%czxz,this%czyz,this%sqrtg, &
         this%chx,this%chy,this%chz,this%weights)
    ! call basic mesh deconstructor
    CALL CloseMesh(this)
!    CALL CloseMesh_midpoint(this)
  END SUBROUTINE CloseMesh_trapezoidal


END MODULE mesh_trapezoidal
