!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: mesh_midpoint.f90                                                 #
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
! mesh module for midpoint quadrature rule
!----------------------------------------------------------------------------!
MODULE mesh_midpoint
  USE mesh_common, InitMesh_common => InitMesh
  USE boundary_common, ONLY : EAST,WEST,NORTH,SOUTH
  USE geometry_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  INTERFACE Divergence
     MODULE PROCEDURE VectorDivergence2D, TensorDivergence2D, TensorDivergence3D
  END INTERFACE 
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: mesh_name = "midpoint"
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
       ! methods
       InitMesh, &
       CloseMesh, &
       InitMesh_midpoint, &
       Divergence, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error, &
       CloseMesh_midpoint
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitMesh(this,meshtype,meshname,geometry,inum,jnum,xmin,xmax, &
                      ymin,ymax,gparam)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    INTEGER           :: meshtype
    CHARACTER(LEN=32) :: meshname
    INTEGER           :: geometry
    INTEGER           :: inum,jnum
    REAL              :: xmin,xmax,ymin,ymax
    REAL, OPTIONAL    :: gparam
    !------------------------------------------------------------------------!
    INTENT(IN)        :: meshtype,meshname,geometry,inum,jnum,xmin,xmax,ymin,ymax, &
                         gparam
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! basic mesh initialization
    CALL InitMesh_common(this,meshtype,meshname,inum,jnum,xmin,xmax,ymin,ymax)

    ! initialize geometry
    CALL InitGeometry(this%geometry,geometry,gparam)
  END SUBROUTINE InitMesh


  SUBROUTINE InitMesh_midpoint(this,meshtype,geometry,inum,jnum,xmin,xmax, &
                               ymin,ymax,gparam)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    INTEGER           :: meshtype
    INTEGER           :: geometry
    INTEGER           :: inum,jnum
    REAL              :: xmin,xmax,ymin,ymax    
    REAL              :: gparam
    !------------------------------------------------------------------------!
    INTEGER           :: err
    !------------------------------------------------------------------------!
    INTENT(IN)        :: meshtype,geometry,inum,jnum,xmin,xmax,ymin,ymax, &
                         gparam
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    
    ! basic mesh and geometry initialization
    CALL InitMesh(this,meshtype,mesh_name,geometry,inum,jnum,xmin,xmax, &
                  ymin,ymax,gparam)

    ! allocate memory for pointers that are specific for midpoint fluxes
    ALLOCATE(this%dAx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%dAy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%dAxdy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%dAydx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%cyxy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%cxyx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%czxz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%czyz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%sqrtg(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         STAT=err)
    IF (err.NE.0) THEN
       CALL Error(this,"InitMesh_midpoint", "Unable to allocate memory.")
    END IF

    ! get geometrical scale factors
    CALL ScaleFactors(this%geometry,this%fpos,this%fhx,this%fhy,this%fhz)   ! faces
    CALL ScaleFactors(this%geometry,this%center,this%bhx,this%bhy,this%bhz) ! centers

    ! surface elements divided by dx or dy
    this%dAxdy(:,:,1) = this%fhz(:,:,1)*this%fhy(:,:,1)       ! perpendicular to x-direction
    this%dAydx(:,:,1) = this%fhz(:,:,3)*this%fhx(:,:,3)       ! perpendicular to y-direction

    ! surface elements
    this%dAx(:,:,1) = this%dAxdy(:,:,1)*this%dy               ! perpendicular to x-direction
    this%dAy(:,:,1) = this%dAydx(:,:,1)*this%dx               ! perpendicular to y-direction

    ! volume elements
    this%volume(:,:) = this%bhx(:,:)*this%bhy(:,:)*this%bhz(:,:)*this%dx*this%dy

    ! inverse volume elements multiplied by dx or dy
    this%dxdV(:,:) = 1./(this%bhx(:,:)*this%bhy(:,:)*this%bhz(:,:)*this%dy + TINY(1.0))
    this%dydV(:,:) = 1./(this%bhx(:,:)*this%bhy(:,:)*this%bhz(:,:)*this%dx + TINY(1.0))

    ! cell bary centers
    this%bcenter(:,:,:)  = this%center(:,:,:)

    ! commutator coefficients 
    this%cyxy(:,:,1) = 0.5*(this%fhz(:,:,2)+this%fhz(:,:,1)) &
         * (this%fhy(:,:,2)-this%fhy(:,:,1)) * this%dydV(:,:)
    this%cxyx(:,:,1) = 0.5*(this%fhz(:,:,4)+this%fhz(:,:,3)) &
         * (this%fhx(:,:,4)-this%fhx(:,:,3)) * this%dxdV(:,:)
    this%czxz(:,:,1) = 0.5*(this%fhy(:,:,2)+this%fhy(:,:,1)) &
         * (this%fhz(:,:,2)-this%fhz(:,:,1)) * this%dydV(:,:)
    this%czyz(:,:,1) = 0.5*(this%fhx(:,:,4)+this%fhx(:,:,3)) &
         * (this%fhz(:,:,4)-this%fhz(:,:,3)) * this%dxdV(:,:)

    ! square root of determinant of metric 
    this%sqrtg(:,:) = this%bhx(:,:)*this%bhy(:,:)*this%bhz(:,:)

    ! center line elements
    this%dlx(:,:) = this%bhx(:,:)*this%dx
    this%dly(:,:) = this%bhy(:,:)*this%dy
  END SUBROUTINE InitMesh_midpoint


  PURE SUBROUTINE VectorDivergence2D(this,vx,vy,divv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX) &
                      :: vx,vy,divv
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,vx,vy
    INTENT(OUT)       :: divv
    !------------------------------------------------------------------------!
!CDIR NODEP
    DO j=this%JGMIN+1,this%JGMAX-1
!CDIR NODEP
      DO i=this%IGMIN+1,this%IGMAX-1
!CDIR IEXPAND
         divv(i,j) = Divergence3D(this%fhx(i,j,NORTH),this%fhx(i,j,SOUTH), &
                                 this%fhy(i,j,EAST),this%fhy(i,j,WEST), &
                                 this%fhz(i,j,NORTH),this%fhz(i,j,SOUTH), &
                                 this%fhz(i,j,EAST),this%fhz(i,j,WEST), &
                                 0.0,0.0,0.0, &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 0.5*(vx(i+1,j)+vx(i,j)),0.5*(vx(i-1,j)+vx(i,j)), &
                                 0.5*(vy(i,j+1)+vy(i,j)),0.5*(vy(i,j-1)+vy(i,j)), &
                                 0.0,0.0,0.0)
! we cannot use the 2D Divergence, because hz .ne. 1 in general
! for 3D rotationally symmetric meshes
!Divergence2D(this%fhx(i,j,NORTH),this%fhx(i,j,SOUTH), &
!                                this%fhy(i,j,EAST),this%fhy(i,j,WEST), &
!                                0.0,0.0,this%dxdV(i,j),this%dydV(i,j), &
!                                0.5*(vx(i+1,j)+vx(i,j)),0.5*(vx(i-1,j)+vx(i,j)), &
!                                0.5*(vy(i,j+1)+vy(i,j)),0.5*(vy(i,j-1)+vy(i,j)), &
!                                0.0,0.0)
      END DO
    END DO
  END SUBROUTINE VectorDivergence2D


  PURE SUBROUTINE TensorDivergence2D(this,Txx,Txy,Tyx,Tyy,divTx,divTy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX) &
                      :: Txx,Txy,Tyx,Tyy,divTx,divTy
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Txx,Txy,Tyx,Tyy
    INTENT(OUT)       :: divTx,divTy
    !------------------------------------------------------------------------!
!CDIR NODEP
    DO j=this%JGMIN+1,this%JGMAX-1
!CDIR NODEP
      DO i=this%IGMIN+1,this%IGMAX-1
         ! x component of tensor divergence
!CDIR IEXPAND
         divTx(i,j) = Divergence2D(this%fhx(i,j,NORTH),this%fhx(i,j,SOUTH), &
                                 this%fhy(i,j,EAST),this%fhy(i,j,WEST), &
                                 this%cxyx(i,j,1),this%cyxy(i,j,1), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 0.5*(Txx(i+1,j)+Txx(i,j)),0.5*(Txx(i-1,j)+Txx(i,j)), &
                                 0.5*(Txy(i,j+1)+Txy(i,j)),0.5*(Txy(i,j-1)+Txy(i,j)), &
                                 Tyx(i,j),Tyy(i,j))
         ! y component of tensor divergence
!CDIR IEXPAND
         divTy(i,j) = Divergence2D(this%fhx(i,j,NORTH),this%fhx(i,j,SOUTH), &
                                 this%fhy(i,j,EAST),this%fhy(i,j,WEST), &
                                 -this%cxyx(i,j,1),-this%cyxy(i,j,1), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 0.5*(Tyx(i+1,j)+Tyx(i,j)),0.5*(Tyx(i-1,j)+Tyx(i,j)), &
                                 0.5*(Tyy(i,j+1)+Tyy(i,j)),0.5*(Tyy(i,j-1)+Tyy(i,j)), &
                                 Txx(i,j),Txy(i,j))
      END DO
    END DO
  END SUBROUTINE TensorDivergence2D

  ! 3D tensor divergence with coordinate symmetry in z-direction
  PURE SUBROUTINE TensorDivergence3D(this,Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz, &
                                     divTx,divTy,divTz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX) &
                      :: Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz,divTx,divTy,divTz
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz
    INTENT(OUT)       :: divTx,divTy,divTz
    !------------------------------------------------------------------------!
!CDIR NODEP
    DO j=this%JGMIN+1,this%JGMAX-1
!CDIR NODEP
      DO i=this%IGMIN+1,this%IGMAX-1
         ! x component of tensor divergence
!CDIR IEXPAND
         divTx(i,j) = Divergence3D(this%fhx(i,j,NORTH),this%fhx(i,j,SOUTH), &
                                 this%fhy(i,j,EAST),this%fhy(i,j,WEST), &
                                 this%fhz(i,j,NORTH),this%fhz(i,j,SOUTH), &
                                 this%fhz(i,j,EAST),this%fhz(i,j,WEST), &
                                 this%cxyx(i,j,1),this%cyxy(i,j,1),this%czxz(i,j,1), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 0.5*(Txx(i+1,j)+Txx(i,j)),0.5*(Txx(i-1,j)+Txx(i,j)), &
                                 0.5*(Txy(i,j+1)+Txy(i,j)),0.5*(Txy(i,j-1)+Txy(i,j)), &
                                 Tyx(i,j),Tyy(i,j),Tzz(i,j))
         ! y component of tensor divergence
!CDIR IEXPAND
         divTy(i,j) = Divergence3D(this%fhx(i,j,NORTH),this%fhx(i,j,SOUTH), &
                                 this%fhy(i,j,EAST),this%fhy(i,j,WEST), &
                                 this%fhz(i,j,NORTH),this%fhz(i,j,SOUTH), &
                                 this%fhz(i,j,EAST),this%fhz(i,j,WEST), &
                                 -this%cxyx(i,j,1),-this%cyxy(i,j,1),this%czyz(i,j,1), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 0.5*(Tyx(i+1,j)+Tyx(i,j)),0.5*(Tyx(i-1,j)+Tyx(i,j)), &
                                 0.5*(Tyy(i,j+1)+Tyy(i,j)),0.5*(Tyy(i,j-1)+Tyy(i,j)), &
                                 Txx(i,j),Txy(i,j),Tzz(i,j))
!CDIR IEXPAND
         divTz(i,j) = Divergence3D(this%fhx(i,j,NORTH),this%fhx(i,j,SOUTH), &
                                 this%fhy(i,j,EAST),this%fhy(i,j,WEST), &
                                 this%fhz(i,j,NORTH),this%fhz(i,j,SOUTH), &
                                 this%fhz(i,j,EAST),this%fhz(i,j,WEST), &
                                 this%czyz(i,j,1),0.0,this%czxz(i,j,1), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 0.5*(Tzx(i+1,j)+Tzx(i,j)),0.5*(Tzx(i-1,j)+Tzx(i,j)), &
                                 0.5*(Tzy(i,j+1)+Tzy(i,j)),0.5*(Tzy(i,j-1)+Tzy(i,j)), &
                                 Tyz(i,j),0.0,-Txz(i,j))
      END DO
    END DO
  END SUBROUTINE TensorDivergence3D


  ! elemental functions

  ! Computes the divergence of a vector or symmetric tensor, given its components
  ! and some mesh properties.
  !
  ! 2D vector or tensor on 2D mesh
  ELEMENTAL FUNCTION Divergence2D(hxNorth,hxSouth,hyEast,hyWest,cxyx,cyxy,dxdV,dydV, &
                        TxxEast,TxxWest,TxyNorth,TxySouth,TyxCent,TyyCent) RESULT(div)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL,INTENT(IN) :: hxNorth,hxSouth,hyEast,hyWest, & ! scale factors
                       cxyx,cyxy,&                      ! commutator coeffs
                       dxdV,dydV, &                     ! dx/dV and dy/dV
                       TxxEast,TxxWest, &               ! tensor components
                       TxyNorth,TxySouth, &             !   on cell faces
                       TyxCent,TyyCent                  !   central
                       ! to compute the y-component of the tensor divergence,
                       ! one has to modify the input according to
                       !   Txx   <->   Tyx
                       !   Txy   <->   Tyy
                       !   cxyx   ->   -cxyx
                       !   cyxy   ->   -cyxy
    REAL            :: div
    !------------------------------------------------------------------------!
    div = dydV*(hyEast*TxxEast-hyWest*TxxWest) &
        + dxdV*(hxNorth*TxyNorth-hxSouth*TxySouth) &
        + cxyx*TyxCent - cyxy*TyyCent
  END FUNCTION Divergence2D

  ! 3D vector or tensor on 2D mesh with rotational symmetry
  ELEMENTAL FUNCTION Divergence3D(hxNorth,hxSouth,hyEast,hyWest, &
                                  hzNorth,hzSouth,hzEast,hzWest, &
                                  cxyx,cyxy,czxz,dxdV,dydV, &
                                  TxxEast,TxxWest,TxyNorth,TxySouth, &
                                  TyxCent,TyyCent,TzzCent) RESULT(div)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL,INTENT(IN) :: hxNorth,hxSouth,hyEast,hyWest, & ! scale factors
                       hzNorth,hzSouth,hzEast,hzWest, & !
                       cxyx,cyxy,czxz, &                ! commutator coeffs
                       dxdV,dydV, &                     ! dx/dV and dy/dV
                       TxxEast,TxxWest, &               ! tensor components
                       TxyNorth,TxySouth, &             !   on cell faces
                       TyxCent,TyyCent,TzzCent          !   central
                       ! to compute the y-component of the tensor divergence,
                       ! one has to modify the input according to
                       !   Txx   <->   Tyx
                       !   Txy   <->   Tyy
                       !   cxyx   ->   -cxyx
                       !   cyxy   ->   -cyxy
                       !   czxz   ->    czyz
                       ! to compute the z-component of the tensor divergence,
                       ! one has to modify the input according to
                       !   Txx    ->   Tzx
                       !   Txy    ->   Tzy
                       !   Tyx    ->   Tyz
                       !   Tyy    ->   0
                       !   Tzz    ->  -Txz
                       !   cxyx   ->   czyz
                       !   cyxy   ->   0
    REAL            :: div
    !------------------------------------------------------------------------!
    div = dydV*(hyEast*hzEast*TxxEast-hyWest*hzWest*TxxWest) &
        + dxdV*(hxNorth*hzNorth*TxyNorth-hxSouth*hzSouth*TxySouth) &
        + cxyx*TyxCent - cyxy*TyyCent - czxz*TzzCent
  END FUNCTION Divergence3D


  SUBROUTINE CloseMesh_midpoint(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%dAx,this%dAy,this%dAydx,this%dAxdy, &
         this%cyxy,this%cxyx,this%czxz,this%czyz,this%sqrtg)
    ! call basic mesh deconstructor
    CALL CloseMesh(this)
  END SUBROUTINE CloseMesh_midpoint

END MODULE mesh_midpoint
