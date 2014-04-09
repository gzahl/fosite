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
!> \author Tobias Illenseer
!!
!! \brief mesh module for midpoint quadrature rule
!!
!! \extends mesh_common
!! \ingroup mesh
!----------------------------------------------------------------------------!
MODULE mesh_midpoint
  USE mesh_common, InitMesh_common => InitMesh
  USE boundary_common, ONLY : EAST,WEST,NORTH,SOUTH
  USE geometry_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE Divergence
     MODULE PROCEDURE VectorDivergence2D_1, VectorDivergence2D_2, &
                      TensorDivergence2D_1, TensorDivergence2D_2, &
                      TensorDivergence3D
  END INTERFACE 
  !> \endcond
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

  SUBROUTINE InitMesh(this,config,meshname)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)          :: this
    TYPE(Dict_TYP),POINTER  :: config
    CHARACTER(LEN=32)       :: meshname
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: this
    INTENT(IN)              :: meshname
    !------------------------------------------------------------------------!
    ! basic mesh initialization
    CALL InitMesh_common(this,meshname,config)

    ! initialize geometry
    CALL InitGeometry(this%geometry,config)
  END SUBROUTINE InitMesh


  SUBROUTINE InitMesh_midpoint(this,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)          :: this
    TYPE(Dict_TYP),POINTER  :: config
    !------------------------------------------------------------------------!
    INTEGER                 :: err
    !------------------------------------------------------------------------!
    INTENT(INOUT)           :: this
    !------------------------------------------------------------------------!
    
    ! basic mesh and geometry initialization
    CALL InitMesh(this,config, mesh_name)
    CALL GetAttr(config, "dz", this%dz)

    ! allocate memory for pointers that are specific for midpoint fluxes
    ALLOCATE(this%dAx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%dAy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%dAxdy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%dAydx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%cyxy(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%cxyx(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%czxz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%czyz(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,1), &
         this%sqrtg(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%invsqrtg(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         STAT=err)
    IF (err.NE.0) THEN
       CALL Error(this,"InitMesh_midpoint", "Unable to allocate memory.")
    END IF

    ! get geometrical scale factors
    CALL ScaleFactors(this%geometry,this%fpos,this%fhx,this%fhy,this%fhz)   ! faces
    CALL ScaleFactors(this%geometry,this%center,this%bhx,this%bhy,this%bhz) ! centers

    ! surface elements divided by dx or dy
    this%dAxdy(:,:,1:2) = this%fhz(:,:,1:2)*this%fhy(:,:,1:2)*this%dz ! perpendicular to x-direction
    this%dAydx(:,:,1:2) = this%fhz(:,:,3:4)*this%fhx(:,:,3:4)*this%dz ! perpendicular to y-direction

    ! surface elements
    this%dAx(:,:,:) = this%dAxdy(:,:,:)*this%dy                 ! perpendicular to x-direction
    this%dAy(:,:,:) = this%dAydx(:,:,:)*this%dx                 ! perpendicular to y-direction

    ! square root of determinant of metric and its inverse 
    this%sqrtg(:,:) = this%bhx(:,:)*this%bhy(:,:)*this%bhz(:,:)
    this%invsqrtg(:,:) = 1.0 / (this%sqrtg(:,:)+TINY(1.0))

    ! volume elements
    this%volume(:,:) = this%sqrtg(:,:)*this%dx*this%dy*this%dz  ! = hx*hy*hz*dx*dy*dz

    ! inverse volume elements multiplied by dx or dy
    this%dxdV(:,:) = this%invsqrtg(:,:) / (this%dy*this%dz)     ! = dx / dV
    this%dydV(:,:) = this%invsqrtg(:,:) / (this%dx*this%dz)     ! = dy / dV

    ! cell bary centers
    this%bcenter(:,:,:)  = this%center(:,:,:)

    ! commutator coefficients 
    this%cyxy(:,:,1) = 0.5*(this%fhz(:,:,2)+this%fhz(:,:,1)) &
         * (this%fhy(:,:,2)-this%fhy(:,:,1)) * this%dydV(:,:)*this%dz
    this%cxyx(:,:,1) = 0.5*(this%fhz(:,:,4)+this%fhz(:,:,3)) &
         * (this%fhx(:,:,4)-this%fhx(:,:,3)) * this%dxdV(:,:)*this%dz
    this%czxz(:,:,1) = 0.5*(this%fhy(:,:,2)+this%fhy(:,:,1)) &
         * (this%fhz(:,:,2)-this%fhz(:,:,1)) * this%dydV(:,:)*this%dz
    this%czyz(:,:,1) = 0.5*(this%fhx(:,:,4)+this%fhx(:,:,3)) &
         * (this%fhz(:,:,4)-this%fhz(:,:,3)) * this%dxdV(:,:)*this%dz

    ! center line elements
    this%dlx(:,:) = this%bhx(:,:)*this%dx
    this%dly(:,:) = this%bhy(:,:)*this%dy
  END SUBROUTINE InitMesh_midpoint


  ! computes the cell centered curvilinear vector divergence
  ! for cell centered 2D vector components vx,vy on the whole mesh
  ! except for the outermost boundary cells
  PURE SUBROUTINE VectorDivergence2D_1(this,vx,vy,divv)
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
    ! we simply use the 3D tensor divergence and set the commutator coefficients
    ! and the off-diagonal tensor components to 0
!CDIR NODEP
    DO j=this%JGMIN+1,this%JGMAX-1
!CDIR NODEP
      DO i=this%IGMIN+1,this%IGMAX-1
!CDIR IEXPAND
         divv(i,j) = Divergence3D(this%dAxdy(i,j,1),this%dAxdy(i,j,2), &
                                 this%dAydx(i,j,1),this%dAydx(i,j,2), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 0.0,0.0,0.0, & ! no commutator coefficients
                                 0.5*(vx(i-1,j)+vx(i,j)),0.5*(vx(i+1,j)+vx(i,j)), &
                                 0.5*(vy(i,j-1)+vy(i,j)),0.5*(vy(i,j+1)+vy(i,j)), &
                                 0.0,0.0,0.0)   ! no off-diagonal tensor components
      END DO
    END DO
  END SUBROUTINE VectorDivergence2D_1


  ! computes the cell centered curvilinear vector divergence
  ! for 2D vector v given on the 4 face centered positions
  PURE SUBROUTINE VectorDivergence2D_2(this,v,divv)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,2) &
                      :: v
    REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX) &
                      :: divv
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,v
    INTENT(OUT)       :: divv
    !------------------------------------------------------------------------!
    ! we simply use the 3D tensor divergence and set all commutator coefficients
    ! and the tensor components Tyx, Tyy, Tzz to zero
!CDIR NODEP
    DO j=this%JGMIN,this%JGMAX
!CDIR NODEP
      DO i=this%IGMIN,this%IGMAX
!CDIR IEXPAND
         divv(i,j) = Divergence3D(this%dAxdy(i,j,1),this%dAxdy(i,j,2), &
                                 this%dAydx(i,j,1),this%dAydx(i,j,2), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 0.0,0.0,0.0, & ! vanishing commutator coefficients
                                 v(i,j,WEST,1),v(i,j,EAST,1), &
                                 v(i,j,SOUTH,2),v(i,j,NORTH,2), &
                                 0.0,0.0,0.0)   ! vanishing tensor components
      END DO
    END DO
  END SUBROUTINE VectorDivergence2D_2


  ! computes the cell centered curvilinear tensor divergence on the whole mesh
  ! except for the outermost boundary cells
  ! input: 2D rank 2 tensor T with components Txx,Txy,Tyx,Tyy given at cell centers
  ! output: 2D vector vector components divTx,divTy
  PURE SUBROUTINE TensorDivergence2D_1(this,Txx,Txy,Tyx,Tyy,divTx,divTy)
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
    ! we simply use the 3D tensor divergence and set the commutator coefficient
    ! and the tensor component related to the z-direction to zero
!CDIR NODEP
    DO j=this%JGMIN+1,this%JGMAX-1
!CDIR NODEP
      DO i=this%IGMIN+1,this%IGMAX-1
         ! x component of tensor divergence
!CDIR IEXPAND
         divTx(i,j) = Divergence3D(this%dAxdy(i,j,1),this%dAxdy(i,j,2), &
                                 this%dAydx(i,j,1),this%dAydx(i,j,2), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 this%cxyx(i,j,1),this%cyxy(i,j,1), &
                                 0.0, & ! czxz = 0 because of 2D
                                 0.5*(Txx(i-1,j)+Txx(i,j)),0.5*(Txx(i+1,j)+Txx(i,j)), &
                                 0.5*(Txy(i,j-1)+Txy(i,j)),0.5*(Txy(i,j+1)+Txy(i,j)), &
                                 Tyx(i,j),Tyy(i,j),0.0) ! Tzz = 0 because of 2D
         ! y component of tensor divergence
!CDIR IEXPAND
         divTy(i,j) = Divergence3D(this%dAxdy(i,j,1),this%dAxdy(i,j,2), &
                                 this%dAydx(i,j,1),this%dAydx(i,j,2), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 -this%cxyx(i,j,1),-this%cyxy(i,j,1), &
                                 0.0, & ! czyz = 0 because of 2D
                                 0.5*(Tyx(i-1,j)+Tyx(i,j)),0.5*(Tyx(i+1,j)+Tyx(i,j)), &
                                 0.5*(Tyy(i,j-1)+Tyy(i,j)),0.5*(Tyy(i,j+1)+Tyy(i,j)), &
                                 Txx(i,j),Txy(i,j),0.0) ! Tzz = 0 because of 2D
      END DO
    END DO
  END SUBROUTINE TensorDivergence2D_1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ATTENTION: TensorDivergence2D_2 is untested, use with care
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! computes the cell centered curvilinear tensor divergence
  ! input: 2D rank 2 tensor T given on the 4 face centered positions
  ! output: 2D vector divT
  PURE SUBROUTINE TensorDivergence2D_2(this,T,divT)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,2,2) :: T
    REAL, DIMENSION(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2)   :: divT
    !------------------------------------------------------------------------!
    REAL              :: Txx,Txy,Tyx,Tyy
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this,T
    INTENT(OUT)       :: divT
    !------------------------------------------------------------------------!
    ! we simply use the 3D tensor divergence and set the commutator coefficient
    ! and the tensor component related to the z-direction to zero
!CDIR NODEP
    DO j=this%JGMIN,this%JGMAX
!CDIR NODEP
      DO i=this%IGMIN,this%IGMAX
         ! compute mean value of all four face values to obtain cell centered values
         Txx = 0.25*SUM(T(i,j,:,1,1))
         Txy = 0.25*SUM(T(i,j,:,1,2))
         Tyx = 0.25*SUM(T(i,j,:,2,1))
         Tyy = 0.25*SUM(T(i,j,:,2,2))         
         ! x component of tensor divergence
!CDIR IEXPAND
         divT(i,j,1) = Divergence3D(this%dAxdy(i,j,1),this%dAxdy(i,j,2), &
                                 this%dAydx(i,j,1),this%dAydx(i,j,2), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 this%cxyx(i,j,1),this%cyxy(i,j,1), &
                                 0.0, & ! czxz = 0 because of 2D
                                 T(i,j,WEST,1,1),T(i,j,EAST,1,1), &
                                 T(i,j,SOUTH,1,2),T(i,j,NORTH,1,2), &
                                 Tyx,Tyy,0.0) ! Tzz = 0 because of 2D
         ! y component of tensor divergence
!CDIR IEXPAND
         divT(i,j,2) = Divergence3D(this%dAxdy(i,j,1),this%dAxdy(i,j,2), &
                                 this%dAydx(i,j,1),this%dAydx(i,j,2), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 -this%cxyx(i,j,1),-this%cyxy(i,j,1), &
                                 0.0, & ! czyz = 0 because of 2D
                                 T(i,j,WEST,2,1),T(i,j,EAST,2,1), &
                                 T(i,j,SOUTH,2,2),T(i,j,NORTH,2,2), &
                                 Txx,Txy,0.0) ! Tzz = 0 because of 2D
      END DO
    END DO
  END SUBROUTINE TensorDivergence2D_2


  ! computes the cell centered tensor divergence in curvilinear coordinates with
  ! symmetry along the z-direction on the whole mesh except for the outermost 
  ! boundary cells
  ! input: 3D rank 2 tensor T with components Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz
  !        given at cell centers
  ! output: 3D vector vector components divTx,divTy,divTz
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
         divTx(i,j) = Divergence3D(this%dAxdy(i,j,1),this%dAxdy(i,j,2), &
                                 this%dAydx(i,j,1),this%dAydx(i,j,2), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 this%cxyx(i,j,1),this%cyxy(i,j,1),this%czxz(i,j,1), &
                                 0.5*(Txx(i-1,j)+Txx(i,j)),0.5*(Txx(i+1,j)+Txx(i,j)), &
                                 0.5*(Txy(i,j-1)+Txy(i,j)),0.5*(Txy(i,j+1)+Txy(i,j)), &
                                 Tyx(i,j),Tyy(i,j),Tzz(i,j))
         ! y component of tensor divergence
!CDIR IEXPAND
         divTy(i,j) = Divergence3D(this%dAxdy(i,j,1),this%dAxdy(i,j,2), &
                                 this%dAydx(i,j,1),this%dAydx(i,j,2), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 -this%cxyx(i,j,1),-this%cyxy(i,j,1),this%czyz(i,j,1), &
                                 0.5*(Tyx(i-1,j)+Tyx(i,j)),0.5*(Tyx(i+1,j)+Tyx(i,j)), &
                                 0.5*(Tyy(i,j-1)+Tyy(i,j)),0.5*(Tyy(i,j+1)+Tyy(i,j)), &
                                 Txx(i,j),Txy(i,j),Tzz(i,j))
         ! z component of tensor divergence
!CDIR IEXPAND
         divTz(i,j) = Divergence3D(this%dAxdy(i,j,1),this%dAxdy(i,j,2), &
                                 this%dAydx(i,j,1),this%dAydx(i,j,2), &
                                 this%dxdV(i,j),this%dydV(i,j), &
                                 this%czyz(i,j,1),0.0,this%czxz(i,j,1), &
                                 0.5*(Tzx(i-1,j)+Tzx(i,j)),0.5*(Tzx(i+1,j)+Tzx(i,j)), &
                                 0.5*(Tzy(i,j-1)+Tzy(i,j)),0.5*(Tzy(i,j+1)+Tzy(i,j)), &
                                 Tyz(i,j),0.0,-Txz(i,j))
      END DO
    END DO
  END SUBROUTINE TensorDivergence3D


  ! the elemental workhorse: basic function for computing 2D/3D vector or tensor divergences
  ! works with planar geometry as well (czxz = 0)
  ! input: area and volume elements (multiplied and devided by dx or dy,
  !        commutator coefficients, tensor components
  ! output: vector divergence (scalar) or x-component of the tensor divergence
  !         call this function with different input to obtain other components
  ELEMENTAL FUNCTION Divergence3D(dAxdyWest,dAxdyEast,dAydxSouth,dAydxNorth, &
                                  dxdV,dydV,cxyx,cyxy,czxz, &
                                  TxxWest,TxxEast,TxySouth,TxyNorth, &
                                  TyxCent,TyyCent,TzzCent) RESULT(div)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL,INTENT(IN) :: dAxdyWest,dAxdyEast, &           ! surface elements
                       dAydxSouth,dAydxNorth, &         !   dAx/dy, dAy/dx
                       dxdV,dydV, &                     ! dx/dV and dy/dV
                       cxyx,cyxy,czxz, &                ! commutator coeffs
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
    !------------------------------------------------------------------------!
    REAL            :: div
    !------------------------------------------------------------------------!
    div = dydV*(dAxdyEast*TxxEast-dAxdyWest*TxxWest) &
        + dxdV*(dAydxNorth*TxyNorth-dAydxSouth*TxySouth) &
        + cxyx*TyxCent - cyxy*TyyCent - czxz*TzzCent
  END FUNCTION Divergence3D


  SUBROUTINE CloseMesh_midpoint(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%dAx,this%dAy,this%dAydx,this%dAxdy, &
         this%cyxy,this%cxyx,this%czxz,this%czyz,this%sqrtg,this%invsqrtg)
    ! call basic mesh deconstructor
    CALL CloseMesh(this)
  END SUBROUTINE CloseMesh_midpoint

END MODULE mesh_midpoint
