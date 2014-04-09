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
!> \addtogroup mesh
!! - general parameters of mesh group as key-values
!! \key{type,INTEGER,spatial integration method
!!      (see \link mesh_generic \endlink for currently supported mesh types)}
!! \key{geometry,INTEGER,geometry of the mesh
!!      (see \link geometry_generic \endlink for currently supported geometries)}
!! \key{inum,INTEGER,x-resolution}
!! \key{jnum,INTEGER,y-resolution}
!! \key{xmin,REAL,minimum of x-coordinates}
!! \key{xmax,REAL,maximum of x-coordinates}
!! \key{ymin,REAL,minimum of y-coordinates}
!! \key{ymax,REAL,maximum of y-coordinates}
!! \key{gparam,REAL,1st geometry parameter}
!! \key{gparam2,REAL,2nd geometry parameter}
!! \key{gparam3,REAL,3rd geometry parameter}
!! \key{dz,REAL,extent of 3rd dimension}
!! - enable/disable output of certain mesh arrays
!! \key{output/bary,INTEGER,cell bary center in cartesian coordinates,1}
!! \key{output/bary_curv,INTEGER,cell bary center in curvilinear coordinates,1}
!! \key{output/corners,INTEGER,cell corners in cartesian coordinates,1}
!! \key{output/dl,INTEGER,line elements,0}
!! \key{output/bh,INTEGER,scale factors,0}
!! \key{output/volume,INTEGER,volume elements,0}
!! \key{output/dA,INTEGER,surface elements,0}
!! \key{output/radius,INTEGER,distance to the origin,0}
!! \key{output/position_vector,INTEGER,curvilinear position vector components,0}
!! \key{output/rotation,INTEGER,rotation angle of local orthonormal frames,0}

!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief generic mesh module
!!
!! \ingroup mesh
!----------------------------------------------------------------------------!
MODULE mesh_generic
  USE mesh_midpoint, InitMesh_common => InitMesh, CloseMesh_common => CloseMesh
  USE mesh_trapezoidal
  USE geometry_generic
  USE common_dict

  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE remap_bounds
    MODULE PROCEDURE remap_bounds_0, remap_bounds_2, remap_bounds_3, &
        remap_bounds_4
  END INTERFACE
  !> \endcond
  !> \name Public Attributes
  !! #### mesh types
  INTEGER, PARAMETER :: MIDPOINT     = 1 !< use midpoint rule to approximate flux integrals
  INTEGER, PARAMETER :: TRAPEZOIDAL  = 2 !< use trapezoidal rule to approximate flux integrals
  !> \}
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
       CARTESIAN, POLAR, LOGPOLAR, TANPOLAR, SINHPOLAR, SINHTANHPOLAR, &
       CYLINDRICAL, TANCYLINDRICAL, LNCOSHCYLINDRICAL, SPHERICAL, SINHSPHERICAL, &
       BIANGLESPHERICAL, OBLATE_SPHEROIDAL, CHANNEL, ELLIPTIC, SINHCARTESIAN, &
       ! methods
       InitMesh, &
       Convert2Cartesian, &
       Convert2Curvilinear, &
       Divergence, &
       remap_bounds, &
       InternalPoint, &
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

  !> \public Constructor of generic mesh module
  !!
  !! This subroutine reads the necessary config data for setting up the mesh.
  !! It initializes the geometry and various mesh data arrays. Some of those
  !! are marked for output. For a detailed discription of the various geometries
  !! and how to setup those see \link geometry \endlink .
  SUBROUTINE InitMesh(this,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)          :: this   !< \param [in,out] this all mesh data
    TYPE(Dict_TYP),POINTER  :: config !< \param [in,out] config sub-dictionary 
                                      !! with mesh configuration data
    TYPE(Dict_TYP),POINTER  :: IO     !< \param [in,out] IO mesh specific i/o data
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: xres,yres
    INTEGER           :: meshtype
    INTEGER           :: i,j,err
    REAL              :: cart_max
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL RequireKey(config, "meshtype")
    CALL GetAttr(config, "meshtype", meshtype)
    
    SELECT CASE(meshtype)
    CASE(MIDPOINT)
       CALL InitMesh_midpoint(this,config)
    CASE(TRAPEZOIDAL)
       CALL InitMesh_trapezoidal(this,config)
    CASE DEFAULT
       CALL Error(this,"InitMesh", "Unknown mesh type.")
    END SELECT

!> \todo: insert passage in mesh_common.f90 (ll. 264), problem: mixed dimensions in
!! bianglespherical for curvilinear and cartesian coordinates -> no generality
!! see also geometry_generic
    ! quick hack for 3D cartesian coordinates
    IF (GetType(this%geometry).EQ.BIANGLESPHERICAL) THEN
      DEALLOCATE(this%bccart,this%ccart,this%fccart,this%posvec,this%bposvec, &
         this%fposvec)
      ALLOCATE(this%bccart(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,3), &
         this%ccart(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,3), &
         this%fccart(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,3), &
         this%posvec(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,3), &
         this%bposvec(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,3), &
         this%fposvec(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,3), &
         STAT=err)
      IF (err.NE.0) &
         CALL Error(this,"InitMesh_generic","Unable to allocate memory!")
    END IF
    
    ! compute cartesian coordinates for bary center values
    CALL Convert2Cartesian(this%geometry,this%bcenter,this%bccart)
    ! compute cartesian coordinates for corner and face positions
    CALL Convert2Cartesian(this%geometry,this%cpos,this%ccart)
    ! compute cartesian coordinates for face centered positions
    CALL Convert2Cartesian(this%geometry,this%fpos,this%fccart)

    ! compute geo. centered radii
    CALL Radius(this%geometry,this%center,this%radius)
    ! compute bary centered radii
    CALL Radius(this%geometry,this%bcenter,this%bradius)
    ! compute face centered radii
    CALL Radius(this%geometry,this%fpos,this%fradius)

    ! compute position vector components with respect to local orthonormal frame
    ! at geometric centers
    CALL PositionVector(this%geometry,this%center,this%posvec)
    ! at bary centers
    CALL PositionVector(this%geometry,this%bcenter,this%bposvec)
    ! at face centers
    CALL PositionVector(this%geometry,this%fpos,this%fposvec)

    ! print some information
    CALL Info(this, " MESH-----> quadrature rule:   " // TRIM(GetName(this)))
    WRITE (xres, '(I0)') this%INUM    ! this is just for better looking output
    WRITE (yres, '(I0)') this%JNUM
    CALL Info(this, "            resolution:        " // TRIM(xres) // " x " // TRIM(yres))
    WRITE (xres, '(ES9.2,A,ES9.2)') this%xmin, " ..", this%xmax
    WRITE (yres, '(ES9.2,A,ES9.2)') this%ymin, " ..", this%ymax
    CALL Info(this, "            computat. domain:  x=" // TRIM(xres))
    CALL Info(this, "                               y="  // TRIM(yres))
#ifdef PARALLEL
    WRITE (xres, '(I0)') this%dims(1)
    WRITE (yres, '(I0)') this%dims(2)
    CALL Info(this, "            MPI partition:     " // TRIM(xres) // " x " // TRIM(yres))
#endif

    ! \todo implement a correct check of the 
    ! curvilinear range (e.g. xi > 0 in polar)


    CALL SetOutput(this,config,IO)

  END SUBROUTINE InitMesh

  !> \private Setup mesh fields for i/o
  SUBROUTINE SetOutput(this,config,IO)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    TYPE(Dict_TYP),POINTER  :: config,IO
    !------------------------------------------------------------------------!
    TYPE(Dict_TYP),POINTER  :: bary,corners,node
    INTEGER                 :: status
    REAL,DIMENSION(:,:,:),POINTER :: cart,curv
    INTEGER                 :: writefields
    !------------------------------------------------------------------------!
    INTENT(INOUT)        :: this
    !------------------------------------------------------------------------! 
    !Set OutputDict
    CALL RequireKey(config, "output/corners", 1)
    CALL GetAttr(config, "output/corners", writefields)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%ccart)) THEN
        CALL AddField(IO,"corners",this%ccart, &
                      Dict("name" / "coordinates:corners"))
    END IF

    CALL RequireKey(config, "output/bary_curv", 1)
    CALL GetAttr(config, "output/bary_curv", writefields)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%bcenter)) THEN
        CALL AddField(IO,"bary_curv",this%bcenter, &
                      Dict("name" / "coordinates:bary_curv"))
    END IF

    CALL RequireKey(config, "output/bary", 1)
    CALL GetAttr(config, "output/bary", writefields)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%bccart)) THEN
        CALL AddField(IO,"bary_centers",this%bccart, &
                      Dict("name" / "coordinates:bary"))
    END IF

    CALL RequireKey(config, "output/rotation", 0)
    CALL GetAttr(config, "output/rotation", writefields)
    IF(writefields.EQ.1) THEN
        ALLOCATE(cart(this%IGMIN:this%IGMAX, this%JGMIN:this%JGMAX, 2), &
                 curv(this%IGMIN:this%IGMAX, this%JGMIN:this%JGMAX, 2), &
                 this%rotation(this%IGMIN:this%IGMAX, this%JGMIN:this%JGMAX), &
                 STAT=status)
        IF(status.NE.0) &
            CALL Error(this, "SetOutput", "Couldn't allocate memory")
        ! initialize with vectors pointing in x-direction
        cart(:,:,1) = 1.
        cart(:,:,2) = 0.
        CALL Convert2Curvilinear(this%geometry, this%bcenter, cart, curv)
        this%rotation = ATAN2(curv(:,:,2),curv(:,:,1))
        DEALLOCATE(cart, curv)
        CALL AddField(IO,"rotation",this%rotation,Dict("name" / "rotation"))
    END IF

    CALL RequireKey(config, "output/dl", 0)
    CALL GetAttr(config, "output/dl", writefields)
    IF(writefields.EQ.1) THEN
        IF (ASSOCIATED(this%dlx)) &
           CALL AddField(IO,"dlx",this%dlx,Dict("name" / "dlx"))
        IF (ASSOCIATED(this%dly)) &
           CALL AddField(IO,"dly",this%dly,Dict("name" / "dly"))
    END IF

    CALL RequireKey(config, "output/bh", 0)
    CALL GetAttr(config, "output/bh", writefields)
    IF(writefields.EQ.1) THEN
        IF (ASSOCIATED(this%bhx)) &
           CALL AddField(IO,"bhx",this%bhx,Dict("name" / "bhx"))
        IF (ASSOCIATED(this%bhy)) &
           CALL AddField(IO,"bhy",this%bhy,Dict("name" / "bhy"))
        IF (ASSOCIATED(this%bhz)) &
           CALL AddField(IO,"bhz",this%bhy,Dict("name" / "bhz"))
    END IF

    CALL RequireKey(config, "output/volume", 0)
    CALL GetAttr(config, "output/volume", writefields)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%volume)) &
        CALL AddField(IO,"volume",this%volume,Dict("name" / "volume"))
    
    CALL RequireKey(config, "output/dA", 0)
    CALL GetAttr(config, "output/dA", writefields)
    IF(writefields.EQ.1) THEN
        IF (ASSOCIATED(this%dAx)) &
           CALL AddField(IO,"dAx",this%dAx,Dict("name" / "dAx"))
        IF (ASSOCIATED(this%dAy)) &
           CALL AddField(IO,"dAy",this%dAy,Dict("name" / "dAy"))
        IF (ASSOCIATED(this%dAz)) &
           CALL AddField(IO,"dAz",this%dAz,Dict("name" / "dAz"))
    END IF

    CALL RequireKey(config, "output/radius", 0)
    CALL GetAttr(config, "output/radius", writefields)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%radius)) &
        CALL AddField(IO,"radius",this%radius,Dict("name" / "radius"))

    CALL RequireKey(config, "output/position_vector", 0)
    CALL GetAttr(config, "output/position_vector", writefields)
    IF((writefields.EQ.1).AND.ASSOCIATED(this%posvec)) &
        CALL AddField(IO,"posvec",this%posvec,Dict("name" / "position_vector"))

  END SUBROUTINE SetOutput

  FUNCTION remap_bounds_0(array) RESULT(ptr)
    IMPLICIT NONE 
    !------------------------------------------------------------------------!
    REAL, DIMENSION(1:,1:), TARGET &
                    :: array
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:), POINTER &
                    :: ptr
    !------------------------------------------------------------------------!
    INTENT(IN)      :: array
    !------------------------------------------------------------------------!
    ptr => array 
  END FUNCTION
  
  FUNCTION remap_bounds_2(this,array) RESULT(ptr)
    IMPLICIT NONE 
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)  :: this 
    REAL, DIMENSION(this%IGMIN:,this%JGMIN:), TARGET &
                    :: array
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:), POINTER &
                    :: ptr
    !------------------------------------------------------------------------!
    INTENT(IN)      :: array
    !------------------------------------------------------------------------!
    ptr => array 
  END FUNCTION

  FUNCTION remap_bounds_3(this,array) RESULT(ptr)
    IMPLICIT NONE 
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)  :: this 
    REAL, DIMENSION(this%IGMIN:,this%JGMIN:,:), TARGET &
                    :: array
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:), POINTER &
                    :: ptr
    !------------------------------------------------------------------------!
    INTENT(IN)      :: array
    !------------------------------------------------------------------------!
    ptr => array 
  END FUNCTION remap_bounds_3

  FUNCTION remap_bounds_4(this,array) RESULT(ptr)
    IMPLICIT NONE 
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)  :: this 
    REAL, DIMENSION(this%IGMIN:,this%JGMIN:,:,:), TARGET &
                    :: array
    !------------------------------------------------------------------------!
    REAL, DIMENSION(:,:,:,:), POINTER &
                    :: ptr
    !------------------------------------------------------------------------!
    INTENT(IN)      :: array
    !------------------------------------------------------------------------!
    ptr => array 
  END FUNCTION remap_bounds_4


  !> \public Check if a given coordinate pair represents an internal point
  !!
  !! This function checks if a point given by its curvilinear coordinates (x,y)
  !! lies inside the computational domain or (if mask is given) inside a
  !! rectangular (with respect to the curvilinear mesh) region specified by an
  !! index mask. 
  PURE FUNCTION InternalPoint(this,x,y,mask) RESULT(ip)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    INTEGER, DIMENSION(4), OPTIONAL :: mask
    REAL        :: x,y
    LOGICAL     :: ip
    !------------------------------------------------------------------------!
    INTEGER     :: imin,imax,jmin,jmax
    !------------------------------------------------------------------------!    
    INTENT(IN)  :: this,x,y,mask
    !------------------------------------------------------------------------!
    ip = .FALSE.
    ! compare with the curvilinear extend (curvilinear domain is always rectangular)
    !
    ! the function behaves in 2 different ways: global and local
    !
    !   1. if mask is given, we do a local comparison, i.e.
    !      we check if the point lies inside the specified rectangular
    !      domain given by its cell indices: imin=mask(1),imax=mask(2),...
    !   2. otherwise we do a global check, i.e. with the whole computational
    !      domain
    IF (PRESENT(mask)) THEN
       ! restrict indices to the actual domain,
       ! which is different for each MPI process in parallel mode
       imin = MAX(this%IGMIN,mask(1))
       imax = MIN(this%IGMAX,mask(2))
       jmin = MAX(this%JGMIN,mask(3))
       jmax = MIN(this%JGMAX,mask(4))
       ! first check if the masked region is at least a subdomain of the actual domain
       IF (((imin.LE.this%IGMAX).AND.(imax.GE.imin)).AND. &
           ((jmin.LE.this%JGMAX).AND.(jmax.GE.jmin))) THEN
          ! compare the curvilinear coordinates at the boundaries of the masked region
          ! with the transformed coordinates of the given point
          IF ((x.GE.this%fpos(imin,jmin,1,1).AND.x.LE.this%fpos(imax,jmax,2,1)).AND. &
              (y.GE.this%fpos(imin,jmin,1,2).AND.x.LE.this%fpos(imax,jmax,2,2))) THEN
             ip = .TRUE.
          END IF
       END IF
    ELSE
       ! do the global check
       IF (((x.GE.this%xmin).AND.(x.LE.this%xmax)).AND. &
           ((y.GE.this%ymin).AND.(y.LE.this%ymax))) THEN
          ip = .TRUE.
       END IF
    END IF
  END FUNCTION InternalPoint


  !> \public Destructor of generic mesh module
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
