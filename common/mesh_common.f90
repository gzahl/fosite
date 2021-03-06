!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: mesh_common.f90                                                   #
!#                                                                           #
!# Copyright (C) 2006-2014                                                   #
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
!> \defgroup mesh mesh
!! \{
!! \brief Family of mesh modules
!! \}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!!
!! \brief basic mesh module
!!
!! \extends common_types
!! \ingroup mesh
!----------------------------------------------------------------------------!
MODULE mesh_common
  USE common_types, &
       GetType_common => GetType, GetName_common => GetName, &
       GetRank_common => GetRank, GetNumProcs_common => GetNumProcs, &
       Initialized_common => Initialized, Info_common => Info, &
       Warning_common => Warning, Error_common => Error
  USE geometry_common, ONLY : Geometry_TYP
  USE common_dict
#ifdef PARALLEL
#ifdef HAVE_MPI_MOD
  USE mpi
#endif
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef HAVE_MPIF_H
  include 'mpif.h'
#endif
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  !> \name Public Attributes
  INTEGER, PARAMETER :: NDIMS = 2          ! dimensions of cartesian topology
  INTEGER, PARAMETER :: VECLEN = &         ! vector length ..
#if defined(NECSX8) || defined(NECSX9)
  256                                      ! .. of NEC SX8/SX9 CPUs 
#else
  1                                        ! .. of everthing else 
#endif
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE GetType
     MODULE PROCEDURE GetMesh, GetType_common
  END INTERFACE
  INTERFACE GetName
     MODULE PROCEDURE GetMeshName, GetName_common
  END INTERFACE
  INTERFACE GetRank
     MODULE PROCEDURE GetMeshRank, GetRank_common
  END INTERFACE
  INTERFACE GetNumProcs
     MODULE PROCEDURE GetMeshNumProcs, GetNumProcs_common
  END INTERFACE
  INTERFACE Initialized
     MODULE PROCEDURE MeshInitialized, Initialized_common
  END INTERFACE
  INTERFACE Info
     MODULE PROCEDURE MeshInfo, Info_common
  END INTERFACE
  INTERFACE Warning
     MODULE PROCEDURE MeshWarning, Warning_common
  END INTERFACE
  INTERFACE Error
     MODULE PROCEDURE MeshError_rank0, MeshError_rankX, Error_common
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  !> index selection type
  TYPE Selection_TYP
     !> \name Variables
     INTEGER           :: IMIN,IMAX       !< selection in x-direction
     INTEGER           :: JMIN,JMAX       !< selection in y-direction
     LOGICAL, POINTER  :: mask(:,:)       !< optional selection mask
  END type Selection_TYP
  !> mesh data structure
  TYPE Mesh_TYP
     !> \name Variables
     TYPE(Common_TYP)  :: mtype           !< mesh type
     TYPE(Geometry_TYP):: Geometry        !< geometrical properties
     INTEGER           :: GNUM            !< number of ghost cells
     INTEGER           :: INUM,JNUM       !< resolution
     INTEGER           :: IMIN,IMAX       !< minimal & maximal index in x-direction
     INTEGER           :: JMIN,JMAX       !< minimal & maximal index in y-direction
     INTEGER           :: IGMIN,IGMAX     !< minimal & maximal index in x-direction with ghost cells
     INTEGER           :: JGMIN,JGMAX     !< minimal & maximal index in y-direction with ghost cells
     REAL              :: xmin, xmax      !< spatial extent of computational domain in x-direction
     REAL              :: ymin, ymax      !< spatial extent of computational domain in y-direction
     REAL              :: dx,dy,dz        !< curvilinear spatial differences 
     REAL              :: invdx,invdy     !< inverse of curvilinear spatial differences 
     !> \name
     !! #### cell coordinates
     REAL, DIMENSION(:,:,:), POINTER :: &
                          center, &       !< geometrical centers
                          bcenter, &      !< bary centers
                          bccart          !< cartesian bary centers
     REAL, DIMENSION(:,:,:,:), POINTER :: &
                          fccart, &       !< cartesian face centered positions
                          ccart, &        !< cartesian corner positions
                          fpos, &         !< curvilinear face centered positions
                          cpos            !< curvilinear corner positions
     !> \name
     !! #### line, area and volume elements
     REAL, DIMENSION(:,:), POINTER :: &
                          dlx,dly, &      !< cell centered line elements
                          volume, &       !< cell volumes
                          dxdV,dydV       !< dx/volume and dy/volume
     REAL, DIMENSION(:,:,:), POINTER :: &
                          dAx,dAy,dAz, &  !< cell surface area elements
                          dAxdy,dAydx     !< dAx/dy and dAy/dx
     !> \name
     !! #### scale factors and commutator coefficients
     REAL, DIMENSION(:,:), POINTER :: &
                          bhx,bhy,bhz     !< scale factors at bary centers
     REAL, DIMENSION(:,:,:), POINTER :: &
                          fhx,fhy,fhz, &  !< scale factors at face centers
                          chx,chy,chz, &  !< scale factors at cell corners
                          cyxy,cxyx,czxz,czyz !< commutator coefficients
     !> \name
     !! #### radius and curvilinear position vector components
     REAL, DIMENSION(:,:), POINTER :: &
                          radius, &       !< geometrically centered radius
                          bradius         !< bary centered radius
     REAL, DIMENSION(:,:,:), POINTER :: &
                          fradius, &      !< face centered radius
                          posvec, &       !< geometrically centered curvilinear position vector components
                          bposvec         !< bary centered curvilinear position vector components
     REAL, DIMENSION(:,:,:,:), POINTER :: &
                          fposvec         !< face centered curvilinear position vector components
     !> \name
     !! #### other geometrial quantities
     REAL, DIMENSION(:,:), POINTER :: &
                          sqrtg, &        !< square root of (det(g_ij))
                          invsqrtg, &     !< 1/sqrtg
                          rotation        !< rotation angle of local curvilinear orthonormal frames
     REAL, DIMENSION(:,:,:,:), POINTER :: &
                          weights         !< interpolation weights
#ifdef PARALLEL
     !> \name Variables in Parallel Mode
     INTEGER                    :: MAXINUM,MAXJNUM !< max. of local INUM,JNUM
     INTEGER                    :: MININUM,MINJNUM !< min. of local INUM,JNUM
     INTEGER                    :: comm_cart       !< cartesian communicator
     INTEGER                    :: Icomm,Jcomm     !< communicators for cartesian rows and cols
     INTEGER, DIMENSION(4)      :: comm_boundaries !< communicators for boundary processes
     INTEGER, DIMENSION(4)      :: rank0_boundaries!< map rank0 -> world rank
     INTEGER, DIMENSION(4)      :: neighbor        !< ranks of neighbor proc.
     INTEGER, DIMENSION(NDIMS)  :: dims            !< dimensions of cart comm
     INTEGER, DIMENSION(NDIMS)  :: mycoords        !< par. proc coordinates
#endif
  END TYPE Mesh_TYP
  !> \}
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Mesh_TYP, &
       Selection_TYP, &
#ifdef PARALLEL
       DEFAULT_MPI_REAL, &
#endif
       ! methods
       InitMesh, &
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

  !> \public Constructor of common mesh class
  SUBROUTINE InitMesh(this,mname,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)          :: this
    TYPE(Dict_TYP),POINTER  :: config
    INTEGER                 :: mtype
    CHARACTER(LEN=*)        :: mname
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
#ifdef PARALLEL
    INTEGER           :: inum, jnum
#endif    
    INTEGER           :: err
    REAL, DIMENSION(4,2) :: cfaces,ccorners
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL RequireKey(config, "meshtype")
    CALL GetAttr(config, "meshtype", mtype)

    CALL InitCommon(this%mtype,mtype,mname)

    ! total resolution
    CALL RequireKey(config, "inum")
    CALL RequireKey(config, "jnum")
    CALL GetAttr(config, "inum", this%inum)
    CALL GetAttr(config, "jnum", this%jnum)

    ! number of ghost rows/columns
    this%GNUM = 2

    ! coordinate domain
    CALL RequireKey(config, "xmin")
    CALL RequireKey(config, "xmax")
    CALL RequireKey(config, "ymin")
    CALL RequireKey(config, "ymax")
    CALL GetAttr(config, "xmin", this%xmin)
    CALL GetAttr(config, "xmax", this%xmax)
    CALL GetAttr(config, "ymin", this%ymin)
    CALL GetAttr(config, "ymax", this%ymax)

    ! coordinate differences in each direction
    this%dx = (this%xmax - this%xmin) / this%inum
    this%dy = (this%ymax - this%ymin) / this%jnum
    ! normaly 1 in case of rotational sym. 2*PI 
    this%dz = 1.0
    
    ! inverse coordinate differences
    this%invdx = 1./this%dx
    this%invdy = 1./this%dy

    ! set index ranges
#ifdef PARALLEL
    CALL InitMesh_parallel(this, config)
    CALL MPI_Barrier(MPI_COMM_WORLD,err)
    inum = this%IMAX - this%IMIN + 1
    CALL MPI_Allreduce(inum,this%MININUM,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,err)
    CALL MPI_Allreduce(inum,this%MAXINUM,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)
    jnum = this%JMAX - this%JMIN + 1
    CALL MPI_Allreduce(jnum,this%MINJNUM,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,err)
    CALL MPI_Allreduce(jnum,this%MAXJNUM,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)
#else
    this%IMIN   = 1
    this%IMAX   = this%INUM
    this%JMIN   = 1
    this%JMAX   = this%JNUM
#endif

    ! index ranges with ghost cells
    this%IGMIN = this%IMIN - this%GNUM
    this%IGMAX = this%IMAX + this%GNUM
    this%JGMIN = this%JMIN - this%GNUM
    this%JGMAX = this%JMAX + this%GNUM

    ! allocate memory for all pointers that are independent of fluxtype
    ALLOCATE(this%center(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%bcenter(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%bccart(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%fpos(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,2), &
         this%cpos(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,2), &
         this%fccart(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,2), &
         this%ccart(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,2), &
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
         this%radius(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%bradius(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX), &
         this%fradius(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4), &
         this%posvec(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%bposvec(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,2), &
         this%fposvec(this%IGMIN:this%IGMAX,this%JGMIN:this%JGMAX,4,2), &
        STAT=err)
    IF (err.NE.0) THEN
       CALL Error(this,"InitMesh_common","Unable to allocate memory!")
    END IF

    ! nullify remaining mesh arrays
    NULLIFY(this%dAx,this%dAy,this%dAz,this%dAxdy,this%dAydx, &
            this%chx,this%chy,this%chz,this%cyxy,this%cxyx, &
            this%czxz,this%czyz,this%sqrtg,this%invsqrtg, &
            this%rotation,this%weights)

    ! translation vectors for cell faces and cell corners
    ! (with respect to geometrical cell center)
    cfaces(1,1) = -0.5*this%dx  ! western x coordinate
    cfaces(1,2) = 0.0           ! western y coordinate
    cfaces(2,1) = 0.5*this%dx   ! eastern x coordinate
    cfaces(2,2) = 0.0           ! eastern y coordinate
    cfaces(3,1) = 0.0           ! southern x coordinate
    cfaces(3,2) = -0.5*this%dy  ! southern y coordinate
    cfaces(4,1) = 0.0           ! northern x coordinate
    cfaces(4,2) = 0.5*this%dy   ! northern y coordinate
    ccorners(1,1) = cfaces(1,1) ! south-west
    ccorners(1,2) = cfaces(3,2)
    ccorners(2,1) = cfaces(2,1) ! south-east
    ccorners(2,2) = cfaces(3,2)
    ccorners(3,1) = cfaces(1,1) ! north-west
    ccorners(3,2) = cfaces(4,2)
    ccorners(4,1) = cfaces(2,1) ! north-east
    ccorners(4,2) = cfaces(4,2)

    ! calculate coordinates
    DO j=this%JGMIN,this%JGMAX
       DO i=this%IGMIN,this%IGMAX
          ! geometrical cell centers
          this%center(i,j,1) = this%xmin + (2*i-1)*0.5*this%dx    ! x coord
          this%center(i,j,2) = this%ymin + (2*j-1)*0.5*this%dy    ! y coord
          ! cell face centered positions
          this%fpos(i,j,1,:) = this%center(i,j,:) + cfaces(1,:)   ! western
          this%fpos(i,j,2,:) = this%center(i,j,:) + cfaces(2,:)   ! eastern
          this%fpos(i,j,3,:) = this%center(i,j,:) + cfaces(3,:)   ! southern
          this%fpos(i,j,4,:) = this%center(i,j,:) + cfaces(4,:)   ! northern
          ! cell corner positions 
          this%cpos(i,j,1,:) = this%center(i,j,:) + ccorners(1,:) ! south-west
          this%cpos(i,j,2,:) = this%center(i,j,:) + ccorners(2,:) ! south-east
          this%cpos(i,j,3,:) = this%center(i,j,:) + ccorners(3,:) ! north-west
          this%cpos(i,j,4,:) = this%center(i,j,:) + ccorners(4,:) ! north-east
       END DO
    END DO

  END SUBROUTINE InitMesh


#ifdef PARALLEL
  !> \public Initialize MPI (parallel mode only)
  SUBROUTINE InitMesh_parallel(this, config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP) :: this
    TYPE(Dict_TYP),POINTER :: config
    !------------------------------------------------------------------------!
    LOGICAL, DIMENSION(NDIMS) :: periods = .FALSE.
    INTEGER        :: num,rem
    INTEGER        :: ierror
    INTEGER        :: i,j
    INTEGER        :: worldgroup,newgroup
    INTEGER, DIMENSION(1) :: rank0in, rank0out
    INTEGER, DIMENSION(NDIMS) :: coords
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ranks
    !------------------------------------------------------------------------!
    INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!

    ! create a cartesian topology of processes
    ! 1. balance number of processes per direction
    this%dims(1)=GetNumProcs(this)
    this%dims(2)=1
    ! account for vector length of vector CPUs
    CALL CalculateDecomposition(this%INUM,this%JNUM,this%GNUM, &
         this%dims(1),this%dims(2))
    IF (this%dims(2).LE.0) THEN
       CALL Error(this,"InitMesh_parallel","Domain decomposition algorithm failed.")
    END IF

    ! Check if the user set the decomposition dims himself and override the 
    ! automatic settings
    CALL RequireKey(config, "decomposition", this%dims(:))
    CALL GetAttr(config, "decomposition", this%dims)
    
    ! If a dimension equals -1, replace it with the total number of processors
    ! This makes it easy to define pure annular ring decompositions in polar 
    ! coordinates like as "decomposition" / (/ -1, 1 /) 
    WHERE(this%dims.EQ.-1) this%dims=GetNumProcs(this)

    IF((this%dims(1).LE.0).OR.(this%dims(2).LE.0).OR.&
       (this%dims(1)*this%dims(2).NE.GetNumProcs(this))) THEN
        CALL Error(this,"InitMesh_parallel","Invalid user-defined MPI domain "&
            //"decomposition with key='/mesh/decomposition'")
    END IF
    
    ! 2. create the cartesian communicator
    ! IMPORTANT: disable reordering of nodes
    CALL MPI_Cart_create(MPI_COMM_WORLD,NDIMS,this%dims,periods,.TRUE., &
         this%comm_cart,ierror)
 
    ! 3. inquire and save the own position
    CALL MPI_Cart_get(this%comm_cart,NDIMS,this%dims,periods,this%mycoords,ierror)

    ! subdivide the mesh and set mesh indices
    ! x-direction
    rem = MOD(this%INUM,this%dims(1))            ! remainder
    num = (this%INUM-rem) / this%dims(1)         ! fraction
    IF (this%mycoords(1).LT.rem) THEN
       ! the first (rem-1) nodes get one more to account for the remainder
       this%IMIN = this%mycoords(1) * (num + 1) + 1
       this%IMAX = this%IMIN + num
    ELSE
       this%IMIN = rem * (num+1) + (this%mycoords(1) - rem) * num + 1
       this%IMAX = this%IMIN + num - 1
    END IF
    ! y-direction
    rem = MOD(this%JNUM,this%dims(2))            ! remainder
    num = (this%JNUM-rem) / this%dims(2)         ! fraction
    IF (this%mycoords(2).LT.rem) THEN
       ! the first (rem-1) nodes get one more to account for the remainder
       this%JMIN = this%mycoords(2) * (num + 1) + 1
       this%JMAX = this%JMIN + num
    ELSE
       this%JMIN = rem * (num+1) + (this%mycoords(2) - rem) * num + 1
       this%JMAX = this%JMIN + num - 1
    END IF

    ! create communicators for all boundaries
    ALLOCATE(ranks(MAX(this%dims(1),this%dims(2))))
    ranks = 0
    CALL MPI_Comm_group(MPI_COMM_WORLD,worldgroup,ierror)
    ! western boundary
    coords(1) = 0
    DO j=0,this%dims(2)-1
       coords(2) = j
       CALL MPI_Cart_rank(this%comm_cart,coords,i,ierror)
       ranks(j+1)=i
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(2),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(1),ierror)
    rank0in(1) = 0
    CALL MPI_Group_translate_ranks(newgroup,1,rank0in,worldgroup,rank0out,ierror)
    this%rank0_boundaries(1) = rank0out(1)
    CALL MPI_Group_free(newgroup,ierror)
    ! eastern boundary
    coords(1) = this%dims(1)-1
    DO j=0,this%dims(2)-1
       coords(2) = j
       CALL MPI_Cart_rank(this%comm_cart,coords,i,ierror)
       ranks(j+1)=i
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(2),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(2),ierror)
    rank0in(1) = 0
    CALL MPI_Group_translate_ranks(newgroup,1,rank0in,worldgroup,rank0out,ierror)
    this%rank0_boundaries(2) = rank0out(1)
    CALL MPI_Group_free(newgroup,ierror)
    ! southern boundary
    coords(2) = 0
    DO i=0,this%dims(1)-1
       coords(1) = i
       CALL MPI_Cart_rank(this%comm_cart,coords,j,ierror)
       ranks(i+1)=j
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(1),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(3),ierror)
    rank0in(1) = 0
    CALL MPI_Group_translate_ranks(newgroup,1,rank0in,worldgroup,rank0out,ierror)
    this%rank0_boundaries(3) = rank0out(1)
    CALL MPI_Group_free(newgroup,ierror)
    ! northern boundary
    coords(2) = this%dims(2)-1
    DO i=0,this%dims(1)-1
       coords(1) = i
       CALL MPI_Cart_rank(this%comm_cart,coords,j,ierror)
       ranks(i+1)=j
    END DO
    CALL MPI_Group_incl(worldgroup,this%dims(1),ranks,newgroup,ierror)
    CALL MPI_Comm_create(this%comm_cart,newgroup,this%comm_boundaries(4),ierror)
    rank0in(1) = 0
    CALL MPI_Group_translate_ranks(newgroup,1,rank0in,worldgroup,rank0out,ierror)
    this%rank0_boundaries(4) = rank0out(1)
    CALL MPI_Group_free(newgroup,ierror)
    CALL MPI_Group_free(worldgroup,ierror)
    DEALLOCATE(ranks)
  END SUBROUTINE InitMesh_parallel


  ! return the best partitioning of processes
  ! pj x pk for a given mesh with resolution nj x nk with
  ! ghost number of ghost cells GNUM
  ! pj : total number of processes (<MAXNUM  see module "factors")
  ! pk = 1 initially
  ! pk return 0, for erroneous input
  SUBROUTINE CalculateDecomposition(nj,nk,gnum,pj,pk)
    USE factors
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER, INTENT(IN)    :: nj,nk,gnum
    INTEGER, INTENT(INOUT) :: pj,pk
    !------------------------------------------------------------------------!
    INTEGER :: p1,p2
    !------------------------------------------------------------------------!
    ! return immediatly for malformed input
    IF (((pj.LT.2).AND.(pk.NE.1)).OR.(pj.GT.MAXNUM) &
        .OR.(nj*nk.LT.pj)) THEN
       pk = 0
       RETURN
    END IF
    p1=Decompose(nj,nk,pj,pk)
    p2=pj/p1
    pj=p1
    pk=p2

  CONTAINS
    
    ! searches for the best domain decomposition accounting
    ! for the costs due to MPI communication (internal boundaries)
    ! and optimize for a given vector length (VECLEN) on vector computers;
    ! parameters:
    !   nj : number of grid cells in first dimension
    !   nk : number of grid cells in second dimension
    !   pj : number of processes  in first dimension (=NumProcs first call)
    !   pk : number of processes  in second dimension (=1 at first call)
    RECURSIVE FUNCTION Decompose(nj,nk,pj,pk) RESULT(pjres)
      IMPLICIT NONE
      !-------------------------------------------------------------------!
      INTEGER, INTENT(IN)    :: nj,nk,pj,pk
      INTEGER :: pjres
      !-------------------------------------------------------------------!
      INTEGER :: pp,ptot
      INTEGER :: p1,p2,pjnew,pknew,pjold,pkold
      INTEGER :: pfmin,pfnew,pfold
      INTEGER :: bl,vl,blnew,vlnew
      REAL    :: bl_gain,vl_gain
      !-------------------------------------------------------------------!
      ! measure the costs of the given configuration
      CALL GetCosts(nj,nk,pj,pk,bl,vl)
!!$PRINT '(A,2(I7),I12,I4)'," costs: ",pj,pk,bl,vl
      ! save the configuration
      pjold=pj
      pkold=pk
      p1 = pj
      p2 = pk
      ! compute the total number of processes
      ptot = pj*pk
      pfmin = GetFactor(pk)      ! get smallest prime factor of pk
      pfold = 1
      pp = pj
      DO ! loop over all prime factors of pj which are larger than
         ! the smallest prime factor of pk
         IF (pp.LE.1) EXIT       ! if pj has been reduced to 1
         ! get smallest prime factor of pp, i.e. pj
         pfnew = GetFactor(pp)
         IF ((pfnew.GT.pfmin).AND.(pfmin.NE.1)) EXIT
         pp = pp/pfnew
         ! skip multiple prime factors 
         IF (pfnew.NE.pfold) THEN
            ! create new configuration
            p1 = p1*pfold/pfnew
            p2 = p2/pfold*pfnew
            pfold = pfnew
            ! get the best configuration possible with p1 x p2 processes
            ! by reducing p1 => recursion
            pjnew = Decompose(nj,nk,p1,p2)
            pknew = ptot/pjnew  ! compute the second factor using the
                                ! total number of processes
            CALL GetCosts(nj,nk,pjnew,pknew,blnew,vlnew)
            bl_gain = bl*(1.0/blnew)             ! smaller is better
            vl_gain = vlnew*(1.0/vl)             ! larger is better
!!$PRINT '(4I7)',pjold,pkold,bl,vl
!!$PRINT '(4I7)',pjnew,pknew,blnew,vlnew
!!$PRINT *,"--------------------------------------------"
!!$PRINT '(A,3F7.1)',"              ",bl_gain,vl_gain,bl_gain*vl_gain
            ! compare new with old configuration
            IF (vl_gain*bl_gain.GT.1.0) THEN
               ! new configuration is better
               pjold = pjnew
               pkold = pknew
               bl = blnew
               vl = vlnew
            END IF
         END IF
      END DO
      ! return optimized number of processes in the first dimension
      ! (for the initial configuration pj x pk)
      pjres=pjold
    END FUNCTION Decompose
    
    ! computes the sum of the length of all internal boundaries (bl)
    ! and the maximal vector length (vl)
    PURE SUBROUTINE GetCosts(n1,n2,p1,p2,bl,vl)
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      INTEGER, INTENT(IN)  :: n1,n2,p1,p2
      INTEGER, INTENT(OUT) :: bl,vl
      !------------------------------------------------------------------------!
      INTEGER :: num,rem
      !------------------------------------------------------------------------!
      ! length of internal boundaries
      bl = n1*(p2-1) + n2*(p1-1)
      ! maximal possible vector length (first array dimension)
      ! if n1 > VECLEN return the length of the remainder
      rem = MOD(n1,p1)
      num = (n1-rem) / p1 + MIN(rem,1) + 2*gnum  ! account for ghost cells
      IF (num.GT.VECLEN) num=MOD(num,VECLEN)
      vl  = VECLEN-MOD(ABS(num-VECLEN),VECLEN)
    END SUBROUTINE GetCosts

  END SUBROUTINE CalculateDecomposition
#endif


  PURE FUNCTION GetMesh(this) RESULT(mt)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN) :: this
    INTEGER :: mt
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    mt=GetType_common(this%mtype)
  END FUNCTION GetMesh


  PURE FUNCTION GetMeshName(this) RESULT(mn)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN) :: this
    CHARACTER(LEN=32) :: mn
    !------------------------------------------------------------------------!
!CDIR IEXPAND
    mn=GetName_common(this%mtype)
  END FUNCTION GetMeshName


  PURE FUNCTION GetMeshRank(this) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN) :: this
    INTEGER :: r
    !------------------------------------------------------------------------!
    r = GetRank_common(this%mtype)
  END FUNCTION GetMeshRank


  PURE FUNCTION GetMeshNumProcs(this) RESULT(p)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN) :: this
    INTEGER :: p
    !------------------------------------------------------------------------!
    p = GetNumProcs_common(this%mtype)
  END FUNCTION GetMeshNumProcs

  PURE FUNCTION MeshInitialized(this) RESULT(i)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN) :: this
    LOGICAL :: i
    !------------------------------------------------------------------------!
    i = Initialized_common(this%mtype)
  END FUNCTION MeshInitialized

 
  SUBROUTINE MeshInfo(this,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: msg
    !------------------------------------------------------------------------!
    CALL Info_common(this%mtype,msg)
  END SUBROUTINE MeshInfo


  SUBROUTINE MeshWarning(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN) :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Warning_common(this%mtype,modproc,msg)
  END SUBROUTINE MeshWarning


  SUBROUTINE MeshError_rank0(this,modproc,msg)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN)    :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    !------------------------------------------------------------------------!
    CALL Error_common(this%mtype,modproc,msg)
  END SUBROUTINE MeshError_rank0


  SUBROUTINE MeshError_rankX(this,modproc,msg,rank)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP), INTENT(IN)    :: this
    CHARACTER(LEN=*),  INTENT(IN) :: modproc,msg
    INTEGER, INTENT(IN)           :: rank
    !------------------------------------------------------------------------!
    CALL Error_common(this%mtype,modproc,msg,rank)
  END SUBROUTINE MeshError_rankX


  !> \public Destructor of common mesh class
  SUBROUTINE CloseMesh(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Mesh_TYP)    :: this
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER :: i,ierror
    !------------------------------------------------------------------------!
    DO i=1,4
       IF (this%comm_boundaries(i).NE.MPI_COMM_NULL) &
          CALL MPI_Comm_free(this%comm_boundaries(i),ierror)
    END DO
#endif
    DEALLOCATE(this%center,this%bcenter,this%bccart, &
         this%fpos,this%cpos,this%fccart,this%ccart, &
         this%bhx,this%bhy,this%bhz, &
         this%fhx,this%fhy,this%fhz, &
         this%volume,this%dxdV,this%dydV,this%dlx,this%dly,&
         this%radius,this%bradius,this%fradius,&
         this%posvec,this%bposvec,this%fposvec)
    IF(ASSOCIATED(this%rotation)) &
        DEALLOCATE(this%rotation)
    CALL CloseCommon(this%mtype)
  END SUBROUTINE CloseMesh

END MODULE mesh_common
