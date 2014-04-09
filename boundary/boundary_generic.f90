!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: boundary_generic.f90                                              #
!#                                                                           #
!# Copyright (C) 2006-2010                                                   #
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
! Subroutines for boundary conditions
!----------------------------------------------------------------------------!
MODULE boundary_generic
  USE timedisc_common, ONLY : Timedisc_TYP
  USE mesh_common, ONLY : Mesh_TYP, GetRank, Error
  USE fluxes_common, ONLY : Fluxes_TYP
  USE reconstruction_common, ONLY : Reconstruction_TYP, PrimRecon
  USE boundary_nogradients, InitBoundary_common => InitBoundary
  USE boundary_periodic
  USE boundary_reflecting, InitBoundary_common1 => InitBoundary
  USE boundary_axis
  USE boundary_folded
  USE boundary_fixed
  USE boundary_extrapolation
  USE boundary_noh
  USE boundary_noslip
  USE boundary_custom
  USE boundary_farfield
  USE physics_generic, ONLY : Physics_TYP, Convert2Primitive, Convert2Conservative
  IMPLICIT NONE
#ifdef PARALLEL
    include 'mpif.h'
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  INTERFACE CloseBoundary
     MODULE PROCEDURE CloseBoundary_one, CloseBoundary_all
  END INTERFACE
  !--------------------------------------------------------------------------!
#ifdef PARALLEL
  INTEGER, PARAMETER :: NONE          = 0
#endif
  INTEGER, PARAMETER :: NO_GRADIENTS  = 1
  INTEGER, PARAMETER :: PERIODIC      = 2
  INTEGER, PARAMETER :: REFLECTING    = 3
  INTEGER, PARAMETER :: AXIS          = 4
  INTEGER, PARAMETER :: FOLDED        = 5
  INTEGER, PARAMETER :: FIXED         = 6
  INTEGER, PARAMETER :: EXTRAPOLATION = 7
  INTEGER, PARAMETER :: NOH2D         = 8
  INTEGER, PARAMETER :: NOH3D         = 9
  INTEGER, PARAMETER :: NOSLIP        = 10
  INTEGER, PARAMETER :: CUSTOM        = 11
  INTEGER, PARAMETER :: FARFIELD      = 12
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Boundary_TYP, &
       ! constants
       WEST, EAST, SOUTH, NORTH, &
       NO_GRADIENTS, PERIODIC, REFLECTING, AXIS, FOLDED, FIXED, EXTRAPOLATION, &
       NOH2D, NOH3D, NOSLIP, CUSTOM, FARFIELD, &
       CUSTOM_NOGRAD, CUSTOM_PERIOD, CUSTOM_REFLECT, CUSTOM_REFLNEG, &
       CUSTOM_EXTRAPOL, CUSTOM_FIXED, &
       ! methods
       InitBoundary, &
       CenterBoundary, &
       CloseBoundary, &
       GetType, &
       GetName, &
       GetDirection, &
       GetDirectionName, &
       GetRank, &
       GetNumProcs, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitBoundary_one(this,Mesh,Physics,btype,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: btype,dir
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER, PARAMETER :: strlen = 32
    CHARACTER(LEN=strlen) :: sendbuf
    CHARACTER(LEN=strlen) :: recvbuf
    INTEGER            :: status(MPI_STATUS_SIZE)
    INTEGER            :: ierror
    INTEGER            :: i
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,btype,dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    ! set boundary properties
    SELECT CASE(btype)
#ifdef PARALLEL
    CASE(NONE)
       CALL InitBoundary_common(this,btype,"none",dir)
#endif
    CASE(NO_GRADIENTS)
       CALL InitBoundary_nogradients(this,btype,dir)
    CASE(PERIODIC)
       CALL InitBoundary_periodic(this,btype,dir)
    CASE(REFLECTING)
       CALL InitBoundary_reflecting(this,Physics,btype,dir)
    CASE(AXIS)
       CALL InitBoundary_axis(this,Physics,btype,dir)
    CASE(FOLDED)
       CALL InitBoundary_folded(this,Mesh,Physics,btype,dir)
    CASE(FIXED)
       CALL InitBoundary_fixed(this,Mesh,Physics,btype,dir)
    CASE(EXTRAPOLATION)
       CALL InitBoundary_extrapolation(this,Physics,btype,dir)
    CASE(NOH2D)
       CALL InitBoundary_noh(this,Mesh,Physics,btype,dir,2)
    CASE(NOH3D)
       CALL InitBoundary_noh(this,Mesh,Physics,btype,dir,3)
    CASE(NOSLIP)
       CALL InitBoundary_noslip(this,Mesh,Physics,btype,dir)
    CASE(CUSTOM)
       CALL InitBoundary_custom(this,Mesh,Physics,btype,dir)
    CASE(FARFIELD)
       CALL InitBoundary_farfield(this,Mesh,Physics,btype,dir)
    CASE DEFAULT
       CALL Error(this, "InitBoundary_one", "Unknown boundary condition.")
    END SELECT

    ! print some information
#ifdef PARALLEL
    ! send boundary information to the rank 0 process;
    ! we only need this to synchronize the output
    IF (GetRank(this) .EQ. 0 .AND. GetRank(this).EQ.Mesh%rank0_boundaries(dir)) THEN
       ! print output without communication 
#endif
       CALL Info(this, " BOUNDARY-> condition:         " //  TRIM(GetDirectionName(this)) &
            // " " // TRIM(GetName(this)), GetRank(this))
#ifdef PARALLEL
    ELSE IF (GetRank(this).EQ.Mesh%rank0_boundaries(dir)) THEN
       ! send info to root
       sendbuf = TRIM(GetDirectionName(this))//" "//TRIM(GetName(this))
       CALL MPI_SEND(sendbuf,strlen,MPI_CHARACTER,0, &
                     0,MPI_COMM_WORLD,ierror)
    ELSE IF (GetRank(this).EQ.0) THEN
       ! receive input from rank0_boundaries(dir)
       CALL MPI_RECV(recvbuf,strlen,MPI_CHARACTER,Mesh%rank0_boundaries(dir),&
                     MPI_ANY_TAG,MPI_COMM_WORLD,status,ierror)
       CALL Info(this, " BOUNDARY-> condition:         " // TRIM(recvbuf),GetRank(this))
    END IF
#endif
  END SUBROUTINE InitBoundary_one


  SUBROUTINE InitBoundary(this,Mesh,Physics,western,eastern,southern,northern)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), DIMENSION(4) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: western, eastern, southern, northern
    !------------------------------------------------------------------------!
    INTEGER            :: new_western, new_eastern, new_southern, new_northern
#ifdef PARALLEL
    INTEGER            :: comm_old
    INTEGER            :: sizeofreal, ignum, jgnum, twoslices
    INTEGER            :: ierr
    LOGICAL, DIMENSION(SIZE(Mesh%dims)) :: periods = .FALSE.
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Physics,western,eastern,southern,northern
    INTENT(INOUT) :: this,Mesh
    !------------------------------------------------------------------------!

    new_western = western
    new_eastern = eastern
    new_southern = southern
    new_northern = northern

#ifdef PARALLEL
    ! define connections
    IF (Mesh%mycoords(1).NE.0)  new_western = NONE
    IF (Mesh%mycoords(1).NE.Mesh%dims(1)-1)  new_eastern = NONE
    IF (Mesh%mycoords(2).NE.0)  new_southern = NONE
    IF (Mesh%mycoords(2).NE.Mesh%dims(2)-1)  new_northern = NONE
#endif
    
    ! initialize every boundary
    ! IMPORTANT: do this before anything else
    CALL InitBoundary_one(this(WEST),Mesh,Physics,new_western,WEST)
    CALL InitBoundary_one(this(EAST),Mesh,Physics,new_eastern,EAST)
    CALL InitBoundary_one(this(SOUTH),Mesh,Physics,new_southern,SOUTH)
    CALL InitBoundary_one(this(NORTH),Mesh,Physics,new_northern,NORTH)
    
#ifdef PARALLEL
    ! check periodicity
    IF (western.EQ.PERIODIC.AND.eastern.EQ.PERIODIC) THEN
       periods(1) = .TRUE.
    ELSE IF (western.EQ.PERIODIC.NEQV.eastern.EQ.PERIODIC) THEN
       CALL Error(this(WEST),"InitBoundary", &
            "Opposite boundary should be periodic.")
    END IF
    IF (southern.EQ.PERIODIC.AND.northern.EQ.PERIODIC) THEN
       periods(2) = .TRUE.
    ELSE IF (southern.EQ.PERIODIC.NEQV.northern.EQ.PERIODIC) THEN
       CALL Error(this(SOUTH),"InitBoundary", &
            "Opposite boundary should be periodic.")
    END IF

    ! create new cartesian communicator using Mesh%comm_cart
    ! and account for the periodicity
    ! IMPORTANT: disable reordering of nodes
    comm_old = Mesh%comm_cart
    CALL MPI_Cart_create(comm_old,SIZE(Mesh%dims),Mesh%dims,periods,.FALSE., &
         Mesh%comm_cart,ierr)

    ! save ranks of neighbor processes
    CALL MPI_Cart_shift(Mesh%comm_cart,0,1,Mesh%neighbor(WEST),Mesh%neighbor(EAST),ierr)
    CALL MPI_Cart_shift(Mesh%comm_cart,1,1,Mesh%neighbor(SOUTH),Mesh%neighbor(NORTH),ierr)

    ! allocate memory for boundary data buffers
    ALLOCATE(this(WEST)%sendbuf(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%vnum), &
         this(WEST)%recvbuf(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%vnum), &
         this(EAST)%sendbuf(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%vnum), &
         this(EAST)%recvbuf(Mesh%GNUM,Mesh%JMIN:Mesh%JMAX,Physics%vnum), &
         this(SOUTH)%sendbuf(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM,Physics%vnum), &
         this(SOUTH)%recvbuf(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM,Physics%vnum), &
         this(NORTH)%sendbuf(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM,Physics%vnum), &
         this(NORTH)%recvbuf(Mesh%IMIN:Mesh%IMAX,Mesh%GNUM,Physics%vnum), &
         STAT=ierr)
    IF (ierr.NE.0) THEN
       CALL Error(this(WEST),"InitBoundary", &
            "Unable to allocate memory for data buffers.")
    END IF

!!$    ! create a data handle for the boundary arrays
!!$    ! in each direction
!!$    CALL MPI_Type_extent(DEFAULT_MPI_REAL,sizeofreal,ierr)
!!$
!!$    ! total number of cells (min to max with ghosts cells)
!!$    ignum = 1+Mesh%IGMAX-Mesh%IGMIN
!!$    jgnum = 1+Mesh%JGMAX-Mesh%JGMIN
!!$    ! 1. western/eastern boundaries (horizontal direction)
!!$    ! create datatype for a 2D section in the x-y plane (one variable)
!!$    CALL MPI_Type_vector(jgnum-2*Mesh%GNUM,Mesh%GNUM,ignum,DEFAULT_MPI_REAL,twoslices,ierr)
!!$    ! create datatype for the full boundary block (all variables)
!!$    CALL MPI_Type_hvector(Physics%vnum,1,ignum*jgnum*sizeofreal,twoslices,Mesh%weblocks,ierr)
!!$    CALL MPI_Type_commit(Mesh%weblocks,ierr) 
!!$
!!$    ! 2. southern/northern boundaries (vertical direction)
!!$    ! create datatype for a 2D section in the x-y plane (one variable)
!!$    CALL MPI_Type_vector(Mesh%GNUM,ignum-2*Mesh%GNUM,jgnum,DEFAULT_MPI_REAL,twoslices,ierr)
!!$    ! create datatype for the full boundary block (all variables)
!!$    CALL MPI_Type_hvector(Physics%vnum,1,ignum*jgnum*sizeofreal,twoslices,Mesh%snblocks,ierr)
!!$    CALL MPI_Type_commit(Mesh%snblocks,ierr)
#endif

  END SUBROUTINE InitBoundary


  SUBROUTINE CenterBoundary(this,Mesh,Fluxes,Physics,time,pvar,cvar)
    IMPLICIT NONE
#ifdef PARALLEL
    include 'mpif.h'
#endif
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this(4)
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Fluxes_TYP)   :: Fluxes
    TYPE(Physics_TYP)  :: Physics
    REAL               :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM) &
                       :: pvar, cvar
    !------------------------------------------------------------------------!
    INTEGER       :: i,j
    INTEGER       :: ierr
#ifdef PARALLEL
    INTEGER       :: req(4)
    INTEGER       :: status(MPI_STATUS_SIZE)
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Mesh,Physics,Fluxes,time
    INTENT(INOUT) :: this,pvar,cvar
    !------------------------------------------------------------------------!
    CALL Convert2Primitive(Physics,Mesh,Mesh%IMIN,Mesh%IMAX,Mesh%JMIN, &
             Mesh%JMAX,cvar,pvar)
#ifdef PARALLEL
    ! IMPORTANT:
    ! First exchange the boundary data of the _internal_ boundaries before
    ! processing the _real_ boundaries of the computational domain

    ! send boundary data to western neighbor
    this(WEST)%sendbuf(:,:,:) = pvar(Mesh%IMIN:Mesh%IMIN+Mesh%GNUM-1,Mesh%JMIN:Mesh%JMAX,1:Physics%vnum)
    CALL MPI_Issend(this(WEST)%sendbuf,Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+WEST,Mesh%comm_cart,req(WEST),ierr)
    ! receive boundary data from eastern neighbor
    CALL MPI_Recv(this(EAST)%recvbuf,Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+WEST,Mesh%comm_cart,status,ierr)
    CALL MPI_Wait(req(WEST),status,ierr)
    pvar(Mesh%IMAX+1:Mesh%IGMAX,Mesh%JMIN:Mesh%JMAX,1:Physics%vnum) = this(EAST)%recvbuf(:,:,:)

    ! send boundary data to eastern neighbor
    this(EAST)%sendbuf(:,:,:) = pvar(Mesh%IMAX-Mesh%GNUM+1:Mesh%IMAX,Mesh%JMIN:Mesh%JMAX,1:Physics%vnum)
    CALL MPI_Issend(this(EAST)%sendbuf,Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(EAST),10+EAST,Mesh%comm_cart,req(EAST),ierr)
    ! receive boundary data from western neighbor
    CALL MPI_Recv(this(WEST)%recvbuf,Mesh%GNUM*(Mesh%JMAX-Mesh%JMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(WEST),10+EAST,Mesh%comm_cart,status,ierr)
    CALL MPI_Wait(req(EAST),status,ierr)
    pvar(Mesh%IGMIN:Mesh%IMIN-1,Mesh%JMIN:Mesh%JMAX,1:Physics%vnum) = this(WEST)%recvbuf(:,:,:)
    
    ! send boundary data to southern neighbor
    this(SOUTH)%sendbuf(:,:,:) = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMIN:Mesh%JMIN+Mesh%GNUM-1,1:Physics%vnum)
    CALL MPI_Issend(this(SOUTH)%sendbuf,Mesh%GNUM*(Mesh%IMAX-Mesh%IMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+SOUTH,Mesh%comm_cart,req(SOUTH),ierr)
    ! receive boundary data from northern neighbor
    CALL MPI_Recv(this(NORTH)%recvbuf,Mesh%GNUM*(Mesh%IMAX-Mesh%IMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+SOUTH,Mesh%comm_cart,status,ierr)
    CALL MPI_Wait(req(SOUTH),status,ierr)
    pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX+1:Mesh%JGMAX,1:Physics%vnum) = this(NORTH)%recvbuf(:,:,:)

    ! send boundary data to northern neighbor
    this(NORTH)%sendbuf(:,:,:) = pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JMAX-Mesh%GNUM+1:Mesh%JMAX,1:Physics%vnum)
    CALL MPI_Issend(this(NORTH)%sendbuf,Mesh%GNUM*(Mesh%IMAX-Mesh%IMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(NORTH),10+NORTH,Mesh%comm_cart,req(NORTH),ierr)
    ! receive boundary data from southern neighbor
    CALL MPI_Recv(this(SOUTH)%recvbuf,Mesh%GNUM*(Mesh%IMAX-Mesh%IMIN+1)*Physics%vnum, &
         DEFAULT_MPI_REAL,Mesh%neighbor(SOUTH),10+NORTH,Mesh%comm_cart,status,ierr)
    CALL MPI_Wait(req(NORTH),status,ierr)
    pvar(Mesh%IMIN:Mesh%IMAX,Mesh%JGMIN:Mesh%JMIN-1,1:Physics%vnum) = this(SOUTH)%recvbuf(:,:,:)

#endif

    ! Set boundary values depending on the selected boundary condition.
    IF (Mesh%INUM.GT.1) THEN
       DO i=1,2
          SELECT CASE(GetType(this(i)))
          CASE(NO_GRADIENTS)
             CALL CenterBoundary_nogradients(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(PERIODIC)
             ! do nothing in parallel version, because periodicity is
             ! handled via MPI communication
#ifndef PARALLEL
             CALL CenterBoundary_periodic(this(i),Mesh,Physics,pvar)
#endif
          CASE(REFLECTING)
             CALL CenterBoundary_reflecting(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(AXIS)
             CALL CenterBoundary_axis(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(FOLDED)
             CALL CenterBoundary_folded(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(FIXED)
             CALL CenterBoundary_fixed(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(EXTRAPOLATION)
             CALL CenterBoundary_extrapolation(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(NOH2D,NOH3D)
             CALL CenterBoundary_noh(this(i),Mesh,Physics,Fluxes,time,pvar)
          CASE(NOSLIP)
             CALL CenterBoundary_noslip(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(CUSTOM)
             CALL CenterBoundary_custom(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(FARFIELD)
             CALL CenterBoundary_farfield(this(i),Mesh,Physics,Fluxes,pvar)
          END SELECT
       END DO
    END IF
    IF (Mesh%JNUM.GT.1) THEN
       DO i=3,4
          SELECT CASE(GetType(this(i)))
          CASE(NO_GRADIENTS)
             CALL CenterBoundary_nogradients(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(PERIODIC)
             ! do nothing in parallel version, because periodicity is
             ! handled via MPI communication
#ifndef PARALLEL
             CALL CenterBoundary_periodic(this(i),Mesh,Physics,pvar)
#endif
          CASE(REFLECTING)
             CALL CenterBoundary_reflecting(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(AXIS)
             CALL CenterBoundary_axis(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(FOLDED)
             CALL CenterBoundary_folded(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(FIXED)
             CALL CenterBoundary_fixed(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(EXTRAPOLATION)
             CALL CenterBoundary_extrapolation(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(NOH2D,NOH3D)
             CALL CenterBoundary_noh(this(i),Mesh,Physics,Fluxes,time,pvar)
          CASE(NOSLIP)
             CALL CenterBoundary_noslip(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(CUSTOM)
             CALL CenterBoundary_custom(this(i),Mesh,Physics,Fluxes,pvar)
          CASE(FARFIELD)
             CALL CenterBoundary_farfield(this(i),Mesh,Physics,Fluxes,pvar)
          END SELECT
       END DO
    END IF

    ! FIXME: implement MPI communication to exchange corner values
    ! this is a quick hack to set defined boundary values in
    ! the corners outside the computational domain;
    ! this is also necessary, because we need some of these values in the
    ! viscosity module
    IF ((Mesh%INUM > 1) .AND. (Mesh%JNUM > 1)) THEN
       ! south west
       pvar(Mesh%IMIN-1,Mesh%JMIN-1,:) = 0.5 * (pvar(Mesh%IMIN,Mesh%JMIN-1,:) &
         + pvar(Mesh%IMIN-1,Mesh%JMIN,:))
       pvar(Mesh%IMIN-1,Mesh%JMIN-2,:) = pvar(Mesh%IMIN-1,Mesh%JMIN-1,:)
       pvar(Mesh%IMIN-2,Mesh%JMIN-1,:) = pvar(Mesh%IMIN-1,Mesh%JMIN-1,:)
       pvar(Mesh%IMIN-2,Mesh%JMIN-2,:) = pvar(Mesh%IMIN-1,Mesh%JMIN-1,:)
       ! south east
       pvar(Mesh%IMAX+1,Mesh%JMIN-1,:) = 0.5 * (pvar(Mesh%IMAX,Mesh%JMIN-1,:) &
         + pvar(Mesh%IMAX+1,Mesh%JMIN,:))
       pvar(Mesh%IMAX+2,Mesh%JMIN-1,:) = pvar(Mesh%IMAX+1,Mesh%JMIN-1,:)
       pvar(Mesh%IMAX+1,Mesh%JMIN-2,:) = pvar(Mesh%IMAX+1,Mesh%JMIN-1,:)
       pvar(Mesh%IMAX+2,Mesh%JMIN-2,:) = pvar(Mesh%IMAX+1,Mesh%JMIN-1,:)
       ! north west
       pvar(Mesh%IMIN-1,Mesh%JMAX+1,:) = 0.5 * (pvar(Mesh%IMIN-1,Mesh%JMAX,:) &
         + pvar(Mesh%IMIN,Mesh%JMAX+1,:))
       pvar(Mesh%IMIN-2,Mesh%JMAX+1,:) = pvar(Mesh%IMIN-1,Mesh%JMAX+1,:)
       pvar(Mesh%IMIN-1,Mesh%JMAX+2,:) = pvar(Mesh%IMIN-1,Mesh%JMAX+1,:)
       pvar(Mesh%IMIN-2,Mesh%JMAX+2,:) = pvar(Mesh%IMIN-1,Mesh%JMAX+1,:)
       ! north east
       pvar(Mesh%IMAX+1,Mesh%JMAX+1,:) = 0.5 * (pvar(Mesh%IMAX+1,Mesh%JMAX,:) &
         + pvar(Mesh%IMAX,Mesh%JMAX+1,:))
       pvar(Mesh%IMAX+2,Mesh%JMAX+1,:) = pvar(Mesh%IMAX+1,Mesh%JMAX+1,:)
       pvar(Mesh%IMAX+1,Mesh%JMAX+2,:) = pvar(Mesh%IMAX+1,Mesh%JMAX+1,:)
       pvar(Mesh%IMAX+2,Mesh%JMAX+2,:) = pvar(Mesh%IMAX+1,Mesh%JMAX+1,:)
    ELSE IF (Mesh%JNUM.EQ.1) THEN
       ! for 1D simulations along the x-direction copy internal data
       ! into the ghost cells at the y-boundaries
       DO j=1,Mesh%GNUM
          pvar(:,Mesh%JMIN-j,:) = pvar(:,Mesh%JMIN,:)
          pvar(:,Mesh%JMAX+j,:) = pvar(:,Mesh%JMAX,:)
       END DO
    ELSE IF (Mesh%INUM.EQ.1) THEN
       ! for 1D simulations along the y-direction copy internal data
       ! into the ghost cells at the x-boundaries
       DO i=1,Mesh%GNUM
          pvar(Mesh%IMIN-i,:,:) = pvar(Mesh%IMIN,:,:)
          pvar(Mesh%IMAX+i,:,:) = pvar(Mesh%IMAX,:,:)
       END DO
    END IF

    ! convert primitive variables in ghost cells
    CALL Convert2Conservative(Physics,Mesh,Mesh%IGMIN,Mesh%IMIN-1,&
         Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CALL Convert2Conservative(Physics,Mesh,Mesh%IMAX+1,Mesh%IGMAX,&
         Mesh%JGMIN,Mesh%JGMAX,pvar,cvar)
    CALL Convert2Conservative(Physics,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
         Mesh%JGMIN,Mesh%JMIN-1,pvar,cvar)
    CALL Convert2Conservative(Physics,Mesh,Mesh%IGMIN,Mesh%IGMAX,&
         Mesh%JMAX+1,Mesh%JGMAX,pvar,cvar)

  END SUBROUTINE CenterBoundary


  SUBROUTINE CloseBoundary_one(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(NO_GRADIENTS,PERIODIC,EXTRAPOLATION)
       ! do nothing
    CASE(REFLECTING)
       CALL CloseBoundary_reflecting(this)
    CASE(AXIS)
       CALL CloseBoundary_axis(this)
    CASE(FOLDED)
       CALL CloseBoundary_folded(this)
    CASE(FIXED)
       CALL CloseBoundary_fixed(this)
    CASE(NOSLIP)
       CALL CloseBoundary_noslip(this)
    CASE(CUSTOM)
       CALL CloseBoundary_custom(this)
    CASE(FARFIELD)
       CALL CloseBoundary_farfield(this)
    END SELECT
  END SUBROUTINE CloseBoundary_one


  SUBROUTINE CloseBoundary_all(this,dir)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Boundary_TYP), DIMENSION(4) :: this
    INTEGER       :: dir
    !------------------------------------------------------------------------!
    INTENT(IN)    :: dir
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    SELECT CASE(dir)
    CASE(WEST,EAST,SOUTH,NORTH)
       CALL CloseBoundary_one(this(dir))
    END SELECT
  END SUBROUTINE CloseBoundary_all


END MODULE boundary_generic
