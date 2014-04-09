!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: fluxes_midpoint.f90                                               #
!#                                                                           #
!# Copyright (C) 2007-2008                                                   #
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
! numerical fluxes module for midpoint quadrature rule
!----------------------------------------------------------------------------!
MODULE fluxes_midpoint
  USE fluxes_common
  USE mesh_common, ONLY : Mesh_TYP
  USE boundary_common, ONLY : &
#ifdef PARALLEL
       DEFAULT_MPI_REAL, &
#endif
       WEST, EAST, SOUTH, NORTH
  USE physics_generic
  USE reconstruction_generic
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
  ! name for the numerical flux function
  CHARACTER(LEN=32), PARAMETER :: qrule_name = "midpoint"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Fluxes_TYP, &
       ! methods
       InitFluxes, &
       InitFluxes_midpoint, &
       CalculateFaceData, &
       CalculateFluxes_midpoint, &
       GetBoundaryFlux, &
       PrimRecon, &
       GetType, &
       GetName, &
       GetRank, &
       GetNumProcs, &
       Initialized, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitFluxes_midpoint(this,qrule)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    INTEGER           :: qrule
    !------------------------------------------------------------------------!
    INTENT(IN)        :: qrule
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL InitFluxes(this,qrule,qrule_name)
  END SUBROUTINE InitFluxes_midpoint


  PURE SUBROUTINE CalculateFaceData(this,Mesh,Physics,pvar,cvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,cvar
    INTENT(INOUT)     :: this,Physics
    !------------------------------------------------------------------------!

    ! reconstruct data on cell faces
    IF (PrimRecon(this%reconstruction)) THEN
       CALL CalculateSlopes(this%Reconstruction,Mesh,Physics,pvar)
       CALL CalculateStates(this%Reconstruction,Mesh,Physics,4,Mesh%bcenter,&
            Mesh%fpos,pvar,this%rstates)
       CALL Convert2Conservative(Physics,Mesh,this%rstates,this%cons)
    ELSE
       CALL CalculateSlopes(this%Reconstruction,Mesh,Physics,cvar)
       CALL CalculateStates(this%Reconstruction,Mesh,Physics,4,Mesh%bcenter,&
            Mesh%fpos,cvar,this%rstates)
       CALL Convert2Primitive(Physics,Mesh,this%rstates,this%prim)
    END IF

    ! get minimal & maximal wave speeds on cell interfaces
    CALL CalculateWaveSpeeds(Physics,Mesh,this%prim)
  END SUBROUTINE CalculateFaceData


  PURE SUBROUTINE CalculateFluxes_midpoint(this,Mesh,Physics,pvar,cvar,xflux,yflux)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)  :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar,cvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: xflux,yflux
    !------------------------------------------------------------------------!
    INTEGER           :: i,j,k
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,pvar,cvar
    INTENT(INOUT)     :: this,Physics
    INTENT(OUT)       :: xflux,yflux
    !------------------------------------------------------------------------!

    ! execute generic tasks common to all flux types
    CALL CalculateFaceData(this,Mesh,Physics,pvar,cvar)

    ! west and east
    IF (Mesh%INUM.GT.1) THEN
       ! physical fluxes
       CALL CalculateFluxesX(Physics,Mesh,1,2,this%prim,this%cons,this%pfluxes)
       ! numerical fluxes
       DO k=1,Physics%VNUM
!CDIR OUTERUNROLL=8
          DO j=Mesh%JGMIN,Mesh%JGMAX
!CDIR NODEP
             DO i=Mesh%IMIN-1,Mesh%IMAX
                ! x-direction
                xflux(i,j,k) = Mesh%dAxdy(i+1,j,1) / &
                     (Physics%amax(i,j) - Physics%amin(i,j)) * &
                     (Physics%amax(i,j)*this%pfluxes(i,j,2,k) - &
                     Physics%amin(i,j)*this%pfluxes(i+1,j,1,k) + &
                     Physics%amin(i,j)*Physics%amax(i,j) * &
                     (this%cons(i+1,j,1,k) - this%cons(i,j,2,k)))
             END DO
          END DO
       END DO
    ELSE
       xflux(:,:,:) = 0.0
    END IF

    ! south and north
    IF (Mesh%JNUM.GT.1) THEN
       ! physical fluxes
       CALL CalculateFluxesY(Physics,Mesh,3,4,this%prim,this%cons,this%pfluxes)
       ! numerical fluxes
       DO k=1,Physics%VNUM
!CDIR OUTERUNROLL=8
          DO j=Mesh%JMIN-1,Mesh%JMAX
!CDIR NODEP
             DO i=Mesh%IGMIN,Mesh%IGMAX
                ! y-direction
                yflux(i,j,k) = Mesh%dAydx(i,j+1,1) / &
                     (Physics%bmax(i,j) - Physics%bmin(i,j)) * &
                     (Physics%bmax(i,j)*this%pfluxes(i,j,4,k) - &
                     Physics%bmin(i,j)*this%pfluxes(i,j+1,3,k) + &
                     Physics%bmin(i,j)*Physics%bmax(i,j) * &
                     (this%cons(i,j+1,3,k) - this%cons(i,j,4,k)))
             END DO
          END DO
       END DO
    ELSE
       yflux(:,:,:) = 0.0
    END IF
  END SUBROUTINE CalculateFluxes_midpoint


  FUNCTION GetBoundaryFlux(this,Mesh,Physics,direction,comm) RESULT(bflux)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Fluxes_TYP)   :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: direction
    REAL, DIMENSION(Physics%VNUM) :: bflux
    INTEGER, OPTIONAL  :: comm                        ! communicator for MPI !
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER, DIMENSION(2) :: coords
    REAL, DIMENSION(Physics%VNUM) :: bflux_all,bflux_local
    INTEGER :: sender_rank(1),dest_ranks(1),rank0(1)
    INTEGER :: dest_comm,union_comm
    INTEGER :: world_group,dest_group,union_group,sender_group
    INTEGER :: ierror
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics,direction  
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    ! create world group
    CALL MPI_Comm_group(MPI_COMM_WORLD,world_group,ierror)
    ! determine the destination for the result
    IF (PRESENT(comm)) THEN
       dest_comm = comm
       ! set destination group
       CALL MPI_Comm_group(dest_comm,dest_group,ierror)
    ELSE
       ! if no communicator is given, send the result
       ! to the rank0 process
       dest_ranks(1) = 0
       CALL MPI_Group_incl(world_group,1,dest_ranks,dest_group,ierror)
    END IF
#endif
    SELECT CASE(direction)
    CASE(WEST) ! western boundary flux
#ifdef PARALLEL
       IF (Mesh%mycoords(1).EQ.0) THEN
          bflux_local(:) = SUM(this%bxflux(:,1,:),1)
       ELSE
          bflux_local(:) = 0.0
       END IF
#else
       bflux(:) = SUM(this%bxflux(:,1,:),1)
#endif
    CASE(EAST) ! eastern boundary flux
#ifdef PARALLEL
       IF (Mesh%mycoords(1).EQ.Mesh%dims(1)-1) THEN
          bflux_local(:) = SUM(this%bxflux(:,2,:),1)
       ELSE
          bflux_local(:) = 0.0
       END IF
#else
       bflux = SUM(this%bxflux(:,2,:),1)
#endif
    CASE(SOUTH) ! southern boundary flux
#ifdef PARALLEL
       IF (Mesh%mycoords(2).EQ.0) THEN
          bflux_local(:) = SUM(this%byflux(:,1,:),1)
       ELSE
          bflux_local(:) = 0.0
       END IF
#else
       bflux(:) = SUM(this%byflux(:,1,:),1)       
#endif
    CASE (NORTH) ! northern boundary flux
#ifdef PARALLEL
       IF (Mesh%mycoords(2).EQ.Mesh%dims(2)-1) THEN
          bflux_local(:) = SUM(this%byflux(:,2,:),1)
       ELSE
          bflux_local(:) = 0.0
       END IF
#else
       bflux(:) = SUM(this%byflux(:,2,:),1)
#endif
    CASE DEFAULT
       CALL Error(this,"GetBoundaryFlux","wrong direction")
    END SELECT
#ifdef PARALLEL
    ! collect and sum up the result in process with rank 0 with respect to the
    ! subset of MPI processes at the boundary under consideration
    IF (Mesh%comm_boundaries(direction).NE.MPI_COMM_NULL) THEN
       CALL MPI_Reduce(bflux_local,bflux,Physics%VNUM,DEFAULT_MPI_REAL, &
            MPI_SUM,0,Mesh%comm_boundaries(direction),ierror)
    ELSE
       bflux(:) = 0.0
    END IF
    ! get sender group
    rank0(1) = Mesh%rank0_boundaries(direction)
    CALL MPI_Group_incl(world_group,1,rank0,sender_group,ierror)
    ! merge sender with destination
    CALL MPI_Group_union(sender_group,dest_group,union_group,ierror)
    ! create a communicator for the union group
    CALL MPI_Comm_create(MPI_COMM_WORLD,union_group,union_comm,ierror)
    IF (union_comm.NE.MPI_COMM_NULL) THEN
       ! get rank of sender in union group
       rank0(1) = 0
       CALL MPI_Group_translate_ranks(sender_group,1,rank0,union_group,sender_rank,ierror)
       IF (sender_rank(1).EQ.MPI_UNDEFINED) &
           CALL Error(this,"GetBoundaryFlux","sender rank undefined")
       ! send result to all processes in communicator 'union_comm'
       CALL MPI_Bcast(bflux,Physics%VNUM,DEFAULT_MPI_REAL,sender_rank(1),union_comm,ierror)
       ! free union communicator
       CALL MPI_Comm_free(union_comm,ierror)
    END IF
    ! free all groups
    CALL MPI_Group_free(union_group,ierror)
    CALL MPI_Group_free(sender_group,ierror)
    CALL MPI_Group_free(world_group,ierror)
    CALL MPI_Group_free(dest_group,ierror)
#endif
  END FUNCTION GetBoundaryFlux

END MODULE fluxes_midpoint
