!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_modeuler.f90                                             #
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
! subroutines for modified Euler i.e. Runge-Kutta methods
!----------------------------------------------------------------------------!
MODULE timedisc_modeuler
  USE timedisc_common
  USE mesh_generic
  USE fluxes_generic
  USE boundary_generic
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE sources_generic
  IMPLICIT NONE
#ifdef PARALLEL
  include 'mpif.h'
#endif
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "modified Euler"
  REAL, DIMENSION(3) :: eta
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! methods 
       InitTimedisc_modeuler, &
       CloseTimedisc_modeuler, &
       SolveODE_modeuler, &
       SetBoundaries_modeuler, &
       GetBoundaryFlux, &
       GetType, &
       GetName, &
       GetOrder, &
       GetCFL, &
       GetRank, &
       GetNumProcs, &
       Info, &
       Warning, &
       Error
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc_modeuler(this,os,order,stoptime,cfl,dtlimit,maxiter)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    INTEGER            :: os, order,maxiter
    REAL               :: stoptime,cfl,dtlimit
    !------------------------------------------------------------------------!
    INTEGER            :: err
    !------------------------------------------------------------------------!
    INTENT(IN)         :: os,order,stoptime,cfl,dtlimit,maxiter
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    CALL InitTimedisc(this,os,ODEsolver_name,order,stoptime,cfl,dtlimit,maxiter)

    ! initialize parameter vector    
    SELECT CASE(GetOrder(this))
    CASE(1) ! time order=1 -> forward Euler
       eta(:) = (/ 0.0, 0.0, 0.0 /)
    CASE(2) ! time order=2 -> two step mod. Euler
       eta(:) = (/ 0.0, 0.5, 0.0 /)
    CASE(3) ! time order=3 -> three step mod. Euler
       eta(:) = (/ 0.0, 0.75, 1.0/3.0 /)
    CASE DEFAULT
       CALL Error(this,"InitTimedisc_modeuler","time order must be one of 1,2,3")
    END SELECT
  END SUBROUTINE InitTimedisc_modeuler


  SUBROUTINE SetBoundaries_modeuler(this,Mesh,Physics,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics,Fluxes
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    IF (PrimRecon(Fluxes%reconstruction)) THEN
       ! set center boundaries for primitive vars
       CALL Convert2Primitive(Physics,Mesh,this%cvar,this%pvar)
       CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,this%time,this%pvar)
    ELSE
       ! set center boundaries for conservative vars
       CALL CenterBoundary(this%boundary,Mesh,Fluxes,Physics,this%time,this%cvar)
       CALL Convert2Primitive(Physics,Mesh,this%cvar,this%pvar)
    END IF
  END SUBROUTINE SetBoundaries_modeuler


  SUBROUTINE SolveODE_modeuler(this,Mesh,Physics,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    !------------------------------------------------------------------------!
    INTEGER            :: n,i,j,k
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh
    INTENT(INOUT)      :: this,Physics,Fluxes
    !------------------------------------------------------------------------!

    DO n=1,GetOrder(this)
       ! get the numerical fluxes
       CALL CalculateFluxes(Fluxes,Mesh,Physics, &
            this%pvar,this%cvar,this%xflux,this%yflux)
       
       ! get geometrical sources for non-cartesian mesh
       IF (GetType(Mesh%geometry).NE.CARTESIAN) THEN 
          CALL GeometricalSources(Physics,Mesh,Fluxes, &
               this%pvar,this%cvar,this%geo_src)
       END IF
       
        ! get source terms due to external forces if present
       IF (ASSOCIATED(Physics%sources)) THEN
          CALL ExternalSources(Physics%sources,Mesh,Fluxes,Physics, &
               this%pvar,this%cvar,this%src)
       END IF

       DO k=1,Physics%VNUM
          ! flux ballance in x-direction
          IF (Mesh%INUM.GT.1) THEN
!CDIR NODEP
             DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NOUNROLL
                DO i=Mesh%IMIN,Mesh%IMAX
                   this%dxflux(i,j,k) = this%xflux(i,j,k) - this%xflux(i-1,j,k)
                END DO
                ! update western/eastern boundary fluxes (for accounting only)
                this%bxflux(j,1,k) = UpdateTimestep(eta(n),this%dt, this%bxfold(j,1,k), &
                     this%bxflux(j,1,k), Mesh%dy * this%xflux(Mesh%IMIN-1,j,k))
                this%bxflux(j,2,k) = UpdateTimestep(eta(n),this%dt, this%bxfold(j,2,k), &
                     this%bxflux(j,2,k), -Mesh%dy * this%xflux(Mesh%IMAX,j,k))
             END DO
          ELSE
             this%dxflux(:,:,k) = 0.
          END IF

          ! flux ballance in y-direction
          IF (Mesh%JNUM.GT.1) THEN
!CDIR NOUNROLL
             DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
                DO i=Mesh%IMIN,Mesh%IMAX
                   this%dyflux(i,j,k) = this%yflux(i,j,k) - this%yflux(i,j-1,k)
                END DO
             END DO
!CDIR NODEP
             DO i=Mesh%IMIN,Mesh%IMAX
                ! update southern/northern boundary fluxes (for accounting only)
                this%byflux(i,1,k) = UpdateTimestep(eta(n),this%dt,this%byfold(i,1,k), &
                     this%byflux(i,1,k), Mesh%dx * this%yflux(i,Mesh%JMIN-1,k))
                this%byflux(i,2,k) = UpdateTimestep(eta(n),this%dt,this%byfold(i,2,k), &
                     this%byflux(i,2,k), -Mesh%dx * this%yflux(i,Mesh%JMAX,k))
             END DO
          ELSE
             this%dyflux(:,:,k) = 0.
          END IF

          ! time step update of cell mean values
          this%cvar(:,:,k) = UpdateTimestep(eta(n),this%dt,this%cold(:,:,k),this%cvar(:,:,k), &
               Mesh%dydV(:,:)*this%dxflux(:,:,k) + Mesh%dxdV(:,:)*this%dyflux(:,:,k) &
               - this%geo_src(:,:,k) - this%src(:,:,k))
       END DO

       ! set boundary values
       CALL SetBoundaries_modeuler(this,Mesh,Physics,Fluxes)
    END DO

    CONTAINS

      ELEMENTAL FUNCTION UpdateTimestep(eta_n,dt,y0,yn,rhs) RESULT(y)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        REAL, INTENT(IN) :: eta_n,dt,y0,yn,rhs
        REAL :: y
        !------------------------------------------------------------------------!
        y = eta_n * y0 + (1.0 - eta_n) * (yn - dt * rhs)
      END FUNCTION UpdateTimestep

  END SUBROUTINE SolveODE_modeuler

  
  FUNCTION GetBoundaryFlux(this,Mesh,Physics,direction,comm) RESULT(bflux)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: direction
    REAL, DIMENSION(Physics%VNUM) :: bflux
    INTEGER, OPTIONAL  :: comm                        ! communicator for MPI !
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    INTEGER, DIMENSION(2) :: coords
    REAL, DIMENSION(Physics%VNUM) :: bflux_all,bflux_local
    INTEGER :: sender_rank
    INTEGER :: dest_comm,union_comm
    INTEGER :: world_group,dest_group,union_group,sender_group
    INTEGER :: ierror
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)         :: this,Mesh,Physics,direction  
    !------------------------------------------------------------------------!
#ifdef PARALLEL
    ! determine the destination for the result
    IF (PRESENT(comm)) THEN
       dest_comm = comm
    ELSE
       ! default is to send the result to all processes
       dest_comm = MPI_COMM_WORLD
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
    sender_rank = 0              ! with respect to Mesh%comm_boundaries(direction)
    CALL MPI_Reduce(bflux_local,bflux,Physics%VNUM,DEFAULT_MPI_REAL, &
         MPI_SUM,sender_rank,Mesh%comm_boundaries(direction),ierror)
    ! get destination group
    CALL MPI_Comm_group(dest_comm,dest_group,ierror)
    ! get world group
    CALL MPI_Comm_group(MPI_COMM_WORLD,world_group,ierror)
    ! get sender group
    CALL MPI_Group_incl(world_group,1,Mesh%rank0_boundaries(direction),sender_group,ierror)
    ! merge sender with destination
    CALL MPI_Group_union(sender_group,dest_group,union_group,ierror)
    ! create a communicator for the union group
    CALL MPI_Comm_create(MPI_COMM_WORLD,union_group,union_comm,ierror)
    ! get rank of sender in union group
    CALL MPI_Group_translate_ranks(sender_group,1,0,union_group,sender_rank,ierror)
    IF (sender_rank.EQ.MPI_UNDEFINED) &
         CALL Error(this,"GetBoundaryFlux","sender rank undefined")
    ! send result to all processes in communicator 'union_comm'
    CALL MPI_Bcast(bflux,Physics%VNUM,DEFAULT_MPI_REAL,sender_rank,union_comm,ierror)
    ! free all groups
    CALL MPI_Group_free(union_group,ierror)
    CALL MPI_Group_free(sender_group,ierror)
    CALL MPI_Group_free(world_group,ierror)
    CALL MPI_Group_free(dest_group,ierror)
#endif
  END FUNCTION GetBoundaryFlux


  SUBROUTINE CloseTimedisc_modeuler(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    !------------------------------------------------------------------------!
  END SUBROUTINE CloseTimedisc_modeuler

END MODULE timedisc_modeuler
