!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_generic.f90                                              #
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
! generic subroutines for time discretization
!----------------------------------------------------------------------------!
MODULE timedisc_generic
  USE timedisc_modeuler, SetBoundaries => SetBoundaries_modeuler
  USE mesh_generic
  USE physics_generic
  USE fluxes_generic
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: MODIFIED_EULER = 1
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! constants
       MODIFIED_EULER, &
       ! methods 
       InitTimedisc, &
       CalcTimestep, &
       SolveODE, &
       SetBoundaries, &
       CloseTimedisc
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitTimedisc(this,Mesh,Physics,method,order,cfl,stoptime,dtlimit,maxiter)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: method
    INTEGER            :: order,maxiter
    REAL               :: cfl,stoptime,dtlimit
    !------------------------------------------------------------------------!
    INTEGER            :: err
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics,method,order,cfl,stoptime,dtlimit,maxiter
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!

    ! call individual constructors
    SELECT CASE(method)
    CASE(MODIFIED_EULER)
       CALL InitTimedisc_modeuler(this,method,order,cfl,stoptime,dtlimit,maxiter)
    CASE DEFAULT
       PRINT *,"ERROR in InitTimedisc: Unknown ODE solver"
       STOP       
    END SELECT

    ! allocate memory for data structures needed in all timedisc modules
    ALLOCATE(this%pvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%cvar(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%pold(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%cold(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%geo_src(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%src(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%xflux(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%yflux(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%amax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         STAT = err)
    IF (err.NE.0) THEN
       PRINT *, "ERROR in InitTimedisc: Can't allocate memory!"
       STOP
    END IF

    ! initialize all arrays
    this%pvar = 0.
    this%cvar = 0.
    this%pold = 0.
    this%cold = 0.
    this%src = 0.
    this%geo_src = 0.
    this%xflux = 0.
    this%yflux = 0.
    this%amax = 0.

    ! print some information
    PRINT "(A,A)", " TIMEDISC-> ODE solver:        ", TRIM(GetName(this))
    PRINT "(A,I1)", "            order:             ", GetOrder(this)
    PRINT "(A,F4.2)", "            CFL number:        ", GetCFL(this)
  END SUBROUTINE InitTimedisc


  SUBROUTINE CalcTimestep(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    !------------------------------------------------------------------------!
    REAL               :: invdt_x, invdt_y
    REAL               :: dt_cfl, dt_src
    !------------------------------------------------------------------------!    
    INTENT(IN)         :: Mesh
    INTENT(INOUT)      :: this,Physics
    !------------------------------------------------------------------------!   
    ! store old values
    this%cold(:,:,:) = this%cvar(:,:,:)
    this%pold(:,:,:) = this%pvar(:,:,:)

    ! CFL condition:
    ! maximal wave speeds in each direction
    CALL MaxWaveSpeeds(Physics,Mesh,this%pvar,this%amax)
    
    ! inverse of time step in each direction
    invdt_x = MAXVAL(this%amax(:,:,1) / Mesh%dlx(:,:))    
    invdt_y = MAXVAL(this%amax(:,:,2) / Mesh%dly(:,:))
  
    ! largest time step due to CFL condition
    dt_cfl = this%cfl / MAX(invdt_x, invdt_y)
    
    ! time step due to source terms
!!$    dt_src = GetSourcesTimescale(Physics%sources,dt_cfl)
!!$    
!!$    dt = MIN(dt_cfl,dt_src)
!!$
    this%dt = dt_cfl
  END SUBROUTINE CalcTimestep


  PURE SUBROUTINE AdjustTimestep(this,Mesh,Physics)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP):: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! adjust time step
    this%dt = this%dt * 0.5
    ! count adjustments for information
    this%n_adj = this%n_adj + 1
    ! set data to old values
    this%cvar(:,:,:) = this%cold(:,:,:)
    this%pvar(:,:,:) = this%pold(:,:,:)
  END SUBROUTINE AdjustTimestep


  SUBROUTINE SolveODE(this,Mesh,Physics,Fluxes)
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

    DO
       SELECT CASE(GetType(this))
       CASE(MODIFIED_EULER)
          CALL SolveODE_modeuler(this,Mesh,Physics,Fluxes)
       END SELECT

       ! check data and adjust time step if necessary
       IF (CheckData(Physics,Mesh,this%pvar,this%pold)) THEN
          CALL AdjustTimestep(this,Mesh,Physics)
          ! abort if time step is to small
          IF (this%dt.LT.this%dtlimit) THEN
             PRINT "(A)", "ERROR in SolveODE: time step to small, aborting .."
             STOP
          END IF
       ELSE
          ! just for information
          this%dtmin = MIN(this%dt,this%dtmin)
          this%time  = this%time + this%dt
          RETURN
       END IF
    END DO
  END SUBROUTINE SolveODE


  SUBROUTINE CloseTimedisc(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    !------------------------------------------------------------------------!

    ! call individual destructors
    SELECT CASE(GetType(this))
    CASE(MODIFIED_EULER)
       CALL CloseTimedisc_modeuler(this)
    CASE DEFAULT
       PRINT *,"ERROR in InitTimedisc: Unknown ODE solver"
       STOP       
    END SELECT

    DEALLOCATE(this%pvar,this%cvar,this%pold,this%cold, &
         this%geo_src,this%src,this%xflux,this%yflux,this%amax)
  END SUBROUTINE CloseTimedisc

END MODULE timedisc_generic
