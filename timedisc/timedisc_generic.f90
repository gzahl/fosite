!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_generic.f90                                              #
!#                                                                           #
!# Copyright (C) 2007-2010                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
  USE timedisc_modeuler
  USE boundary_generic
  USE mesh_generic
  USE physics_common, ONLY : GetErrorMap
  USE physics_generic
  USE fluxes_generic
  USE sources_generic, CalcTimestep_sources => CalcTimestep
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
       CloseTimedisc, &
       Print_Checkdata, &
       GetOrder, &
       GetCFL, &
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

  SUBROUTINE InitTimedisc(this,Mesh,Physics,method,order,stoptime,cfl,dtlimit,maxiter)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    INTEGER            :: method
    INTEGER            :: order,maxiter
    REAL               :: stoptime
    REAL, OPTIONAL     :: cfl,dtlimit
    !------------------------------------------------------------------------!
    INTEGER            :: err
    CHARACTER(LEN=8)   :: order_str, cfl_str
    REAL               :: cfl_def,dtlimit_def
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics,method,order,stoptime,cfl,dtlimit,maxiter
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    IF (.NOT.Initialized(Physics).OR..NOT.Initialized(Mesh)) &
         CALL Error(this,"InitTimedisc","physics and/or mesh module uninitialized")
    ! set default values
    IF (PRESENT(cfl)) THEN
       cfl_def = cfl
    ELSE
       cfl_def = 0.4
    END IF
    IF (PRESENT(dtlimit)) THEN
       dtlimit_def = dtlimit
    ELSE
       dtlimit_def = 1.0E-16 * stoptime
    END IF

    ! call individual constructors
    SELECT CASE(method)
    CASE(MODIFIED_EULER)
       CALL InitTimedisc_modeuler(this,method,order,stoptime,cfl_def,&
            dtlimit_def,maxiter)
    CASE DEFAULT
       CALL Error(this,"InitTimedisc", "Unknown ODE solver.")
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
         this%dxflux(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%dyflux(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum), &
         this%amax(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         STAT = err)
    IF (err.NE.0) THEN
       CALL Error(this,"InitTimedisc", "Unable to allocate memory.")
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
    WRITE (order_str, '(I0)') GetOrder(this)
    WRITE (cfl_str, '(F4.2)') GetCFL(this)
    CALL Info(this," TIMEDISC-> ODE solver:        " // TRIM(GetName(this)) // ACHAR(10) //&
                   "            order:             " // TRIM(order_str) // ACHAR(10) // &
                   "            CFL number:        " // TRIM(cfl_str))
  END SUBROUTINE InitTimedisc


  SUBROUTINE CalcTimestep(this,Mesh,Physics,Fluxes)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
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
    Fluxes%bxfold(:,:,:) = Fluxes%bxflux(:,:,:)
    Fluxes%byfold(:,:,:) = Fluxes%byflux(:,:,:)

    ! CFL condition:
    ! maximal wave speeds in each direction
    CALL MaxWaveSpeeds(Physics,Mesh,this%pvar,this%amax)

    ! inverse of time step in each direction
    IF (Mesh%INUM.GT.1) THEN
       invdt_x = MAXVAL(this%amax(:,:,1) / Mesh%dlx(:,:))
    ELSE
       ! set to zero, i.e. no CFL limit in x-direction
       invdt_x = 0.0
    END IF
    IF (Mesh%JNUM.GT.1) THEN
       invdt_y = MAXVAL(this%amax(:,:,2) / Mesh%dly(:,:))
    ELSE
       ! set to zero, i.e. no CFL limit in y-direction
       invdt_y = 0.0
    END IF
  
    ! largest time step due to CFL condition
    dt_cfl = this%cfl / MAX(invdt_x, invdt_y)
    
    ! initialize this to be sure dt_src > 0
    dt_src = dt_cfl
    CALL CalcTimestep_sources(Physics%sources,Mesh,Physics,this%pvar,this%cvar,dt_src)
    this%dt = MIN(dt_cfl,dt_src)
  END SUBROUTINE CalcTimestep


  SUBROUTINE SolveODE(this,Mesh,Physics,Fluxes)
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
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Fluxes_TYP)   :: Fluxes
    !------------------------------------------------------------------------!
    INTEGER            :: bad_data
#ifdef PARALLEL
    INTEGER            :: bad_data_all
    INTEGER            :: ierror
#endif
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh
    INTENT(INOUT)      :: this,Physics,Fluxes
    !------------------------------------------------------------------------!

    DO
       SELECT CASE(GetType(this))
       CASE(MODIFIED_EULER)
          CALL SolveODE_modeuler(this,Mesh,Physics,Fluxes)
       END SELECT

       ! check data
       bad_data = CheckData(Physics,Mesh,this%pvar,this%pold)
       IF (bad_data.NE.0) THEN
          IF ((this%dt * 0.5).LT.this%dtlimit) THEN
             PRINT *,"Return value of CheckData: ", bad_data
             CALL Print_Checkdata(this,Mesh,Physics)
             CALL Error(this,"SolveODE", "Time step to small, aborting.",GetRank(this))
          END IF
       END IF
#ifdef PARALLEL
       CALL MPI_Allreduce(bad_data,bad_data_all,1,MPI_LOGICAL,MPI_LOR,Mesh%comm_cart,ierror)
       bad_data = bad_data_all
#endif

       ! adjust time step if necessary
       IF (bad_data.NE.0) THEN
          ! adjust time step
          this%dt = this%dt * 0.5
          ! abort if time step is to small
          ! count adjustments for information
          this%n_adj = this%n_adj + 1
          ! set data to old values
          this%cvar(:,:,:) = this%cold(:,:,:)
          this%pvar(:,:,:) = this%pold(:,:,:)
          Fluxes%bxflux(:,:,:) = Fluxes%bxfold(:,:,:)
          Fluxes%byflux(:,:,:) = Fluxes%byfold(:,:,:)
       ELSE
          ! just for information
          this%dtmin = MIN(this%dt,this%dtmin)
          this%time  = this%time + this%dt
          RETURN
       END IF
    END DO
  END SUBROUTINE SolveODE

 SUBROUTINE Print_Checkdata(this,Mesh,Physics)
   IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    !------------------------------------------------------------------------!
    INTEGER, DIMENSION(2,6) :: mapping 
    INTEGER  :: i, j                
    REAL     :: rhomin, pmin
    INTEGER, DIMENSION(4) :: meshrange
    !------------------------------------------------------------------------!

    IF ((Mesh%JNUM > 3) .AND. (Mesh%INUM > 3)) THEN
       !print boundary mark (North)
       write(*, '(A)', advance='YES') "       N"
       j = 5
       mapping(1,:) = (/ Mesh%IGMIN, Mesh%IMIN, Mesh%IMIN+1, Mesh%IMAX, Mesh%IMAX+1, Mesh%IGMAX+1 /)
       mapping(2,:) = (/ Mesh%JGMIN, Mesh%JMIN, Mesh%JMIN+1, Mesh%JMAX, Mesh%JMAX+1, Mesh%JGMAX+1 /)
    ELSE IF (Mesh%JNUM > 3) THEN
       !print boundary mark (North)
       write(*, '(A)', advance='YES') "    N"
       j = 5
       mapping(1,3:4) = (/ Mesh%IGMIN, Mesh%IGMAX+1 /)
       mapping(2,:) = (/ Mesh%JGMIN, Mesh%JMIN, Mesh%JMIN+1, Mesh%JMAX, Mesh%JMAX+1, Mesh%JGMAX+1 /)
    ELSE
       j = 3
       mapping(1,:) = (/ Mesh%IGMIN, Mesh%IMIN, Mesh%IMIN+1, Mesh%IMAX, Mesh%IMAX+1, Mesh%IGMAX+1 /)
       mapping(2,3:4) = (/ Mesh%JGMIN, Mesh%JGMAX+1 /)
    END IF
    
    DO
      !print seperator
      IF ((Mesh%INUM > 3) .AND. (j == 1 .or. j == 4)) THEN
         write(*, '(A)', advance='YES') "    -------"     
      ELSE IF (j == 1 .or. j == 4) THEN
         write(*, '(A)', advance='YES') "    -" 
      END IF
      !print boundary mark (West)
      IF ((Mesh%INUM > 3) .AND. (j == 3)) THEN
         write(*,'(A)',advance='NO') " W  "
      ELSE 
         write(*,'(A)',advance='NO') "    "
      END IF

      IF (Mesh%INUM > 3) THEN 
         i = 1
      ELSE 
         i = 3
      END IF
      DO 
         !print seperator
         IF (i == 2 .or. i == 5) write(*, '(A)', advance='NO') "|"     
      
         meshrange(1) = mapping(1,i)
         meshrange(2) = mapping(1,i+1)-1
         meshrange(3) = mapping(2,j)
         meshrange(4) = mapping(2,j+1)-1
      
         write(*,'(A)',advance='NO') GetErrorMap(Physics, CheckData(Physics,Mesh,this%pvar,this%pold,meshrange))
         i = i+1
         if ((i > 5) .OR. .NOT.(Mesh%INUM > 3)) exit
      END DO
      !print boundary mark (East)
      IF ((Mesh%INUM > 3) .AND. (j == 3)) THEN 
         write (*,*) " E"
      ELSE 
         write (*,*) "  "
      END IF
      j = j-1
      if ((j < 1) .OR. .NOT.(Mesh%JNUM > 3)) exit
    END DO

    IF ((Mesh%JNUM > 3) .AND. (Mesh%INUM > 3)) THEN
       !print boundary mark (South)
       write(*, '(A)', advance='YES') "       S"
    ELSE IF (Mesh%JNUM > 3) THEN
       write(*, '(A)', advance='YES') "    S"
    END IF

  END SUBROUTINE Print_Checkdata


  SUBROUTINE CloseTimedisc(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP)   :: this
    !------------------------------------------------------------------------!
    ! call boundary destructors
    CALL CloseBoundary(this%Boundary,WEST)
    CALL CloseBoundary(this%Boundary,EAST)
    CALL CloseBoundary(this%Boundary,SOUTH)
    CALL CloseBoundary(this%Boundary,NORTH)

    ! call individual destructors
    SELECT CASE(GetType(this))
    CASE(MODIFIED_EULER)
       CALL CloseTimedisc_modeuler(this)
    END SELECT

    DEALLOCATE(this%pvar,this%cvar,this%pold,this%cold, &
         this%geo_src,this%src,this%xflux,this%yflux,this%dxflux,this%dyflux, &
         this%amax)
  END SUBROUTINE CloseTimedisc

END MODULE timedisc_generic
