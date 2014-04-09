!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# program file: main.f90                                                    #
!#                                                                           #
!# Copyright (C) 2006 Tobias Illenseer <tillense@ita.uni-heidelberg.de>      #
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
! Main program file
!----------------------------------------------------------------------------!
PROGRAM fosite
  USE physics_generic
  USE sources_generic
  USE mesh_generic, ONLY : Mesh_TYP, CloseMesh
  USE fluxes_generic
  USE boundary_generic
  USE output_generic
  USE logio_generic
  USE reconstruction_generic
  USE timedisc_generic
  USE init
  USE integration
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  TYPE(Mesh_TYP)       :: Mesh
  TYPE(Fluxes_TYP)     :: Fluxes
  TYPE(Physics_TYP)    :: Physics
  TYPE(Output_TYP)     :: Output
  TYPE(Timedisc_TYP)   :: Timedisc
  TYPE(Logio_TYP)      :: Logio
  
  INTEGER              :: n
  INTEGER              :: rtc_count,rtc_rate  ! real time clock count/rate   ! 
  INTEGER              :: rtc_max             ! rtc max count                ! 
  REAL                 :: log_count           ! time for next log output     !
  REAL                 :: start_time          ! system clock start time      !
  REAL                 :: end_time            ! system clock end time        !
  !--------------------------------------------------------------------------!

  CALL InitIntegration

  ! print some information
  PRINT "(A)", "+---------------------------------------------------------+"
  PRINT "(A)", "|          Solution of 2D advection problems              |"
  PRINT "(A)", "+---------------------------------------------------------+"
  PRINT *, "Initializing program setup:"
  
  ! setup simulation
  CALL InitProgram(Mesh,Physics,Fluxes,Timedisc,Output,Logio)

  ! allocate memory for physics and fluxes modules
  CALL MallocPhysics(Physics,Mesh)
  CALL MallocFluxes(Fluxes,Mesh,Physics)

  ! set boundary values
  CALL SetBoundaries(Timedisc,Mesh,Physics,Fluxes)

  PRINT *, "==================================================================="
  PRINT *, "Starting calculation..."


  CALL SYSTEM_CLOCK(rtc_count,rtc_rate,rtc_max)
  log_count = rtc_count

  CALL CPU_TIME(start_time)

  ! main loop
  DO n=1,Timedisc%maxiter
     ! calculate timestep
     CALL CalcTimestep(Timedisc,Mesh,Physics)

     ! adjust timestep for output
     IF ((Timedisc%time+Timedisc%dt)/Output%time.GT.1.0) THEN
        Timedisc%dt = Output%time - Timedisc%time
     ELSE IF((Timedisc%time+1.5*Timedisc%dt)/Output%time.GT.1.0) THEN
        Timedisc%dt = 0.5*(Output%time - Timedisc%time)
     END IF

     ! advance the solution in time
     CALL SolveODE(Timedisc,Mesh,Physics,Fluxes)

     ! log output
     CALL SYSTEM_CLOCK(COUNT=rtc_count)
     IF (rtc_count.GE.log_count) THEN
        CALL WriteLogdata(Logio,Mesh,Physics,Timedisc)
        log_count = MODULO(rtc_count + rtc_rate*GetLogstep(Logio),rtc_max)
     END IF

     ! write output
     IF (ABS(1.0-Timedisc%time/Output%time).LT.1.0E-5) THEN
        CALL WriteOutput(Output,Mesh,Physics,Timedisc%pvar)
        CALL WriteLogdata(Logio,Mesh,Physics,Timedisc)
        CALL PrintInfo(n,Timedisc%time,Timedisc%dtmin,Timedisc%n_adj)
        ! reset dt_min and n_adj
        Timedisc%dtmin = Timedisc%stoptime
        Timedisc%n_adj = 0
        IF (Output%time.GT.GetEnd(Output)) THEN
           Output%time = Timedisc%stoptime
        END IF
     END IF

     ! finish simulation if stoptime is reached
     IF (ABS(1.0-Timedisc%time/Timedisc%stoptime).LT.1.0E-05) EXIT
  END DO

  PRINT *, "==================================================================="

  CALL CPU_TIME(end_time)
  PRINT "(A,F10.2,A)", " main loop runtime: ", end_time - start_time, " sec."

  IF (n.LT.Timedisc%maxiter) THEN
     PRINT *, "calculation finished correctly."
  ELSE
     PRINT *, "too many iterations, aborting!"
  END IF


  CALL CloseTimedisc(Timedisc)
  CALL CloseLogio
  CALL CloseOutput(Output)
  IF (ASSOCIATED(Physics%sources)) CALL CloseSources(Physics%sources,Fluxes)
  CALL CloseFluxes(Fluxes)
  CALL ClosePhysics(Physics)
  CALL CloseBoundary(Mesh%Boundary,WEST)
  CALL CloseBoundary(Mesh%Boundary,EAST)
  CALL CloseBoundary(Mesh%Boundary,SOUTH)
  CALL CloseBoundary(Mesh%Boundary,NORTH)
  CALL CloseMesh(Mesh,Fluxes)
  CALL CloseIntegration


CONTAINS

  SUBROUTINE PrintInfo(i,t,d,na)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER     :: i, na
    REAL        :: t, d
    !------------------------------------------------------------------------!
    INTENT(IN)  :: i,t,d
    !------------------------------------------------------------------------!
    PRINT "(A,I8,A,ES11.3,A,ES11.3,A,I5)", " Iteration ", i, &
         "  time ", t, "  min dt ", d, "  adj ", na
  END SUBROUTINE PrintInfo


END PROGRAM fosite
