!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_cashkarp.f90                                             #
!#                                                                           #
!# Copyright (C) 2011,2013                                                   #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
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
!> \author Björn Sperling
!! \author Tobias Illenseer
!!
!! \brief subroutines for embedded Runge-Kutta method
!!
!! Reference: G.Engeln-Müllges & F.Reutter; .....
!!
!! \extends timedisc_common
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_cashkarp
  
USE timedisc_common
  USE mesh_generic
  USE fluxes_generic
  USE boundary_generic
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE sources_generic
  USE timedisc_rkfehlberg, SolveODE_cashkarp => SolveODE_rkfehlberg, &
       CloseTimedisc_cashkarp => CloseTimedisc_rkfehlberg, &
       CalcTimestep_cashkarp => CalcTimestep_rkfehlberg
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "Cash-Karp method"

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! methods 
       InitTimedisc_cashkarp, &
       CloseTimedisc_cashkarp, &
       SolveODE_cashkarp, &
       CalcTimestep_cashkarp, &
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

  SUBROUTINE InitTimedisc_cashkarp(this,Mesh,Physics,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Timedisc_TYP) :: this
    TYPE(Mesh_TYP)     :: Mesh
    TYPE(Physics_TYP)  :: Physics
    TYPE(Dict_TYP), POINTER &
                       :: config
    !------------------------------------------------------------------------!
    INTEGER            :: err,method
    !------------------------------------------------------------------------!
    INTENT(IN)         :: Mesh,Physics
    INTENT(INOUT)      :: this
    !------------------------------------------------------------------------!
    ! set default order 
    CALL RequireKey(config, "order", 5)
    CALL GetAttr(config, "order", this%order)

    CALL GetAttr(config, "method", method)
    CALL InitTimedisc(this,method,ODEsolver_name)

!CDIR IEXPAND
    SELECT CASE(GetOrder(this))    
    CASE(5)
       !set number of coefficients
       this%m = 6 
       ! allocate memory 
       ALLOCATE(this%coeff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,this%m), &
                this%A1(this%m),this%A2(this%m),this%a(this%m),this%b(this%m,this%m), &
         STAT = err)
       IF (err.NE.0) THEN
          CALL Error(this,"timedisc_cashkarp", "Unable to allocate memory.")
       END IF
       !set coefficient scheme of Cash-Karp
       this%A1 = (/ 37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0 /)
       this%A2 = (/ 2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 0.25 /)
       this%a  = (/ 0.0, 0.2, 0.3, 0.6, 1.0, 7.0/8.0  /)
       this%b  = RESHAPE((/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                            0.2, 0.0, 0.0, 0.0, 0.0, 0.0, &
                            0.075, 0.225, 0.0, 0.0, 0.0, 0.0, &
                            0.3, -0.9, 1.2, 0.0, 0.0, 0.0, &
                           -11.0/54.0 , 2.5, -70.0/27.0, 35.0/27.0, 0.0, 0.0, &
                            1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0, 0.0/),(/this%m,this%m/))
    CASE DEFAULT
       CALL Error(this,"timedisc_cashkarp","time order must be 5")
    END SELECT
    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL Error(this,"timedisc_cashkarp", &
         "error tolerance levels must be greater than 0")
    IF (this%tol_rel.GT.1.0) THEN
         CALL Warning(this,"timedisc_cashkarp", &
            "adaptive step size control disabled (tol_rel>1)")
    ELSE IF(this%tol_rel.GE.0.01) THEN
         CALL Warning(this,"timedisc_cashkarp", &
             "You chose a relatively high tol_rel (in comparison to order)")
    END IF
  END SUBROUTINE InitTimedisc_cashkarp


END MODULE timedisc_cashkarp
