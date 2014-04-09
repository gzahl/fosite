!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_ssprk.f90                                                #
!#                                                                           #
!# Copyright (C) 2013                                                        #
!# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
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
!> \author Manuel Jung
!!
!! \brief subroutines for strong stability preserving (SSP) Runge Kutta methods
!!
!! Reference: Gottlieb et. al (2011): Strong stability preserving runge-kutta
!! and multistep time discretization (book)
!!
!! \extends timedisc_common
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_ssprk
  
USE timedisc_common
  USE mesh_generic
  USE fluxes_generic
  USE boundary_generic
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE sources_generic
  USE timedisc_rkfehlberg, SolveODE_ssprk => SolveODE_rkfehlberg, &
       CloseTimedisc_ssprk => CloseTimedisc_rkfehlberg, &
       CalcTimestep_ssprk => CalcTimestep_rkfehlberg
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "SSP Runge-Kutta method"

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! methods 
       InitTimedisc_ssprk, &
       CloseTimedisc_ssprk, &
       SolveODE_ssprk, &
       CalcTimestep_ssprk, &
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

  SUBROUTINE InitTimedisc_ssprk(this,Mesh,Physics,config)
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
    CASE(3)
       CALL Warning(this, "timedisc_ssprk", &
         "This 3rd order embedded SSP RK scheme has been constructed from a second " // &
         "third order SSP RK scheme by hand! It seems to work, but i am not sure, that one " // &
         "is allowed to do so. An optimal embedded third order SSP RK scheme is described " // &
         "in chapter 6.3 of the main reference (see above), but still needs to be translated into " // &
         "a butchers tableau. => Better use the 5th order scheme!")
       !set number of coefficients
       this%m = 3
       ! allocate memory 
       ALLOCATE(this%coeff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,this%m), &
                this%A1(this%m),this%A2(this%m),this%a(this%m),this%b(this%m,this%m), &
         STAT = err)
       IF (err.NE.0) THEN
          CALL Error(this,"timedisc_ssprk", "Unable to allocate memory.")
       END IF
       this%A1 = (/ 1./6., 1./6., 2./3. /)
       this%A2 = (/ 0.5,   0.5,   0.    /)
       this%a  = (/ 0.,    1.,    0.5   /)
       this%b  = RESHAPE((/ 0.0,  0.0,  0.0, &
                            1.0,  0.0,  0.0, &
                            0.25, 0.25, 0.0 /),(/this%m,this%m/))
    CASE(5) 
       ! RK(5,4) SSP(5,3) 
       ! Reference: Colin Barr Macdonald: Constructing High-Order Runge_lutta Methods with Embedded 
       ! Strong-Stability-PReserving Pairs (2003)
       !set number of coefficients
       this%m = 5
       ! allocate memory 
       ALLOCATE(this%coeff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,this%m), &
                this%A1(this%m),this%A2(this%m),this%a(this%m),this%b(this%m,this%m), &
         STAT = err)
       IF (err.NE.0) THEN
          CALL Error(this,"timedisc_ssprk", "Unable to allocate memory.")
       END IF
       this%A1 = (/ 0.17279, 0.094505, 0.12947, 0.29899, 0.30424 /)
       this%A2 = (/ 0.12293, 0.31981, -0.15316, 0.31887, 0.39155 /)
       this%a  = (/ 0., 0.36717, 0.58522, 0.44156, 0.8464 /)
       this%b  = RESHAPE((/ 0., 0., 0., 0., 0., &
                            0.36717, 0., 0., 0., 0., &
                            0.26802, 0.3172, 0., 0., 0., &
                            0.11606, 0.13735, 0.18816, 0., 0., &
                            0.11212, 0.13269, 0.18178, 0.4198, 0. /),(/this%m,this%m/))
    CASE DEFAULT
       CALL Error(this,"timedisc_ssprk","time order must be 3 or 5")
    END SELECT
    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL Error(this,"timedisc_ssprk", &
         "error tolerance levels must be greater than 0")
    IF (this%tol_rel.GT.1.0) THEN
         CALL Warning(this,"timedisc_ssprk", &
            "adaptive step size control disabled (tol_rel>1)")
    ELSE IF(this%tol_rel.GE.0.01) THEN
         CALL Warning(this,"timedisc_ssprk", &
             "You chose a relatively high tol_rel (in comparison to order)")
    END IF
  END SUBROUTINE InitTimedisc_ssprk


END MODULE timedisc_ssprk
