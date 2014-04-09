!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: timedisc_dormand_prince.f90                                       #
!#                                                                           #
!# Copyright (C) 2013                                                        #
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
!> \author Björn Sperling
!!
!! \brief subroutines for Dormand-Prince method
!!
!! Reference: Dormand, J. R.; Prince, P. J. (1980),
!!            "A family of embedded Runge-Kutta formulae",
!!            Journal of Computational and Applied Mathematics 6 (1): 19–26,
!!            doi:10.1016/0771-050X(80)90013-3
!!
!! \extends timedisc_common
!! \ingroup timedisc
!----------------------------------------------------------------------------!
MODULE timedisc_dormand_prince
  USE timedisc_common
  USE mesh_generic
  USE fluxes_generic
  USE boundary_generic
  USE physics_generic, GeometricalSources_Physics => GeometricalSources, &
       ExternalSources_Physics => ExternalSources
  USE sources_generic
  USE timedisc_rkfehlberg, SolveODE_dormand_prince => SolveODE_rkfehlberg, &
       CloseTimedisc_dormand_prince => CloseTimedisc_rkfehlberg, &
       CalcTimestep_dormand_prince => CalcTimestep_rkfehlberg
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: ODEsolver_name = "Dormand-Prince method"

  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Timedisc_TYP, &
       ! methods 
       InitTimedisc_dormand_prince, &
       CloseTimedisc_dormand_prince, &
       SolveODE_dormand_prince, &
       CalcTimestep_dormand_prince, &
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

  SUBROUTINE InitTimedisc_dormand_prince(this,Mesh,Physics,config)
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
       this%m = 7 
       ! allocate memory 
       ALLOCATE(this%coeff(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%VNUM,this%m), &
                this%A1(this%m),this%A2(this%m),this%a(this%m),this%b(this%m,this%m), &
         STAT = err)
       IF (err.NE.0) THEN
          CALL Error(this,"timedisc_dormand_prince", "Unable to allocate memory.")
       END IF
       !set coefficient scheme of dormand_prince 
       this%A1 = (/ 35./384.,0.,500./1113.,125./192.,-2187./6784.,11./84.,0. /)
       this%A2 = (/ 5179./57600.,0.,7571./16695.,393./640.,-92097./339200.,187./2100.,1./40. /)
       this%a  = (/ 0.0, 0.2, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0 /)
       this%b  = &
            RESHAPE((/ 0.,0.,0.,0.,0.,0.,0., &
                       0.2,0.,0.,0.,0.,0.,0., &
                       3.0/40.0,9.0/40.0,0.,0.,0.,0.,0., &
                       44./45.,-56./15.,32./9.,0.,0.,0.,0., &
                       19372./6561.,-25360./2187.,64448./6561.,-212./729.,0.,0.,0.,&
                       9017./3168.,-355./33.,46732./5247.,49./176.,-5103./18656,0.,0., &
                       35./384.,0.,500./1113.,125./192.,-2187./6784.,11./84.,0. /),(/this%m,this%m/))
    CASE DEFAULT
       CALL Error(this,"timedisc_dormand_prince","time order must be 5")
    END SELECT
    IF ((this%tol_rel.LT.0.0).OR.MINVAL(this%tol_abs(:)).LT.0.0) &
         CALL Error(this,"timedisc_dormand_prince", &
         "error tolerance levels must be greater than 0")
    IF (this%tol_rel.GT.1.0) THEN
         CALL Warning(this,"timedisc_dormand_prince", &
            "adaptive step size control disabled (tol_rel>1)")
    ELSE IF(this%tol_rel.GE.0.01) THEN
         CALL Warning(this,"timedisc_dormand_prince", &
             "You chose a relatively high tol_rel (in comparison to order)")
    END IF
  END SUBROUTINE InitTimedisc_dormand_prince

END MODULE timedisc_dormand_prince
