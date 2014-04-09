!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: multipole_generic.f90                                            #
!#                                                                           #
!# Copyright (C) 2009-2011                                                   #
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
! generic module for multipole expansion
!----------------------------------------------------------------------------!
MODULE multipole_generic
  USE mesh_common, ONLY : Selection_TYP
  USE physics_common, ONLY : Physics_TYP
  USE multipole_spherical, InitMultipole_basic => InitMultipole, &
       CloseMultipole_common => CloseMultipole
  USE multipole_cylindrical
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  INTEGER, PARAMETER :: SPHERMULTEXPAN       = 1
  INTEGER, PARAMETER :: CYLINMULTEXPAN       = 2
  !--------------------------------------------------------------------------!
  INTERFACE CalculatePotential
     MODULE PROCEDURE CalculatePotential_1, CalculatePotential_2
  END INTERFACE
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Multipole_TYP, &
       ! constants
       SPHERMULTEXPAN, CYLINMULTEXPAN, &
       ! methods
       InitMultipole, &
       CalculatePotential, &
       CloseMultipole, &
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

  SUBROUTINE InitMultipole(this,etype,cart_coords,volume,imin,imax,jmin,jmax, &
       ireg,oreg,order)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP) :: this
    TYPE(Selection_TYP) :: ireg,oreg(:)
    INTEGER             :: etype,imin,imax,jmin,jmax
    REAL, DIMENSION(imin:imax,jmin:jmax,2) :: cart_coords
    REAL, DIMENSION(imin:imax,jmin:jmax)   :: volume
    INTEGER, OPTIONAL   :: order
    !------------------------------------------------------------------------!
    INTEGER             :: order_def
    CHARACTER(LEN=32)   :: info_str
    !------------------------------------------------------------------------!
    INTENT(IN)    :: etype,cart_coords,volume,imin,imax,jmin,jmax,ireg,oreg,order
    INTENT(INOUT) :: this
    !------------------------------------------------------------------------!
    ! set default values for optional arguments
    IF (PRESENT(order)) THEN
       order_def = order
    ELSE
       order_def = 1
    END IF
    SELECT CASE(etype)
    CASE(SPHERMULTEXPAN)
       CALL InitMultipole_spherical(this,etype,cart_coords,volume,imin,imax, &
            jmin,jmax,ireg,oreg,order)
    CASE(CYLINMULTEXPAN)
       ! order of cylindrical multipole expansion is allways 1
       CALL InitMultipole_cylindrical(this,etype,cart_coords,volume,imin,imax, &
            jmin,jmax,ireg,oreg,1)
    CASE DEFAULT
       CALL Error(this,"InitMultipole","unknown expansion type")
    END SELECT
    ! print some information
    WRITE (info_str,'(I0)') this%ORDER
    CALL Info(this, " MULTIPOLE> expansion type:    " // TRIM(GetName(this)) // &
         ACHAR(10)//"            order:             " // TRIM(info_str))
  END SUBROUTINE InitMultipole


  PURE SUBROUTINE CalculatePotential_1(this,Physics,Selection,rho,Phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP)  :: this
    TYPE(Selection_TYP)  :: Selection(:)
    TYPE(Physics_TYP)    :: Physics
    REAL, DIMENSION(this%IMIN:this%IMAX,this%JMIN:this%JMAX) :: rho,Phi
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Physics,Selection,rho
    INTENT(INOUT) :: this,Phi
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(SPHERMULTEXPAN)
       CALL CalculatePotential_spherical(this,Physics,Selection,rho,Phi)
    CASE(CYLINMULTEXPAN)
       CALL CalculatePotential_cylindrical(this,Physics,Selection,rho,Phi)
    END SELECT
  END SUBROUTINE CalculatePotential_1


  PURE SUBROUTINE CalculatePotential_2(this,Physics,rho,Phi)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP)  :: this
    TYPE(Physics_TYP)    :: Physics
    REAL, DIMENSION(this%IMIN:this%IMAX,this%JMIN:this%JMAX) :: rho,Phi
    !------------------------------------------------------------------------!
    INTENT(IN)    :: Physics,rho
    INTENT(INOUT) :: this,Phi
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(SPHERMULTEXPAN)
       CALL CalculatePotential_spherical(this,Physics,this%oregion,rho,Phi)
    CASE(CYLINMULTEXPAN)
       CALL CalculatePotential_cylindrical(this,Physics,rho,Phi)
    END SELECT
  END SUBROUTINE CalculatePotential_2


  SUBROUTINE CloseMultipole(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Multipole_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)  :: this
    !------------------------------------------------------------------------!
    SELECT CASE(GetType(this))
    CASE(SPHERMULTEXPAN)
       CALL CloseMultipole_spherical(this)
    CASE(CYLINMULTEXPAN)
       CALL CloseMultipole_cylindrical(this)
    END SELECT
  END SUBROUTINE CloseMultipole

END MODULE multipole_generic
