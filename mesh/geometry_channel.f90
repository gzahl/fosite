!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: geometry_channel.f90                                              #
!#                                                                           #
!# Copyright (C) 2011,2014                                                   #
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
!! \brief define properties of a 2.5D channel mesh for direct numerical 
!! simulations of a channel flow
!!
!! (Kim and Moin, J. Fluid Mech. 118. pp 341-377, 1982)
!!
!! \extends geometry_cartesian
!! \ingroup geometry
!----------------------------------------------------------------------------!
MODULE geometry_channel
  USE functions, ONLY : Atanh
  USE geometry_cartesian
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! exclude interface block from doxygen processing
  !> \cond InterfaceBlock
  INTERFACE Convert2Cartesian_channel
     MODULE PROCEDURE Channel2Cartesian_coords,Channel2Cartesian_vectors
  END INTERFACE
  INTERFACE Convert2Curvilinear_channel
     MODULE PROCEDURE Cartesian2Channel_coords,Cartesian2Channel_vectors
  END INTERFACE
  !> \endcond
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: geometry_name = "channel geometry"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
      InitGeometry_channel, &
      ScaleFactors_channel, &
      Radius_channel, &
      PositionVector_channel, &
      Convert2Cartesian_channel, &
      Convert2Curvilinear_channel, &
      Channel2Cartesian_coords, &
      Cartesian2Channel_coords, &
      Channel2Cartesian_vectors, &
      Cartesian2Channel_vectors
  !--------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE InitGeometry_channel(this,gt,gp)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Geometry_TYP), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: gt
    REAL, INTENT(IN)    :: gp
    !------------------------------------------------------------------------!
    CALL InitGeometry(this,gt,geometry_name)
    CALL SetScale(this,gp)
  END SUBROUTINE InitGeometry_channel
    

  ELEMENTAL SUBROUTINE ScaleFactors_channel(a,eta,hxi,heta,hz)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: a,eta
    REAL, INTENT(OUT) :: hxi,heta,hz
    !------------------------------------------------------------------------!
    REAL              :: aa
    !------------------------------------------------------------------------!
    aa  = Atanh(a)
    hxi = 1.
    heta= aa/(a*COSH(eta*aa)**2)
    hz  = 1.
  END SUBROUTINE ScaleFactors_channel

  ELEMENTAL FUNCTION Radius_channel(a,xi,eta) RESULT(radius)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: a,xi,eta
    REAL :: radius
    REAL :: tmp
    !------------------------------------------------------------------------!
    tmp = COSH(Atanh(a)*eta)
    radius = SQRT((a*a+xi*xi+1.0)*tmp*tmp-1.0) / (a*tmp)
  END FUNCTION Radius_channel

  ELEMENTAL SUBROUTINE PositionVector_channel(a,xi,eta,rx,ry)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: a,xi,eta
    REAL, INTENT(OUT) :: rx,ry
    !------------------------------------------------------------------------!
    rx = xi
    ry = TANH(Atanh(a)*eta)/a
  END SUBROUTINE PositionVector_channel


  ! coordinate transformation
  ! channel -> cartesian
  ELEMENTAL SUBROUTINE Channel2Cartesian_coords(a,xi,eta,x,y)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: a,xi,eta
    REAL, INTENT(OUT) :: x,y
    !------------------------------------------------------------------------!
    x = xi
    y = TANH(eta*Atanh(a))/a
  END SUBROUTINE Channel2Cartesian_coords


  ! cartesian -> channel
  ELEMENTAL SUBROUTINE Cartesian2Channel_coords(a,x,y,xi,eta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x,y,a
    REAL, INTENT(OUT) :: xi,eta
    !------------------------------------------------------------------------!
    xi = x
    eta= Atanh(a*y) / Atanh(a)
  END SUBROUTINE Cartesian2Channel_coords

! vector transformations
  ELEMENTAL SUBROUTINE Channel2Cartesian_vectors(vxi,veta,vx,vy)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: vxi,veta
    REAL, INTENT(OUT) :: vx,vy
    !------------------------------------------------------------------------!
    vx = vxi
    vy = veta
  END SUBROUTINE Channel2Cartesian_vectors

  ELEMENTAL SUBROUTINE Cartesian2Channel_vectors(vx,vy,vxi,veta)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: vx,vy
    REAL, INTENT(OUT) :: vxi,veta
    !------------------------------------------------------------------------!
    vxi = vx
    veta = vy
  END SUBROUTINE Cartesian2Channel_vectors

END MODULE geometry_channel
