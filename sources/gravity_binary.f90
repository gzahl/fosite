!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: gravity_binary.f90                                                #
!#                                                                           #
!# Copyright (C) 2010-2014                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!# Björn Sperling   <sperling@astrophysik.uni-kiel.de>                       #
!# Anna Feiler      <afeiler@astrophysik.uni-kiel.de>                        #
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
!> \addtogroup gravity
!! - parameters of \link gravity_binary \endlink as key-values
!! \key{mass1,REAL,mass of primary component,1.0}
!! \key{mass2,REAL,mass of secondary component,1.0}
!! \key{excentricity,REAL,excentricity,0.0}
!! \key{semimayoraxis,REAL,semi major axis,1.0}
!! \key{softening1,REAL,softening parameter of primary component,0.0}
!! \key{softening2,REAL,softening parameter of secondary component,0.0}
!! \key{switchon1,REAL,soft switch on,-1.0}
!----------------------------------------------------------------------------!
!> \author Tobias Illenseer
!! \author Björn Sperling
!! \author Anna Feiler
!!
!! \brief source terms module for gravitational acceleration due to two
!! pointmasses
!!
!! \extends gravity_common
!! \ingroup gravity
!----------------------------------------------------------------------------!
MODULE gravity_binary
  USE gravity_common
  USE gravity_pointmass, ONLY: CalcDiskHeight_binary => CalcDiskHeight_pointmass
  USE physics_generic
  USE mesh_generic
  USE common_dict
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  PRIVATE
  CHARACTER(LEN=32), PARAMETER :: gravity_name = "binary point masses"
  !--------------------------------------------------------------------------!
  PUBLIC :: &
       ! types
       Gravity_TYP, &
       ! methods
       InitGravity_binary, &
       InfoGravity_binary, &
       GetAccelGravity_binary, &
       GetDiskHeight_binary, &
       GetinvDistanceCO_binary, &
       GetDistVector_binary,&
       CloseGravity_binary
  !--------------------------------------------------------------------------!

CONTAINS
 !> \public Constructor of gravity binary module
  !!
  !! This subroutine reads the necessary config data for the binary.
  !! It initializes the gravity type and various mesh data arrays. Some of those
  !! are marked for output.
  SUBROUTINE InitGravity_binary(this,Mesh,Physics,config)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this !< \param [in,out] this all gravity data
    TYPE(Mesh_TYP)    :: Mesh          !< \param [in] Mesh mesh type
    TYPE(Physics_TYP) :: Physics       !< \param [in] Physics physics type
    TYPE(Dict_TYP),POINTER :: config  !< \param [in,out] config sub-dictionary 
                                      !! with binary configuration data       
    !------------------------------------------------------------------------!
    INTEGER           :: gtype 
    INTEGER           :: err
    REAL              :: omega_min
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    !------------------------------------------------------------------------!
    CALL GetAttr(config, "gtype", gtype)
    CALL InitGravity(this,gtype,gravity_name)

    SELECT CASE(GetType(Physics))
    CASE(EULER2D, EULER2D_ISOTHERM, EULER2D_SGS, &
        EULER2D_ISOIAMT, EULER2D_IAMT, EULER2D_IAMROT, EULER2D_ISOIAMROT)
       ! do nothing
    CASE DEFAULT
       CALL Error(this,"InitGravity_binary", &
            "Physics not supported for binary gravity term")
    END SELECT

    ALLOCATE(this%accel(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%DIM), &
         this%height(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%omega2(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2), &
         this%omega(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX), &
         this%invr(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX),&
         this%invr_sec(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX),&
         this%r_prim(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2),&
         this%r_sec(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2),&
         STAT = err)
    IF (err.NE.0) CALL Error(this,"InitGravity_pointmass", "Unable allocate memory!")

    ! mass of primary component
    CALL RequireKey(config, "mass1", 1.0)
    CALL GetAttr(config, "mass1", this%mass)
    
    ! mass of secondary component
    CALL RequireKey(config, "mass2", 1.0)
    CALL GetAttr(config, "mass2", this%mass2)
    
    ! excentricity
    CALL RequireKey(config, "excentricity", 0.0)
    CALL GetAttr(config, "excentricity", this%excent)
    
    ! semi major axis
    CALL RequireKey(config, "semimayoraxis", 1.0)
    CALL GetAttr(config, "semimayoraxis", this%semaaxis)
    
    ! Softening parameter
    CALL RequireKey(config, "softening1", 0.0)
    CALL GetAttr(config, "softening1", this%eps1)

    CALL RequireKey(config, "softening2", 0.0)
    CALL GetAttr(config, "softening2", this%eps2)

    ! soft switch on
    CALL RequireKey(config, "switchon1", -1.0)
    CALL GetAttr(config, "switchon1", this%switchon)

    ! reset mass flux
    this%mdot = 0.0

    ! set period
    this%period = SQRT(4.*PI**2*this%semaaxis**3/(this%mass+this%mass2))

    ! initialize acceleration
    this%accel  = 0.0

    ! initialize speed of sound pointer
    this%cs => GetSoundSpeeds(Physics)

    ! reset omega and omega**2
    this%omega  = 0.0
    this%omega2 = 0.0

    ! set ghost cell data to one to avoid zero division
    this%omega(Mesh%IGMIN:Mesh%IMIN-1,:) = 1.0
    this%omega(Mesh%IMAX+1:Mesh%IGMAX,:) = 1.0
    this%omega(:,Mesh%JGMIN:Mesh%JMIN-1) = 1.0
    this%omega(:,Mesh%JMAX+1:Mesh%JGMAX) = 1.0

    ! reset start value for time variable, to guarantee omega^2 and
    ! disk height will be set to new values when the calculation starts
    this%time = -1.0
  END SUBROUTINE InitGravity_binary

!> \public write binary parameters to screen 
  SUBROUTINE InfoGravity_binary(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    !------------------------------------------------------------------------!
    CHARACTER(LEN=32) :: mass1_str,mass2_str,excent_str,sema_str,omega_str
    !------------------------------------------------------------------------!
    WRITE (mass1_str,'(ES8.2)') this%mass
    WRITE (mass2_str,'(ES8.2)') this%mass2
    WRITE (excent_str,'(ES8.2)') this%excent
    WRITE (sema_str,'(ES8.2)') this%semaaxis
    CALL Info(this,"            primary mass:      " // TRIM(mass1_str) // &
        ACHAR(10)//"            secondary mass:    " // TRIM(mass2_str) // &
        ACHAR(10)//"            excentricity:      " // TRIM(excent_str) // &
        ACHAR(10)//"            semi major axis:   " // TRIM(sema_str))
    IF(this%period.NE.0.0) THEN
        WRITE(omega_str,'(ES8.2)') 2.0*PI/this%period
        CALL Info(this,"            pattern speed:     " // TRIM(omega_str))
    END IF
  END SUBROUTINE InfoGravity_binary


  SUBROUTINE UpdateBinary(this,Mesh,Physics,pvar,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL              :: time
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,pvar,time
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! compute square of Keplerian velocities for updated time value
    CALL UpdateInverseDistances_binary(this,Mesh,Physics,time)
!CDIR IEXPAND
    ! local Keplerian velocity squared 
    ! with respect to PRIMARY star
    this%omega2(:,:,1) = Physics%Constants%GN*GetMass_binary(this,time) &
                          *this%invr(:,:)**3
    ! and SECONDARY star
    this%omega2(:,:,2) = Physics%Constants%GN*this%mass2 &
                          *this%invr_sec(:,:)**3
    
    ! there may be more than two components
    this%accel(:,:,:) = 0.

    ! initialize gravitational acceleration an Keplerian angular velocity
    ! curvilinear components of the gravitational acceleration
    ! -d Phi / dr = -r * omega^2 * e_r = -omega^2 * r_prim
!CDIR COLLAPSE
    DO j=Mesh%JMIN,Mesh%JMAX
!CDIR NODEP
      DO i=Mesh%IGMIN,Mesh%IGMAX
        this%accel(i,j,1:2) = -this%omega2(i,j,1)* this%r_prim(i,j,1:2)&
                              -this%omega2(i,j,2)* this%r_sec(i,j,1:2)
      END DO
    END DO

  END SUBROUTINE UpdateBinary

!> Update the inverse of the distances to the central stars
  SUBROUTINE UpdateInverseDistances_binary(this,Mesh,Physics,time)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    CALL UpdatePositions_binary(this,time)

    ! define position vector and inverse radius
    ! shifted point mass position:
    ! compute curvilinear components of shift vector
    this%r_prim(:,:,1) = this%binpos(1,1)
    this%r_prim(:,:,2) = this%binpos(2,1)
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,this%r_prim,this%r_prim)
    this%r_sec(:,:,1) = this%binpos(1,2)
    this%r_sec(:,:,2) = this%binpos(2,2)
    CALL Convert2Curvilinear(Mesh%geometry,Mesh%bcenter,this%r_sec,this%r_sec)
    ! subtract the result from the position vector:
    ! this gives you the curvilinear components of all vectors pointing 
    ! from the point mass to the bary center of any cell on the mesh
    this%r_prim(:,:,:) = Mesh%bposvec(:,:,:) - this%r_prim(:,:,:)
    this%r_sec(:,:,:) = Mesh%bposvec(:,:,:) - this%r_sec(:,:,:)
    ! compute the inverse of its absolute value
    this%invr(:,:) = 1.0 / SQRT(this%r_prim(:,:,1)**2+this%r_prim(:,:,2)**2)
    this%invr_sec(:,:) = 1.0 / SQRT(this%r_sec(:,:,1)**2+this%r_sec(:,:,2)**2)

  END SUBROUTINE UpdateInverseDistances_binary

!> Solves the Kepler equation to get the positions of the stars for a given time. 
!! A detailed explanation and a derivation of the relevant equations can 
!! be found in :
!!  Taff, L. G. :Celestial Mechanics  A computational guide for the 
!!    practitioner, John Wiley & Sons (1985), Chapter 2
  SUBROUTINE UpdatePositions_binary(this,time)
    USE roots
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    REAL              :: time
    !------------------------------------------------------------------------!
    REAL              :: E,cosE,r,phi,r1
    REAL              :: tau,excent
    REAL, DIMENSION(2):: plist
    !------------------------------------------------------------------------!
    INTENT(IN)        :: time
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    !> First we map the dimensionless orbital time into the interval \f$[0,2\pi]\f$
    !! to calculate the mean anomaly \f$ \tau \f$.
    !! Then we have to solve the Kepler equation:
    !! \f[
    !!    \tau= E-\epsilon\sin{E}
    !! \f]
    !! to calculate the eccentric anomaly \f$ E \f$. 
    !! \f$ \epsilon \f$ is the eccentricity of the orbit which is given as an 
    !! input parameter of the binary.
    !! We can solve the equation using the function \link roots::GetRoot_regfalsi\endlink.
    !! To do this, we have to know an Intervall that includes the correct solution 
    !! for \f$ E \f$. This should be as small as possible to speed up the calculation.
    !! Looking at the Kepler equation we can see, that if \f$ \tau\f$ is in 
    !! the interval \f$[0,\pi]\f$ \f$E\f$ also has a value in  \f$[0,\pi]\f$.
    !! The same is true for the other possible intervall of values for \f$ \tau \f$
    !! and \f$ E \f$ \f$]\pi,2\pi[\f$. 
    !! We can further constrain the Intervall for the correct solution for \f$ E \f$
    !! by looking at these separate cases, keeping in mind that the eccentricity 
    !! \f$ \epsilon \f$ can only have values of:
    !! \f[
    !!    0 \leq \epsilon \leq 1
    !! \f]
    !! Case 1 : \f$ E,\tau \in [0,\pi] \f$
    !! \f{eqnarray*}{
    !!    E &=& \tau + \epsilon\sin{E} \leq \tau + \epsilon \\
    !!    E &=& \tau + \epsilon\sin{E} \geq \tau  \\
    !! \f}
    !! Case 2 : \f$ E,\tau \in ]\pi,2\pi[ \f$
    !! \f{eqnarray*}{
    !!    E &=& \tau + \epsilon\sin{E} < \tau  \\
    !!    E &=& \tau + \epsilon\sin{E} > \tau - \epsilon \\
    !! \f}
    !! If we know the mean anomaly \f$ E \f$ we can calculate the distance
    !! between the stars
    !! \f[
    !!    r = a (1-\epsilon \cos{E}).
    !! \f]
    !! \f$ a \f$ is the semimajor axis of the binary. 
    !! We can also calculate the angle \f$ \phi \f$ between the semimajoraxis of 
    !! the ellipse and the line connecting the primary component to the center of mass
    !! \f[
    !!    \phi = \arctan\left(\sqrt{\frac{1+\epsilon}{1-\epsilon}} \tan{\frac{E}{2}} \right)
    !! \f] 
    !! This equation can be written as 
    !! \f[
    !!    \phi = \arctan\left(\pm\sqrt{\frac{1+\epsilon}{1-\epsilon} \frac{1+\cos{\frac{E}{2}}}{1-\cos{\frac{E}{2}}}} \right)
    !! \f] 
    !! to avoid one evaluation of \f$ \tan \f$.
    !! To choose the coorect sign in this expression we have to look at the sign of 
    !! \f$ \tan(\frac{E}{2}) \f$ which is positive for  \f$ E\in[0,\pi] \f$ and negative for \f$ E\in[\pi,2\pi]\f$.
    !! When we know the angle \f$ \phi \f$ and the distance between the stars \f$ r \f$
    !! we can calculate the positions of the stars \f$ \vec{r}_1,\vec{r}_2 \f$ using
    !! the equations:
    !! \f{eqnarray*}{
    !!    r_1 &=& \frac{m_2}{m_1+m_2} r\\
    !!    \vec{r}_1&=& -\vec{r}_2 \\
    !! \f}
    !! with the masses of the binary components \f$ m_1, m_2 \f$ and \f$ r_1=|\vec{r_1}|\f$.

    ! mean anomaly
    tau = 2*PI*MODULO(time,this%period)/this%period
    excent = this%excent

    plist(1)=excent
    plist(2)=tau
    ! solve the Kepler equation: E-tau-excent*sin(E) = 0
    !   i.e. compute the intersection of f(E) = E-tau and g(E) = excent*sin(E)
    IF ((tau.GE.0.0) .AND. (tau.LE.PI)) THEN
       ! intersection lies ABOVE abscissa
       E = GetRoot(func,MAX(0.0,tau),MIN(PI,excent+tau),plist)
    ELSE
       ! intersection lies BELOW abscissa      
       E = GetRoot(func,MAX(PI,tau-excent),MIN(2*PI,tau),plist)       
    END IF

    ! compute positions
    cosE = COS(E)
    ! distance to the origin of primary star position
    r = this%semaaxis * (1.0 - excent*cosE) 
    ! position angle of primary star
    !   phi1 = 2.0*ATAN(SQRT((1.0+excent)/(1.0-excent))*TAN(0.5*E))
    !   A little bit more complicated but this avoids one evaluation of TAN;
    !   uses TAN(x/2) = +/- SQRT((1-cos(x))/(1+cos(x)))
    phi = 2.0*ATAN(SIGN(SQRT(((1.+excent)*(1.-cosE)) &
         /((1.-excent)*(1.+cosE))),PI-E))
       
    ! cartesian coordinates of PRIMARY component
    r1               = this%mass2/(this%mass+this%mass2)*r
    this%binpos(1,1) = r1*COS(phi)
    this%binpos(2,1) = r1*SIN(phi)
    ! and SECONDARY COMPONENT
    this%binpos(:,2) = -this%binpos(:,1) * this%mass/this%mass2
  END SUBROUTINE UpdatePositions_binary

!> Update the pressure scale height for binary potential
  SUBROUTINE UpdateDiskHeight_binary(this,Mesh,Physics,time,pvar)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    REAL              :: sqrtgamma
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    ! compute disk height
    sqrtgamma = SQRT(Physics%gamma)
    this%omega(:,:) = SQRT(this%omega2(:,:,1) + this%omega2(:,:,2))

    DO j=Mesh%JMIN,Mesh%JMAX
      DO i=Mesh%IGMIN,Mesh%IGMAX
        this%height(i,j) &
          = CalcDiskHeight_binary(sqrtgamma,this%cs(i,j),this%omega(i,j))
      END DO
    END DO
  END SUBROUTINE UpdateDiskHeight_binary

!> \public Return the gravitational acceleration due to the binary
  FUNCTION GetAccelGravity_binary(this,Mesh,Physics,time,pvar) RESULT(ac)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL,DIMENSION(:,:,:), POINTER &
                      :: ac
    !------------------------------------------------------------------------!
    INTEGER           :: i,j
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar
    !------------------------------------------------------------------------!
    ! update Binary
    CALL UpdateBinary(this,Mesh,Physics,pvar,time)
    !> \todo Pointer?
    ac => this%accel
  END FUNCTION GetAccelGravity_binary

!> \public Returns the Pressure scale height for the gravitational potential of 
!! the binary
  FUNCTION GetDiskHeight_binary(this,Mesh,Physics,time,pvar) RESULT(height)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: height
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar
    !------------------------------------------------------------------------!
    ! update DiskHeight -- other things are up to date due to gravity_generic
    CALL UpdateDiskHeight_binary(this,Mesh,Physics,time,pvar)
    !> \todo Pointer?
    height(:,:) = this%height(:,:)
  END FUNCTION GetDiskHeight_binary

!> \public Returns the inverse of the distances between the stars of the binary and 
!! the center of mass
  FUNCTION GetinvDistanceCO_binary(this,Mesh,Physics,time,pvar,component) RESULT(invdis)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX) &
                      :: invdis
    INTEGER, OPTIONAL :: component
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar,component
    !------------------------------------------------------------------------!
    ! update Binary
    CALL UpdateBinary(this,Mesh,Physics,pvar,time)
    !> \todo Pointer?
    IF (present(component) .AND. (component .EQ. 2)) THEN
      invdis(:,:) = this%invr_sec(:,:)
    ELSE 
      invdis(:,:) = this%invr(:,:)
    END IF
  END FUNCTION GetinvDistanceCO_binary

!> \public Returns an array of vectors between the stars and the cell bary 
!! centers 
  FUNCTION GetDistVector_binary(this,Mesh,Physics,time,pvar,component) RESULT(r)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP), POINTER :: this
    TYPE(Mesh_TYP)    :: Mesh
    TYPE(Physics_TYP) :: Physics
    REAL              :: time
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,Physics%vnum) &
                      :: pvar
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX,2) &
                      :: r
    INTEGER, OPTIONAL :: component
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,time,pvar,component
    !------------------------------------------------------------------------!
    ! update Binary
    CALL UpdateBinary(this,Mesh,Physics,pvar,time)
    !> \todo Pointer?
    IF (present(component) .AND. (component .EQ. 2)) THEN
      r(:,:,:) = this%r_sec(:,:,:)
    ELSE 
      r(:,:,:) = this%r_prim(:,:,:)
    END IF
  END FUNCTION GetDistVector_binary

  !> get mass of the primary component, eventually scaled by switch on factor
  FUNCTION GetMass_binary(this,time) RESULT(mass)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    REAL              :: time
    !------------------------------------------------------------------------!
    REAL              :: mass
    !------------------------------------------------------------------------!
    INTENT(IN)        :: this, time
    !------------------------------------------------------------------------!
    IF(time.LE.this%switchon) THEN
        mass = this%mass * SIN(0.5*PI*time/this%switchon)**2
    ELSE
        mass = this%mass
    END IF
  END FUNCTION GetMass_binary

!> \public Closes the binary source term
  SUBROUTINE CloseGravity_binary(this)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    TYPE(Gravity_TYP) :: this
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: this
    !------------------------------------------------------------------------!
    DEALLOCATE(this%accel,this%omega2,this%omega,&
               this%invr,this%invr_sec,this%r_prim,this%r_sec)
    CALL CloseGravity(this)
  END SUBROUTINE CloseGravity_binary


!> find root of the Kepler equation 
  SUBROUTINE func(x,fx,plist)
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL, INTENT(IN)  :: x
    REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: plist 
    !plist[1]=excent, plist[2]=tau
    REAL, INTENT(OUT) :: fx
    !------------------------------------------------------------------------!
    fx  = x - plist(1)*SIN(x) - plist(2)
  END SUBROUTINE func

END MODULE gravity_binary
