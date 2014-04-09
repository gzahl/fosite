# flags: Have to be redefined for export to python

class physics:
	class problem:
		# flags from physics/physics_generic.f90
		EULER2D          = 1
		EULER2D_ISOTHERM = 2
		EULER3D_ROTSYM   = 3
		EULER3D_ROTAMT   = 4
		EULER3D_ROTSYMSGS= 5
	
	class units:
		# flags from physics/constants_generic.f90
		SI          = 1
		CGS         = 2
		GEOMETRICAL = 3

class mesh:
	class geometry:
		# flags from mesh/geometry_generic.f90
		CARTESIAN         = 1
		POLAR             = 20
		LOGPOLAR          = 21
		TANPOLAR          = 22
		SINHPOLAR         = 23
		SINHTANHPOLAR	  = 24
		BIPOLAR           = 29
		CYLINDRICAL       = 30
		TANCYLINDRICAL    = 31
		SPHERICAL         = 40
		SINHSPHERICAL     = 41
		OBLATE_SPHEROIDAL = 50
	class meshtype: 
		# flags from fluxes/fluxes_generic.f90
		MIDPOINT     = 1
		TRAPEZOIDAL  = 2
  
class fluxes:
	class order:
		# flags from fluxes/reconstruction_generic.f90
		CONSTANT = 1
		LINEAR   = 2

	class limiter:
		# flags from fluxes/reconstruction_linear.f90
		MINMOD   = 1
		MONOCENT = 2
		SWEBY    = 3
		SUPERBEE = 4
		OSPRE    = 5
		PP       = 6

	class variables:
		# flags from common/reconstruction_common.f90
		PRIMITIVE    = True
		CONSERVATIVE = False

class boundary:
	# flags from boundary/boundary_generic.f90
	NO_GRADIENTS  = 1
	PERIODIC      = 2
	REFLECTING    = 3
	AXIS          = 4
	FOLDED        = 5
	FIXED         = 6
	EXTRAPOLATION = 7
	NOH2D         = 8
	NOH3D         = 9
	NOSLIP        = 10
	CUSTOM        = 11
	FARFIELD      = 12
	ABSORBING     = 13

	WEST		  = 1
	EAST		  = 2
	SOUTH		  = 3
	NORTH		  = 4

	CUSTOM_NOGRAD   = 1
	CUSTOM_PERIOD   = 2
	CUSTOM_REFLECT  = 3
	CUSTOM_REFLNEG  = 4
	CUSTOM_EXTRAPOL = 5
	CUSTOM_FIXED    = 6

class file:
	class fileformat:
		# flags from io/fileio_generic.f90
		BINARY  = 1
		GNUPLOT = 2
		NETCDF  = 3
		VTK     = 4
		NPY     = 5

class sources:
	class stype:
		# flags from sources/sources_generic.f90
		POINTMASS        = 1
		DISK_THOMSON     = 2
		VISCOSITY        = 3
		C_ACCEL          = 4
		COOLING          = 5
		POISSON          = 6
		ROTATING_FRAME   = 20
		POINTMASS_BINARY = 21
		GRAVMONOPOL      = 22
		SGS              = 23
		DISK_COOLING     = 24

	class solver:
		# flags from sources/poisson_generic.f90
		MULTIGRID    = 1
		#MULTIPOL     = 2,
		#FFT          = 3,
		SPECTRAL     = 4

	class potential:
		# flags from sources/sources_pointmass.f90
		NEWTON = 1
		WIITA  = 2

	class vismodel:
		# flags from sources/sources_viscosity.f90
		MOLECULAR = 1     # constant viscosity
		ALPHA     = 2     # Shakura-Sunyaev prescription
		BETA      = 3     # Duschl prescription
		PRINGLE   = 4     # constant kinematic viscosity
		ALPHA_ALT = 5     # alternative Shakura-Sunyaev

class timedisc:
	class method:
		# flags from timedisc/timedisc_generic.f90
		MODIFIED_EULER = 1

