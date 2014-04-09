#!/usr/bin/env python
#############################################################################
#                                                                           #
# fosite - 2D hydrodynamical simulation program                             #
# module: gauss2d.py                                                        #
#                                                                           #
# Copyright (C) 2006-2012                                                   #
# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
# Manuel Jung <mjung@astrophysik.uni-kiel.de>                               #
#                                                                           #
# This program is free software; you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 2 of the License, or (at     #
# your option) any later version.                                           #
#                                                                           #
# This program is distributed in the hope that it will be useful, but       #
# WITHOUT ANY WARRANTY; without even the implied warranty of                #
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, GOOD TITLE or        #
# NON INFRINGEMENT.  See the GNU General Public License for more            #
# details.                                                                  #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with this program; if not, write to the Free Software               #
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
#                                                                           #
#############################################################################

#----------------------------------------------------------------------------#
# Program and data initialization for 2D Gaussian pressure pulse
#----------------------------------------------------------------------------#

import numpy as np
from baldr import *
import baldr_flags as f

b = baldr()

# simulation parameter
TSIM   = 0.3         # simulation time
GAMMA  = 1.4         # ratio of specific heats
CSISO  = 0.0         # if .ne. 0.0 -> isothermal simulation
                     #   with CSISO as sound speed
# initial condition (dimensionless units)
RHO0   = 1.0         # ambient density
P0     = 1.0         # ambient pressure
AMP    = 1.0         # amplitude of the pulse
PWIDTH = 0.06        # half width of the Gaussian
ETA    = 0.0         # dynamic viscosity (0.0 disables)
# mesh settings
MGEO = f.mesh.geometry.CARTESIAN    # geometry of the mesh
#MGEO = POLAR
#MGEO = LOGPOLAR
#MGEO = TANPOLAR
#MGEO = SINHPOLAR
XRES = 100          # resolution
YRES = 100
RMIN = 1.0E-2       # inner radius for polar grids
RMAX = 0.3          # outer radius
GPAR = 0.2          # geometry scaling parameter
# output file parameter
ONUM = 10           # number of output time steps
ODIR = "./"         # output directory
OFNAME = 'gauss2d'  # output data file name


#! mesh settings
if(MGEO == f.mesh.geometry.CARTESIAN):
	mesh = dict(
		meshtype = f.mesh.meshtype.MIDPOINT,
		geometry = MGEO,
		inum = XRES,        # resolution in x and            #
		jnum = YRES,        #   y direction                  #
		gparam = GPAR,
		xmin = -0.5,
    	xmax = 0.5,
	    ymin = -0.5,
    	ymax = 0.5
	)
	boundary = dict(
  		western = f.boundary.NO_GRADIENTS,
	    eastern = f.boundary.NO_GRADIENTS,
       	southern = f.boundary.NO_GRADIENTS,
       	northern = f.boundary.NO_GRADIENTS
	)
elif(MGEO == POLAR):
	mesh = dict(
		meshtype = f.mesh.meshtype.MIDPOINT,
		geometry = MGEO,
		inum = XRES,        # resolution in x and            #
		jnum = YRES,        #   y direction                  #
		gparam = GPAR,
		xmin = RMIN,
		xmax = 0.5*SQRT(2.0),
		ymin = 0.0,
		ymax = 2*PI
	)
	boundary = dict(
  		western = f.boundary.NO_GRADIENTS,
	    eastern = f.boundary.NO_GRADIENTS,
		southern = f.boundary.PERIODIC,
		northern = f.boundary.PERIODIC,
	)
elif(MGEO == LOGPOLAR):
	mesh = dict(
		meshtype = f.mesh.meshtype.MIDPOINT,
		geometry = MGEO,
		inum = XRES,        # resolution in x and            #
		jnum = YRES,        #   y direction                  #
		gparam = GPAR,
		xmin = LOG(RMIN/GPAR),
		xmax = LOG(0.5*SQRT(2.0)/GPAR),
		ymin = 0.0,
		ymax = 2*PI
	)
	boundary = dict(
  		western = f.boundary.NO_GRADIENTS,
	    eastern = f.boundary.NO_GRADIENTS,
		southern = f.boundary.PERIODIC,
		northern = f.boundary.PERIODIC,
	)
elif(MGEO == TANPOLAR):
	mesh = dict(
		meshtype = f.mesh.meshtype.MIDPOINT,
		geometry = MGEO,
		inum = XRES,        # resolution in x and            #
		jnum = YRES,        #   y direction                  #
		gparam = GPAR,
		xmin = ATAN(RMIN/GPAR),
		xmax = ATAN(0.5*SQRT(2.0)/GPAR),
		ymin = 0.0,
		ymax = 2*PI
	)
	boundary = dict(
  		western = f.boundary.NO_GRADIENTS,
	    eastern = f.boundary.NO_GRADIENTS,
		southern = f.boundary.PERIODIC,
		northern = f.boundary.PERIODIC,
	)
elif(MGEO == SINHPOLAR):
	x1 = RMIN/GPAR                # temporary
	x2 = 0.5*SQRT(2.0)/GPAR       # temporary
	mesh = dict(
		meshtype = f.mesh.meshtype.MIDPOINT,
		geometry = MGEO,
		inum = XRES,        # resolution in x and            #
		jnum = YRES,        #   y direction                  #
		gparam = GPAR,
		xmin = LOG(x1+SQRT(1.0+x1*x1)),  # = ASINH(RMIN/sc)
		xmax = LOG(x2+SQRT(1.0+x2*x2)),  # = ASINH(RMAX/sc)
		ymin = 0.0,
		ymax = 2*PI
	)
	boundary = dict(
  		western = f.boundary.NO_GRADIENTS,
	    eastern = f.boundary.NO_GRADIENTS,
		southern = f.boundary.PERIODIC,
		northern = f.boundary.PERIODIC,
	)
else:
	print "geometry should be one of cartesian,polar,logpolar,tanpolar,sinhpolar or bipolar"

# physics settings
if(CSISO > 1.0E-10): # 1E-10==TINY
	physics = dict(
		problem = f.physics.problem.EULER2D_ISOTHERM,
		cs      = CSISO)                       # isothermal sound speed  #
else:
	physics = dict(
		problem   = f.physics.problem.EULER2D, 
		gamma     = GAMMA,             # ratio of specific heats        #
		dpmax     = 1.0)               # for advanced time step control #

# flux calculation and reconstruction method
fluxes = dict(
	order     = f.fluxes.order.LINEAR, 
	variables = f.fluxes.variables.CONSERVATIVE,  # vars. to use for reconstruction#
	limiter   = f.fluxes.limiter.MONOCENT,        # one of: minmod, monocent,...   #
	theta     = 1.2)                              # optional parameter for limiter #

# viscosity source term
viscosity = dict(
	stype    = f.sources.stype.VISCOSITY, 
	vismodel = f.sources.vismodel.MOLECULAR, 
	dynconst = ETA)

if(ETA > 1.0E-10):
	sources = [viscosity]
else:
	sources = list()

# time discretization settings
timedisc = dict(
	method   = f.timedisc.method.MODIFIED_EULER, 
	order    = 3, 
	cfl      = 0.4, 
	stoptime = TSIM, 
	tol_rel  = 0.01, 
	dtlimit  = 1.0E-5, 
	maxiter  = 1000000)

# initialize log input/output
logfile = dict(
	fileformat = f.file.fileformat.BINARY,
	filename   = ODIR + OFNAME + 'log', 
	filecycles = 1)

# initialize data input/output
datafile = dict(
#	fileformat = f.file.fileformat.VTK,
	fileformat = f.file.fileformat.GNUPLOT, filecycles = 0,
	filename   = ODIR + OFNAME,
	count      = ONUM)


config = dict(
	physics = physics,
	fluxes = fluxes,
	mesh = mesh,
	boundary = boundary,
	sources = sources,
	timedisc = timedisc,
	datafile = datafile,
#	logfile = logfile,
	logfile = None
)

b.setup(config)

# center of the pressure pulse (give in cartesian coordinates)
x0, y0 = 0.0, 0.0

# velocities
b.pvar[:,:,b.XVELOCITY] = 0.0
b.pvar[:,:,b.YVELOCITY] = 0.0

if(physics['problem']==f.physics.problem.EULER2D):
	# non-isothermal setup with constant density and pressure pulse
	b.pvar[:,:,b.DENSITY] = RHO0
	b.pvar[:,:,b.PRESSURE] = P0 + AMP*np.exp(-np.log(2.0) * \
		((b.cart_coords[:,:,0]-x0)**2+(b.cart_coords[:,:,1]-y0)**2)/PWIDTH**2)
elif(physics['problem']==f.physics.problem.EULER2D):
	# in isothermal configurations the pressure is proportional to
	# the density; thus the pulse is applied to the density
	b.pvar[:,:,b.DENSITY] = RHO0 + AMP*np.exp(-np.log(2.0) * \
		((b.cart_coords[:,:,0]-x0)**2+(b.cart_coords[:,:,1]-y0)**2)/PWIDTH**2)

b.pvar.save()
b.convert2Conservative()
print " DATA-----> initial condition: 2D gaussian pressure pulse"

b.run()

b.close()
