#!/bin/env python
#############################################################################
#                                                                           #
# fosite - 2D hydrodynamical simulation program                             #
# program file: baldr.py                                                    #
#                                                                           #
# Copyright (C) 2011-2012                                                   #
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

# some import magic, because the other modules are not in PYTHONPATH
# look here: http://stackoverflow.com/questions/279237/python-import-a-module-from-a-folder
import os, sys, inspect
# cmd_folder = os.path.dirname(os.path.abspath(__file__)) # DO NOT USE __file__ !!!
# __file__ fails if script is called in different ways on Windows
# __file__ fails if someone does os.chdir() before
# sys.argv[0] also fails because it doesn't not always contains the path
cmd_folder = os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0])
if cmd_folder not in sys.path:
	sys.path.insert(0, cmd_folder)

import baldr_impl as impl
from baldr_data import data
import signal, sys
import numpy as np


class baldr(object):
	def __init__(self):
		self.impl = impl
		impl.baldr.initbaldr()
		signal.signal(signal.SIGINT, self._handle_sigint)

	def _initPhysics(self, physics):
		impl.baldr.initphysicsbaldr(**physics)
		self.VNUM = impl.baldr.vnum
		self.DENSITY = impl.baldr.density
		self.PRESSURE = impl.baldr.pressure
		self.ENERGY = impl.baldr.energy
		self.XVELOCITY = impl.baldr.xvelocity
		self.XMOMENTUM = impl.baldr.xmomentum
		self.YVELOCITY = impl.baldr.yvelocity
		self.YMOMENTUM = impl.baldr.ymomentum
		self.ZVELOCITY = impl.baldr.zvelocity
		self.ZMOMENTUM = impl.baldr.zmomentum
		self.GN = data(self.impl, 'gn')

	def _initFluxes(self, fluxes):
		impl.baldr.initfluxesbaldr(**fluxes)
	
	def _initMesh(self, mesh):
		impl.baldr.initmeshbaldr(**mesh)
		self.imin = impl.baldr.imin
		self.imax = impl.baldr.imax
		self.jmin = impl.baldr.jmin
		self.jmax = impl.baldr.jmax
		self.igmin = impl.baldr.igmin
		self.igmax = impl.baldr.igmax
		self.jgmin = impl.baldr.jgmin
		self.jgmax = impl.baldr.jgmax
		#self.inum = impl.baldr.inum
		#self.jnum = impl.baldr.jnum
		self.inum = data(self.impl,'inum')
		self.jnum = data(self.impl,'jnum')
		self.gnum = data(self.impl,'gnum')
		self.bhx = data(self.impl,'bhx')
		self.bhy = data(self.impl,'bhy')
		self.bhz = data(self.impl,'bhz')
		self.fhx = data(self.impl,'fhx')
		self.fhy = data(self.impl,'fhy')
		self.fhz = data(self.impl,'fhz')
		self.radius = impl.baldr.radius
		self.curv_coords = impl.baldr.curv_coords
		self.cart_coords = data(self.impl,'cart_coords')
		self.dlx = impl.baldr.getdlx(mesh['inum']+4,mesh['jnum']+4)
		self.dly = impl.baldr.getdlx(mesh['inum']+4,mesh['jnum']+4)
		self.volume = data(self.impl,'volume')
		#self.dlx = impl.baldr.dlx
		#self.dly = impl.baldr.dly

	def _initBoundary(self, boundary):
		impl.baldr.initboundarybaldr(**boundary)

	def _initSources(self, sources):
		outbound = sources.pop('outbound',None)
		impl.baldr.initsourcesbaldr(**sources)
		if outbound==False:
			print "No Outbound!"
			impl.baldr.nooutboundbaldr()

	def _initTimedisc(self, timedisc):
		impl.baldr.inittimediscbaldr(**timedisc)
		self.pvar = data(self.impl,'pvar')
		self.pvar[...] = 0.0
		self.time = data(self.impl,'time')
		self.dt = data(self.impl,'dt')
		self.dtmin = data(self.impl,'dtmin')
		self.iter = data(self.impl,'iter')

	def _initDatafile(self, datafile):
		if(datafile!=None):
			impl.baldr.initdatafilebaldr(**datafile)

	def _initLogfile(self, logfile):
		if(logfile!=None):
			impl.baldr.initlogfilebaldr(**logfile)

	def setup(self, config):
		self._initMesh(config['mesh'])
		self._initPhysics(config['physics'])
		self._initFluxes(config['fluxes'])
		self._initBoundary(config['boundary'])
		if config['sources']!=None:
			[ self._initSources(src) for src in config['sources'] ]
		self._initTimedisc(config['timedisc'])
		self._initDatafile(config['datafile'])
		self._initLogfile(config['logfile'])

	def step(self):
		return impl.baldr.stepbaldr()

	def run(self, callback=None):
		#impl.baldr.run()
		self.running = True
		while (not self.step()) and self.running:
			if(callback!=None):
				callback(self)
		self.running = False

	def convert2Conservative(self):
		impl.baldr.convert2conservativebaldr()

	def adddata(self, str):
		return data(self.impl, str)
	
	def close(self):
		impl.baldr.closebaldr()

	def _handle_sigint(self, signum, frame):
		self.running = False
		#sys.exit(0)

	def __del__(self):
		self.close()


