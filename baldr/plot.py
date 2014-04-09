#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.animation as animation
import os.path
import os
import glob
import sys

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']))
#rc('text', usetex=True)

class baldr_plot:
	def __init__(self, filename, var):
		filename = os.path.abspath(args.filename)
		self.filename = filename
		self.filetemplate = filename.replace('0000','%04d')
		self.var = var
		self.allvars = []
		self.frame = -1
		self.logscale = False
		self.InitCoordinates()
		self.InitVars()
		
		self.fig = plt.figure()

		self.fig.canvas.mpl_connect('key_press_event', self.press)

	def FtoC(self, a):
		shape = a.shape
		a = a.reshape(shape).T
		return a

	def getFilename(self,n):
		if(self.filename.find('0000')==-1):
			f = self.filename
		else:
			f = self.filetemplate % n
		return f
		
	def InitCoordinates(self):
		f = self.getFilename(0)
		root, ext = os.path.splitext(f)
		if(ext=='.nc'):
			from netCDF4 import Dataset
			d = Dataset(f,'r')
			c = self.FtoC(d.variables['corners'][:,:,:,:])

		elif(ext=='.h5'):
			import h5py
			d = h5py.File(f, 'r')
			c = d['corners'][:,:,:,:]
		
		self.dim = c[:,:,0,0].shape
		cdim = tuple((i+1 for i in self.dim))
		cx = np.empty(cdim)
		cy = np.empty(cdim)

		cx[:-1,:-1] = c[:,:,0,0]
		cy[:-1,:-1] = c[:,:,0,1]
		cx[-1,:-1] = c[-1,:,1,0]
		cy[-1,:-1] = c[-1,:,1,1]
		cx[:-1,-1] = c[:,-1,2,0]
		cy[:-1,-1] = c[:,-1,2,1]
		cx[-1,-1] = c[-1,-1,3,0]
		cy[-1,-1] = c[-1,-1,3,1]

		self.cx, self.cy = cx, cy
		if (ext=='.nc'):
			if(self.dim[0]==1 or self.dim[1]==1):
				self.bc = self.FtoC(d.variables['bary_centers'][:,:])
			else:
				self.bc = self.FtoC(d.variables['bary_centers'][:,:,:])
			d.close()

	def InitVars(self):
		f = self.getFilename(0)
		#print "file= ",f
		if(os.path.exists(f)):
			root, ext = os.path.splitext(f)
			if(ext=='.nc'):
				from netCDF4 import Dataset
				d = Dataset(f,'r')
				for varname in d.variables.keys():
					val = d.variables[varname]
					for attr in val.ncattrs():
						if(attr=='field'):
							self.allvars.append(varname)
				d.close()
			elif(ext=='.h5'):
				import h5py
				d = h5py.File(f, 'r')
				print "Not supported yet!"
		#print "variables: "self.allvars
		
	def emptyVar(self):
		return 0.0, np.ma.masked_all(self.dim)


	def vrange(self, val):
		vmin, vmax = val.min(), val.max()
		#if(np.abs(np.log10(vmax)-np.log10(vmin)) < 1.0):
		#	vmax = 10**(np.ceil(np.log10(vmin)+1.0))
		#print vmin, vmax
		return vmin,vmax


	def meshPlot(self, frame, v=(None,None)):
		if frame==None:
			frame = -1
			self.time, self.val = self.emptyVar()
		else:
			self.time, self.val = self.getVar(frame)
		#self.mesh = plt.pcolormesh(self.cx,self.cy,self.val,cmap=plt.cm.jet,norm=LogNorm(*self.vrange(self.val)))
		self.mesh = plt.pcolormesh(self.cx,self.cy,self.val,cmap=plt.cm.jet)
		#mesh = plt.pcolor(self.cx,self.cy,val)
		ax = self.mesh.get_axes()
		ax.set_aspect('equal')
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		#ax.set_title('n = %i, t = %.3e' % (frame, time))
		#self.fig.canvas.draw()
		return self.mesh

	def averagePlot(self, frame, axis=0):
		if frame==None:
			frame = -1
			self.time, self.val = self.emptyVar()
		else:
			self.time, self.val = self.getVar(frame)
		
		tinitial, initial = self.getVar(0)
		self.val = np.mean(self.val/initial, axis)
		if(axis==0):
			xaxis=1
		else:
			xaxis=0
		x = np.mean(self.bc[:,:,xaxis], axis)
		line = plt.plot(x,self.val)[0]
		ax = line.get_axes()
		ax.set_xlabel('axis %i' % axis)
		ax.set_ylabel('averaged over axis %i' % axis)
		return line

	def getVar(self, frame):
		self.frame = frame
		#print "frame: %i" % self.frame
		#f = self.filetemplate % frame
		f = self.getFilename(frame)
		if(os.path.exists(f)):
			#print "file=",f
			root, ext = os.path.splitext(f)
			if(ext=='.nc'):
				from netCDF4 import Dataset
				d = Dataset(f,'r')
				if(self.filename.find('0000')==-1):
					#multiple frames in one file
					time = d.variables['time'][frame]
					if(self.dim[0]==1 or self.dim[1]==1):
						#1D data
						data = d.variables[self.var][:][frame]
						data = data.reshape(self.dim)
					else:
						data = self.FtoC(d.variables[self.var][:,:][frame])
				else:
					#each frame in a different filename
					time = d.variables['time'][0]
					if(self.dim[0]==1 or self.dim[1]==1):
						#1D data
						data = d.variables[self.var][:]
						data = data.reshape(self.dim)
					else:
						#2D data
						data = self.FtoC(d.variables[self.var][:,:])
				#time, data = d.variables['time'][0], d.variables[self.var][:,:]
				d.close()
			elif(ext=='.h5'):
				import h5py
				d = h5py.File(f, 'r')
				timedisc = d['timedisc']
				time, data = 0.0, timedisc[self.var]
			if(self.logscale):
				return time, np.log10(data)
			else:
				return time, data
		else:
			print "file '%s' does not exist!" % f
			return self.emptyVar()
			

	def update(self, i):
		self.time, self.val = self.getVar(i)
		self.mesh.set_array(self.val.ravel())
		ax = self.mesh.get_axes()
		ax.set_title('n = %i, t = %.3e' % (self.frame, self.time))
		#disabled, because slow for animate()
		#self.fig.canvas.draw()
		return self.mesh,

	def plot(self, frame):
		self.fig.clf()
		frame = self.checkRange(frame)
		self.mesh = self.meshPlot(frame)
		self.mesh.set_figure(self.fig)
		ax = self.mesh.get_axes()
		ax.set_title('n = %i, t = %.3e' % (frame, self.time))
		self.makeCbar()
		return self.mesh


	def makeCbar(self):
		self.cbar = self.fig.colorbar(self.mesh)
		self.cbar.set_label(self.var)
		#self.cbar.ax.minorticks_on()
		# do not use offset
		self.cbar.formatter.set_useOffset(False)
		self.cbar.update_ticks()

		# to use scaling factor
		self.cbar.formatter.set_scientific(True)
		self.cbar.formatter.set_powerlimits((0,0))

		self.cbar.update_ticks()
		return self.cbar

	def updateCbar(self):
		self.mesh.set_clim(*self.vrange(self.mesh.get_array()))

	def plot1d(self, frame, axis=0):
		self.fig.clf()
		frame = self.checkRange(frame)
		line = self.averagePlot(frame, axis)
		line.set_figure(self.fig)
		ax = line.get_axes()
		ax.set_title('n = %i, t = %.3e' % (frame, self.time))


	def checkRange(self, frames):
		def exists(n):
			f = self.getFilename(n)
			return os.path.exists(f)
		
		def ReplaceNegative(i):
			if (i<0):
				if(self.filename.find('0000')!=-1):
					#multiple data files
					return sorted([ int(j[-7:-3]) for j in glob.iglob(self.filename.replace('0000','*'))])[i]
				else:
					#single data file
					root, ext = os.path.splitext(self.filename)
					if(ext=='.nc'):
						from netCDF4 import Dataset
						d = Dataset(self.filename,'r')
						#time, data = d.variables['time'][0], d.variables[self.var][:,:]
						times = d.variables['time']
						k = range(times.shape[0])[i]
						d.close()
					elif(ext=='.h5'):
						import h5py
						d = h5py.File(self.filename, 'r')
						#timedisc = d['timedisc']
						#time, data = 0.0, timedisc[self.var]
						k=i
					return k
			else:
				return i
		
		try:
			#tuple:
			a, b = (ReplaceNegative(i) for i in frames)
			frames = np.arange(a, b+1)

		except TypeError:
			# integer
			frames = np.array([ReplaceNegative(frames)])
		
		e = np.vectorize(exists)
		mask = e(frames)
		return frames[mask]

	def press(self, event):
		#print 'press', event.key
		maxframe = self.checkRange(-1)
		newframe = self.frame
		if(event.key=='right'):
			newframe += 1
		elif(event.key=='up'):
			newframe += 10
		elif(event.key=='left'):
			newframe -= 1
		elif(event.key=='down'):
			newframe -= 10
		elif(event.key=='home'):
			newframe = 0
		elif(event.key=='end'):
			newframe = maxframe
		elif(event.key=='pageup'):
			i = self.allvars.index(self.var)+1
			self.var = self.allvars[i % len(self.allvars)]
			self.cbar.set_label(self.var)
		elif(event.key=='pagedown'):
			i = self.allvars.index(self.var)-1
			self.var = self.allvars[i % len(self.allvars)]
			self.cbar.set_label(self.var)
		elif(event.key=='o'):
			self.logscale = not self.logscale
		elif(event.key=='backspace'):
			self.updateCbar()
		elif(event.key=='q'):
			sys.exit(0)

		if(newframe>maxframe):
			newframe = maxframe
		elif(newframe<0):
			newframe = 0

		self.update(newframe)
		if(event.key==';'):
			self.updateCbar()

		self.fig.canvas.draw()


	def animate(self, i):
		self.plot(i[-1])
		ani = animation.FuncAnimation(self.fig, 
			 						  #lambda x: (self.meshPlot(x),),
									  self.update,
									  self.checkRange(i),
									  interval=10,
									  #init_func=lambda: (self.meshPlot(None),),
									  blit=True)
		return ani


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='baldr plot routines')
	parser.add_argument('action', metavar='action', help="[plot, animate]")
	parser.add_argument('filename', metavar='filename', help="first output file")
	parser.add_argument('--var','-v',default='density')
	parser.add_argument('--range','-r',default='0:-1')
	parser.add_argument('--axis','-a',type=int,default=1)

	args = parser.parse_args()
	b = baldr_plot(args.filename,args.var)
	argrange = [int(i) for i in args.range.split(':')]
	valid_action = True
	if(args.action=='plot'):
		b.plot(argrange[len(argrange)-1])
	elif(args.action=='animate'):
		ani = b.animate(argrange)
	elif(args.action=='plot1d'):
		b.plot1d(argrange[len(argrange)-1],args.axis)
	else:
		valid_action = False
		print "Enter a valid action!"
		parser.print_help()

	if(valid_action):
		print "Basic keyboard control:"
		print "  left/right   previous/next frame"
		print "  up/down      move 10 frames forward/back"
		print "  PgUp/PgDn    next/previous variable"
		print "  backspace    rescale color bar"
		print "  o            toggle logarithmic scale"
		print "  q            exit program"
		plt.show()


