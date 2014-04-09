import numpy as np

class data(np.ndarray):
	def __new__(cls, impl, name):
		# Input array is an already formed ndarray instance
		# We first cast to be our class type
		sh = getattr(impl.baldr, 'shape'+name)()
		if(np.all(sh==np.array([1]))): sh=()
		obj = np.asarray(getattr(impl.baldr, 'get'+name)(*sh)).view(cls)
		# add the new attribute to the created instance
		obj._impl = impl
		obj._name = name
		# Finally, we must return the newly created object:
		return obj

	def __array_finalize__(self, obj):
		# see InfoArray.__array_finalize__ for comments
		if obj is None: return
		self._name = getattr(obj, '_name', None)
		self._impl = getattr(obj, '_impl', None)

	def load(this):
		this[...] = data(this._impl,this._name)
		return this
	
	def save(this):
		getattr(this._impl.baldr, 'set'+this._name)(this)
		return this
