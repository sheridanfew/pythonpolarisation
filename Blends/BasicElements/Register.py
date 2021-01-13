from Singleton import Singleton
import inspect
import re

"""RegisteredObject is a base class for an object which only ever exists as a reference 
(the Index) to the instance itself, which is stored in a dict (the Register).
Whenever any index is deleted, so is the underlying instance.  Each Register is a Singleton, 
and is itself registered by the Singleton RegisterOfRegisters. 

An easy way to access all objects of a given type with the GetRegister('class type') method,
which returns an iterable list for the whole Register.

This is a useful base class for when an object needs to be stored in several collections, 
and modification or deletion of an object in any collection must automatically lead to 
updating or deleting of all objects.

The class is also an implementation of the Register design pattern, allowing access of all 
objects of a given type globally, but without using a global variable.

"""
class Index():
	"""The index is an object that allows to control the deletion of objects from the Register.

	Atom.__new__ is going to load an atoms onto the register and return an index.
	Index() will return the instance of the Atom from the register!
		
	"""
	def __init__(self, parent,n):
		self._parent = parent
		self._n = n

	def __del__(self): 
		"""the delete method removes the instance from register """
		if not self._parent is None:
			if self._n in self._parent:
				del self._parent[self._n]

	def __call__(self):
		""" calling and index a() will return the actual instance """
		if self._n in self._parent:
			return self._parent[self._n]
		else:
			raise TypeError('instance has been deleted')


class Register(dict):
	"""The register forms a list of objects and labels them incrementally from 0
		
	the idea is that all the atoms will be stored in here only, it will therefore be very
	easy to loop over all atoms. Indexes will be used to access the instances	
	
	"""

	def __init__(self, cls):
		self._n =0
		self._type =cls
		if not isinstance(cls,str):
			raise TypeError('a register must be initialised with string naming it\'s type. cls had class:' + cls.__class__)

        def append(self, element, key=None):
		if key == None:
			self._n+=1
			self[self._n] = element
			return Index(self, self._n)
		else:
			self[key] = element

@Singleton
class RegisterOfRegisters(Register):
	"""Stores all registers"""
	def __init__(self):
		Register.__init__(self, 'RegisterOfRegisters')
	
def make_a_register(cls):
	"""Factory for registers"""
	if not str(cls) in RegisterOfRegisters.Instance():
		@Singleton
		class aRegister(Register):
			def __init__(self):
				Register.__init__(self, str(cls))
		RegisterOfRegisters.Instance().append(aRegister, str(cls))
	return RegisterOfRegisters.Instance()[str(cls)]
		
class RegisteredObject(object):
	"""a base class ensuring that an object is created as a reference to an instance in a register"""
	def __new__(cls,*args, **kwargs):
		bases = inspect.getmro(cls)
		for i in range(len(bases)):
			b = bases[-i-1]
			try :	
				it = b.__new__(cls) 
				break
			except TypeError: 
				pass 
		reg = make_a_register(cls)
		index = reg.Instance().append(it)
		it._index = Index(reg.Instance(), index._n)
		it.__init__(*args, **kwargs)
		return index

def GetRegister(cls):
	""" Usage: GetRegister('Atom') will return an iterable list of all the atoms"""
	pattern  = '\.'+cls+'\'>'
	for k,reg in RegisterOfRegisters.Instance().iteritems():
		if re.search(pattern, k) != None:
			return reg.Instance().iteritems()
	raise TypeError('No key found that matches: ', cls)
