from Register import Index as Index
from Register import RegisteredObject
from Rotation import Rotation as rot
from Position import Position
from Polarizability import Polarizability
import numpy as np

class Atom(RegisteredObject, object):
	"""an Atom returns an Index object
	
	the universal singleton register of all atoms is stored in RegisterAtoms

	"""
	def __init__(self, parent=None, **kwargs):
		"""an atom is initialised with a position, a polarizability, a charge and an element name

		example: a=Atom(pos=Position([0.,0.,0.]), pol=Polarizability(iso=10.), crg=-1., elname='H')
		if an atom belongs to a molecule mol, the best initilisation is to set the "parent" kewyword:
		a=Atom(mol._index,pos=Position([0.,0.,0.]), pol=Polarizability(iso=10.), crg=-1., elname='H')
		mol._index is the integer that identifies the position of molecule mol in the RegisterMolecules dictionary

		"""
#		print "kwargs:" , kwargs
		self._pos = kwargs.get('pos', Position([0.,0.,0.]))
		self._pol = kwargs.get('pol', Polarizability(iso=0.))
		self._crg = kwargs.get('crg',0.)
		self._elname = kwargs.get('elname', 'noname')
		self._parent  = parent
		self._bondparnters = kwargs.get('bondpartners',0.)
		self._cut = kwargs.get('cut', 0.)
		
		convert = kwargs.get('convert_to_bohr' , False)
		if convert :
			self._pos = self._pos *0.521
    
		if np.linalg.det(self._pol) == 0.0:
			self._haspol = False
		else:
			self._haspol = True
	
	def __repr__(self):
		return self._elname + ' pos: ' + str(self._pos) + 'parent: ' + str(self._parent) + ' haspol? ' + str(self._haspol) + '; charge: ' + str(self._crg) + ' ||pol|| : ' + str(np.linalg.det(self._pol)) + '|pol|' + str(np.sqrt( (self._pol[0,0]*self._pol[0,0]) + (self._pol[1,1]*self._pol[1,1]) + (self._pol[2,2]*self._pol[2,2]) ) )

	def rotate(self, r):	
		"""rotates the atom by r

		this function will change the position of the molecule 
		- so long as its position != 0. 0. 0.
		
		"""
		self._pos = self._pos.rotate(r)
		self._pol = self._pol.rotate(r) 	
					

	def move(self, d):
		"""moves the atom by d


		"""
		self._pos+=d

	def setcharge(self, q):
		"""sets atom charge to q
		"""
		self._crg = q 

	def setpol(self, pol):
		"""sets atom polarisability to q
		"""
		self._pol = pol 

	def field_at(self, pos):
		"""finds the field from the atom at position pos



		"""
		d = pos -self._pos
		#print 'pos, self._pos, d, d.abs(), self._crg',pos, self._pos, d, d.abs(), self._crg
		if self._crg=='0.0':
			return Position([0.,0.,0.])
		else:
			return (self._crg/(d.abs())**3 ) * d
	
	def potential_at(self, pos):
		"""finds the potential from the atom at position pos



		"""
		d = pos -self._pos
		#print 'self._crg,d.abs()',self._crg,d.abs()
		if self._crg=='0.0':
			return np.float('0.0')
		else:
			return (self._crg/d.abs())
	


