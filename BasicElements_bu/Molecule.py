from Atom import Atom
from Rotation import Rotation as rot
from Position import Position
#from RegisterMolecules import RegisterMolecules as _rm
from Register import RegisteredObject
import copy

class Molecule ( RegisteredObject, list):
	"""A base molecule class

		
	"""
	#	def __new__(cls, *args, **kwargs):
#		"""Molecule, like Atom returns an "Index" object. 	
#		
#		each molecule contains an _index integer which corresponds to the unique Index identifier
#
#		"""
##		print "entered the baseclass new"
#		it = list.__new__(cls)
#		list.__init__(it)
#		index =  _rm.Instance().append(it)
#		it._index = index._n
#		print "registerd instance: ", it, " with index: ", it._index
#		it.__init__(*args, **kwargs)
#		return index

	def __init__(self, **kwargs):
		list.__init__(self)
		self._type = "BaseClass"
	
	def append(self, ind):
		if not isinstance(ind(), Atom):
			raise TypeError('only atoms can be added to a molecule')
		list.append(self, ind)
 
	
	def add_atom(self, **kwargs):
		"""this is the main interface for adding atoms.

		It adds an atom defined by the dictionary kwarg to the register,
		puts the index onto the molecule and registers the atom as belonging to this molecule
	
		"""
		index = Atom(self._index,**kwargs)
	#	print "molecule index: ", self._index, " atom index: ", index, " atom instance: ", index()
		Molecule.append(self, index)

	def __del__(self):
		#print "entered  molecule delete here"
		pass


	def __repr__(self):
		res = ''
		for i in self:
			res += i().__repr__() + '\n'
		return res	

	def copy(self):
		""" deep copies a molecule. it calls the initialiser of the subtype, 
		then overwrites the position and orientation, copying the values from the atoms
		"""
#		print "enter copy"
		subtype = self.__class__
		copied_molecule = subtype.__new__(subtype)
		copied_molecule._type  = self._type
#		print "now start copying atom information: "
		for i in range(len(copied_molecule())): 
#			print "copying informatio for atom: i=", i
			copied_molecule()[i]()._pos   = self[i]()._pos.copy()
			copied_molecule()[i]()._pol   = self[i]()._pol.copy()
			copied_molecule()[i]()._haspol= self[i]()._haspol
		return copied_molecule

	def rotate(self, r):
		"""rotates the molecule by rotation r


		"""
		if not isinstance(r,rot):
			raise TypeError('must rotate molecules by a rotation')
		for atom in self:
			atom().rotate(r)

	def get_com(self):
		""" returns the com of the molecule.

		I am assuming that the mass of each atom type is the same, in principle
		atoms could have a characteristic called "mass"
		
		"""
		com = Position ([0.,0.,0.])
		for i in self:
			com += i()._pos
		return com/len(self)

	def place_at_d(self, d):
		"""places the com of the molecule at d


		"""
#		print "com: ", self.get_com()
		dx = d - self.get_com()
#		print "dx: ", dx
		for i in self:
			i().move(dx)	

	def move(self, d):
		"""moves molecule by d


		"""
		for i in self:
			i().move(d)		

	def RedCoordMove(self, TVs, movevec):
		"""moves molecule by movevec in reduced coords (requires TVs)
		"""
	   #TVs should be 3*3 matrix containing translation vectors. TVs[0,1,2] are the three translation vectors.
		#movevec specifies number of each TVs to move.
		#print 'movevec', movevec
		#print 'TVs', TVs
		d=(movevec[0]*TVs[0]+movevec[1]*TVs[1]+movevec[2]*TVs[2])
		for i in self:
			i().move(d)

