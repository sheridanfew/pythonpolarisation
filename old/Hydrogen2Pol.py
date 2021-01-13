from BasicElements import Molecule,Position, Polarizability

class Hydrogen2Pol(Molecule):
	"""an example with 2 polarizable sites
	"""
	def __init__(self, polariz=5.):
		self._type = "h2 with 1 polarizable site"
		self.add_atom(pos=Position([-0.5, 0.,0.]), pol=Polarizability(iso=polariz), elname='H', crg=0.) 
		self.add_atom(pos=Position([ 0.5, 0.,0.]), pol=Polarizability(iso=polariz), elname='H', crg=0.) 
