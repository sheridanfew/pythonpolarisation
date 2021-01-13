from BasicElements import Molecule,Position, Polarizability

class Hydrogen1Pol(Molecule):
	def __init__(self, polariz=10.):
		self._type = "h2 with 1 polarizable site"
		self.add_atom(pos=Position([-0.5, 0.,0.]), elname='H', crg=0.) 
		self.add_atom(pos=Position([ 0.5, 0.,0.]), elname='H', crg=0.) 
		self.add_atom(pos=Position([ 0. , 0.,0.]), pol=Polarizability(iso=polariz), elname='dipole', crg=0.) 