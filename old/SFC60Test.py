import sys
from BasicElements import Molecule,Position, Polarizability
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule

ReadMoleculeType('BasicElements/tests/H2.xyz')
h2 = GetMolecule('BasicElements/tests/H2.xyz')
print h2()

ReadMoleculeType('extras/c60_neut_B3LYP_631GSTAR_mul.xyz')
c60 = GetMolecule('extras/c60_neut_B3LYP_631GSTAR_mul.xyz')
print c60()

"""class H1Pol(Molecule):
	def __init__(self, polariz=10.):
		self._type = "h2 with 1 polarizable site"
		self.add_atom(pos=Position([-0.5, 0.,0.]), elname='H', crg=0.) 
		self.add_atom(pos=Position([ 0.5, 0.,0.]), elname='H', crg=0.) 
		self.add_atom(pos=Position([ 0. , 0.,0.]), pol=Polarizability(iso=polariz), elname='dipole', crg=0.) """
