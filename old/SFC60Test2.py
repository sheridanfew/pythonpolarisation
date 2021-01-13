import sys
from BasicElements import Molecule,Position, Polarizability
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from BasicElements import *
import numpy as np

def ModifyPolarizability(molecule, isopol):
	""" takes a c60 and changes the isotropic polarizability
	for each atom
	"""
	pol = Polarizability(iso=isopol)
	for atom in molecule: 
		atom._pol = pol

ReadMoleculeType('c60.xyz')
c60 = GetMolecule('c60.xyz')

polrange    = np.arange(4., 8., 4.)
cutoffrange = np.arange(3., 6., 1.)

for p in polrange 
	ModifyPolarizability(c60(), p)
	print c60()
d = get_dipoles(E0=np.matrix([0.,0.,1.]),cutoff=5)

#print "SF Dipoles:"
#print d


split_d = split_dipoles_onto_atoms(d)
tot = np.matrix([0.,0.,0.])
for dd in split_d:
#     print dd
     tot += dd	

print "total dipole: " , tot
    
"""print "SF c60 type"
print c60()._type
print "SF c60:"
print c60()"""

"""based on James simpletest.py"""

