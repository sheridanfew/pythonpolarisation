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
		atom()._pol= pol

ReadMoleculeType('c60.xyz')
c60 = GetMolecule('c60.xyz')
#ModifyPolarizability(c60(), 4.)
#print c60()

polrange    = np.arange(4., 9., .5)
cutoffrange = np.arange(3., 7., .5)

print "# polarizability(au^3) cutoff(au) total dipole  "
for p in polrange:
	ModifyPolarizability(c60(), p)
	for  c in cutoffrange:
		d = get_dipoles(E0=np.matrix([0.,0.,10.]),cutoff=c)
		split_d = split_dipoles_onto_atoms(d)
		tot = np.matrix([0.,0.,0.])
		for dd in split_d:
		     tot += dd	
		print p, c, tot[0,2] 


