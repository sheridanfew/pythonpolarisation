import sys
from BasicElements import Molecule,Position, Polarizability
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from BasicElements import *
import numpy as np

def ModifyPolarizability(molecule, C, H, S):
	""" takes a mol and changes the isotropic polarizability
	for each atom
	"""
	C = Polarizability(iso=C)
	H = Polarizability(iso=H)
	S = Polarizability(iso=S)
	for atom in molecule:
		if (atom()._elname=="C"):
			atom()._pol= C
		elif (atom()._elname=="H"):
			atom()._pol= H
		elif (atom()._elname=="S"):
			atom()._pol= S

ReadMoleculeType('12T_no_charge.xyz')
mol = GetMolecule('12T_no_charge.xyz')

#ModifyPolarizability(mol(), C=4.,H=5.,S=5.)
#print mol()
#ModifyPolarizability(mol(), C=5.,H=5.,S=4.)
#print mol()


polrangeC    = np.arange(1., 10., .5)
polrangeS    = np.arange(1., 10., .5)
polrangeH    = np.arange(1., 10., .5)
cutoffrange = np.arange(4., 12., .5)

print "# polarizability(au^3) cutoff(au) total dipole  "
for p in polrange:
	ModifyPolarizability(mol(), C=p,S=p,H=0.1 )
	for  c in cutoffrange:
		d = get_dipoles(E0=np.matrix([1.,0.,0.]),cutoff=c)
		split_d = split_dipoles_onto_atoms(d)
		tot = np.matrix([0.,0.,0.])
		for dd in split_d:
		     tot += dd	
		print p, c, tot


"""
	switch = {
	    'C': atom()._pol= C
	    'H': atom()._pol= H
	    'S': atom()._pol= S
	    }

	try:
	    result = switch[choice]
	except KeyError:
	    print 'I didn\'t understand your choice.'
	else:
	    result()
"""
