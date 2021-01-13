import sys
sys.path.append('../../')
from BasicElements import *
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

ReadMoleculeType('../../8T_no_charge.xyz')
octoT = GetMolecule('../../8T_no_charge.xyz')

ModifyPolarizability(octoT(), C=4.,H=5.,S=5.)

print "8T"

print octoT()

ReadMoleculeType('../../c60.xyz')
C60 = GetMolecule('../../c60.xyz')

ModifyPolarizability(C60(), C=4.,H=5.,S=5.)

print "C60"

print C60()

Molecule.place_at_d(C60(), [0.,0.,10.])

print "C60 moved"

print C60()


cut=4.

polrangeC    = np.arange(1., 10., .5)
polrangeS    = np.arange(1., 10., .5)
polrangeH    = np.arange(1., 10., .5)
cutoffrange = np.arange(4., 12., .5)

#polrangeC    = np.arange(1., 2., .5)
#polrangeS    = np.arange(1., 2., .5)
#polrangeH    = np.arange(1., 2., .5)
#cutoffrange = np.arange(4., 5., .5)

#print "polarizabilityC(au^3) polarizabilityS(au^3) polarizabilityH(au^3) cutoff(au) total dipole  "


E0 = np.matrix([0.,0.,1.])
#print "E0"
#print E0
d = get_dipoles(E0=E0,cutoff=cut)
split_d = split_dipoles_onto_atoms(d)
tot = np.matrix([0.,0.,0.])
for dd in split_d:
	tot += dd	
	print cut, tot

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
