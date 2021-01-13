import sys
sys.path.append('../../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability.GetEnergy import get_energy
from Polarizability import *
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

ReadMoleculeType('../../extras/8T_mul.xyz')
octoT = GetMolecule('../../extras/8T_mul.xyz')
print "8T"
#print octoT()

ReadMoleculeType('../../c60.xyz')
C60 = GetMolecule('../../c60.xyz')
print "C60"
#print C60()

ModifyPolarizability(octoT(), C=8.,H=7.5,S=5.5)
ModifyPolarizability(C60(), C=8.,H=7.5,S=5.5)
cut=4.

distancerangeA = np.arange(0.5, 4., 0.5)

for dist in distancerangeA:

	C60 = GetMolecule('../../c60.xyz')
	distbohr=np.divide(dist,0.529177)
	Molecule.place_at_d(C60(), [0.,0.,distbohr])

	E0 = np.matrix([0.,0.,0.])
	#print "E0"
	#print E0
	d = get_dipoles(E0=E0,cutoff=cut)
	split_d = split_dipoles_onto_atoms(d)
	energy = get_energy(E0=E0,cutoff=cut)
	tot = np.matrix([0.,0.,0.])
	for dd in split_d:
		tot += dd
	print "Hello!"	
	print dist, distbohr, tot, energy
	

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
