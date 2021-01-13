import sys
sys.path.append('../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import *
from Polarizability.GetEnergy import get_energy
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

ReadMoleculeType('../Molecules/neutrontest.xyz')
proton = GetMolecule('../Molecules/neutrontest.xyz')
proton2 =  GetMolecule('../Molecules/neutrontest.xyz')

#ReadMoleculeType('../Molecules/H2test.xyz')
#proton=GetMolecule('../Molecules/H2test.xyz')

cut= 4.00007106634
#ModifyPolarizability(proton(), C=9.072879660383,H=0.100004083159,S=13.2644932888)
#	print "8T"
#	print proton()
#	print "C60"

distancerange = np.arange(1., 5.1, 1.0)

for dist in distancerange:
	print "Hello!"
	print proton()
	print "Cutoff:"
	print cut
	Molecule.place_at_d(proton2(),np.matrix([0.,0.,dist]) )
	E0 = np.matrix([1.,0.,0.])
	d = get_dipoles(E0=E0,cutoff=cut)
	split_d = split_dipoles_onto_atoms(d)
	tot = np.matrix([0.,0.,0.])
	for dd in split_d:
		tot += dd
	energy = get_energy(E0=E0,cutoff=cut)
	energyev = np.multiply(energy,27.211)

	print "dist, tot, energy, energyev"
	print dist, tot, energy, energyev

