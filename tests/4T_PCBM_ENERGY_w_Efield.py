import sys
sys.path.append('../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import *
from Polarizability.GetEnergy import get_energy
import numpy as np

inputs=np.array([[8.0,9.5,2.0,8.5]])
#[8.0, 9.5, 2.5, 8.5],
#[8.0, 9.0, 2.5, 8.5],
#[7.5, 9.5, 3.0, 8.5],
#[7.5, 9.5, 3.5, 8.5]])
# format polC, polH, polS, cut, eta
steps=np.arange(0,inputs.shape[0],1)
distancerangeA = np.arange(2., 11., 1.)

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

ReadMoleculeType('../Molecules/4T_cation_mul.xyz')
octoT = GetMolecule('../Molecules/4T_cation_mul.xyz')
ReadMoleculeType('../Molecules/PCBM_EDITEDO_anion_mul.xyz')
C60 = GetMolecule('../Molecules/PCBM_EDITEDO_anion_mul.xyz')

for i in steps:
	cut=inputs[i,3]
	ModifyPolarizability(octoT(), C=inputs[i,0],H=inputs[i,1],S=inputs[i,2])
	ModifyPolarizability(C60(), C=inputs[i,0],H=inputs[i,1],S=inputs[i,2])
#	print "4T"
#	print octoT()
#	print "C60"
#	print C60()
	print "Cutoff:"
	print cut
	print "cut, dist, distbohr, tot, energy, energyev"
	for dist in distancerangeA:
		distbohr=np.divide(dist+3.3,0.529177)
		Molecule.place_at_d(C60(), [0.,0.,distbohr])
		#print C60()

		E0 = np.matrix([1000.,0.,0.])
		#print "E0"
		#print E0
		d = get_dipoles(E0=E0,cutoff=cut)
		split_d = split_dipoles_onto_atoms(d)
		tot = np.matrix([0.,0.,0.])
		for dd in split_d:
			tot += dd
		energy = get_energy(E0=E0,cutoff=cut)
		energyev = np.multiply(energy,27.211)



		print cut, dist, distbohr, tot, energy, energyev

