import sys
sys.path.append('../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import *
from Polarizability.GetEnergy import *
#from Polarizability.GetElectricField import get_electric_field_multirow
from Polarizability.JMatrix import * 

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

ReadMoleculeType('../Molecules/protontest.xyz')
proton = GetMolecule('../Molecules/protontest.xyz')
proton2 =  GetMolecule('../Molecules/protontest.xyz')

#ReadMoleculeType('../Molecules/H2test.xyz')
#proton=GetMolecule('../Molecules/H2test.xyz')

cut= 1.00007106634
#ModifyPolarizability(proton(), C=9.072879660383,H=0.100004083159,S=13.2644932888)
#	print "8T"
#	print proton()
#	print "C60"

distancerange = np.arange(2., 2.1, 1.0)

for dist in distancerange:
	print "Hello!"
	print proton()
	print "Cutoff:"
	print cut
	Molecule.place_at_d(proton2(),np.matrix([0.,0.,dist]) )
	E0 = np.matrix([0.,0.,0.])
	d = get_dipoles(E0=E0,cutoff=cut)
	split_d = split_dipoles_onto_atoms(d)
	tot = np.matrix([0.,0.,0.])
	for dd in split_d:
		tot += dd
#	energy = get_energy(E0=E0,cutoff=cut)
#	energyev = np.multiply(energy,27.211)
	jm = get_jmat(E0=E0,cutoff=cut)
	Efield = get_electric_field(E0)
	potential = get_potential()

	print "dist, tot, "
#energy, energyev"
	print dist, tot, 
#energy, energyev

	print "jm"
	print jm

	print "jm._m * Efield.T"
	print jm * Efield.T

	print "d * Efield.T"
	print d * Efield.T

	print "d * jm._m * d.T"
	print d * jm * d.T

	print "Efield.T"
	print Efield.T

	print "jm._m * d.T"
	print jm * d.T


	print "jm*d"
	print jm * d.T


	print "E"
	print Efield

	natom   = get_atoms_with_polarizability()
	Efield_multi = np.reshape(get_electric_field(E0),(natom,3))

	print Efield_multi

	print "Potential"
	
	print potential

	Uqqfromfunct=get_U_qq()

	print 'Uqq'
	print Uqqfromfunct

	Udipsdips=get_U_dipdip(E=E0,cutoff=cut)

	print 'Udipsdips'
	print Udipsdips

	Uqdips=get_U_qdip(E=E0,cutoff=cut)

	print 'Uqdips'
	print Uqdips

	utot=get_energy(E0=E0,cutoff=cut)
	print 'Energy tot'
	print utot

	print "dipoles:"
	print d
	dipoles_multi = np.reshape(d,(natom,3))
	print dipoles_multi

	atoms=np.arange(0,natom,1)
	Uqdiptot=0

	print 'atoms:'
	print proton()
	print proton2()
#		Uqdip = dipoles_multi[atom] * Efield_multi[atom].T 
#		Uqdiptot += Uqdip
