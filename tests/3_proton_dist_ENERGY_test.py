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
proton3 =  GetMolecule('../Molecules/protontest.xyz')

#ReadMoleculeType('../Molecules/H2test.xyz')
#proton=GetMolecule('../Molecules/H2test.xyz')

cut= 4.00007106634
#ModifyPolarizability(proton(), C=9.072879660383,H=0.100004083159,S=13.2644932888)
#	print "8T"
#	print proton()
#	print "C60"

Molecule.place_at_d(proton3(),np.matrix([0.,0.,-8]) )

distancerange = np.arange(8., 8.1, 1.0)

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

	print "E"
	print Efield

	natom   = get_atoms_with_polarizability()
	Efield_multi = np.reshape(get_electric_field(E0),(natom,3))

	print Efield_multi

	print "Potential"
	
	print potential


	Uqqtot=0.0

	for k, atom in GetRegister('Atom'):
#		Uqq = atom._crg * potential[atom].T 
		print 'charge'
		print atom._crg
#		print 'pot at:'
#		print potential[atom].T
#		print 'pot at with k:'
#		print potential[0,k-1].T
		print 'Uqq'
		print atom._crg * potential[0,k-1].T
		Uqq = atom._crg*float(potential[0,k-1].T)
#		print type(Uqq)
#		print type(Uqqtot)
		Uqqtot += Uqq/2


#	Uqq=get_U_q_q()

	print 'Uqqtot'

	print Uqqtot

	Uqqfromfunct=get_U_qq()

	print 'Uqqfromfunct'
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

	print 'Jmat'
	print JMatrix(cutoff=cut)._m

	print 'Jmat.dips'
	print d * JMatrix(cutoff=cut)._m *d.T

	atoms=np.arange(0,natom,1)
	Uqdiptot=0
#	for atom in atoms:
#		print 'atom:'
#		print atom
#		Uqdip = dipoles_multi[atom] * Efield_multi[atom].T 
#		Uqdiptot += Uqdip
