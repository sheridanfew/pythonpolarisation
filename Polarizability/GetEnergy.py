import sys
sys.path.append('../')

from BasicElements import *
from GetElectricField import get_electric_field
from GetPotential import get_potential
from JMatrix import *
import numpy as np

DEBUG = False

def get_atoms():
	n=0
	for i, atom in GetRegister('Atom'):
			n+=1
	return n

def get_U_dipdip(**kwargs):
	"""SF edited from GetDipoles.py
	computes the dipole on each atom due to an external field E0 and due to the charges on each atom in the box, then uses to calc energy of system.

	the input can be E0=XXX, cutoff=YYY where XXX is an np.matrix object containing the external field
	and YYY is a float with the cutoff for the J_BB' elements	
	
	"""
	E0 = kwargs.get('E0', np.matrix([0.,0.,0.]))
	cutoff = kwargs.get('cutoff', 0.)
	if DEBUG:
		print "cutoff used:", cutoff
	jm = JMatrix(cutoff=cutoff)
	Efield = get_electric_field(E0)
	dips = np.linalg.solve(jm._m, Efield.T)
#	return 0.5 * dips.T * jm._m * dips
	U_dd=0.5 * dips.T * jm._m * dips
	print 'U_dd:'
	print U_dd
	return U_dd[0][0]

def get_U_qdip(**kwargs):
	"""SF edited from GetDipoles.py
	computes the dipole on each atom due to an external field E0 and due to the charges on each atom in the box, then uses to calc energy of system.

	the input can be E0=XXX, cutoff=YYY where XXX is an np.matrix object containing the external field
	and YYY is a float with the cutoff for the J_BB' elements	
	
	"""
	E0 = kwargs.get('E0', np.matrix([0.,0.,0.]))
	cutoff = kwargs.get('cutoff', 0.)
	if DEBUG:
		print "cutoff used:", cutoff
	jm = JMatrix(cutoff=cutoff)
	Efield = get_electric_field(E0)
	natom = get_atoms_with_polarizability()
	Efield_multi = np.reshape(get_electric_field(E0),(natom,3))

	dips = np.linalg.solve(jm._m, Efield.T)
	dipoles_multi = np.reshape(dips,(natom,3))

	atoms=np.arange(0,natom,1)
	Uqdiptot=0.0
	for atom in atoms:
		Uqdip = -(dipoles_multi[atom] * Efield_multi[atom].T)
		Uqdiptot += Uqdip

	return Uqdiptot[0]
#	return dips.T
#	return 0.5 * dips.T * jm._m * dips

def get_U_qq(**kwargs):
	potential = get_potential()
#	print potential
#	print type(potential)
#	print potential.shape
#	print potential.ndim

	#matrix with entries for potential
	Uqqtot=0.0
	for k, atom in GetRegister('Atom'):
#		print k
		if atom._crg != 0.0:
			Uqq = float(atom._crg)*float(potential[0,k-1].T)
			Uqqtot += Uqq/2
	#print Uqqtot
	return Uqqtot

def get_energy(**kwargs):
	E0 = kwargs.get('E0', np.matrix([0.,0.,0.]))
	cutoff = kwargs.get('cutoff', 0.)
	energy=get_U_dipdip(E0=E0,cutoff=cutoff)+get_U_qdip(E0=E0,cutoff=cutoff)+get_U_qq()
	Uqq=get_U_qq()
	Uqd=get_U_qdip(E0=E0,cutoff=cutoff)
	Udd=get_U_dipdip(E0=E0,cutoff=cutoff)

	print 'Utot, Uqq, Udd, Uqd'
	print energy, Uqq, Udd, Uqd
	return energy


def get_jmat(**kwargs):
	"""SF edited from GetDipoles.py
	computes the dipole on each atom due to an external field E0 and due to the charges on each atom in the box, then uses to calc energy of system.

	the input can be E0=XXX, cutoff=YYY where XXX is an np.matrix object containing the external field
	and YYY is a float with the cutoff for the J_BB' elements	
	
	"""
	E0 = kwargs.get('E0', np.matrix([0.,0.,0.]))
	cutoff = kwargs.get('cutoff', 0.)
	if DEBUG:
		print "cutoff used:", cutoff
	jm = JMatrix(cutoff=cutoff)
	Efield = get_electric_field(E0)
	dips = np.linalg.solve(jm._m, Efield.T)
#	return dips.T
#	return -0.5 * dips.T * jm._m * dips
	return jm._m


