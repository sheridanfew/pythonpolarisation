import sys
sys.path.append('../')

#from BasicElements import *
from BasicElements.Register import GetRegister
from GetElectricField import get_electric_field
from GetPotential import get_potential
from JMatrix import get_atoms_with_polarizability
import numpy as np

def get_atoms():
	n=0
	for i, atom in GetRegister('Atom'):
			n+=1
	return n

def get_U_dipdip(**kwargs):
	"""SF edited from GetDipoles.py
	Uses Jmat and dips to get Udd
	
	"""
	Efield = kwargs.get('Efield')
	jm = kwargs.get('jm')
	dips = kwargs.get('dips')
	U_dd=0.5 * dips.T * jm * dips
	#print 'U_dd:'
	#print U_dd
	return float(U_dd[0][0])

def get_U_qdip(Efield,dips):
	"""SF edited from GetDipoles.py
	SF edited from GetDipoles.py
	Uses E-field and dips to get Uqd
	
	"""

	natom = get_atoms_with_polarizability()
	Efield_multi = np.reshape(Efield,(natom,3))

	dipoles_multi = np.reshape(dips,(natom,3))

	atoms=np.arange(0,natom,1)

	Uqdiptot=0.0
	for atom in atoms:
		Uqdip = -(dipoles_multi[atom] * Efield_multi[atom].T)
		Uqdiptot += Uqdip

	return float(Uqdiptot[0][0])

def get_U_qq(**kwargs):
	'''Uses potential (from get_potential) for Uqq'''
	potential = kwargs.get('potential')
	#matrix with entries for potential

	Uqqtot=0.0
	for k, atom in GetRegister('Atom'):
#		print k
		if atom._crg != 0.0:
			Uqq = float(atom._crg)*float(potential[0,k-1].T)
			Uqqtot += Uqq/2
	#print Uqqtot
	return float(Uqqtot)

def get_energy(**kwargs):
	E0 = kwargs.get('E0', np.matrix([0.,0.,0.]))
	cutoff = kwargs.get('cutoff', 0.)
	Udd=get_U_dipdip(jm=jm,dips=dips)
	Uqd=get_U_qdip(dips=dips,Efield=Efield)
	Uqq=get_U_qq(potential=potential)
	Utot=Udd+Uqd+Udd

	print 'Utot, Uqq, Udd, Uqd'
	print energy, Uqq, Udd, Uqd
	return float(energy)


def get_jmat(**kwargs):
	"""SF edited from GetDipoles.py
	computes the dipole on each atom due to an external field E0 and due to the charges on each atom in the box, then uses to calc energy of system.

	the input can be E0=XXX, cutoff=YYY where XXX is an np.matrix object containing the external field
	and YYY is a float with the cutoff for the J_BB' elements	
	
	"""
	E0 = kwargs.get('E0', np.matrix([0.,0.,0.]))
	cutoff = kwargs.get('cutoff', 0.)
	#if DEBUG:
	#	print "cutoff used:", cutoff
	jm = JMatrix(cutoff=cutoff)
	Efield = get_electric_field(E0)
	dips = np.linalg.solve(jm._m, Efield.T)
#	return dips.T
#	return -0.5 * dips.T * jm._m * dips
	return jm._m


