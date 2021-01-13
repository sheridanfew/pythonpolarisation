import sys
sys.path.append('../')

from BasicElements import *
from GetElectricField import get_electric_field
#from JMatrix import JMatrix
import numpy as np


def get_dipoles(**kwargs):
	"""computes the dipole on each atom due to an external field E0 and due to the charges on each atom in the box

	the input can be E0=XXX, cutoff=YYY where XXX is an np.matrix object containing the external field
	and YYY is a float with the cutoff for the J_BB' elements	
	
	"""
	for key in kwargs:
			if key == 'Efield':
				Efield = kwargs[key]
			elif key == 'E0':
				Efield = get_electric_field(kwargs[key])

	jm = kwargs.get('jm')

	dips = np.linalg.solve(jm, Efield.T)
	return dips.T

def split_dipoles_onto_atoms(dips):
	res = []
	i=0
	while i < dips.shape[1]:
		dipat = np.matrix( [dips[0,i], dips[0,i+1] , dips[0,i+2]] )
		res.append(dipat)
		i+=3
	return res
