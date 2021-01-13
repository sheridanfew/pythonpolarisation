import sys
sys.path.append('../')

from BasicElements import *
from BasicElements.Register import GetRegister
import numpy as np

def get_potential():
	""" computes the electric field at each polarizable atom due to an external field E0  and due to the charges on all atoms

	E0 is a constant electric field, the other contribution is from the atoms. I do not allow atoms on the same molecule to
	add an electric field onto each other

	"""
        print GetRegister('Atom')
        type(GetRegister('Atom'))
	res = [0]*10000
	for k1, atom1 in GetRegister('Atom'):
		v = 0
		for k2, atom2 in GetRegister('Atom'):
			if atom2._parent != atom1._parent:
				v += atom2.potential_at(atom1._pos) # add the potential at atom2 from atom1
		res[k1]=v
#	print res
	return np.matrix(res) # return the result as a numpy matrix object 
		
'''
OLDWAY:

        print GetRegister('Atom')
        type(GetRegister('Atom'))
	res = []
	for k1, atom1 in GetRegister('Atom'):
		v = 0
		for k2, atom2 in GetRegister('Atom'):
			if atom2._parent != atom1._parent:
				v += atom2.potential_at(atom1._pos) # add the potential at atom2 from atom1
		res.append(v) 
'''
