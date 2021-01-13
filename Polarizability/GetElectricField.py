import sys
sys.path.append('../')

from BasicElements import *
from BasicElements.Register import GetRegister
import numpy as np

def get_electric_field(E0=np.matrix([0.,0.,0.])):
	""" computes the electric field at each polarizable atom due to an external field E0  and due to the charges on all atoms

	E0 is a constant electric field, the other contribution is from the atoms. I do not allow atoms on the same molecule to
	add an electric field onto each other

	"""

	res = []
	for k1, atom1 in GetRegister('Atom'):
		if not atom1._haspol:  # skip centers with no polarisability
			continue
		e = E0.copy()
		for k2, atom2 in GetRegister('Atom'):
			if atom2._parent != atom1._parent:
				#print 'atom2.field_at(atom1._pos)', atom2.field_at(atom1._pos)
				#print 'type(atom2.field_at(atom1._pos))',type(atom2.field_at(atom1._pos))
				e += atom2.field_at(atom1._pos) # add the field at atom2 from atom1
		for i in range(3):	
			res.append(e[0,i]) 
#	print res
	return np.matrix(res) # return the result as a numpy matrix object 
		
""" New way to avoid problems after molecule deletion, doesn't yet work:

	res = [0]*100000
	# make list size 100000 to begin with rather than appending to avoid atom number confusions upon deleting atoms
	for k1, atom1 in GetRegister('Atom'):
		if k1 >= 33333:
			print 'Need to make potential list larger'
			exit
		if not atom1._haspol:  # skip centers with no polarisability
			continue
		e = E0.copy()
		for k2, atom2 in GetRegister('Atom'):
			if atom2._parent != atom1._parent:
				e += atom2.field_at(atom1._pos) # add the field at atom2 from atom1
		for i in range(3):	
			#res.append(e[0,i])
			print 'k1', k1
			print 'i', i
			print 'k1+i', (k1+i) 
			res[k1+i]=e[0,i]
#	print res
	return np.matrix(res) # return the result as a numpy matrix object 
		
"""
