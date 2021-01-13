import sys
sys.path.append('../')

from BasicElements import *
from GetElectricField import get_electric_field
from JMatrix import JMatrix
import numpy as np

DEBUG = False
def get_energy(**kwargs):
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

	print 'E field:'
	print Efield	

	energyEfield= np.dot(Efield,dips)

	print 'Energy E field:'
	print energyEfield

	print 'Energy E field ev approx:'
	print energyEfield*27	

#	print dips.T
#	print dips.T * jm._m * dips 
	
	energydipoles=dips.T * jm._m * dips

	print 'Energy dipoles ev approx:'
	print energydipoles*27			

	Energytotev=(energyEfield*27) + (energydipoles*27)
	print 'Energy tot, ev:'
	print Energytotev

