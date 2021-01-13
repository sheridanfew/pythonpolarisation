import sys
sys.path.append('../')
from BasicElements import *
from BasicElements.Register import GetRegister
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from BasicElements.Crystal import *
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import *
from Polarizability.GetEnergyFromDips import *
from Polarizability.JMatrix import JMatrix
import numpy as np
from math import *
from time import gmtime, strftime
import os

print strftime("%a, %d %b %Y %X +0000", gmtime())

name='Ethane'

for n in range(1,2,1):
	exec( "from Molecules.pol_C2H6 import pol")
	calpol=pol
	namefile= str( 'C2H6_w_connectivity_Tholeexp.xyz' )
	ReadMoleculeType('../Molecules/' + namefile)
	mol = GetMolecule('../Molecules/' + namefile)

	jm=JMatrix(jmtype='TholeExpIso')

	pypol=np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
	etamat=np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])  

	for i in np.arange(0. ,2.1 ,1. ):
		E0 = np.matrix([0.,0.,0.])
		E0[0,i]=1.
		d = get_dipoles(E0=E0,jm=jm._m)
		split_d = split_dipoles_onto_atoms(d)
		tot = np.matrix([0.,0.,0.])
       		for dd in split_d:
        		tot += dd
		print 'tot'
		print tot
       		pypol.T[i] = tot

	etamat=np.multiply((pypol-calpol),(pypol-calpol))/np.multiply(calpol,calpol)

#	If calibration polarisability is 0, eta not appropriate fit here (and previous funct will have made undefined), set value in etamat to 0
	for i in np.arange(0,3,1):
		for j in np.arange(0,3,1):
                	if calpol[i,j] == 0:
                    		etamat[i,j] = '0'

	eta=0
	eta = etamat[0,0]+etamat[1,0]+etamat[1,1]+etamat[2,0]+etamat[2,1]+etamat[2,2]

	print 'namefile: ', namefile
	print 'eta', eta

	print '\ncalpol:\n'
	print calpol
	print '\npypol\n'
	print pypol

'''
f = open('properies_%s.dat' % namefile, 'w')
f.write
f.write(TVstr)
f.flush()
f.close()
'''

print 'Job Completed Successfully.'
