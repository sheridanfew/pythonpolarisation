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

basename='nTs'
g = open('Polarisabilities_%s.csv' % basename, 'w')
g.write(basename)
g.write('\nLength\tMethod\tFit\tLong_calpol\tLong_pypol\tLong_ratio\tShort_calpol\tShort_pypol\tShort_ratio\tFace_calpol\tFace_pypol\tFace_ratio')

for method in ['Lin','Exp']:
	for fit in ['components', 'mean', 'empirical']:
		name=basename + '_' + method + '_' + fit
		f = open('properies_%s.dat' % name, 'w')
		f.write(name)
		for n in range(1,9,1):
			print 'N', n, ' Namefile: ', name
			exec( "from Molecules.pol_thio_" + str(n) + "T_neut import pol")
			calpol=pol
			calpoldiag=list(np.diag(calpol))

			calpol_max = max(calpoldiag)
			calpol_max_index = calpoldiag.index(calpol_max)

			calpol_min = min(calpoldiag)
			calpol_min_index = calpoldiag.index(calpol_min)

			calpol_thirdindex=[v for v in [0,1,2] if not ( v == calpol_max_index or v == calpol_min_index)][0]

			namefile= str( 'thio_' + str(n) + 'T_neut_aniso_chelpg_thole_' + method + '_' + fit + '.xyz' )
			ReadMoleculeType('../Molecules/' + namefile)
			mol = GetMolecule('../Molecules/' + namefile)

			jm=JMatrix(jmtype='Thole' + method + 'Iso')

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

			pypoldiag=list(np.diag(pypol))

			etamat=np.multiply((pypol-calpol),(pypol-calpol))/np.multiply(calpol,calpol)
			ratios=np.divide(pypol,calpol)
			ratiosdiag=list(np.diag(ratios))

			g.write('\n' + str(n) + '\t' + method + '\t' + fit + '\t' + str(calpoldiag[calpol_max_index]) + '\t' + str(pypoldiag[calpol_max_index]) + '\t' + str(ratiosdiag[calpol_max_index]) + '\t' + str(calpoldiag[calpol_thirdindex]) + '\t' + str(pypoldiag[calpol_thirdindex]) + '\t' + str(ratiosdiag[calpol_thirdindex]) + '\t' + str(calpoldiag[calpol_min_index]) + '\t' + str(pypoldiag[calpol_min_index]) + '\t' + str(ratiosdiag[calpol_min_index]))

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
			f.write(str(n) + '\nNamefile:' + namefile + '\n\nRatios:\n' + str(ratios) + '\n\nCalpol:\n' + str(calpol) + '\n\nPypol:\n' + str(pypol) + '\n\n')

		f.flush()
		f.close()
g.flush()
g.close()
print 'Job Completed Successfully.'
