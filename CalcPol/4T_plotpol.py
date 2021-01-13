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

for n in range(1,9,1):
	name= str( 'thio_' + str(n) + 'T_E_100' )
	print 'N', n, ' Namefile: ', name
	namefile= str( 'thio_' + str(n) + 'T_neut_theta_aniso_chelpg.xyz' )
	ReadMoleculeType('../Molecules/' + namefile)
	mol = GetMolecule('../Molecules/' + namefile)

	jm=JMatrix(jmtype='Stern',cutoff=8.0)

	E0 = np.matrix([1.,0.,0.])
	d = get_dipoles(E0=E0,jm=jm._m)
	split_d = split_dipoles_onto_atoms(d)

	# print dipoles
	if not os.path.exists('Dips_Posns_TVs'): os.makedirs('Dips_Posns_TVs')
	f = open('Dips_Posns_TVs/%s_dipoles.dat' % name, 'w')
	for dd in split_d:
		dstr=str(dd)
		f.write(dstr)
		f.write('\n')
	f.flush()
	f.close()

	f = open('Dips_Posns_TVs/%s_posns.dat' % name, 'w')
	f.write('Molecule Starts:\n')
	for atom in GetRegister('Atom'):
		astr=str(atom)
		f.write(astr)
		f.write('\n')
	f.write('Molecule Ends.')
	f.flush()
	f.close()

print 'Job Completed Successfully.'
