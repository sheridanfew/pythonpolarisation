# 12Tneutral:PC71BM_neutral_B3LYP_6_31gSTAR_opt_mul polarisation model. from files:
#12T_neutral_B3LYP_6_31gSTAR_opt_mul.xyz
#PC71BM_neutral_B3LYP_6_31gSTAR_opt_mul.xyz
#Produced with makeblends_pol_model.sh

import sys
sys.path.append('../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import *
from Polarizability.GetEnergy import *
import numpy as np
from math import *
from time import gmtime, strftime

print strftime("%a, %d %b %Y %X +0000", gmtime())

ReadMoleculeType('../Molecules/./TIPS_Pc_from_cif_rm_extrachains___TIPS_Pc_neutral_cifcoords_HF_6_31PLUSgSTAR_cifstruct_mul.xyz')

centre='TIPS_Pc_from_cif_rm_extrachains___TIPS_Pc_neutral_cifcoords_HF_6_31PLUSgSTAR_cifstruct_mul.xyz'
surroundings='TIPS_Pc_from_cif_rm_extrachains___TIPS_Pc_neutral_cifcoords_HF_6_31PLUSgSTAR_cifstruct_mul'

#Translation vectors:
TV=np.matrix([[7.56500000, 0.00000000, 0.00000000],
[0.85985150, 7.70215265, 0.00000000],
[3.37939010, -0.12597319, 16.49184884]])

cut=5.0

f = open('../Crystals/TIPS_Pc_neutral_neutral_2013_06_20_new.csv', 'w')
f.write ("Crystal Struct\tCentre\tSurroundings\tSize\tEnergy\tUqq\tUqd\tUdd\tDipole Moment")
f.flush()
f.close()


# size given number of translation vectors ion each direction from central Pc

#mols=[]
#for size in np.arange(0,6,1):
#	for y in np.arange(-size,size+0.1,1):
#		for x in np.arange(-size,size+0.1,1):
#			for z in np.arange(-size,size+0.1,1):
#				mols.apped[Pc()]

for size in np.arange(1,2,1):
	mols = [[[] for i in range((2*size)+1)] for i in range((2*size)+1)]
	print mols
	for x in np.arange(-size,size+1,1):
		xlist=x+size
		print 'x', x
		print 'xlist', xlist
		print xlist
		for y in np.arange(-size,size+1,1):
			ylist=y+size
			for z in np.arange(-size,size+1,1):
				zlist=z+size
				if x == 0 and y == 0 and z == 0:
					#central mol
					mols[xlist][ylist].append(GetMolecule('../Molecules/./TIPS_Pc_from_cif_rm_extrachains___TIPS_Pc_neutral_cifcoords_HF_6_31PLUSgSTAR_cifstruct_mul.xyz'))
				else:
# Places molecule only if not in centre (centre mol places elsewhere)
					mols[xlist][ylist].append(GetMolecule('../Molecules/./TIPS_Pc_from_cif_rm_extrachains___TIPS_Pc_neutral_cifcoords_HF_6_31PLUSgSTAR_cifstruct_mul.xyz'))
					#print mols
					for atom in mols[xlist][ylist][zlist]():
						atom().move(x*TV[0]+y*TV[1]+z*TV[2])
	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'mols', mols
	print mols[size][size][size]()
#	for atom in mols[size][size][size]():
#		print atom
#Remove duplicated centre molecule
#	del (mols[size][size][size])
	print 'size', size
	print 'atoms:'
	for atom in Register.GetRegister('Atom'):
		print atom
	E0 = np.matrix([0.,0.,0.])
#	jm = get_jmat(E0=E0,cutoff=cut)
#	print jm
	d = get_dipoles(E0=E0,cutoff=cut)
	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'dips', d
	split_d = split_dipoles_onto_atoms(d)
	tot = np.matrix([0.,0.,0.])
	for dd in split_d:
		tot += dd
	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'tot', tot
	Uqq = np.multiply(get_U_qq(),27.211)
	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'Uqq', Uqq
	Uqd = np.multiply(get_U_qdip(),27.211)
	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'Uqd', Uqd
	Udd = np.multiply(get_U_dipdip(),27.211)
	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'Udd', Udd
	energyev = Udd+Uqd+Uqq
	print 'energyev', energyev
	energy=energyev/27.211

#csv file



	f = open('../Crystals/TIPS_Pc_neutral_neutral_2013_06_20_new.csv', 'a')
	f.write ('\nTIPS-Pc\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (centre,surroundings,size,energyev,Uqq,Uqd,Udd,tot))
	f.flush()
	f.close()

#_new.dat for gnuplot

	f = open('TIPS_neutral_neutral_size_%s_new.dat' % size, 'w')
	f.write('Molecule Starts:\n')
	for atom in Register.GetRegister('Atom'):
		astr=str(atom)
		f.write(astr)
		f.write('\n')
	f.write('Molecule Ends.\n\n')
	f.write ("Crystal Struct\tCentre\tSurroundings\tSize\tEnergy\tUqq\tUqd\tUdd\tDipole Moment")
	f.write ('\nTIPS-Pc\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (centre,surroundings,size,energyev,Uqq,Uqd,Udd,tot))
	f.write('\n\ndipoles start:\n')
	for dd in split_d:
		dstr=str(dd)
		f.write(dstr)
		f.write('\n')
	f.flush()
	f.close()
exit()

