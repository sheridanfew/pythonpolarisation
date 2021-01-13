# 12Tanion:PC71BM_anion_B3LYP_6_31gSTAR_opt_mul polarisation model. from files:
#12T_anion_B3LYP_6_31gSTAR_opt_mul.xyz
#PC71BM_anion_B3LYP_6_31gSTAR_opt_mul.xyz
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

ReadMoleculeType('../Molecules/./PENCEN_neut_HF631PLUSGSTAR_aniso_w_connectivity_cifstruct_chelpg.xyz')
ReadMoleculeType('../Molecules/./PENCEN_neut_HF631PLUSGSTAR_aniso_w_connectivity_cifstruct_chelpg.xyz')

centre='PENCEN_neut_HF631PLUSGSTAR_aniso_w_connectivity_cifstruct_chelpg.xyz'
surroundings='PENCEN_neut_HF631PLUSGSTAR_aniso_w_connectivity_cifstruct_chelpg.xyz'


#From cif:
'''
_symmetry_equiv_pos_as_xyz
1 x,y,z
2 -x,-y,-z
_cell_length_a                   7.900
_cell_length_b                   6.060
_cell_length_c                   16.010
_cell_angle_alpha                101.90
_cell_angle_beta                 112.60
_cell_angle_gamma                85.80
_cell_volume                     692.384
'''
#Get translation vectors:

a=7.900
b=6.060
c=16.010

alpha=101.90*(pi/180)
beta=112.60*(pi/180)
gamma=85.80*(pi/180)

cell_volume=sqrt(1 - (cos(alpha)**2) - (cos(beta)**2) - (cos(gamma)**2) + (2*cos(alpha)*cos(beta)*cos(gamma)))
cif_cell_volume=692.384

print 'calc cell_vol, cif cell_vol:',cell_volume, ', ' cif_cell_volume

matrix_to_cartesian=np.matrix( [[a, b*cos(gamma), c*cos(beta)],
[0, b*sin(gamma), c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)],
[0, 0, c*cell_volume/sin(gamma)]])

TV=matrix_to_cartesian.T

#Reflection Properties (see http://mathworld.wolfram.com/Reflection.html):

'''No reflection for Pc
#Getting this sort of info v. tricky for some molecules, have glide planes, screw axes, etc. See Display,Symmetry Elements in mercury for individual structs
a_refl=-1
b_refl=1
c_refl=-1
d_refl=1.5

n_refl=np.matrix( [[a_refl, b_refl, c_refl]] )/ np.sqrt(np.dot(np.matrix( [[a_refl, b_refl, c_refl]] ),np.matrix( [[a_refl, b_refl, c_refl]] ).T))

'''

#TVs, TV[0,1,2] are the three translation vectors.

print 'TV:'
print TV

print 'TV[0,1,2] are the three translation vectors.'

cut=5.0

f = open('../Crystals/TIPS_Pc_anion_neutral_2013_06_20_new.csv', 'w')
f.write ("Crystal Struct\tCentre\tSurroundings\tReorg\tSize\tEnergy\tUqq\tUqd\tUdd\tDipole Moment")
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
					mols[xlist][ylist].append(GetMolecule('../Molecules/./sexithiophene_wH_anion_aniso_w_connectivity_cifstruct_chelpg.xyz'))
					print 'xlist,ylist,zlist',xlist,ylist,zlist
				else:
# Places molecule only if not in centre (centre mol places elsewhere)
					mols[xlist][ylist].append(GetMolecule('../Molecules/./sexithiophene_wH_neut_aniso_w_connectivity_cifstruct_chelpg.xyz'))
					#print mols
					for atom in mols[xlist][ylist][zlist]():
						atom().move(x*TV[0]+y*TV[1]+z*TV[2])
	print strftime("%a, %d %b %Y %X +0000", gmtime())
#	print 'mols', mols
#	print mols[size][size][size]()
#	for atom in mols[size][size][size]():
#		print atom
#Remove duplicated centre molecule
#	del (mols[size][size][size])
	print 'size', size
	print 'atoms:'
	for atom in Register.GetRegister('Atom'):
		print atom
	E0 = np.matrix([0.,0.,0.])
#	Reorg energy part
	for cen_atom, x_atom in zip(mols[size][size][size](), mols[size][size][size+1]()):
    		print(cen_atom, x_atom)
		print 'cen q:', cen_atom()._crg
		print 'x q:', x_atom()._crg
		cen_initial=cen_atom()._crg
		x_initial=x_atom()._crg
		print 'cenini, xini:', cen_initial, x_initial
		cen_atom()._crg=x_initial
		x_atom()._crg=cen_initial
		print 'cen q final:', cen_atom()._crg
		print 'x q final:', x_atom()._crg


exit()


"""	d = get_dipoles(E0=E0,cutoff=cut)
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



	f = open('../Crystals/TIPS_Pc_anion_neutral_2013_06_20_new.csv', 'a')
	f.write ('\nTIPS-Pc\t%s\t%s\tN/A\t%s\t%s\t%s\t%s\t%s\t%s' % (centre,surroundings,size,energyev,Uqq,Uqd,Udd,tot))
	f.flush()
	f.close()

#_new.dat for gnuplot

	f = open('TIPS_anion_neutral_size_%s_new.dat' % size, 'w')
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
	f.close() """

