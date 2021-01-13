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

name='neutprottest'

ReadMoleculeType('../Molecules/./TEST_neutron_pos000.xyz')
ReadMoleculeType('../Molecules/./TEST_neutron_pos000.xyz')

ReadMoleculeType('../Molecules/./TEST_neutron_pos000.xyz')
ReadMoleculeType('../Molecules/./TEST_neutron_pos000.xyz')

centre='TEST_neutron_pos000.xyz'
surroundings='TEST_neutron_pos000.xyz'

centreb='TEST_neutron_pos000.xyz'
surroundingsb='TEST_neutron_pos000.xyz'

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

print 'calc cell_vol, cif cell_vol:',cell_volume, ', ', cif_cell_volume

#Converts frac coords to carts
matrix_to_cartesian=np.matrix( [[a, b*cos(gamma), c*cos(beta)],
[0, b*sin(gamma), c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)],
[0, 0, c*cell_volume/sin(gamma)]])

#carts to frac
matrix_to_fractional=matrix_to_cartesian.I

#TVs, TV[0,1,2] are the three translation vectors.
TV=matrix_to_cartesian.T

print 'TV:'
print TV

print 'TV[0,1,2] are the three translation vectors.'



#Reflection Properties (see http://mathworld.wolfram.com/Reflection.html):
''' No reflection in Pc
a_refl=-1
b_refl=1
c_refl=-1
d_refl=1.5

#in fractional coords:
#perp vector to refl plane in fractional coords
n_unnorm_frac=np.matrix( [[a_refl, b_refl, c_refl]] )
#perp vector to refl plane in fractional coords, normalised
n_refl_frac=np.matrix( [[a_refl, b_refl, c_refl]] )/ np.sqrt(np.dot(np.matrix( [[a_refl, b_refl, c_refl]] ),np.matrix( [[a_refl, b_refl, c_refl]] ).T))
#dist from origin to refl plane in fractional coords, will be in direction of normal to plane
dist_refl_norm_frac=d_refl/np.sqrt(np.dot(np.matrix( [[a_refl, b_refl, c_refl]] ),np.matrix( [[a_refl, b_refl, c_refl]] ).T))

print 'n_refl_frac', n_refl_frac

print 'd_refl', d_refl

#in cartesians
# Perpendicular cartesian vector, 
n_refl=n_refl_frac*TV
n_refl_norm=n_refl/sqrt(np.dot(n_refl,n_refl.T))

print 'n_refl', n_refl
print 'n_refl_norm', n_refl_norm

'''



cut=8.0

# size given number of translation vectors ion each direction from central Pc

#mols=[]
#for size in np.arange(0,6,1):
#	for y in np.arange(-size,size+0.1,1):
#		for x in np.arange(-size,size+0.1,1):
#			for z in np.arange(-size,size+0.1,1):
#				mols.apped[Pc()]

for size in np.arange(1,2,1):
	mols = [[[] for i in range((2*size)+1)] for i in range((2*size)+1)]
	molsb = [[[] for i in range((2*size)+1)] for i in range((2*size)+1)]
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
					mols[xlist][ylist].append(GetMolecule('../Molecules/./TEST_neutron_pos000.xyz'))
					print 'xlist,ylist,zlist',xlist,ylist,zlist
					#flippedmol
					molsb[xlist][ylist].append(GetMolecule('../Molecules/./TEST_neutron_pos000.xyz'))
					#print mols
				else:
					# Places molecule of this type only if not in centre
					mols[xlist][ylist].append(GetMolecule('../Molecules/./TEST_neutron_pos000.xyz'))
					#print mols
					for atom in mols[xlist][ylist][zlist]():
						atom().move(x*TV[0]+y*TV[1]+z*TV[2])
						print 'atom pos:', atom()._pos
					#FLIPPED EXTRA MOL IN UNIT CELL
					molsb[xlist][ylist].append(GetMolecule('../Molecules/./TEST_neutron_pos000.xyz'))
					#print mols

	print strftime("%a, %d %b %Y %X +0000", gmtime())
#	print 'mols', mols
#	print mols[size][size][size]()
#	for atom in mols[size][size][size]():
#		print atom
#Remove duplicated centre molecule
#	del (mols[size][size][size])
	print 'size', size
	print 'atoms:'
	#for atom in Register.GetRegister('Atom'):
	#	print atom
	E0 = np.matrix([0.,0.,0.])
	
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

	#Reorg energy part
	'''#iterates over two together, cen_atom with charge, and x_atom moving to. Swaps charges.
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
		print 'x q final:', x_atom()._crg\
	'''

	#with other mol in unit cells
	print 'Moving charges...'
	for cen_atom, x_atom in zip(mols[size][size][size](), molsb[size][size][size]()):
    		#print(cen_atom, x_atom)
		#print 'cen q:', cen_atom()._crg
		#print 'x q:', x_atom()._crg
		cen_initial=cen_atom()._crg
		x_initial=x_atom()._crg
		#print 'cenini, xini:', cen_initial, x_initial
		cen_atom()._crg=x_initial
		x_atom()._crg=cen_initial
		#print 'cen q final:', cen_atom()._crg
		#print 'x q final:', x_atom()._crg\

	#Calc energy following reorg

	natom = get_atoms_with_polarizability()
	dipoles_multi = np.reshape(d,(natom,3))
	Efield_multi = np.reshape(get_electric_field(E0),(natom,3))

	atoms=np.arange(0,natom,1)
	Uqdiptot=0.0
	for atom in atoms:
		Uqdip = -(dipoles_multi[atom] * Efield_multi[atom].T)
		Uqdiptot += Uqdip

	newUqd=27.211*float(Uqdiptot[0])
	print 'New Uqd (eV): ', newUqd
	NewEnergyev = Udd+newUqd+Uqq
	print 'New Utot (eV): ',NewEnergyev
	reorgeV=newUqd-Uqd
	print 'Reorginisation Energy (eV): ',reorgeV

	time=strftime("%Y %b %d %X", gmtime())

#csv file



	f = open('../Crystals/crystals.csv', 'a')
	f.write ('\n%s\t%s\t%s\tN/A\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,centre,surroundings,size,energyev,Uqq,Uqd,Udd,tot[0,0],tot[0,1],tot[0,2]))
	f.flush()
	f.close()

#_new.dat for gnuplot
	print 'making .dat'


	f = open('pencen_neutral_neutral_size_%s_new.dat' % size, 'w')
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

#_new.dats for gnuplot with cross sections
	for x in np.arange(-size,size+1,1):
		xlist=x+size
		for y in np.arange(-size,size+1,1):
			ylist=y+size
			for z in np.arange(-size,size+1,1):
				zlist=z+size

				f = open('%s_size_%s_new_00%s.dat' % (name,size,xlist+1), 'w')
				f.write('Molecule Starts:\n')
				for atom in molsa[xlist][ylist][zlist]():
					astr=str(atom)
					f.write(astr)
					f.write('\n')		
				f.write('Molecule Ends.\n\n')
				f.write('\n\ndipoles start:\n')
				#Currently prints all dips, pol plotting script can tell which belong to which atoms based on numbers from pos'ns.
				for dd in split_d:
					dstr=str(dd)
					f.write(dstr)
					f.write('\n')
				f.flush()
				f.close()

				f = open('%s_size_%s_new_0%s0.dat' % (name,size,xlist+1), 'w')
				f.write('Molecule Starts:\n')
				for atom in molsa[ylist][xlist][zlist]():
					astr=str(atom)
					f.write(astr)
					f.write('\n')		
				f.write('Molecule Ends.\n\n')
				f.write('\n\ndipoles start:\n')
				#Currently prints all dips, pol plotting script can tell which belong to which atoms based on numbers from pos'ns.
				for dd in split_d:
					dstr=str(dd)
					f.write(dstr)
					f.write('\n')
				f.flush()
				f.close()

				f = open('%s_size_%s_new_%s0.dat' % (name,size,xlist+1), 'w')
				f.write('Molecule Starts:\n')
				for atom in molsa[zlist][ylist][xlist]():
					astr=str(atom)
					f.write(astr)
					f.write('\n')		
				f.write('Molecule Ends.\n\n')
				f.write('\n\ndipoles start:\n')
				#Currently prints all dips, pol plotting script can tell which belong to which atoms based on numbers from pos'ns.
				for dd in split_d:
					dstr=str(dd)
					f.write(dstr)
					f.write('\n')
				f.flush()
				f.close()

