import sys
sys.path.append('../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import *
from Polarizability.GetEnergyFromDips import *
from Polarizability.JMatrix import JMatrix
import numpy as np
from math import *
from time import gmtime, strftime

print strftime("%a, %d %b %Y %X +0000", gmtime())

name='PENCEN_cat_neut'

centre='PENCEN_cation_aniso_w_connectivity_cifstruct_chelpg.xyz'
surroundings='PENCEN_aniso_w_connectivity_cifstruct_chelpg.xyz'

centreb='PENCEN_2ndincell_aniso_w_connectivity_cifstruct_chelpg.xyz'
surroundingsb='PENCEN_2ndincell_aniso_w_connectivity_cifstruct_chelpg.xyz'

ReadMoleculeType('../Molecules/./%s' % centre)
ReadMoleculeType('../Molecules/./%s' % surroundings)

ReadMoleculeType('../Molecules/./%s' % centreb)
ReadMoleculeType('../Molecules/./%s' % surroundingsb)

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

# Place Molecules
for size in np.arange(0,1,1):
	mols = [[[] for i in range((2*size)+1)] for i in range((2*size)+1)]
	molsb = [[[] for i in range((2*size)+1)] for i in range((2*size)+1)]
	print mols
	for x in np.arange(-size,size+1,1):
		xlist=x+size
		for y in np.arange(-size,size+1,1):
			ylist=y+size
			for z in np.arange(-size,size+1,1):
				zlist=z+size
				if x == 0 and y == 0 and z == 0:
					#central mol
					mols[xlist][ylist].append(GetMolecule('../Molecules/./%s' % centre))
					print 'xlist,ylist,zlist',xlist,ylist,zlist
					#flippedmol
					molsb[xlist][ylist].append(GetMolecule('../Molecules/./%s' % centreb))
					#print mols
				else:
					# Places molecule of this type if not in centre
					mols[xlist][ylist].append(GetMolecule('../Molecules/./%s' % surroundings))
					#print mols
					for atom in mols[xlist][ylist][zlist]():
						atom().move(x*TV[0]+y*TV[1]+z*TV[2])
						print 'atom pos:', atom()._pos
					#FLIPPED EXTRA MOL IN UNIT CELL
					molsb[xlist][ylist].append(GetMolecule('../Molecules/./%s' % surroundingsb))

					for atom in molsb[xlist][ylist][zlist]():
						atom().move(x*TV[0]+y*TV[1]+z*TV[2])
						print 'atom pos:', atom()._pos


	#Calculate Properties:

	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'size', size
	E0 = np.matrix([0.,0.,0.])

	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'Calc jm'
	jm = JMatrix(cutoff=cut)
	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'Calc dips:'
	d = get_dipoles(E0=E0,jm=jm._m,cutoff=cut)
	print strftime("%a, %d %b %Y %X +0000", gmtime())
	Efield = get_electric_field(E0)
	potential = get_potential()

	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'dips', d
	split_d = split_dipoles_onto_atoms(d)
	tot = np.matrix([0.,0.,0.])
	for dd in split_d:
		tot += dd

	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'tot', tot
	Uqq = np.multiply(get_U_qq(potential=potential),27.211)
	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'Uqq', Uqq
	Uqd = np.multiply(get_U_qdip(jm=jm._m,dips=d,Efield=Efield),27.211)
	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'Uqd', Uqd
	Udd = np.multiply(get_U_dipdip(jm=jm._m,dips=d.T),27.211)
	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'Udd', Udd
	energyev = Udd+Uqd+Uqq
	print 'energyev', energyev
	energy=energyev/27.211

	print strftime("%a, %d %b %Y %X +0000", gmtime())
	print 'Making .dat cross sections for gnuplot'

#_new.dats for gnuplot with cross sections
	for x in np.arange(-size,size+1,1):
		xlist=x+size
		f = open('%s_size_%s_new_00%s.dat' % (name,size,xlist+1), 'a')
		g = open('%s_size_%s_new_0%s0.dat' % (name,size,xlist+1), 'a')
		h = open('%s_size_%s_new_%s00.dat' % (name,size,xlist+1), 'a')
		f.write('Molecule Starts:\n')
		g.write('Molecule Starts:\n')
		h.write('Molecule Starts:\n')
		f.flush()
		f.close()
		g.flush()
		g.close()
		h.flush()
		h.close()

		for y in np.arange(-size,size+1,1):
			ylist=y+size
			for z in np.arange(-size,size+1,1):
				zlist=z+size

				f = open('%s_size_%s_new_00%s.dat' % (name,size,xlist+1), 'a')
				for molecule in mols[xlist][ylist][zlist]():
					#this is a bit clunky, but puts pos'ns and atom nos of desied molecules into .dat file. Checks each atom to see if parent molecule is one we are interested in.
					molstr=str(molecule())
					for atom in Register.GetRegister('Atom'):
						astr=str(atom)
						if molstr in astr:
							f.write(astr)
							f.write('\n')

				for molecule in molsb[xlist][ylist][zlist]():
					molstr=str(molecule())
					for atom in Register.GetRegister('Atom'):
						astr=str(atom)
						if molstr in astr:
							f.write(astr)
							f.write('\n')	
				f.flush()
				f.close()

				f = open('%s_size_%s_new_0%s0.dat' % (name,size,xlist+1), 'a')
				for molecule in mols[ylist][xlist][zlist]():
					molstr=str(molecule())
					for atom in Register.GetRegister('Atom'):
						astr=str(atom)
						if molstr in astr:
							f.write(astr)
							f.write('\n')
				for molecule in molsb[ylist][xlist][zlist]():
					molstr=str(molecule())
					for atom in Register.GetRegister('Atom'):
						astr=str(atom)
						if molstr in astr:
							f.write(astr)
							f.write('\n')		
				f.flush()
				f.close()

				f = open('%s_size_%s_new_%s00.dat' % (name,size,xlist+1), 'a')
				for molecule in mols[zlist][ylist][xlist]():
					molstr=str(molecule())
					for atom in Register.GetRegister('Atom'):
						astr=str(atom)
						if molstr in astr:
							f.write(astr)
							f.write('\n')	
				for molecule in molsb[zlist][ylist][xlist]():
					molstr=str(molecule())
					for atom in Register.GetRegister('Atom'):
						astr=str(atom)
						if molstr in astr:
							f.write(astr)
							f.write('\n')		
				f.flush()
				f.close()

		#Currently prints all dips, pol plotting script can tell which belong to which atoms based on numbers from pos'ns.
		f = open('%s_size_%s_new_00%s.dat' % (name,size,xlist+1), 'a')
		g = open('%s_size_%s_new_0%s0.dat' % (name,size,xlist+1), 'a')
		h = open('%s_size_%s_new_%s00.dat' % (name,size,xlist+1), 'a')
		f.write('Molecule Ends.\n\n\n\ndipoles start:\n')
		g.write('Molecule Ends.\n\n\n\ndipoles start:\n')
		h.write('Molecule Ends.\n\n\n\ndipoles start:\n')
		for dd in split_d:
			dstr=str(dd)
			f.write(dstr)
			f.write('\n')
			g.write(dstr)
			g.write('\n')
			h.write(dstr)
			h.write('\n')
		f.flush()
		f.close()
		g.flush()
		g.close()
		h.flush()
		h.close()

	#REORG ENERGY with other mol in unit cells
	print 'Moving charges...'
	for x in np.arange(-size,size+1,1):
		print 'entered reorg loop'
		xlist=x+size
		print 'x', x
		print 'xlist', xlist
		print xlist
		for y in np.arange(-size,size+1,1):
			ylist=y+size
			for z in np.arange(-size,size+1,1):
				zlist=z+size
				print 'swapping charges of molecule at 0,0,0, (mols([',size,'][',size,'][',size,']) and molecule at ', x, ',', y, ',', z, '.'
				for cen_atom, x_atom in zip(mols[size][size][size](), mols[xlist][ylist][zlist]()):
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

				Efield=get_electric_field(E0)
				newUqd = np.multiply(get_U_qdip(jm=jm._m,dips=d,Efield=Efield),27.211)

				print 'New Uqd (eV): ', newUqd
				NewEnergyev = Udd+newUqd+Uqq
				print 'New Utot (eV): ',NewEnergyev
				reorgeV=newUqd-Uqd
				print 'Reorginisation Energy (', x, ',', y, ',', z, ') (eV): ',reorgeV

				time=strftime("%Y %b %d %X", gmtime())
				print time

				#Put charge back
				for cen_atom, x_atom in zip(mols[size][size][size](), mols[xlist][ylist][zlist]()):
					cen_initial=cen_atom()._crg
					x_initial=x_atom()._crg
					cen_atom()._crg=x_initial
					x_atom()._crg=cen_initial


				#Repeat for molsb
				print 'swapping charges of molecule at 0,0,0, (mols([',size,'][',size,'][',size,']) and molecule b at ', x, ',', y, ',', z, '.'
				for cen_atom, x_atom in zip(mols[size][size][size](), molsb[xlist][ylist][zlist]()):
					cen_initial=cen_atom()._crg
					x_initial=x_atom()._crg
					cen_atom()._crg=x_initial
					x_atom()._crg=cen_initial

					#Calc energy following reorg

				#Calc energy following reorg

				Efield=get_electric_field(E0)
				newUqd = np.multiply(get_U_qdip(jm=jm._m,dips=d,Efield=Efield),27.211)

				print 'New Uqd (eV): ', newUqd
				NewEnergyev = Udd+newUqd+Uqq
				print 'New Utot (eV): ',NewEnergyev
				reorgeV=newUqd-Uqd
				print 'Reorginisation Energy (', x, ',', y, ',', z, ') (eV): ',reorgeV

				time=strftime("%Y %b %d %X", gmtime())
				print time

				#put charge back
				for cen_atom, x_atom in zip(mols[size][size][size](), molsb[xlist][ylist][zlist]()):
					cen_initial=cen_atom()._crg
					x_initial=x_atom()._crg
					cen_atom()._crg=x_initial
					x_atom()._crg=cen_initial

#csv file



	f = open('../Crystals/crystals.csv', 'a')
	f.write ('\n%s\t%s\t%s\tN/A\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,centre,surroundings,size,energyev,Uqq,Uqd,Udd,tot[0,0],tot[0,1],tot[0,2]))
	f.flush()
	f.close()

#_new.dat for gnuplot
	print 'making .dat'


	f = open('%s_size_%s_new.dat' % (name, size), 'w')
	f.write('Molecule Starts:\n')
	for atom in Register.GetRegister('Atom'):
		astr=str(atom)
		f.write(astr)
		f.write('\n')
	f.write('Molecule Ends.\n\n')
	f.write ("Crystal Struct\tCentre\tSurroundings\tSize\tEnergy\tUqq\tUqd\tUdd\tDipole Moment")
	f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (name,centre,surroundings,size,energyev,Uqq,Uqd,Udd,tot))
	f.write('\n\ndipoles start:\n')
	for dd in split_d:
		dstr=str(dd)
		f.write(dstr)
		f.write('\n')
	f.flush()
	f.close()

