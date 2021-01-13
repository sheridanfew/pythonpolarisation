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

ReadMoleculeType('../Molecules/./sexithiophene_wH_neut_aniso_w_connectivity_cifstruct_chelpg.xyz')
ReadMoleculeType('../Molecules/./sexithiophene_wH_anion_aniso_w_connectivity_cifstruct_chelpg.xyz')

centre='sexithiophene_wH_anion_aniso_w_connectivity_cifstruct_chelpg.xyz'
surroundings='sexithiophene_wH_neut_aniso_w_connectivity_cifstruct_chelpg.xyz'

#Get translation vectors:

a=44.7086
b=7.8513
c=6.0292

alpha=90*(pi/180)
beta=90.762*(pi/180)
gamma=90*(pi/180)

cell_volume=sqrt(1 - (cos(alpha)**2) - (cos(beta)**2) - (cos(gamma)**2) + (2*cos(alpha)*cos(beta)*cos(gamma)))

#Converts frac coords to carts
matrix_to_cartesian=np.matrix( [[a, b*cos(gamma), c*cos(beta)],
[0, b*sin(gamma), c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)],
[0, 0, c*cell_volume/sin(gamma)]])

#carts to frac
matrix_to_fractional=matrix_to_cartesian.I

TV=matrix_to_cartesian.T

#Reflection Properties (see http://mathworld.wolfram.com/Reflection.html):

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
#TVs, TV[0,1,2] are the three translation vectors.

print 'TV:'
print TV

print 'TV[0,1,2] are the three translation vectors.'

#in cartesians
# Perpendicular cartesian vector, 
n_refl=n_refl_frac*TV
n_refl_norm=n_refl/sqrt(np.dot(n_refl,n_refl.T))

print 'n_refl', n_refl
print 'n_refl_norm', n_refl_norm

cut=8.0

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

for size in np.arange(0,1,1):
	mols = [[[] for i in range((2*size)+1)] for i in range((2*size)+1)]
	molsflip = [[[] for i in range((2*size)+1)] for i in range((2*size)+1)]
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
					#flippedmol
					molsflip[xlist][ylist].append(GetMolecule('../Molecules/./sexithiophene_wH_neut_aniso_w_connectivity_cifstruct_chelpg.xyz'))
					#print mols
					for atom in molsflip[xlist][ylist][zlist]():
						#Flip molecule coords
						print 'ATOM POSITION FLIP:'
						print 'initial atom pos', atom()._pos
						# Convert coords of atom (preflip) to fractional coords
						atom_pos_fract=(matrix_to_fractional*atom()._pos.T).T
						print 'Convert coords of atom (preflip) to fractional coords. atom_pos_fract', atom_pos_fract
						# Dist to plane, see http://mathworld.wolfram.com/Point-PlaneDistance.html
						D = ( (np.dot(n_unnorm_frac,atom_pos_fract.T)) + d_refl)/sqrt(np.dot(n_unnorm_frac,n_unnorm_frac.T))
						#print 'dist to plane, D = ( (np.dot(n_unnorm_frac,atom_pos_fract.T)) + d_refl)/sqrt(np.dot(n_unnorm_frac,n_unnorm_frac.T))'
						#print 'n_unnorm_frac, atom_pos_fract , d_refl', n_unnorm_frac, atom_pos_fract , d_refl
						#print 'np.dot(n_unnorm_frac,n_unnorm_frac.T)',np.dot(n_unnorm_frac,n_unnorm_frac.T)
					
						print 'Dist to plane', D
						movevec_frac=-2*D*n_refl_frac
						print 'Vector by which atom must move in fractional coords for reflection, movevec_frac', movevec_frac
						movevec=matrix_to_cartesian*movevec_frac.T
						print 'Vector by which atom must move in cartesians coords for reflection, movevecc', movevec
						print 'old type(atom()._pos)', type(atom()._pos)
						atom().move(movevec.T)
						print 'new type(atom()._pos)', type(atom()._pos)
						print 'New position of atom (Cartesian)',atom()._pos

						# Dold = projection of position of molecule onto vector running perp. to the plane. Gives distance from point to plane if plane ran through origin
						# Dold=(np.dot(n_refl_frac,atom_pos_fract.T))

						#Polarisation tensor flip. NB. This method will only work for a diagonal tensor! Not thought about nondiag!
						#Convert diagonal elements to vector, then rotate in plane as before, except in this case move plane to cross origin. (ie set d_refl=0)

						print 'POLARISABILITY FLIP:'
						print 'atom()._pol', atom()._pol


						#Vector containing diagonal elements of polarisability tensor
						atom_pol_diags=np.matrix( [[ atom()._pol[0,0], atom()._pol[1,1], atom()._pol[2,2] ]])
						print 'Vector containing diagonal elements of polarisability tensor, atom_pol_diags', atom_pol_diags

						# Dist to plane through origin, see http://mathworld.wolfram.com/Point-PlaneDistance.html
						print 'n_refl_norm, atom_pol_diags', n_refl_norm, atom_pol_diags
						Dpol=(np.dot(n_refl_norm,atom_pol_diags.T))
						print 'Dist to plane through origin, Dpol', Dpol
						# Reflect in plane -> New diagonals
						atom_pol_diags_flipped=(atom_pol_diags-(2*Dpol*n_refl_norm)).T
						print 'Reflect in plane -> New diagonals in cartesian coords, atom_pol_diags_flipped (N.B -ve of vector has been taken for +ve polarisabilities)', atom_pol_diags_flipped

						print 'OLD type(atom()._pol)', type(atom()._pol)
						print 'type(atom()._pol[0,0])', type(atom()._pol[0,0])	
						print 'type(atom()._pol[0,1])', type(atom()._pol[0,1])	


						atom()._pol=Polarizability(noniso =[ [ np.sqrt(atom_pol_diags_flipped[0,0]**2),0.0,0.0], [0.0,np.sqrt(atom_pol_diags_flipped[1,0]**2),0.0], [0.0,0.0,np.sqrt(atom_pol_diags_flipped[2,0]**2)] ])			
						print 'NEW type(atom()._pol)', type(atom()._pol)
						print 'type(atom()._pol[0,0])', type(atom()._pol[0,0])	
						print 'type(atom()._pol[0,1])', type(atom()._pol[0,1])	

						'''ALTAPPROACH
						#Convert this to fractional coords
						atom_pol_diags_frac=(matrix_to_fractional* atom_pol_diags.T)
						print 'Convert this to fractional coords, atom_pol_diags_frac', atom_pol_diags_frac
						# Dist to plane through origin, see http://mathworld.wolfram.com/Point-PlaneDistance.html
						print 'n_refl_frac, atom_pol_diags_frac', n_refl_frac, atom_pol_diags_frac
						Dpol_frac=(np.dot(n_refl_frac,atom_pol_diags_frac))
						print 'Dist to plane through origin, Dpol_frac', Dpol_frac
						# Reflect in plane -> New diagonals in frac coords
						atom_pol_diags_frac_flipped=-(atom_pol_diags_frac.T-(2*Dpol_frac*n_refl_frac)).T
						print 'Reflect in plane -> New diagonals in frac coords, atom_pol_diags_frac_flipped (N.B -ve of vector has been taken for +ve polarisabilities)', atom_pol_diags_frac_flipped
						# Convert back to Cartesians
						atom_pol_diags_flipped=matrix_to_cartesian*atom_pol_diags_frac_flipped
						print 'Convert back to Cartesians, atom_pol_diags_flipped', atom_pol_diags_flipped
			
						newpol='Polarizability(noniso =[ [', str(atom_pol_diags_flipped[0]),',0.0,0.0], [0.0,',atom_pol_diags_flipped[1],',0.0], [0.0,0.0,',atom_pol_diags_flipped[2],'] ])}'
						print 'newpol'
						atom()._pol=Polarizability(noniso =[ [ atom_pol_diags_flipped[0,0],0.0,0.0], [0.0,atom_pol_diags_flipped[1,0],0.0], [0.0,0.0,atom_pol_diags_flipped[2,0]] ])
						'''


						print 'new atom()._pol', atom()._pol
				else:
					# Places molecule of this type only if not in centre
					mols[xlist][ylist].append(GetMolecule('../Molecules/./sexithiophene_wH_neut_aniso_w_connectivity_cifstruct_chelpg.xyz'))
					#print mols
					for atom in mols[xlist][ylist][zlist]():
						atom().move(x*TV[0]+y*TV[1]+z*TV[2])
						print 'atom pos:', atom()._pos
					#FLIPPED EXTRA MOL IN UNIT CELL
					molsflip[xlist][ylist].append(GetMolecule('../Molecules/./sexithiophene_wH_neut_aniso_w_connectivity_cifstruct_chelpg.xyz'))
					#print mols
					for atom in molsflip[xlist][ylist][zlist]():
						#Flip molecule coords
						print 'ATOM POSITION FLIP:'
						#print 'initial atom pos', atom()._pos
						# Convert coords of atom (preflip) to fractional coords
						atom_pos_fract=(matrix_to_fractional*atom()._pos.T).T
						#print 'Convert coords of atom (preflip) to fractional coords. atom_pos_fract', atom_pos_fract
						# Dist to plane, see http://mathworld.wolfram.com/Point-PlaneDistance.html
						D = ( (np.dot(n_unnorm_frac,atom_pos_fract.T)) + d_refl)/sqrt(np.dot(n_unnorm_frac,n_unnorm_frac.T))
						#print 'Dist to plane', D
						movevec_frac=-2*D*n_refl_frac
						#print 'Vector by which atom must move in fractional coords for reflection, movevec_frac', movevec_frac
						movevec=matrix_to_cartesian*movevec_frac.T
						#print 'Vector by which atom must move in cartesians coords for reflection, movevecc', movevec
						atom().move(movevec.T)
						#print 'New position of atom (Cartesian)',atom()._pos

						# Dold = projection of position of molecule onto vector running perp. to the plane. Gives distance from point to plane if plane ran through origin
						# Dold=(np.dot(n_refl_frac,atom_pos_fract.T))

						#Polarisation tensor flip. NB. This method will only work for a diagonal tensor! Not thought about nondiag!
						#Convert diagonal elements to vector, then rotate in plane as before, except in this case move plane to cross origin. (ie set d_refl=0)

						print 'POLARISABILITY FLIP:'
						#print 'atom()._pol', atom()._pol

						
						#Vector containing diagonal elements of polarisability tensor
						atom_pol_diags=np.matrix( [[ atom()._pol[0,0], atom()._pol[1,1], atom()._pol[2,2] ]])
						print 'Vector containing diagonal elements of polarisability tensor, atom_pol_diags', atom_pol_diags
						#Convert this to fractional coords
						atom_pol_diags_frac=(matrix_to_fractional* atom_pol_diags.T)
						print 'Convert this to fractional coords, atom_pol_diags_frac', atom_pol_diags_frac
						# Dist to plane through origin, see http://mathworld.wolfram.com/Point-PlaneDistance.html
						print 'n_refl_frac, atom_pol_diags_frac', n_refl_frac, atom_pol_diags_frac
						Dpol_frac=(np.dot(n_refl_frac,atom_pol_diags_frac))
						print 'Dist to plane through origin, Dpol_frac', Dpol_frac
						# Reflect in plane -> New diagonals in frac coords
						atom_pol_diags_frac_flipped=-(atom_pol_diags_frac.T-(2*Dpol_frac*n_refl_frac)).T
						print 'Reflect in plane -> New diagonals in frac coords, atom_pol_diags_frac_flipped (N.B -ve of vector has been taken for +ve polarisabilities)', atom_pol_diags_frac_flipped
						# Convert back to Cartesians
						atom_pol_diags_flipped=matrix_to_cartesian*atom_pol_diags_frac_flipped
						print 'Convert back to Cartesians, atom_pol_diags_flipped', atom_pol_diags_flipped
			

						atom()._pol=Polarizability(noniso =[ [ atom_pol_diags_flipped[0,0],0.0,0.0], [0.0,atom_pol_diags_flipped[1,0],0.0], [0.0,0.0,atom_pol_diags_flipped[2,0]] ])
						print 'newpol',atom()._pol

						atom().move(x*TV[0]+y*TV[1]+z*TV[2])

						#print 'new atom()._pol', atom()._pol
#atom pos:', atom()._pos



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
	'''	#Reorg energy part
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

	print 'N.B. Second molecule in unit cell is in wrong orientation, eed to change n_refl to n_refl_norm and fraco distance (currently calculated in fractional coordinated, need to get to cartesians, (multiply frac_to_cart_matrix*dist_refl_ frac*n_refl_frac should give vector of correct length!)'

#flip plane is in Cartesian coords, should be same but in fractional cell coords' '''
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



#	f = open('../Crystals/TIPS_Pc_anion_neutral_2013_06_20_new.csv', 'a')
#	f.write ('\nTIPS-Pc\t%s\t%s\tN/A\t%s\t%s\t%s\t%s\t%s\t%s' % (centre,surroundings,size,energyev,Uqq,Uqd,Udd,tot))
#	f.flush()
#	f.close()

#_new.dat for gnuplot
	print 'making .dat'


	f = open('sexithiophene_anion_neutral_size_%s_new.dat' % size, 'w')
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

