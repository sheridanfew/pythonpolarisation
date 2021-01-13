for state in cation neut anion
do
  for size in 0 1 2
  do
    cat << EOF > C8_BTBT_${state}_neut_size${size}.py
import sys
sys.path.append('../../../')
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

name='C8_BTBT_${state}_neut'

centrea='C8_BTBT_mola_${state}_aniso_cifstruct_chelpg.xyz'
surroundingsa='C8_BTBT_mola_neut_aniso_cifstruct_chelpg.xyz'

centreb='C8_BTBT_molb_neut_aniso_cifstruct_chelpg.xyz'
surroundingsb='C8_BTBT_molb_neut_aniso_cifstruct_chelpg.xyz'

nmolsincell=2

ReadMoleculeType('../../../Molecules/./%s' % centrea)
ReadMoleculeType('../../../Molecules/./%s' % surroundingsa)

ReadMoleculeType('../../../Molecules/./%s' % centreb)
ReadMoleculeType('../../../Molecules/./%s' % surroundingsb)

#From cif:
'''
C8_BTBT
_cell_length_a                   5.927(7)
_cell_length_b                   7.88(1)
_cell_length_c                   29.18(4)
_cell_angle_alpha                90
_cell_angle_beta                 92.443(4)
_cell_angle_gamma                90
_cell_volume                     1361.61

'''
#Get translation vectors:

a=5.9277/0.5291772109217
b=7.881/0.5291772109217
c=29.184/0.5291772109217

alpha=90*(pi/180)
beta=92.4434*(pi/180)
gamma=90*(pi/180)

cell_volume=sqrt(1 - (cos(alpha)**2) - (cos(beta)**2) - (cos(gamma)**2) + (2*cos(alpha)*cos(beta)*cos(gamma)))
cif_cell_volume=1361.61/(0.5291772109217**3)


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

# Place Molecules
for size in np.arange($size,$(( $size + 1 )),1):
	mols=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(nmolsincell)]
	#mols contains all molecules.
	#mols[0] contains a list of all molecules in position a, mols[1] all mols in pos'n b, etc.
	#mols[0][x,y,z] contains molecule a in position x,y,z
	#mols may as such be iterated over in a number of ways to consider different molecules.

	print mols[0]
	for x in np.arange(-size,size+1,1):
		xlist=x+size
		for y in np.arange(-size,size+1,1):
			ylist=y+size
			for z in np.arange(-size,size+1,1):
				zlist=z+size
				if x == 0 and y == 0 and z == 0:
					#central mol
					mols[0][xlist][ylist].append(GetMolecule('../../../Molecules/./%s' % centrea))
					print 'xlist,ylist,zlist',xlist,ylist,zlist
					#flippedmol
					mols[1][xlist][ylist].append(GetMolecule('../../../Molecules/./%s' % centreb))
					#print mols[0]
				else:
					# Places molecule of this type if not in centrea
					mols[0][xlist][ylist].append(GetMolecule('../../../Molecules/./%s' % surroundingsa))
					mols[1][xlist][ylist].append(GetMolecule('../../../Molecules/./%s' % surroundingsb))

					for molincell in np.arange(0,nmolsincell,1):
						for atom in mols[molincell][xlist][ylist][zlist]():
							atom().move(x*TV[0]+y*TV[1]+z*TV[2])


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

#_new.dats for gnuplot with individual mols
	for x in np.arange(-size,size+1,1):
		xlist=x+size
		for y in np.arange(-size,size+1,1):
			ylist=y+size
			for z in np.arange(-size,size+1,1):
				zlist=z+size
				for molincell in np.arange(0,nmolsincell,1):
					f = open('%s_size_%s_%s%s%s_mol%s.dat' % (name,size,xlist,ylist,zlist,molincell), 'w')
					f.write('Molecule Starts:\n')
					for molecule in mols[molincell][zlist][ylist][xlist]():
						molstr=str(molecule())
						for atom in Register.GetRegister('Atom'):
							astr=str(atom)
							if molstr in astr:
								f.write(astr)
								f.write('\n')

					f.write('Molecule Ends.\n\n\n\ndipoles start:\n')

					for dd in split_d:
						dstr=str(dd)
						f.write(dstr)
						f.write('\n')
					f.write('Dipoles End.\n')
					f.flush()
					f.close()

	#REORG ENERGY with other mol in unit cells
	reorg=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(nmolsincell)]
	Uqdlist=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(nmolsincell)]
	Uqqlist=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(nmolsincell)]
	#reorg[i][x][y][z] contains this property for charge of central molecule moved to molecule in poition [i] in unit cell at [x][y][z].
	#Other matrices same for each peroperty

	print 'Moving charges...'
	for x in np.arange(-size,size+1,1):
		xlist=x+size
		for y in np.arange(-size,size+1,1):
			ylist=y+size
			for z in np.arange(-size,size+1,1):
				zlist=z+size
				for molincell in np.arange(0,nmolsincell,1):
					print 'swapping charges of molecule at 0,0,0, (mols[0]([',size,'][',size,'][',size,']) and molecule',molincell,' at ', x, ',', y, ',', z, '.'
					for cen_atom, x_atom in zip(mols[0][size][size][size](), mols[molincell][xlist][ylist][zlist]()):
						cen_initial=cen_atom()._crg
						x_initial=x_atom()._crg
						cen_atom()._crg=x_initial
						x_atom()._crg=cen_initial

					#Calc energy following reorg

					Efield=get_electric_field(E0)
					newUqd = np.multiply(get_U_qdip(jm=jm._m,dips=d,Efield=Efield),27.211)

					print 'New Uqd (eV): ', newUqd
					NewEnergyev = Udd+newUqd+Uqq
					print 'New Utot (eV): ',NewEnergyev
					reorgeV=newUqd-Uqd
					print 'Reorginisation Energy (', x, ',', y, ',', z, ') (eV): ',reorgeV
					newpotential = get_potential()
					newUqq = np.multiply(get_U_qq(potential=potential),27.211)
					print 'New Uqq (eV) (not used for reorg energy calc as should be same in perfect cryst): ',newUqq
					print 'Change in Uqq (eV): ',(newUqq-Uqq)

					reorg[molincell][xlist][ylist].append(newUqd-Uqd)
					Uqdlist[molincell][xlist][ylist].append(newUqd)
					Uqqlist[molincell][xlist][ylist].append(newUqd-Uqd)

					time=strftime("%Y %b %d %X", gmtime())
					print time

					f = open('%s_size_%s_%s%s%s_mol%s.dat' % (name,size,xlist,ylist,zlist,molincell), 'a')
					f.write('Reorganisation energy:\t%s eV\n' % str(reorgeV) )
					f.write('Polaron Binding\t%s eV\n' % str(Uqd+Uqq) )
					f.write('Other Properties:\n')
					f.write ('Date\tname\tcentrea\tsurroundingsa\tsize\tpos_a\tpos_b\tpos_c\tNewEnergyev\treorgeV\tUqd\tUqq\tUdd')
					f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,centrea,surroundingsa,size,x,y,z,molincell,NewEnergyev,reorgeV,newUqd,newUqq,Udd))
					f.flush()
					f.close()

					f = open('../../../Crystals/reorg_energies.csv', 'a')
					f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,centrea,surroundingsa,size,x,y,z,molincell,NewEnergyev,reorgeV,newUqd,newUqq,Udd))
					f.flush()
					f.close()

					f = open('./%s_size_%s_new.csv' % (name,size), 'w')
					f.write ('Date\tname\tcentrea\tsurroundingsa\tsize\tpos_a\tpos_b\tpos_c\tNewEnergyev\treorgeV\tUqd\tUqq\tUdd')
					f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (time,name,centrea,surroundingsa,size,x,y,z,molincell,energyev,reorgeV,Uqd,newUqq,Udd))
					f.flush()
					f.close()

					#Put charge back
					for cen_atom, x_atom in zip(mols[0][size][size][size](), mols[molincell][xlist][ylist][zlist]()):
						cen_initial=cen_atom()._crg
						x_initial=x_atom()._crg
						cen_atom()._crg=x_initial
						x_atom()._crg=cen_initial


#csv file



	f = open('../../../Crystals/crystals.csv', 'a')
	f.write ('\n%s\t%s\t%s\tN/A\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,centrea,surroundingsa,size,energyev,Uqq,Uqd,Udd,tot[0,0],tot[0,1],tot[0,2]))
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
	f.write ("Crystal Struct\tcentrea\tsurroundingsa\tSize\tEnergy\tUqq\tUqd\tUdd\tDipole Moment")
	f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (name,centrea,surroundingsa,size,energyev,Uqq,Uqd,Udd,tot))
	f.write('\n\ndipoles start:\n')
	for dd in split_d:
		dstr=str(dd)
		f.write(dstr)
		f.write('\n')
	f.flush()
	f.close()
EOF

  done
done
