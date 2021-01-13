import sys
sys.path.append('../')
from BasicElements import *
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

print strftime("%a, %d %b %Y %X +0000", gmtime())

###################################
#START OF MOLECULE SPECIFIC SECTION
###################################
name='Sexithiophene_anion_neut'

size=0

mols_cen=['sexithiophene_mola_anion_aniso_cifstruct_chelpg.xyz','sexithiophene_molb_neut_aniso_cifstruct_chelpg.xyz','sexithiophene_molc_neut_aniso_cifstruct_chelpg.xyz','sexithiophene_mold_neut_aniso_cifstruct_chelpg.xyz']

mols_sur=['sexithiophene_mola_neut_aniso_cifstruct_chelpg.xyz','sexithiophene_molb_neut_aniso_cifstruct_chelpg.xyz','sexithiophene_molc_neut_aniso_cifstruct_chelpg.xyz','sexithiophene_mold_neut_aniso_cifstruct_chelpg.xyz']

#From cif:
'''
_symmetry_equiv_pos_as_xyz
1 x,y,z
2 1/2-x,1/2+y,1/2-z
3 -x,-y,-z
4 1/2+x,1/2-y,1/2+z
_cell_length_a                   44.708(6)
_cell_length_b                   7.851(3)
_cell_length_c                   6.029(2)
_cell_angle_alpha                90
_cell_angle_beta                 90.76(2)
_cell_angle_gamma                90
_cell_volume                     2116.01
'''
#Get translation vectors:

a=44.7086/0.5291772109217
b=7.8513/0.5291772109217
c=6.0292/0.5291772109217

alpha=90*(pi/180)
beta=90.762*(pi/180)
gamma=90*(pi/180)

cif_unit_cell_volume=2116.01/(a*b*c*(0.5291772109217**3))

###################################
#END OF MOLECULE SPECIFIC SECTION
###################################


cell_volume=sqrt(1 - (cos(alpha)**2) - (cos(beta)**2) - (cos(gamma)**2) + (2*cos(alpha)*cos(beta)*cos(gamma)))

print 'calc cell_volume, cif cell_volume:',cell_volume, ', ', cif_unit_cell_volume

if cell_volume > (cif_unit_cell_volume*1.1) or cell_volume < (cif_unit_cell_volume*0.9):
	print 'angles wrong?'
	exit (2)

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

cut=8.0

# Place Molecules

sexi_cry=Crystal(name='sexi_an_neut',mols_cen=mols_cen,mols_sur=mols_sur,TVs=TV,cenpos=[1,1,1],maxTVs,lena=2,lenb=2,lenc=2)

#sexi_cry._mols contains all molecules.
#mols[0] contains a list of all molecules in position a, mols[1] all mols in pos'n b, etc.
#mols[0][x,y,z] contains molecule a in position x,y,z
#mols may as such be iterated over in a number of ways to consider different molecules.


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

# print dipoles
f = open('%s_size_%s_dipoles.dat' % (name,size), 'w')
for dd in split_d:
	dstr=str(dd)
	f.write(dstr)
	f.write('\n')
f.flush()
f.close()


# new.dats for gnuplot with individual molecules
for x in np.arange(-size,size+1,1):
	xlist=x+size
	for y in np.arange(-size,size+1,1):
		ylist=y+size
		for z in np.arange(-size,size+1,1):
			zlist=z+size
			for molincell in np.arange(0,len(centres),1):
				if mols[molincell][zlist][ylist][xlist] != 'empty':
		                        f = open('%s_size_%s_%s%s%s_mol%s.dat' % (name,size,xlist,ylist,zlist,molincell), 'w')
		                        f.write('Molecule Starts:\n')
		                        for molecule in mols[molincell][zlist][ylist][xlist]():
		                                molstr=str(molecule())
		                                for atom in Register.GetRegister('Atom'):
		                                        astr=str(atom)
		                                        if molstr in astr:
		                                                f.write(astr)
		                                                f.write('\n')
		                        f.write('Molecule Ends.')

		                        f.flush()
		                        f.close()

'''
#REORG ENERGY with other mol in unit cells
reorg=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(len(centres))]
Uqdlist=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(len(centres))]
Uqqlist=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(len(centres))]
#reorg[i][x][y][z] contains this property for charge of central molecule moved to molecule in poition [i] in unit cell at [x][y][z].
#Other matrices same for each peroperty

print 'Moving charges...'
for x in np.arange(-size,size+1,1):
	xlist=x+size
	for y in np.arange(-size,size+1,1):
		ylist=y+size
		for z in np.arange(-size,size+1,1):
			zlist=z+size
			if mols[molincell][zlist][ylist][xlist] != 'empty':
        	                for molincell in np.arange(0,len(centres),1):
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
		                        f.write('\n\nReorganisation energy:\t%s eV\n' % str(reorgeV) )
		                        f.write('Polaron Binding\t%s eV\n' % str(Uqd+Udd) )
		                        f.write('Other Properties:\n')
		                        f.write ('Date\tname\tcentres[0]\tsurroundings[0]\tsize\tpos_a\tpos_b\tpos_c\tNewEnergyev\treorgeV\tUqd\tUqq\tUdd')
		                        f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,centres[0],surroundings[0],size,x,y,z,molincell,NewEnergyev,reorgeV,newUqd,newUqq,Udd))
		                        f.flush()
		                        f.close()

		                        f = open('../../../Crystals/reorg_energies.csv', 'a')
		                        f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,centres[0],surroundings[0],size,x,y,z,molincell,NewEnergyev,reorgeV,newUqd,newUqq,Udd))
		                        f.flush()
		                        f.close()

		                        f = open('./%s_size_%s_new.csv' % (name,size), 'w')
		                        f.write ('Date\tname\tcentres[0]\tsurroundings[0]\tsize\tpos_a\tpos_b\tpos_c\tNewEnergyev\treorgeV\tUqd\tUqq\tUdd')
		                        f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (time,name,centres[0],surroundings[0],size,x,y,z,molincell,energyev,reorgeV,Uqd,newUqq,Udd))
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
f.write ('\n%s\t%s\t%s\tN/A\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,centres[0],surroundings[0],size,energyev,Uqq,Uqd,Udd,tot[0,0],tot[0,1],tot[0,2]))
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
f.write ("Crystal Struct\tcentres[0]\tsurroundings[0]\tSize\tEnergy\tUqq\tUqd\tUdd\tDipole Moment")
f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (name,centres[0],surroundings[0],size,energyev,Uqq,Uqd,Udd,tot))
f.flush()
f.close()

'''
