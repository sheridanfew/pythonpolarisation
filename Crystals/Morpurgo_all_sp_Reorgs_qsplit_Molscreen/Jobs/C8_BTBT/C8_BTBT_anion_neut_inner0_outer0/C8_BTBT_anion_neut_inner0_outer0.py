import sys
sys.path.append('../../../../../')
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

qdict={"anion": -1.0, "neut": 0.0, "cation": 1.0}


name='C8_BTBT_anion_neut_inner0_outer0'
#For crystals here, all cubic and centred at centre
insize=0
#number of TVs in each dir central mol is from edge of inner region
outsize=0
state='anion'
mols_cen=['sp_C8_BTBT_mola_neut.xyz','sp_C8_BTBT_molb_neut.xyz']
mols_sur=['sp_C8_BTBT_mola_neut.xyz','sp_C8_BTBT_molb_neut.xyz']
mols_outer=['sp_C8_BTBT_mola_neut.xyz','sp_C8_BTBT_molb_neut.xyz']

Natoms=72

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

cif_unit_cell_volume=1361.61/(a*b*c*(0.5291772109217**3))
cell_volume=sqrt(1 - (cos(alpha)**2) - (cos(beta)**2) - (cos(gamma)**2) + (2*cos(alpha)*cos(beta)*cos(gamma)))

#Converts frac coords to carts
matrix_to_cartesian=np.matrix( [[a, b*cos(gamma), c*cos(beta)],
[0, b*sin(gamma), c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)],
[0, 0, c*cell_volume/sin(gamma)]])

#carts to frac
matrix_to_fractional=matrix_to_cartesian.I

#TVs, TV[0,1,2] are the three translation vectors.
TV=matrix_to_cartesian.T

cut=8.0

totsize=insize+outsize
#number of TVs in each dir nearest c inner mol is from edge of outer region
cenpos=[totsize,totsize,totsize]
length=[2*totsize+1,2*totsize+1,2*totsize+1]
maxTVs=insize
outer_maxTVs=insize+outsize
#for diamond outer, don't specify for cube and will fill to cube edges.

print 'name: ',name,'mols_cen: ', mols_cen,' mols_sur: ',mols_sur,' TVs: ', TV 

# Place Molecules

prot_neut_cry=Crystal(name=name,mols_cen=mols_cen,mols_sur=mols_sur,cenpos=cenpos,length=length,TVs=TV,maxTVs=maxTVs,mols_outer=mols_outer,outer_maxTVs=outer_maxTVs)

#prot_neut_cry._mols contains all molecules.
#mols[0] contains a list of all molecules in position a, mols[1] all mols in pos'n b, etc.
#mols[0][x,y,z] contains molecule a in position x,y,z
#mols may as such be iterated over in a number of ways to consider different molecules.

print 'state',state
print 'q: ', qdict[state]

for atom in prot_neut_cry()._mols[0][prot_neut_cry()._cenpos[0]][prot_neut_cry()._cenpos[1]][prot_neut_cry()._cenpos[2]]():
	atom()._crg=qdict[state]


prot_neut_cry().print_posns()


#Calculate Properties:

print strftime("%a, %d %b %Y %X +0000", gmtime())
E0 = np.matrix([0.,0.,0.])

print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'Calc jm'
screenradius=1.6623/(Natoms)
# Thole paper screenradius value for fit to components of pol. tensor divided by no. atoms in mol. We choose this screenradius value for smearing of charge as, with near planar mols, in some dirs we have molecule-like polarisabilities with near atom-like separations. 
#This form of screenradius will result in charge being smeared along the separation axis of molecules by NAtoms*(Thole's value for a single atom)
jm = JMatrix(jmtype='TholeExp',screenradius=screenradius)

print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'Calc dips:'
d = get_dipoles(E0=E0,jm=jm._m,cutoff=cut)
print strftime("%a, %d %b %Y %X +0000", gmtime())
Efield = get_electric_field(E0)
potential = get_potential()

print strftime("%a, %d %b %Y %X +0000", gmtime())
#print 'dips', d
print 'splitting dips onto atoms'
split_d = split_dipoles_onto_atoms(d)
print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'summing dips:'
tot = np.matrix([0.,0.,0.])
for dd in split_d:
	tot += dd

print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'total dip moment', tot
Uqq = np.multiply(get_U_qq(potential=potential),27.211)
print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'Uqq', Uqq
Uqd = np.multiply(get_U_qdip(dips=d,Efield=Efield),27.211)
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

# print TVs
if not os.path.exists('Dips_Posns_TVs'): os.makedirs('Dips_Posns_TVs')
f = open('Dips_Posns_TVs/%s_TVs.dat' % name, 'w')
TVstr=str(str(TV[0,0]) + ' ' + str(TV[0,1]) + ' ' +  str(TV[0,2]) +  '\n' +  str(TV[1,0]) + ' ' +  str(TV[1,1]) + ' ' +  str(TV[1,2]) +  '\n' +  str(TV[2,0]) + ' ' +  str(TV[2,1]) + ' ' +  str(TV[2,2])+  '\n')
f.write(TVstr)
f.flush()
f.close()


# print dipoles
if not os.path.exists('Dips_Posns_TVs'): os.makedirs('Dips_Posns_TVs')
f = open('Dips_Posns_TVs/%s_dipoles.dat' % name, 'w')
for dd in split_d:
	dstr=str(dd)
	f.write(dstr)
	f.write('\n')
f.flush()
f.close()

# print properties for charge in centrepos
time=strftime("%a, %d %b %Y %X +0000", gmtime())
f = open('%s_properties.csv' % name, 'w')
f.write ('time\tname\tmols_cen\tmols_sur\tmols_outer\tinsize\toutsize\tenergyev\tUqq\tUqd\tUdd\tTotdip_x\tTotdip_y\tTotdip_z')
f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,mols_cen,mols_sur,mols_outer,insize,outsize,energyev,Uqq,Uqd,Udd,tot[0,0],tot[0,1],tot[0,2]))
f.flush()
f.close()
# print header for reorgs
f = open('reorg_energies_%s_properties.csv' % name, 'w')
f.write ('time\tname\tmols_cen\tmols_sur\tmols_outer\tinsize\toutsize\ta\tb\tc\tmolincell\tReorg(eV)')
f.flush()
f.close()

# REORGANISATION ENERGIES
#Note that this assumes a cube, and values for which 

for dist in range(0,(length[0]/2)+1,1):
	print '\n\nDIST: ', dist, '\n'
	for a in range(prot_neut_cry()._cenpos[0]-dist,prot_neut_cry()._cenpos[0]+dist+1,1):
			for b in range(prot_neut_cry()._cenpos[1]-dist,prot_neut_cry()._cenpos[1]+dist+1,1):
				for c in range(prot_neut_cry()._cenpos[2]-dist,prot_neut_cry()._cenpos[2]+dist+1,1):
					print strftime("%a, %d %b %Y %X +0000", gmtime())
					print 'a,b,c',a,b,c
					for molincell in range(0,len(prot_neut_cry()._mols),1):
						prot_neut_cry().calc_reorg_shareq(a1=prot_neut_cry()._cenpos[0],b1=prot_neut_cry()._cenpos[1],c1=prot_neut_cry()._cenpos[2],molincell1=0,a2=a,b2=b,c2=c,molincell2=molincell,jm=jm._m,oldUqd=Uqd)								
						print 'Reorg: ', prot_neut_cry()._reorgs_shareq[molincell][a][b][c]

						f = open('reorg_energies_%s_properties.csv' % name, 'a')
						f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,mols_cen,mols_sur,mols_outer,insize,outsize,a,b,c,molincell,prot_neut_cry()._reorgs_shareq[molincell][a][b][c]))
						f.flush()
						f.close()
						# Redo this and overwrite after each set to ensure we have some even if not all reorgs complete
						prot_neut_cry().print_reorgs_shareq()

print 'Job Completed Successfully.'
