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

###################################
#START OF MOLECULE SPECIFIC SECTION
###################################
name='prot_neut'

mols_cen=['protontest.xyz','TEST_neutron_pos100.xyz']

mols_sur=['TEST_neutron_pos000.xyz','TEST_neutron_pos100.xyz']

mols_outer=['TEST_neutron_pos000.xyz','TEST_neutron_pos100.xyz']

#mols_cen=['protontest.xyz']

#mols_sur=['TEST_neutron_pos000.xyz']


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

a=1.0/0.5291772109217
b=1.0/0.5291772109217
c=1.0/0.5291772109217

alpha=90*(pi/180)
beta=90*(pi/180)
gamma=90*(pi/180)

#cif_unit_cell_volume=2116.01/(a*b*c*(0.5291772109217**3))

###################################
#END OF MOLECULE SPECIFIC SECTION
###################################

cell_volume=sqrt(1 - (cos(alpha)**2) - (cos(beta)**2) - (cos(gamma)**2) + (2*cos(alpha)*cos(beta)*cos(gamma)))

#print 'calc cell_volume, cif cell_volume:',cell_volume, ', ', cif_unit_cell_volume

#if cell_volume > (cif_unit_cell_volume*1.1) or cell_volume < (cif_unit_cell_volume*0.9):
#	print 'angles wrong?'
#	exit (2)

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

#For crystals here, all cubic and centred at centre
'''
cen=2
out=cen+1
cenpos=[out,out,out]
#cenpos=[cen,cen,cen]
#length=[2*cen+1,2*cen+1,2*cen+1]
length=[2*out+1,2*out+1,2*out+1]
'''

insize=1
#number of TVs in each dir central mol is from edge of inner region
outsize=0
out=insize+outsize
#number of TVs in each dir nearest c inner mol is from edge of outer region
cenpos=[insize+outsize,insize+outsize,insize+outsize]
length=[2*out+1,2*out+1,2*out+1]
maxTVs=insize
outer_maxTVs=insize+outsize
#for diamond outer, don't specify for cube and will fill to cube edges.


'''
cenpos=[2,2,2]
length[0]=5
length[1]=5
length[2]=5
'''

print 'name: ',name,'mols_cen: ', mols_cen,' mols_sur: ',mols_sur,' TVs: ', TV 

#prot_neut_cry=Crystal(name=name,mols_cen=mols_cen,mols_sur=mols_sur,TVs=TV,cenpos=[1,1,1],length[0]=2,length[1]=2,length[2]=2)
prot_neut_cry=Crystal(name=name,mols_cen=mols_cen,mols_sur=mols_sur,cenpos=cenpos,length=length,TVs=TV,maxTVs=maxTVs,mols_outer=mols_outer,outer_maxTVs=outer_maxTVs)

#prot_neut_cry._mols contains all molecules.
#mols[0] contains a list of all molecules in position a, mols[1] all mols in pos'n b, etc.
#mols[0][x,y,z] contains molecule a in position x,y,z
#mols may as such be iterated over in a number of ways to consider different molecules.

prot_neut_cry().print_posns()


#Calculate Properties:

print strftime("%a, %d %b %Y %X +0000", gmtime())
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
#print 'dips', d
split_d = split_dipoles_onto_atoms(d)
tot = np.matrix([0.,0.,0.])
for dd in split_d:
	tot += dd

print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'tot', tot
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
f = open('Dips_Posns_TVs/%s_size_%s_TVs.dat' % (name,(length[0]+1)*(length[1]+1)*(length[2]+1)), 'w')
TVstr=str(str(TV[0,0]) + ' ' + str(TV[0,1]) + ' ' +  str(TV[0,2]) +  '\n' +  str(TV[1,0]) + ' ' +  str(TV[1,1]) + ' ' +  str(TV[1,2]) +  '\n' +  str(TV[2,0]) + ' ' +  str(TV[2,1]) + ' ' +  str(TV[2,2])+  '\n')
f.write(TVstr)
f.flush()
f.close()


# print dipoles
if not os.path.exists('Dips_Posns_TVs'): os.makedirs('Dips_Posns_TVs')
f = open('Dips_Posns_TVs/%s_size_%s_dipoles.dat' % (name,(length[0]+1)*(length[1]+1)*(length[2]+1)), 'w')
for dd in split_d:
	dstr=str(dd)
	f.write(dstr)
	f.write('\n')
f.flush()
f.close()

'''
print 'prot_neut_cry() central qs'

for atom in prot_neut_cry()._mols[0][prot_neut_cry()._cenpos[0]][prot_neut_cry()._cenpos[1]][prot_neut_cry()._cenpos[2]]():
	print atom()._crg

for molincell in range(0,len(prot_neut_cry()._mols),1):
	prot_neut_cry().calc_reorg(a1=prot_neut_cry()._cenpos[0],b1=prot_neut_cry()._cenpos[1],c1=prot_neut_cry()._cenpos[2],molincell1=0,a2=prot_neut_cry()._cenpos[0],b2=prot_neut_cry()._cenpos[1],c2=prot_neut_cry()._cenpos[2],molincell2=molincell,dips=d,oldUqd=Uqd)
	print 'Reorg: ', prot_neut_cry()._reorgs[molincell][prot_neut_cry()._cenpos[0]][prot_neut_cry()._cenpos[1]][prot_neut_cry()._cenpos[2]]

#Note that this assumes a cube, and values for which 

	print 'prot_neut_cry()._reorgs[0][0][0][0]: ', prot_neut_cry()._reorgs[0][0][0][0]

'''

mols_cen=['protontest.xyz','TEST_neutron_pos100.xyz']

mols_sur=['TEST_neutron_pos000.xyz','TEST_neutron_pos100.xyz']

mols_outer=['TEST_neutron_pos000.xyz','TEST_neutron_pos100.xyz']

time=strftime("%a, %d %b %Y %X +0000", gmtime())
f = open('%s_insize_%s_outsize_%s_properties.csv' % (name,insize,outsize), 'w')
f.write ('time\tname\tmols_cen\tmols_sur\tmols_outer\tinsize\toutsize\tenergyev\tUqq\tUqd\tUdd\tTotdip_x\tTotdip_y\tTotdip_z')
f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,mols_cen,mols_sur,mols_outer,insize,outsize,energyev,Uqq,Uqd,Udd,tot[0,0],tot[0,1],tot[0,2]))
f.flush()
f.close()

f = open('reorg_energies_%s_insize_%s_outsize_%s_properties.csv' % (name,insize,outsize), 'w')
f.write ('time\tname\tmols_cen\tmols_sur\tmols_outer\tinsize\toutsize\ta\tb\tc\tmolincell\tReorg(eV)')
f.flush()
f.close()


#Note that this assumes a cube, and values for which 

for dist in range(0,(length[0]/2)+1,1):
	print '\n\nDIST: ', dist, '\n'
	for a in range(prot_neut_cry()._cenpos[0]-dist,prot_neut_cry()._cenpos[0]+dist+1,1):
			for b in range(prot_neut_cry()._cenpos[1]-dist,prot_neut_cry()._cenpos[1]+dist+1,1):
				for c in range(prot_neut_cry()._cenpos[2]-dist,prot_neut_cry()._cenpos[2]+dist+1,1):
					print 'a,b,c',a,b,c
					for molincell in range(0,len(prot_neut_cry()._mols),1):
#						if prot_neut_cry()._reorgs[molincell][a][b][c] == []:
						prot_neut_cry().calc_reorg(a1=prot_neut_cry()._cenpos[0],b1=prot_neut_cry()._cenpos[1],c1=prot_neut_cry()._cenpos[2],molincell1=0,a2=a,b2=b,c2=c,molincell2=molincell,dips=d,oldUqd=Uqd)								
						#prot_neut_cry().calc_reorg(a1=prot_neut_cry()._cenpos[0],b1=prot_neut_cry()._cenpos[1],c1=prot_neut_cry()._cenpos[2],molincell1=0,a2=prot_neut_cry()._cenpos[0]+(a*posneg),b2=prot_neut_cry()._cenpos[1]+(b*posneg),c2=prot_neut_cry()._cenpos[2]+(c*posneg),molincell2=molincell,dips=d,oldUqd=Uqd)	
						print 'Reorg: ', prot_neut_cry()._reorgs[molincell][a][b][c]

						f = open('reorg_energies_%s_insize_%s_outsize_%s_properties.csv' % (name,insize,outsize), 'a')
						f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,mols_cen,mols_sur,mols_outer,insize,outsize,a,b,c,molincell,prot_neut_cry()._reorgs[molincell][a][b][c]))
						f.flush()
						f.close()

#						else:
#							print 'skipping, already done for previous size.'
	prot_neut_cry().print_reorgs()

print 'Done.'
