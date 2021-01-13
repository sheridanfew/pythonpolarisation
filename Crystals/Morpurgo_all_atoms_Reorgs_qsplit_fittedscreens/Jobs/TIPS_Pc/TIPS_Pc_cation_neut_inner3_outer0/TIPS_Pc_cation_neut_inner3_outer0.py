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


name='TIPS_Pc_cation_neut_inner3_outer0'
#For crystals here, all cubic and centred at centre
insize=3
#number of TVs in each dir central mol is from edge of inner region
outsize=0
state='cation'
mols_cen=['TIPS_Pc_mola_cation_aniso_cifstruct_chelpg.xyz']
mols_sur=['TIPS_Pc_mola_neut_aniso_cifstruct_chelpg.xyz']
mols_outer=['TIPS_Pc_mola_aniso_cifstruct_chelpg.xyz']

screenradius=2.4704022921

#From cif:
'''
TIPS

data_k01029 

_cell_length_a                    7.5650(15) 
_cell_length_b                    7.7500(15) 
_cell_length_c                    16.835(3) 
_cell_angle_alpha                 89.15(3) 
_cell_angle_beta                  78.42(3) 
_cell_angle_gamma                 83.63(3) 

_cell_volume                      960.9(3) 

'''
#Get translation vectors:

a=7.565015/0.5291772109217
b=7.750015/0.5291772109217
c=16.8353/0.5291772109217

alpha=89.153*(pi/180)
beta=78.423*(pi/180)
gamma=83.633*(pi/180)

cif_unit_cell_volume=960.9/(a*b*c*(0.5291772109217**3))
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

crystal=Crystal(name=name,mols_cen=mols_cen,mols_sur=mols_sur,cenpos=cenpos,length=length,TVs=TV,maxTVs=maxTVs,mols_outer=mols_outer,outer_maxTVs=outer_maxTVs)

crystal().ModifyPolarizabilityCry(jmtype='TholeExp',fittype='empirical')

#crystal._mols contains all molecules.
#mols[0] contains a list of all molecules in position a, mols[1] all mols in pos'n b, etc.
#mols[0][x,y,z] contains molecule a in position x,y,z
#mols may as such be iterated over in a number of ways to consider different molecules.

print 'state',state
#print 'q: ', qdict[state]

#for atom in crystal()._mols[0][crystal()._cenpos[0]][crystal()._cenpos[1]][crystal()._cenpos[2]]():
#	atom()._crg=qdict[state]


crystal().print_posns()


#Calculate Properties:

print strftime("%a, %d %b %Y %X +0000", gmtime())
E0 = np.matrix([0.,0.,0.])

print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'Calc jm'
#screenradius=1.6623/(Natoms**2)
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
	for a in range(crystal()._cenpos[0]-dist,crystal()._cenpos[0]+dist+1,1):
			for b in range(crystal()._cenpos[1]-dist,crystal()._cenpos[1]+dist+1,1):
				for c in range(crystal()._cenpos[2]-dist,crystal()._cenpos[2]+dist+1,1):
					print strftime("%a, %d %b %Y %X +0000", gmtime())
					print 'a,b,c',a,b,c
					for molincell in range(0,len(crystal()._mols),1):
						crystal().calc_reorg_shareq(a1=crystal()._cenpos[0],b1=crystal()._cenpos[1],c1=crystal()._cenpos[2],molincell1=0,a2=a,b2=b,c2=c,molincell2=molincell,jm=jm._m,oldUqd=Uqd)								
						print 'Reorg: ', crystal()._reorgs_shareq[molincell][a][b][c]

						f = open('reorg_energies_%s_properties.csv' % name, 'a')
						f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,mols_cen,mols_sur,mols_outer,insize,outsize,a,b,c,molincell,crystal()._reorgs_shareq[molincell][a][b][c]))
						f.flush()
						f.close()
						# Redo this and overwrite after each set to ensure we have some even if not all reorgs complete
						crystal().print_reorgs_shareq()

print 'Job Completed Successfully.'
