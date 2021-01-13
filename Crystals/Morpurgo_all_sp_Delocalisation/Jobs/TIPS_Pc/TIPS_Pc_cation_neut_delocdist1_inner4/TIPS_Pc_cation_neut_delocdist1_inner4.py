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

name='TIPS_Pc_cation_neut_delocdist1_inner4'
state='cation'
# Max no of TVs from centre to delocalise the charge
delocdist=1
#For crystals here, all diamond shape and centred at centre
insize=4
outsize=0
#number of TVs in each dir central mol is from edge of inner region
mols_cen=['sp_TIPS_Pc_neut.xyz']
mols_sur=['sp_TIPS_Pc_neut.xyz']
mols_outer=['sp_TIPS_Pc_neut.xyz']

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
charge={"cation": 1.0, "neut": 0.0, "anion": -1.0}
# Number of unit cells charge is delocalised over as a function of delocdist (distance delocalising from centre)
ndeloc=[1.0,7.0,25.0]

# Number of sites over which charge is delocalised (could probably calculate this from delocdist, or make array of values)

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

prot_neut_cry().print_posns()

# DELOCALISE CHARGE

print 'prot_neut_cry()._cenpos[0]:', prot_neut_cry()._cenpos[0]
print 'len(prot_neut_cry()._mols[0]):', len(prot_neut_cry()._mols[0])


print '\n\ndelocdist: ', delocdist, '\n'
for a in range(prot_neut_cry()._cenpos[0]-delocdist,prot_neut_cry()._cenpos[0]+delocdist+1,1):
		for b in range(prot_neut_cry()._cenpos[1]-delocdist,prot_neut_cry()._cenpos[1]+delocdist+1,1):
			for c in range(prot_neut_cry()._cenpos[2]-delocdist,prot_neut_cry()._cenpos[2]+delocdist+1,1):
				print 'a,b,c',a,b,c
				if (abs(a-prot_neut_cry()._cenpos[0]) + abs(b-prot_neut_cry()._cenpos[1]) + abs(c-prot_neut_cry()._cenpos[2]))  <= delocdist:
					for molincell in range(0,len(prot_neut_cry()._mols),1):
						'''
						print 'prot_neut_cry'
						print prot_neut_cry()
						print prot_neut_cry()._mols
						print 'type: prot_neut_cry()._mols[molincell][a][b][c]:'
						print type(prot_neut_cry()._mols[molincell][a][b][c])
						print prot_neut_cry()._mols[molincell][a][b][c]
						'''
						for atom in prot_neut_cry()._mols[molincell][a][b][c]():
							print atom()._crg
							atom()._crg = charge[state]/(ndeloc[delocdist]*len(prot_neut_cry()._mols))
							#atom()._pol = charge[state]/(ndeloc[delocdist]*len(prot_neut_cry()._mols))
totq=0.0

print '\n\nCharges on each deloc molecule:\n'
print '\n\ndelocdist: ', delocdist, '\n'
for a in range(prot_neut_cry()._cenpos[0]-delocdist,prot_neut_cry()._cenpos[0]+delocdist+1,1):
		for b in range(prot_neut_cry()._cenpos[1]-delocdist,prot_neut_cry()._cenpos[1]+delocdist+1,1):
			for c in range(prot_neut_cry()._cenpos[2]-delocdist,prot_neut_cry()._cenpos[2]+delocdist+1,1):
				if (abs(a-prot_neut_cry()._cenpos[0]) + abs(b-prot_neut_cry()._cenpos[1]) + abs(c-prot_neut_cry()._cenpos[2]))  <= delocdist:
					for molincell in range(0,len(prot_neut_cry()._mols),1):
						for atom in prot_neut_cry()._mols[molincell][a][b][c]():
							print atom()._crg
							totq=totq+atom()._crg

print 'TOT q: ', totq
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

# print properties
time=strftime("%a, %d %b %Y %X +0000", gmtime())
f = open('%s_properties.csv' % name, 'w')
f.write ('time\tname\tmols_cen\tmols_sur\tmols_outer\tinsize\toutsize\tenergyev\tUqq\tUqd\tUdd\tTotdip_x\tTotdip_y\tTotdip_z')
f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,mols_cen,mols_sur,mols_outer,insize,outsize,energyev,Uqq,Uqd,Udd,tot[0,0],tot[0,1],tot[0,2]))
f.flush()
f.close()

print 'Job Completed Successfully.'
