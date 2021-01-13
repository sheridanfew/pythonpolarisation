import sys
sys.path.append('../')
from math import *
from BasicElements.Crystal import *

name='prot_neut'

mols_cen=['protontest.xyz','TEST_neutron_pos100.xyz']

mols_sur=['TEST_neutron_pos000.xyz','TEST_neutron_pos100.xyz']

#Get translation vectors:

a=1.0/0.5291772109217
b=1.0/0.5291772109217
c=1.0/0.5291772109217

alpha=90*(pi/180)
beta=90*(pi/180)
gamma=90*(pi/180)

cell_volume=sqrt(1 - (cos(alpha)**2) - (cos(beta)**2) - (cos(gamma)**2) + (2*cos(alpha)*cos(beta)*cos(gamma)))

matrix_to_cartesian=np.matrix( [[a, b*cos(gamma), c*cos(beta)],
[0, b*sin(gamma), c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)],
[0, 0, c*cell_volume/sin(gamma)]])

TV=matrix_to_cartesian.T

prot_neut_cry=Crystal(name=name,mols_cen=centres,mols_sur=surroundings,cenpos=[1,1,1],lena=2,lenb=2,lenc=2,TVs=TV)

