# 4TCONSTRAINED:PCBM_EDITEDO_anion_mul polarisation model. from files:
#4T_CONSTRAINED_cation_mul.xyz
#PCBM_EDITEDO_anion_mul.xyz
#Produced with makeblends_pol_model.sh

import sys
sys.path.append('../../../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import *
from Polarizability.GetEnergy import *
import numpy as np
from math import *

ReadMoleculeType('../../../Molecules/./4T_CONSTRAINED_cation_mul.xyz')
don = GetMolecule('../../../Molecules/./4T_CONSTRAINED_cation_mul.xyz')
donatomlist=[]
for atom in don():
	donatomlist.append(atom())

distrangeA=np.arange(2.0, 8.0, 0.5)
rotrange=np.arange(0.0, 330.1, 30.0)

cut=1.1

f = open('../../../Blends/PT_PCBM_EDITEDO_anion_mul/4T_PCBM_2013_05_29/4T_PCBM_2013_05_29.csv', 'w')
f.write ("System\tAlignment\tAngle\tDistance\tEnergy\tUqq\tUqd\tUdd\tDipole Moment")
f.flush()
f.close()

regionnames=['2th_ring']
atomnos=np.array([[8, 9, 10, 11, 12, 13, 14]])

regioniteration=np.arange(0,len(atomnos),1)

#iterate over regions
for i in regioniteration:
#find centres
	sumpos=np.matrix([[0., 0., 0.]])
	for j in atomnos[i]:
		sumpos=sumpos + donatomlist[j]._pos
#??? Check syntax above
	centre=np.zeros(shape=(2,3))
	centre[i]=sumpos/len(atomnos[i])

#find perpendicular versor to face of region
	v1pos=donatomlist[atomnos[i][0]]._pos-donatomlist[atomnos[i][1]]._pos
	v2pos=donatomlist[atomnos[i][0]]._pos-donatomlist[atomnos[i][2]]._pos
#clunky, otherwise v1,2 class is pos'n, and python has trouble w. cross. want matrix.

	v1=np.matrix([[v1pos[0,0],v1pos[0,1],v1pos[0,2]]])
	v2=np.matrix([[v2pos[0,0],v2pos[0,1],v2pos[0,2]]])
	v1=v1
	v2=v2

	perp=np.zeros(shape=(2,3))
	perp[i]=np.cross(v1,v2)/np.sqrt(np.dot(np.cross(v1,v2),np.cross(v1,v2).T))
	print 'perp'
	print perp

	if perp[i][2] < 0:
		perp=-perp
#		if debug == 'TRUE':
#			print 'negative z perp reversed'

#get rotation matrix to this versor from z (may not be most efficient way. Also, psi taken as 0, affects rotational orientation of fullerene face aliogned w. polymer, don't think important)

	if perp[i][2] == 1:
		phi=0.
	else:
		phi=np.arctan(perp[i][1]/perp[i][0])
	psi=0

	for rotangle in rotrange:
		print 'type perp, rotangle'
		tp=type(perp[i][0])
		tr=type(rotangle)
		print tp, tr

