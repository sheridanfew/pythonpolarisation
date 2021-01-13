# 8Tcatmul.xyz:PCBM_EDITEDO_anion_mul polarisation model. from files:
#8Tcat_mul.xyz
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

ReadMoleculeType('../../../Molecules/./8Tcat_mul.xyz')
don = GetMolecule('../../../Molecules/./8Tcat_mul.xyz')
donatomlist=[]
for atom in don():
	donatomlist.append(atom())

distrangeA=np.arange(2.0, 8.0, 0.5)
rotrange=np.arange(0.0, 330.1, 30.0)

f = open('../../../Blends/PTcat_PCBM_EDITEDO_anion_mul/8T_PCBMtest_2013_03_13/8T_PCBMtest_2013_03_13.csv', 'w')
f.write ("System\tAlignment\tAngle\tDistance\tEnergy\tUqq\tUqd\tUdd\t")
f.flush()
f.close()

regionnames=['1th_ring', '2th_ring', '3th_ring', '4th_ring', '5th_ring', '6th_ring', '7th_ring', '8th_ring']
atomnos=np.array([[1, 2, 3, 4, 5, 6, 7, 8],
[9, 10, 11, 12, 13, 14, 15],
[16, 17, 18, 19, 20, 21, 22],
[23, 24, 25, 26, 27, 28, 29],
[30, 31, 32, 33, 34, 35, 36],
[37, 38, 39, 40, 41, 42, 43],
[44, 45, 46, 47, 48, 49, 50],
[51, 52, 53, 54, 55, 56, 57, 58]])

regioniteration=np.arange(0,len(atomnos),1)

#iterate over regions
for i in regioniteration:
#find centres
	sumpos=np.matrix([[0., 0., 0.]])
	for j in atomnos[i]:
		sumpos=sumpos + donatomlist[j]._pos
#??? Check syntax above
	centre=np.zeros(shape=(9,3))
	centre[i]=sumpos/len(atomnos[i])

#find perpendicular versor to face of region
	v1pos=donatomlist[atomnos[i][0]]._pos-donatomlist[atomnos[i][1]]._pos
	v2pos=donatomlist[atomnos[i][0]]._pos-donatomlist[atomnos[i][2]]._pos
#clunky, otherwise v1,2 class is pos'n, and python has trouble w. cross. want matrix.

	v1=np.matrix([[v1pos[0,0],v1pos[0,1],v1pos[0,2]]])
	v2=np.matrix([[v2pos[0,0],v2pos[0,1],v2pos[0,2]]])
	v1=v1
	v2=v2


	t=type(v1)
	print t
	print v1, v2
	print np.cross(v1,v2)
	print np.dot(np.cross(v1,v2).T,np.cross(v1,v2))
	perp=np.zeros(shape=(9,3))
	perp[i]=np.cross(v1,v2)/np.sqrt(np.dot(np.cross(v1,v2),np.cross(v1,v2).T))
	print 'perp'
	print perp

	if perp[i][2] < 0:
		perp=-perp
		if debug == 'TRUE':
			print 'negative z perp reversed'

#get rotation matrix to this versor from z (may not be most efficient way. Also, psi taken as 0, affects rotational orientation of fullerene face aliogned w. polymer, don't think important)

	if perp[i][2] == 1:
		phi=0.
	else:
		phi=np.arctan(perp[i][1]/perp[i][0])
	psi=0

	for rotangle in rotrange:
		if perp[i][2] == 1:
			theta=0.
		else:
			theta=np.arctan((((perp[i][0]^2)+(perp[i][1]^2))/(perp[i][2]^2))^0.5)+(rotangle*pi/180)

#From JK Rotation 
		D = np.matrix([[cos(phi), sin(phi),0.],[-sin(phi), cos(phi), 0.],[0.,0.,1.]])
		B = np.matrix([[cos(psi), sin(psi),0.],[-sin(psi), cos(psi), 0.],[0.,0.,1.]])
		C = np.matrix([[1.,0.,0.],[0., cos(theta), sin(theta)], [0, -sin(theta),cos(theta)]])
		ROTMAT=(B*C*D)

		print ROTMAT

		dispvec=np.zeros(shape=(9,3))
		dispvec[i]=ROTMAT*np.matrix([[0., 0., 1. ]]).T

		ReadMoleculeType('../../../Molecules/./PCBM_EDITEDO_anion_mul.xyz')
		acc = GetMolecule('../../../Molecules/./PCBM_EDITEDO_anion_mul.xyz')

#Crude rotation, couldn't figure out Rotation module in basic elements.

		for atom in acc():
			atom().move(ROTMAT*atom()._pos-atom()._pos)

#Place acc. at centre of unit

		for atom in acc():
			atom().move(centre[i])

#Place at correct Distance, places too close then pushes away until correct dist.

		for distA in distrangeA:
			distbohr=(np.divide(distA,0.529177))
			disttest=distbohr - 0.3

			for atom in acc():
				atom().move(disttest*dispvec[i])

			#closest approach must be within 0.05 Bohr of desired distance
			thresholddist = distbohr - 0.05
			thresholdmet = 'FALSE'

			while thresholdmet == 'FALSE':
				#move acceptor by 0.1 bohr
				for atom in acc():
					atom().move(0.1*dispvec[i])

				#get list of interatomic distances
				interatomic_distances=[]
				for atoma in acc():
					for atomb in don():
						distvec = atoma()._pos - atomb()._pos
						distabs = m.sqrt(np.dot(distvec, distvec.T)[0,0])
						interatomic_distances.append(distabs)

				#check if smallest distance is above threshold dist
				if min(distances) > thresholddist:
					thresholdmet='TRUE'

			E0 = np.matrix([0.,0.,0.])
			d = get_dipoles(E0=E0,cutoff=cut)
			split_d = split_dipoles_onto_atoms(d)
			tot = np.matrix([0.,0.,0.])
			for dd in split_d:
				tot += dd
			Uqq = np.multiply(get_U_qq(),27.211)
			Uqd = np.multiply(get_U_qdip(),27.211)
			Udd = np.multiply(get_U_dipdip(),27.211)
			energy = get_energy(E0=E0,cutoff=cut)
			energyev = np.multiply(energy,27.211)
#			jm = get_jmat(E0=E0,cutoff=cut)

#csv file

			f.open('../../../Blends/PTcat_PCBM_EDITEDO_anion_mul/8T_PCBMtest_2013_03_13/8T_PCBMtest_2013_03_13.csv', 'a')
			f.write ('8Tcatmul.xyz:PCBM_EDITEDO_anion_mul\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (regionname[i],rotangle,distA,energyev,Uqq,Uqd,Udd))
			f.flush()
			f.close()

#.dat for gnuplot

			f = open('8Tcatmul.xyz_PCBM_EDITEDO_anion_mul_%s_%gdeg_%gA.dat' % (regionname[i],rotangle,distA), 'w')
			f.write('Molecule Starts:\n')
			for atom in Register.GetRegister('Atom'):
				astr=str(atom)
				f.write(astr)
				f.write('\n')
			f.write('Molecule Ends.\n\n')
			f.write ('dist, %s \n distbohr, %s \n tot, %s \n energy, %s \n energyev, %s \n' % (dist,distbohr,tot,energy,energyev))
			f.write('\n\ndipoles start:\n')
			for dd in split_d:
				dstr=str(dd)
				f.write(dstr)
				f.write('\n')
			f.flush()
			f.close()
exit()
