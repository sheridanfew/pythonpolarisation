# 33TIcation:PC71BM_anion_B3LYP_6_31gSTAR_opt_mul polarisation model. from files:
#33TI_cation_B3LYP_6_31gSTAR_opt_mul.xyz
#PC71BM_anion_B3LYP_6_31gSTAR_opt_mul.xyz
#Produced with makeblends_pol_model_coulonly.sh

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

ReadMoleculeType('../../../Molecules/./33TI_cation_B3LYP_6_31gSTAR_opt_mul.xyz')
don = GetMolecule('../../../Molecules/./33TI_cation_B3LYP_6_31gSTAR_opt_mul.xyz')
donatomlist=[]
for atom in don():
	donatomlist.append(atom())

ReadMoleculeType('../../../Molecules/./PC71BM_anion_B3LYP_6_31gSTAR_opt_mul.xyz')
acc = GetMolecule('../../../Molecules/./PC71BM_anion_B3LYP_6_31gSTAR_opt_mul.xyz')


distrangeA=np.arange(3.5, 3.8, 0.5)
rotrange=np.arange(0.0, 15.0, 30.0)

cut=5.0

f = open('../../../Blends/P3TI_PC71BM_anion_B3LYP_6_31gSTAR_opt_mul/33TI_PC71BM_coulomb_2013_05_30/33TI_PC71BM_coulomb_2013_05_30_2.csv', 'w')
f.write ("System\tAlignment\tAngle\tDistance\tUqq")
f.flush()
f.close()

regionnames=['1I', '2T1', '2T2', '2T3', '2I']
atomnos=np.array([[22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49],
[50, 51, 52, 53, 54, 55, 56],
[57, 58, 59, 60, 61, 62, 63],
[64, 65, 66, 67, 68, 69, 70],
[71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98]])
regioniteration=np.arange(((len(atomnos) - 4 )/2),(((len(atomnos)+4) )/2),1)

#iterate over regions
for i in regioniteration:
#find centres
        print 'doing unit %s, %s' % (i , regionnames[i])
	sumpos=np.matrix([[0., 0., 0.]])
	for j in atomnos[i]:
		sumpos=sumpos + donatomlist[j]._pos
#??? Check syntax above
	centre=np.zeros(shape=(6,3))
	centre[i]=sumpos/len(atomnos[i])

#find perpendicular versor to face of region
	v1pos=donatomlist[atomnos[i][0]]._pos-donatomlist[atomnos[i][1]]._pos
	v2pos=donatomlist[atomnos[i][0]]._pos-donatomlist[atomnos[i][2]]._pos
#clunky, otherwise v1,2 class is pos'n, and python has trouble w. cross. want matrix.

	v1=np.matrix([[v1pos[0,0],v1pos[0,1],v1pos[0,2]]])
	v2=np.matrix([[v2pos[0,0],v2pos[0,1],v2pos[0,2]]])
	v1=v1
	v2=v2

	perp=np.zeros(shape=(6,3))
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
		if perp[i][2] == 1:
			theta=0. +(rotangle*pi/180)
		else:
			theta=np.arctan((((perp[i][0]**2)+(perp[i][1]**2))/(perp[i][2]**2))**0.5)+(rotangle*pi/180)

#From JK Rotation 
		D = Rotation([[cos(phi), sin(phi),0.],[-sin(phi), cos(phi), 0.],[0.,0.,1.]])
		B = Rotation([[cos(psi), sin(psi),0.],[-sin(psi), cos(psi), 0.],[0.,0.,1.]])
		C = Rotation([[1.,0.,0.],[0., cos(theta), sin(theta)], [0, -sin(theta),cos(theta)]])
		ROTMAT=(B*C*D)

		print B,C,D
		print 'ROTMAT'
		print ROTMAT

		ROTMATINV=np.transpose(ROTMAT)
		print ROTMATINV

		dispversor=np.zeros(shape=(6,3))
		dispversor[i]=(ROTMAT*np.matrix([[0., 0., 1. ]]).T).T

#Crude rotation, couldn't figure out Rotation module in basic elements.
		print 'rotating acc:'
		for atom in acc():
			atom().move((ROTMAT*atom()._pos.T).T-atom()._pos)

#Place acc. at centre of unit

                print 'moving acc to centre of unit:'
		for atom in acc():
			atom().move(centre[i])

#Place at correct Distance, places too close then pushes away until correct dist.

                print 'Placing acc. at correct distance:'
		for distA in distrangeA:
			distbohr=(np.divide(distA,0.529177))
                        print 'going for %s A, %s Bohr' % (distA,distbohr)
#			disttest=distbohr - 0.3
#
#			for atom in acc():
 #                               atom().move(disttest*dispversor[i])

			#closest approach must be within 0.05 Bohr of desired distance
			thresholddist = distbohr - 0.05
			thresholdmet = 'FALSE'
			for atom in acc():
				atom().move((distbohr-0.3)*dispversor[i])

			distmoved=distbohr-0.3
			while thresholdmet == 'FALSE':
				#move acceptor by 0.1 bohr
				for atom in acc():
					atom().move(0.1*dispversor[i])
				
				distmoved=distmoved+0.1

				#get list of interatomic distances
				interatomic_distances=[]
				for atoma in acc():
					for atomb in don():
						distvec = atoma()._pos - atomb()._pos
						distabs = sqrt(np.dot(distvec, distvec.T)[0,0])
						interatomic_distances.append(distabs)
				#check if smallest distance is above threshold dist
				print 'min atom dist = %s Bohr' % min(interatomic_distances)
						                                                
				if min(interatomic_distances) > thresholddist:
					thresholdmet='TRUE'
                                        print 'Threshold met.'

			E0 = np.matrix([0.,0.,0.])
			d = get_dipoles(E0=E0,cutoff=cut)
			split_d = split_dipoles_onto_atoms(d)
			tot = np.matrix([0.,0.,0.])
			for dd in split_d:
				tot += dd
                        print 'Getting Uqq'
			Uqq = np.multiply(get_U_qq(),27.211)
			Uqd = np.multiply(get_U_qdip(),27.211)
			Udd = np.multiply(get_U_dipdip(),27.211)
#			energy = get_energy(E0=E0,cutoff=cut)
			energyev = Uqq + Uqd + Udd
#			jm = get_jmat(E0=E0,cutoff=cut)

#csv file
			f = open('../../../Blends/P3TI_PC71BM_anion_B3LYP_6_31gSTAR_opt_mul/33TI_PC71BM_coulomb_2013_05_30/33TI_PC71BM_coulomb_2013_05_30_2.csv', 'a')
			f.write ('\n33TIcation:PC71BM_anion_B3LYP_6_31gSTAR_opt_mul\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (regionnames[i],rotangle,distA,energyev,Uqq,Uqd,Udd,tot))
			f.flush()
			f.close()

#.dat for gnuplot

			f = open('33TIcation_PC71BM_anion_B3LYP_6_31gSTAR_opt_mul_%s_%gdeg_%gA.dat' % (regionnames[i],rotangle,distA), 'w')
			f.write('Molecule Starts:\n')
			for atom in Register.GetRegister('Atom'):
				astr=str(atom)
				f.write(astr)
				f.write('\n')
			f.write('Molecule Ends.\n\n')
			f.write ('distA, %s \n distbohr, %s \n tot, %s \n energyev, %s \n' % (distA,distbohr,tot,energyev))
			f.write('\n\ndipoles start:\n')
			for dd in split_d:
				dstr=str(dd)
				f.write(dstr)
				f.write('\n')
			f.flush()
			f.close()

			for atom in acc():
				atom().move((-distmoved)*dispversor[i])
				atom().move(-centre[i])
				atom().move((ROTMATINV*atom()._pos.T).T-atom()._pos)

exit()

