'''
- Cartesian 0,0,0 to CoM
- Rotate molecule such that longaxis oriented correctly
- Set direction (perp to longaxis dir, eg. 'longaxis dir * (1,0,0)' ) to define a zero point for angle about longaxis dir)
- Get angle by dot of this and pos'n of perpatom.
- Rotate by angle alpha which is difference between, about longaxis dir, using Rodriguez formula.

Multiply these two rotation matrices to get overall rotation, and apply this to the pol
'''

import sys
sys.path.append('../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
import numpy as np
from math import *
import os.path

qdict={"anion": -1.0, "neut": 0.0, "cation": 1.0}

atomsforlongaxis={"PDIFCN2": 21}
perpatoms={"PDIFCN2": 1}
nmolsincell={"PDIFCN2": 1}

'''
atomsforlongaxis={"TIPS_Pc": 8, "C8_BTBT": 3, "sexithiophene": 6, "Rubrene": 4, "Pc":1, "PTCDA":12}
perpatoms={"TIPS_Pc": 46, "C8_BTBT": 16, "sexithiophene": 30, "Rubrene": 51, "Pc":11, "PTCDA":14}
nmolsincell={"TIPS_Pc": 1, "C8_BTBT": 2, "sexithiophene": 4, "Rubrene": 4, "Pc":2, "PTCDA":1}
'''

#for molecule in ['Pc','C8_BTBT','sexithiophene','TIPS_Pc','Rubrene','PTCDA']:
for molecule in ['PDIFCN2']:
	for molincell in ['a','b','c','d']:
		for charge in ['anion','neut','cation']:
			Importmol_cif=None
			if os.path.isfile(molecule + '_mol' + molincell + '_' + charge + '_aniso_cifstruct_mul.xyz'):
				Importmol_cif=(molecule + '_mol' + molincell + '_' + charge + '_aniso_cifstruct_mul.xyz')
				nameout=str('sp_' + molecule + '_mol' + molincell + '_' + charge + '.xyz')
			if os.path.isfile(molecule + '_mol' + molincell + '_' + charge + '_fit_aniso_cifstruct_mul.xyz'):
				Importmol_cif=(molecule + '_mol' + molincell + '_' + charge + '_fit_aniso_cifstruct_mul.xyz')
				nameout=str('sp_' + molecule + '_mol' + molincell + '_' + charge + '.xyz')
			if (nmolsincell[molecule] == 1) and (os.path.isfile(molecule + '_' + charge + '_aniso_cifstruct_mul.xyz')):
				Importmol_cif=(molecule + '_' + charge + '_aniso_cifstruct_mul.xyz')
				nameout=str('sp_' + molecule + '_' + charge + '.xyz')
			if (nmolsincell[molecule] == 1) and (os.path.isfile(molecule + '_' + charge + '_fit_aniso_cifstruct_mul.xyz')):
				Importmol_cif=(molecule + '_' + charge + '_fit_aniso_cifstruct_mul.xyz')
				nameout=str('sp_' + molecule + '_' + charge + '.xyz')

			if Importmol_cif and not os.path.isfile(nameout):
				print Importmol_cif
				Importmol_gauss=(molecule + '_' + charge + '_fit_aniso_mul.xyz')
				exec( "from Molecules.pol_" + molecule + "_" + charge + " import *")
				pol_gauss=pol
				q=qdict[charge]

				#Initial positions as defined by gaussian, and in cif file.
				atomforlongaxis=atomsforlongaxis[molecule]
				#This should be an atom far from the CoM. This will be used to orient polarisability tensors long axis
				perpatom=perpatoms[molecule]
				#This should be an atom far from the vector from CoM to aotmforlongaxis, this is used to rotate polarisability tensor to correct orientation about long axis

				ReadMoleculeType('../Molecules/./%s' % Importmol_gauss)
				ReadMoleculeType('../Molecules/./%s' % Importmol_cif)

				Mol_gauss=GetMolecule('../Molecules/./%s' % Importmol_gauss)
				Mol_cif=GetMolecule('../Molecules/./%s' % Importmol_cif)

				#Get list of atom posns:

				posns_gauss=[]
				for atom in Mol_gauss():
					print atom()._pos
					posns_gauss.append(np.array(atom()._pos))

				posns_cif=[]
				for atom in Mol_cif():
					print atom()._pos
					posns_cif.append(np.array(atom()._pos))


				#Centres in each case:

				centre_gauss=sum(posns_gauss)/len(posns_gauss)
				centre_cif=sum(posns_cif)/len(posns_cif)

				#posns w. CoM at 0,0,0

				posns_gauss_centre=[]
				posns_cif_centre=[]

				for i in posns_gauss:
					posns_gauss_centre.append(i - centre_gauss)

				for i in posns_cif:
					posns_cif_centre.append(i - centre_cif)


				#Rotation vector to align gauss's molecule with that in the cif file

				print 'doublealign=[', posns_gauss_centre[atomforlongaxis], ',', posns_cif_centre[atomforlongaxis], ',', posns_gauss_centre[perpatom], ',', posns_cif_centre[perpatom], ']'

				rot=Rotation(doublealign=[posns_gauss_centre[atomforlongaxis],posns_cif_centre[atomforlongaxis],posns_gauss_centre[perpatom],posns_cif_centre[perpatom]])

				print 'rot'
				print rot

				#Apply this rotation to gauss moleule:

				posns_gauss_rot=[]
				for i in posns_gauss_centre:
					pos=rot*i.T
					posns_gauss_rot.append(pos.T)


				print '\n\nAtom atomforlongaxis Posn:'
				print 'gauss original atomforlongaxis posn (CoM at 000):'
				print posns_gauss_centre[atomforlongaxis]
				print 'gauss rot (should be same as cif centred):'
				print posns_gauss_rot[atomforlongaxis]
				print 'Cif centred:'
				print posns_cif_centre[atomforlongaxis]

				print '\n\nAtom perpatom Posn:'
				print 'gauss original perpatom posn (CoM at 000):'
				print posns_gauss_centre[perpatom]
				print 'gauss rot (should be same as cif centred):'
				print posns_gauss_rot[perpatom]
				print 'Cif centred:'
				print posns_cif_centre[perpatom]

				stdout=sys.stdout			
				sys.stdout = open('pol_rot_test.txt', 'a')
				print '\n\nMolname: ', Importmol_cif
				print '\nAtom atomforlongaxis Posn:'
				print 'gauss original atomforlongaxis posn (CoM at 000):'
				print posns_gauss_centre[atomforlongaxis]
				print 'gauss one rot:'
				print posns_gauss_rot[atomforlongaxis]
				print 'gauss two rot:'
				print posns_gauss_rot[atomforlongaxis]
				print 'Cif centred:'
				print posns_cif_centre[atomforlongaxis]

				print '\nAtom perpatom Posn:'
				print 'gauss original perpatom posn (CoM at 000):'
				print posns_gauss_centre[perpatom]
				print 'gauss one rot:'
				print posns_gauss_rot[perpatom]
				print 'gauss two rot:'
				print posns_gauss_rot[perpatom]
				print 'Cif centred:'
				print posns_cif_centre[perpatom]
				sys.stdout = stdout

				pol_cif=rot*pol_gauss*rot.T
				#Maybe check T is on correct rotation

				print 'pol gauss:'
				print pol_gauss

				print 'pol cif:'
				print pol_cif



				'''
				print '\n\nALLPOSNS'
				for i in range(0,len(posns_cif_centre),1):
					#print 'Atom ', i
					print 'C ',posns_cif_centre[i][0][0],' ',posns_cif_centre[i][0][1],' ',posns_cif_centre[i][0][2]
					print 'C ',posns_gauss_rot[i][0,0],' ',posns_gauss_rot[i][0,1],' ',posns_gauss_rot[i,0,2]
				'''

				print 'centre_cif', centre_cif

				stdout=sys.stdout			
				sys.stdout = open(nameout, 'w')
				print "{'elname':'", molecule, "','pos':Position([", centre_cif[0,0], ",", centre_cif[0,1], ",", centre_cif[0,2],"]), 'crg':", q, ", 'pol':Polarizability(noniso =[ [", pol_cif[0,0] ,",", pol_cif[0,1] ,",", pol_cif[0,2] ,"], [", pol_cif[1,0] ,",", pol_cif[1,1] ,",", pol_cif[1,2] ,"], [", pol_cif[2,0] ,",", pol_cif[2,1] ,",", pol_cif[2,2] ,"] ])}"
				sys.stdout = stdout

