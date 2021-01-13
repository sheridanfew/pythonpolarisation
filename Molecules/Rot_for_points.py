'''
- Cartesian 0,0,0 to CoM
- Rotate molecule such that longaxis oriented correctly
- Set direction (perp to longaxis dir, eg. 'longaxis dir * (1,0,0)' ) to define a zero point for angle about longaxis dir)
- Get angle by dot of this and pos'n of perpatom.
- Rotate by angle alpha which is difference between, about longaxis dir, using Rodriguez formula.

Multiply these two rotation matrices to get overall rotation, and apply this to the polarigaussility
'''

import sys
sys.path.append('../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
import numpy as np

from Molecules.pol_TIPS_Pc_anion import *
pol_gauss=pol

#Initial positions as difined by gauss, and in cif file.

atomforlongaxis=11
#This should be an atom far from the CoM. This will be used to orient polarisability tensors long axis
perpatom=84
#This should be an atom far from the vector from CoM to aotmforlongaxis, this is used to rotate polarisability tensor to correct orientation about long axis


Importmol_gauss=('TIPS_Pc_anion_aniso_no_charge.xyz')
Importmol_cif=('TIPS_Pc_anion_aniso_chelpg.xyz')

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

#Versor in direction of long axis from centre

longaxis_vers_gauss=posns_gauss_centre[atomforlongaxis]/np.linalg.norm(posns_gauss_centre[atomforlongaxis])
longaxis_vers_cif=posns_cif_centre[atomforlongaxis]/np.linalg.norm(posns_cif_centre[atomforlongaxis])

print 'longaxis_vers_gauss', longaxis_vers_gauss
print 'longaxis_vers_cif', longaxis_vers_cif

print 'length longaxis_vers_gauss', (longaxis_vers_gauss[0,0]**2 + longaxis_vers_gauss[0,1]**2 + longaxis_vers_gauss[0,2]**2)
print 'length longaxis_vers_cif', (longaxis_vers_cif[0,0]**2 + longaxis_vers_cif[0,1]**2 + longaxis_vers_cif[0,2]**2)

#Versor perpendicular to both of these

perpvers=np.matrix(np.cross(longaxis_vers_gauss,longaxis_vers_cif))

#Rotation angle required about this vector

longaxisangle=float(np.arccos(np.dot(longaxis_vers_gauss,longaxis_vers_cif.T)))

print 'longaxisangle',longaxisangle

#Rotation vector to align gauss's molecule with that in the cif file

rot1=Rotation.Rotation(rodrigues=[perpvers[0,0],perpvers[0,1],perpvers[0,2],longaxisangle])

#Apply this rotation to gauss moleule:

posns_gauss_rot=[]
for i in posns_gauss_centre:
	pos=rot1*i.T
	posns_gauss_rot.append(pos.T)

#Rotate about longaxis direction
#Define versor between closest point to line running through longaxis dir and perpatom, use procedure as above to rotate gauss molecule.
#Defining versor, see http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation

perpvect_gauss=posns_gauss_rot[perpatom]-(np.dot(posns_gauss_rot[perpatom],longaxis_vers_cif.T)*longaxis_vers_cif)
perpvers_gauss=perpvect_gauss/np.linalg.norm(perpvect_gauss)

perpvect_cif=posns_cif[perpatom]-(np.dot(posns_cif[perpatom],longaxis_vers_cif.T))*longaxis_vers_cif
perpvers_cif=perpvect_cif/np.linalg.norm(perpvect_cif)


print 'perpvers_gauss', perpvers_gauss
print 'perpvers_cif', perpvers_cif

#Versor perpendicular to both of these

perpvers=np.matrix(np.cross(perpvers_gauss,perpvers_cif))

#Rotation angle required about this vector

rotangle=float(np.arccos(np.dot(perpvers_gauss,perpvers_cif.T)))

rot2=Rotation.Rotation(rodrigues=[perpvers[0,0],perpvers[0,1],perpvers[0,2],rotangle])

posns_gauss_fin=[]
for i in posns_gauss_rot:
	pos=rot2*i.T
	posns_gauss_fin.append(pos.T)

print 'Atom atomforlongaxis Posn:'
print 'gauss original follow centre:'
print posns_gauss_centre[atomforlongaxis]
print 'gauss one rot:'
print posns_gauss_rot[atomforlongaxis]
print 'gauss two rot:'
print posns_gauss_fin[atomforlongaxis]
print 'Cif:'
print posns_cif_centre[atomforlongaxis]

print 'Atom perpatom Posn:'
print 'gauss original follow centre:'
print posns_gauss_centre[perpatom]
print 'gauss one rot:'
print posns_gauss_rot[perpatom]
print 'gauss two rot:'
print posns_gauss_fin[perpatom]
print 'Cif:'
print posns_cif_centre[perpatom]


rot_tot=rot1*rot2

polarigaussility_cif=rot_tot*polarigaussility_gauss*rot_tot.T
#Maybe check T is on correct rotation

print 'Polarigaussility gauss:'
print polarigaussility_gauss

print 'Polarigaussility cif:'
print polarigaussility_cif
