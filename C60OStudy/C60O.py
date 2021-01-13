import sys
sys.path.append('../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
#from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import *
from Polarizability.GetEnergyFromDips import *
#from Polarizability.JMatrix import JMatrix
import numpy as np
from math import *
from time import gmtime, strftime

print strftime("%a, %d %b %Y %X +0000", gmtime())

###################################
#START OF MOLECULE SPECIFIC SECTION
###################################

moltypea='C60O_CAM_B3LYP_6_31PLUSgSTAR_2014_02_10_ESP_ESP.xyz'
moltypeb='C60O_CAM_B3LYP_6_31PLUSgSTAR_2014_02_10_ESP_ESP.xyz'

#moltypea='C60O_HF_6_31PLUSgSTAR_2014_02_10_ESP_ESP.xyz'
#moltypeb='C60O_HF_6_31PLUSgSTAR_2014_02_10_ESP_ESP.xyz'

cut=8.0

ReadMoleculeType('Molecules/./%s' % moltypea)
ReadMoleculeType('Molecules/./%s' % moltypeb)

mola=GetMolecule('Molecules/./%s' % moltypea)
molb=GetMolecule('Molecules/./%s' % moltypeb)

#Furthest in + x 6.786025
#{'elname':'C', 'pos':Position([ 6.78602555060951227679, 1.31847895185204502385, .00010960411748654325]), 'crg':0.0004380, 'pol':Polarizability(iso=9.0799)}
#Furthest in -x -8.903505
#{'elname':'O', 'pos':Position([ -8.90350529085531273834, .00010771439132298216, .00001133835698136654]), 'crg':-0.4767920, 'pol':Polarizability(iso=5.160)}


# Move to 0 disp point where highest x atom of mola is in same x pos. as lowest x mol in molb.
for atom in molb():
	atom().move(np.matrix([[ (6.78602555060951227679 + 8.90350529085531273834), 0., 0. ]]))

print 'same orient'
print 'dist	Energy(eV)'

for distA in np.arange(0.0,6.2,0.2):
#	print 'same orient, dist', distA, 'A'
	dist=distA/0.5291772109217
	for atom in molb():
		atom().move(np.matrix([ dist, 0., 0. ]))

	potential = get_potential()
	Uqq = np.multiply(get_U_qq(potential=potential),27.211)
	print distA, Uqq
	molacom=mola().get_com()
	molbcom=molb().get_com()

#	print 'dist ', dist, ', molacom ', molacom, ', molbcom', molbcom

	for atom in molb():
		atom().move(np.matrix([ -dist, 0., 0. ]))

for atom in molb():
	atom().move(np.matrix([[ (-6.78602555060951227679 - 8.90350529085531273834), 0., 0. ]]))
	atom()._pos=-atom()._pos
	atom().move(np.matrix([[ (-2 * 8.90350529085531273834), 0., 0. ]]))

for distA in np.arange(0.0,4.0,0.5):
	print 'Os facing, dist', distA, 'A'
	dist=distA/0.5291772109217
	for atom in molb():
		atom().move(np.matrix([ -dist, 0., 0. ]))

	potential = get_potential()
	Uqq = np.multiply(get_U_qq(potential=potential),27.211)
	print 'Uqq', Uqq
	molacom=mola().get_com()
	molbcom=molb().get_com()

	print 'dist ', dist, ', molacom ', molacom, ', molbcom', molbcom

	for atom in molb():
		atom().move(np.matrix([ dist, 0., 0. ]))

print 'Align O of mol b with Cs nearest O in mol a in x. Move up in z'

for atom in molb():
	atom().move(np.matrix([[ (8.90350529085531273834 -6.28540874398972593679), 0., 0. ]])) 

for distA in np.arange(0.0,3.0,0.1):
	dist=distA/0.5291772109217
	for atom in molb():
		atom().move(np.matrix([ 0., 0., dist ]))
	potential = get_potential()
	Uqq = np.multiply(get_U_qq(potential=potential),27.211)
	print distA, Uqq
	molacom=mola().get_com()
	molbcom=molb().get_com()

	#print 'dist ', dist, ', molacom ', molacom, ', molbcom', molbcom


#{'elname':'C', 'pos':Position([ -1.38855188772305382771, -1.88071028003474072579, -6.28540874398972593679]), 'crg':0.0000630, 'pol':Polarizability(iso=9.0799)}

'''
f = open('%s_size_%s_new.dat' % (name, size), 'w')
f.write('Molecule Starts:\n')
for atom in Register.GetRegister('Atom'):
	astr=str(atom)
	f.write(astr)
	f.write('\n')
f.write('Molecule Ends.\n\n')
f.write ("Crystal Struct\tcentres[0]\tsurroundings[0]\tSize\tEnergy\tUqq\tUqd\tUdd\tDipole Moment")
f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (name,centres[0],surroundings[0],size,energyev,Uqq,Uqd,Udd,tot))
f.flush()
f.close()
'''

