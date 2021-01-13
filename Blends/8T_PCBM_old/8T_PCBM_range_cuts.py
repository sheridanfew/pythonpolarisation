import sys
sys.path.append('../../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import *
from Polarizability.GetEnergy import get_energy
import numpy as np

inputs=np.array([[3.0,9.0,7.5,4.0,1.05941444760194798725],
[4.0,1.0,8.0,4.5,1.42233774433311802557],
[4.5,6.5,9.5,5.0,1.30409941807122658767],
[5.5,9.5,7.5,5.5,1.30223385333305501688],
[4.5,1.0,1.5,6.0,1.20533292771491485967],
[5.5,1.5,1.0,6.5,.85757871464495679800],
[6.5,1.0,1.5,7.0,.60523347988331955086],
[6.5,1.0,4.5,7.5,.39327635436282564453],
[7.0,1.0,6.5,8.0,.22848396756942164825],
[7.0,1.0,9.5,8.5,.10394912751790246790],
[7.5,4.5,9.5,9.0,.03925800265012225692],
[7.5,9.5,9.5,9.5,.02144662554554283826],
[8.5,9.5,9.5,10.0,.10201904078902975420],
[8.5,9.5,9.5,10.5,.22501199177461277630],
[8.5,9.5,9.5,11.0,.32148915550187491390],
[8.0,9.5,9.5,11.5,.38881548168799920433]])
# format polC, polH, polS, cut, eta
steps=np.arange(0,inputs.shape[0],1)
distancerangeA = np.arange(2., 11., 1.)

def ModifyPolarizability(molecule, C, H, S):
	""" takes a mol and changes the isotropic polarizability
	for each atom
	"""
	C = Polarizability(iso=C)
	H = Polarizability(iso=H)
	S = Polarizability(iso=S)
	for atom in molecule:
		if (atom()._elname=="C"):
			atom()._pol= C
		elif (atom()._elname=="H"):
			atom()._pol= H
		elif (atom()._elname=="S"):
			atom()._pol= S

ReadMoleculeType('../../8T_no_charge.xyz')
octoT = GetMolecule('../../8T_no_charge.xyz')
ReadMoleculeType('../../c60.xyz')
C60 = GetMolecule('../../c60.xyz')

for i in steps:
	cut=inputs[i,3]
	ModifyPolarizability(octoT(), C=inputs[i,0],H=inputs[i,1],S=inputs[i,2])
	ModifyPolarizability(C60(), C=inputs[i,0],H=inputs[i,1],S=inputs[i,2])
#	print "8T"
#	print octoT()
#	print "C60"
#	print C60()
	print "Hello!"
	print "Cutoff:"
	print cut
	print "cut, dist, distbohr, tot, energy, energyev"
	for dist in distancerangeA:
		distbohr=np.divide(dist,0.529177)
		Molecule.place_at_d(C60(), [0.,0.,distbohr])
		#print C60()

		E0 = np.matrix([0.,0.,0.])
		#print "E0"
		#print E0
		d = get_dipoles(E0=E0,cutoff=cut)
		split_d = split_dipoles_onto_atoms(d)
		tot = np.matrix([0.,0.,0.])
		for dd in split_d:
			tot += dd
		energy = get_energy(E0=E0,cutoff=cut)
		energyev = np.divide(energy,27.211)

		print cut, dist, distbohr, tot, energy, energyev
print "C60 moved"

print C60()



#print "polarizabilityC(au^3) polarizabilityS(au^3) polarizabilityH(au^3) cutoff(au) total dipole  "


E0 = np.matrix([0.,0.,1.])
#print "E0"
#print E0
d = get_dipoles(E0=E0,cutoff=cut)
split_d = split_dipoles_onto_atoms(d)
tot = np.matrix([0.,0.,0.])
for dd in split_d:
	tot += dd	
	print cut, tot

"""
	switch = {
	    'C': atom()._pol= C
	    'H': atom()._pol= H
	    'S': atom()._pol= S
	    }

	try:
	    result = switch[choice]
	except KeyError:
	    print 'I didn\'t understand your choice.'
	else:
	    result()
"""
