import sys
sys.path.append('../../')
from BasicElements import *
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import *
from Polarizability.GetEnergy import get_energy
import numpy as np

#inputs=np.array([[8.0,9.5,2.0,8.5]])
#[8.0, 9.5, 2.5, 8.5],
#[8.0, 9.0, 2.5, 8.5],
#[7.5, 9.5, 3.0, 8.5],
#[7.5, 9.5, 3.5, 8.5]])
# format polC, polH, polS, cut, eta
#steps=np.arange(0,inputs.shape[0],1)
distancerangeA = np.arange(0., 8.0, 0.5)

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

ReadMoleculeType('../../Molecules/8Tcat_mul.xyz')
octoT = GetMolecule('../../Molecules/8Tcat_mul.xyz')
ReadMoleculeType('../../Molecules/PCBM_EDITEDO_anion_mul.xyz')
C60 = GetMolecule('../../Molecules/PCBM_EDITEDO_anion_mul.xyz')

cut= 4.00007106634
ModifyPolarizability(octoT(), C=9.072879660383,H=0.100004083159,S=13.2644932888)
ModifyPolarizability(C60(), C=6.27522923483,H=0.100004083159,S=0)
#	print "8T"
#	print octoT()
#	print "C60"
#print C60()
print "Hello!"
print "Cutoff:"
print cut
print "dist, distbohr, tot, energy, energyev"

f = open('8T_PCBM_dist_polaxis_data.csv', 'w')
f.write('dist, distbohr, pol, energy, energyev\n')
f.flush()
f.close()

for dist in distancerangeA:
#		Molecule.get_com(C60())
#		matrix([[ -1.07281416,  -0.02480016,  10.45065994]])
	distbohr=(np.divide(dist+3.3,0.529177))
	for atom in C60():
		moveto=np.matrix([30.13+distbohr,0.,-np.divide(3.3,0.529177)])
		atom().move(moveto)
#		atom().move([30.13+distbohr,0.,-np.divide(3.3,0.529177)])
#	Molecule.place_at_d(C60(), moveto)
#	print C60()
	E0 = np.matrix([0.,0.,0.])
	#print "E0"
	#print E0
	d = get_dipoles(E0=E0,cutoff=cut)
	split_d = split_dipoles_onto_atoms(d)
	tot = np.matrix([0.,0.,0.])
	for dd in split_d:
		tot += dd
	energy = get_energy(E0=E0,cutoff=cut)
	energyev = np.multiply(energy,27.211)

	print dist, distbohr, tot, energy, energyev

	f = open('8T_PCBM_dist_polaxis_data.csv', 'a')
	propstr='%s 	%s	 %s 	%s 	%s' % (dist, distbohr, tot, energy, energyev)
#	propstr=str(dist, distbohr, tot, energy, energyev)
	f.write(propstr)
	f.write('\n')
	f.flush()
	f.close()


	f = open('8T_PCBM_polaxis_%gA.dat' % dist, 'w')
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
	for atom in C60():
		atom().move(-moveto)
#		atom().move([0.,0.,-distbohr])
#print "C60 moved"

#print C60()



#print "polarizabilityC(au^3) polarizabilityS(au^3) polarizabilityH(au^3) cutoff(au) total dipole  "


#E0 = np.matrix([0.,0.,1.])
#print "E0"
#print E0
#d = get_dipoles(E0=E0,cutoff=cut)
#split_d = split_dipoles_onto_atoms(d)
#tot = np.matrix([0.,0.,0.])
#for dd in split_d:
#	tot += dd	
#	print cut, tot

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
