import sys
sys.path.append('../')
from BasicElements import Molecule,Position, Polarizability
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
#from Polarizability.GetDipoles import *
from BasicElements import *
import numpy as np

def ModifyPolarizability(molecule, **kwargs):
	""" takes a mol and changes the isotropic polarizability
	for each atom
	"""
	jmtype = kwargs.get('jmtype', '')
	fittype = kwargs.get('fittype', '')

	print "keywords understood:"
	for key in kwargs:
       		print "%s: %s" % (key, kwargs[key])

	if jmtype=='TholeLin' and fittype=='empirical':
		H=3.5020
		C=10.1756
		N=7.6048
		O=6.3940
		F=2.9413
		S=19.7422
		Cl=16.1145
		Br=22.6671
		I=34.2434

	if jmtype=='TholeExp' and fittype=='empirical':
		H=2.7929
		C=8.6959
		N=6.5565
		O=5.7494
		F=3.0013
		S=16.6984
		Cl=16.1979
		Br=23.5714
		I=36.9880

	if jmtype=='TholeLin' and fittype=='mean':
		H=1.5350
		C=8.7622
		N=5.2663
		O=2.6175
		F=0.5819
		S=15.0931
		Cl=7.0741
		Br=12.1193
		I=19.6286

	if jmtype=='TholeExp' and fittype=='mean':
		H=1.3368
		C=7.9421
		N=4.5428
		O=2.5694
		F=1.1109
		S=13.9500
		Cl=9.0475
		Br=15.7347
		I=23.4488

	if jmtype=='TholeLin' and fittype=='components':
		H=0.5974
		C=7.7764
		N=4.7575
		O=2.6442
		F=0.8149
		S=11.0208
		Cl=7.5682
		Br=12.7339
		I=20.2867

	if jmtype=='TholeExp' and fittype=='components':
		H=1.3849
		C=6.5780
		N=2.9855
		O=2.3184
		F=1.4906
		S=10.8651
		Cl=6.7129
		Br=12.3384
		I=20.0081

	# Set bonds to (close to) zero unless told otherwise
	single=0.0001
	double=0.0001
	trile=0.0001
	pi=0.0001
	dpi=0.0001

	# Sets values for atom types if specified
	#for key in kwargs:
	#	exec ( "" + key + "=" kwargs[key] )

	for atom in molecule:
		for atomtype in ['H','C','N','O','F','S','Cl','Br','I','single','double','triple','pi','dpi']:
			if (atom()._elname.upper()==atomtype.upper()):
				if atom()._pol[0,0]==atom()._pol[1,1] and atom()._pol[0,0]==atom()._pol[2,2]:
					#If isotropic, then make iso to this value
					exec( "atom().setpol(3.0*" + atomtype + "*atom()._pol/np.trace(atom()._pol))")
				else:
					#If anisotropic, then make trace = this value
					exec( "atom().setpol(" + atomtype + "*atom()._pol/np.trace(atom()._pol))")
					#Will conserve np.trace at 1 * aniso value



"""
Usage example:

mol = GetMolecule('12T_no_charge.xyz')

#ModifyPolarizability(mol(), 4.)
#print mol()

polrange    = np.arange(4., 9., .5)
cutoffrange = np.arange(3., 10., .5)

print "# polarizability(au^3) cutoff(au) total dipole  "
for p in polrange:
	ModifyPolarizability(mol(), p)
	for  c in cutoffrange:
		d = get_dipoles(E0=np.matrix([10.,0.,0.]),cutoff=c)
		split_d = split_dipoles_onto_atoms(d)
		tot = np.matrix([0.,0.,0.])
		for dd in split_d:
		     tot += dd	
		print p, c, tot[0,0]
"""
