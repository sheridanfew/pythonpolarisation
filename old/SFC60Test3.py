from Hydrogen2Pol import Hydrogen2Pol
from Polarizability.GetDipoles import get_dipoles as get_dips
from Polarizability.GetDipoles import split_dipoles_onto_atoms as split_dips
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from BasicElements import *
import numpy as np

""" here I want to show how smart the code is 

Inside your main loop you create a new hydrogen molecule at each step. the code can keep track of whether the molecules exist or not 
and thus the get_dips code will still work. Notice that no matter what you do this code will NEVER reproduce the result with 1 single
polarizable atom. The total dipole never goes over ~20 or so.

"""

"""ReadMoleculeType('extras/c60_neut_B3LYP_631GSTAR_chelpg.xyz')
c60 = GetMolecule('extras/c60_neut_B3LYP_631GSTAR_chelpg.xyz')
print c60()"""

for polariz in np.arange(5., 100., 1.):
	c60 = GetMolecule('extras/c60_neut_B3LYP_631GSTAR_mul.xyz')
	c60(polariz)
	dips = get_dips(E0=np.matrix([0.,0.,10]))
	res = split_dips(dips)
	total = np.matrix([0.,0.,0.])
	for dip in res:
		total += dip
	print polariz, total	,len(res)

