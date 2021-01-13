from SFHydrogen1TestDefPol import Hydrogen1Pol
from Polarizability.GetDipoles import get_dipoles as get_dips
from Polarizability.GetDipoles import split_dipoles_onto_atoms as split_dips
from BasicElements import *
import numpy as np

""" cutoff should make no difference & dipole should always be 10*E

"""
for cutoff in np.arange(0., 2., 0.02):
	h = Hydrogen1Pol(100.)
	dips = get_dips(E0=np.matrix([0.,0.,10]), cutoff=cutoff)
	res = split_dips(dips)
	total = np.matrix([0.,0.,0.])
	for dip in res:
		total += dip
		#adds dipoles on each polarisable atom and calcs total
	print cutoff, total	,len(res)

