from Hydrogen2Pol import Hydrogen2Pol
from Polarizability.GetDipoles import get_dipoles as get_dips
from Polarizability.GetDipoles import split_dipoles_onto_atoms as split_dips
from BasicElements import *
import numpy as np

""" make sure linear polarisability

"""
for cutoff in np.arange(0., 2., 0.02):
        for E in np.arange(2., 10., 2.):
                h = Hydrogen2Pol(100.)
                dips = get_dips(E0=np.matrix([0.,0.,E]), cutoff=cutoff)
                res = split_dips(dips)
                total = np.matrix([0.,0.,0.])
                for dip in res:
                        total += dip
                        #adds dipoles on each polarisable atom and calcs total
                print E, cutoff, total/E , len(res)

