from Hydrogen1Pol import Hydrogen1Pol
from Polarizability.GetDipoles import get_dipoles
from BasicElements import *
import numpy as np

h1 = Hydrogen1Pol()

""" I expect that since there is one polarizable dipole with polarizability 10 and a constant electric field 10.,
the total dipole moment should be 100. 100. 100.
I also demonstrate that you can pring the molecules simply enough. If  you want to customise this change the Atom.__repr__(self) method
"""
d = get_dipoles(E0=np.matrix([0.,0.,10.]))
print d
print h1()._type
print h1()
