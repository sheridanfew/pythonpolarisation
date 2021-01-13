import sys
from BasicElements import Molecule,Position, Polarizability
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import GetElectricField
from BasicElements import *
import numpy as np


ReadMoleculeType('extras/c60_neut_B3LYP_631GSTAR_mul.xyz')
c60 = GetMolecule('extras/c60_neut_B3LYP_631GSTAR_mul.xyz')
print c60()

d = get_dipoles(E0=np.matrix([0.,0.,10.]))

print "SF Dipoles:"
print d

print "SF Split Dipoles:"
split_d = split_dipoles_onto_atoms(d)
for dd in split_d:
     print dd
    
"""print "SF c60 type"
print c60()._type
print "SF c60:"
print c60()"""

"""based on James simpletest.py"""

