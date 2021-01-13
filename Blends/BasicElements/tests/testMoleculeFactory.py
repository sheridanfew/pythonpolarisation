import sys
sys.path.append('../')

from MoleculeFactory import ReadMoleculeType
from MoleculeFactory import GetMolecule

ReadMoleculeType('H2.xyz')
h2 = GetMolecule('H2.xyz')
h3 = GetMolecule('H2.xyz')	
print h2()


