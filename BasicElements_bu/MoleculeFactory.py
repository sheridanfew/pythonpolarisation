from Singleton import Singleton
from Molecule  import Molecule
from Position import Position
from Polarizability import Polarizability
import re

class MoleculeFactory(dict):

        def register(self, cls, registration):
                if registration not in self:
                        self[registration] = cls

        def unregister(self, registration):
                del self[registration]

MF = MoleculeFactory()

def ReadMoleculeType(namefile):
	"""used to read in a molecule and register it in the molecule factory.

	usage: 
		>>ReadMoleculeType('H2.xyz')   # read in the molecule in H2.xyz and register it with that name
		>>h2 = GetMoleculeType('H2.xyz') # returns an instance of molecule H2.xyz


	 """
        class aMolecule(Molecule):
                def __init__(self):
                        self._type = namefile
			filein = open(namefile, 'r')
			if not filein:
				raise IOError('cannot open ' + namefile )
			for line in filein: 
				self.add_atom(**dict(eval(line)));

	
        MF.register(aMolecule,namefile)
 
def GetMolecule(namefile):
	return MF[namefile]()

