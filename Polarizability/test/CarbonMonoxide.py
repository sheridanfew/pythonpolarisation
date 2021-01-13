import sys
sys.path.append('../../')
from BasicElements import *

class CarbonMonoxide(Molecule):
	""" an example molecule
	
	

	"""	
	def __init__(self, **kwargs):
		Molecule.__init__(self)
		self._type = "CarbonMonoxide" # when you make your own molecules, make sure that you change the "_type"  string 
					# this will overwrite the "BaseClass" type of Molecule!
		
		#now I want to allow for hydrogen to have user entered polarizability
		# in principle I could also have the charges entered by the user
		polC = kwargs.get('polC', Polarizability(iso=0.) ) 
		polO = kwargs.get('polO', Polarizability(iso=0.) ) 
		self.add_atom(pos=Position([-0.75, 0.,0.]), pol=polC, crg=.2, elname='C')
		self.add_atom(pos=Position([ 0.75, 0.,0.]), pol=polO, crg=-.2, elname='O')

