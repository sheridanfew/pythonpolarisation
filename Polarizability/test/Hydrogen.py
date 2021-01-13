import sys
sys.path.append('../../')
from BasicElements import *

class Hydrogen(Molecule):
	""" an example molecule
	
	

	"""	
	def __init__(self, **kwargs):
		Molecule.__init__(self)
		self._type = "Hydrogen" # when you make your own molecules, make sure that you change the "_type"  string 
					# this will overwrite the "BaseClass" type of Molecule!
		
		#now I want to allow for hydrogen to have user entered polarizability
		# in principle I could also have the charges entered by the user
		polH = kwargs.get('pol', Polarizability(iso=10.) ) 
		self.add_atom(pos=Position([-0.5, 0.,0.]), pol=polH, crg=0., elname='H')
		self.add_atom(pos=Position([ 0.5, 0.,0.]), pol=polH, crg=0., elname='H')

