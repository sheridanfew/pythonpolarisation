from Hydrogen import Hydrogen
from CarbonMonoxide import CarbonMonoxide 
import numpy as np

import sys
sys.path.append('../')
sys.path.append('../../')
from BasicElements.Register import GetRegister
from BasicElements.Register import RegisterOfRegisters
from BasicElements.Rotation import Rotation
from BasicElements.Position import Position
from JMatrix import JMatrix
from GetElectricField import get_electric_field
from GetDipoles import get_dipoles


#define the position of all the molecules
co = CarbonMonoxide()
h1 = Hydrogen()
h2 = Hydrogen()
rot = Rotation(euler=[np.pi/4, np.pi/2, -np.pi/3])
h2().rotate(rot) 
co().place_at_d(Position([4.,0.,0.]))
h2().place_at_d(Position([-3.5, 0.,2.]))

np.set_printoptions(precision=2, suppress=True)
jm = JMatrix()
print jm._m


for i,atom in GetRegister('Atom'):
	print atom._elname , atom._pos

print "all registered class types:"
for k in RegisterOfRegisters.Instance():
	print k
reg = RegisterOfRegisters.Instance()['<class \'BasicElements.Atom.Atom\'>']
example = jm.get_element(reg.Instance()[2], reg.Instance()[4])
print example

example2 = jm.get_element(reg.Instance()[1], reg.Instance()[3])
print example2
print example2.__class__

print "now test the get electric field thing"
field = get_electric_field()
print field

field2 = get_electric_field(np.matrix([10.,10., 10.]))
print field2

print "now test the dipoles thingy"
dips = get_dipoles(E0=np.matrix([10.,20., 30.]))
print dips 
