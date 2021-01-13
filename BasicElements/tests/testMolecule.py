import sys
sys.path.append('../')
from Hydrogen import Hydrogen
from CarbonMonoxide import CarbonMonoxide 
from Position import Position
from Polarizability import Polarizability
from Register import GetRegister

h1 = Hydrogen()
h2 = Hydrogen(pol=Polarizability(iso=23.))
h2().place_at_d(Position([10.,0.,0.]) )
h3 = h2().copy()
h3().place_at_d(Position([3.4, 0.,0.]) )
co = CarbonMonoxide()

print "print the class of h1 and the number of atoms in it"
print h1().__class__
print len(h1())
for i in h1():
	print i()._elname, i()._pos
print "print all the registered atoms"
for i, atom in GetRegister('Atom'):
	print atom
print "del the first molecule and pring all the registered atoms"
del(h1)
for i, atom in GetRegister('Atom'):
	print atom
print h2()

print "print all hydrogens:"
for i, mol in GetRegister('Hydrogen'):
	print mol
