import sys
sys.path.append('../')
from Atom import Atom
from Position import Position
from Polarizability import Polarizability

a=Atom(pos=Position([0.,0.,0.]), pol=Polarizability(iso=10.), crg=-1., elname='H')
print a()
print "here"
print a().__class__
print a.__class__
print a()._pos
print a()._pol
print a()._crg
print a()._elname
print a().__class__
print a()._pos.__class__
print a()._pol.__class__
print a()._crg.__class__
print a()._elname.__class__
