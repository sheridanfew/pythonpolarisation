import sys
sys.path.append('../')
from Polarizability import Polarizability
from numpy import pi
from Rotation import Rotation 

a = Polarizability (iso=10.)
b = Polarizability (noniso=[[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])

print "This should give the identity"
print a.T*a.I
print "this should give the matrix: "
print b.T

print "test a+=a"
a+=a
print a
print "test *scalar a=3*a"
a=a*3
print a
print "test c=b b=a a=c c*=2"
c=b.copy()
print c is b
b=a.copy()
a=c.copy()
c*=2
print "a: ", id(a)
print "b: ", id(b)
print "c: ", id(c)

print "rotate the isotropic pol!"
r = Rotation(euler=[pi/4.,0.,0.])
print "rot: ", r
print "b: ", b
print type(b)
c = b.rotate( r)
print "c = b.rotate(r):", c
print type(c)


