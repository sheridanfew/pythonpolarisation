import sys
import numpy as np
sys.path.append('../')
from Position import Position as pos
from Rotation import Rotation 

p = pos([1.,2.,3.])
r = Rotation(euler=[0.1, 0.,-0.1])

print type(r)
print "before rot:"
print type(p)
print p
print "after rot:"
p=p.rotate(r)
print p
print type(p)


a = pos([1.,0.,0.])
b = pos([2.,0.,0.])

print np.dot(a,b.T)
print np.dot(b,b.T)

print a.distance(a)
print a.distance(b)
