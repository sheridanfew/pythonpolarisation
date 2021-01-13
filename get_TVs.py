import math as m
import numpy as np

a=7.5650
b=7.7500
c=16.835

alpha=89.15*(m.pi/180)
beta=78.42*(m.pi/180)
gamma=83.63*(m.pi/180)

cell_volume=m.sqrt(1 - (m.cos(alpha)**2) - (m.cos(beta)**2) - (m.cos(gamma)**2) + (2*m.cos(alpha)*m.cos(beta)*m.cos(gamma)))

matrix_to_cartesian=np.matrix( [[a, b*m.cos(gamma), c*m.cos(beta)],
[0, b*m.sin(gamma), c*(m.cos(alpha) - m.cos(beta)*m.cos(gamma))/m.sin(gamma)],
[0, 0, c*cell_volume/m.sin(gamma)]])

TV=matrix_to_cartesian.T

#TVs, TV[0,1,2] are the three translation vectors.

print 'TV:'
print TV

print 'TV[0,1,2] are the three translation vectors.'

