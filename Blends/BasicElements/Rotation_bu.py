import numpy as np
from numpy import cos, sin, matrix
from math import *

def parse_input( **kwargs):
	for key in kwargs:
		if key == 'euler':
			#Specify by rot=Rotation.Rotation(euler=[theta,phi,psi])
			p = kwargs[key]
			theta = p[0]
			phi   = p[1]
			psi   = p[2]
			D = Rotation([[cos(phi), sin(phi),0.],[-sin(phi), cos(phi), 0.],[0.,0.,1.]])
			B = Rotation([[cos(psi), sin(psi),0.],[-sin(psi), cos(psi), 0.],[0.,0.,1.]])
			C = Rotation([[1.,0.,0.],[0., cos(theta), sin(theta)], [0, -sin(theta),cos(theta)]])
			return  (B*C*D)
		elif key == 'rodrigues':
			#See http://math.stackexchange.com/questions/142821/matrix-for-rotation-around-a-vector
			p = kwargs[key]
			ux =  p[0]/sqrt(p[0]**2+p[1]**2+p[2]**2)
			uy =  p[1]/sqrt(p[0]**2+p[1]**2+p[2]**2)
			uz =  p[2]/sqrt(p[0]**2+p[1]**2+p[2]**2)

			rotangle =  p[3]

			W=Rotation([[0., -uz, uy],[uz ,0., -ux],[-uy, ux, 0.]])
			I=Rotation(np.eye(3))

			return (I + (sin(rotangle)*W) + ( 2*((sin(rotangle/2)**2))*(W**2) ))
		elif key == 'align':
			#Rot to align vector a with b using Rodrigues. NOTE: There is still a free axis of rotation about the vector from the origin to b.
			#Specify by rot=Rotation.Rotation(align=[np.matrix([a0,a1,a2]),np.matrix([b0,b1,b2])])
			p = kwargs[key]
			startvec = p[0]/np.linalg.norm(p[0])
			finvec   = p[1]/np.linalg.norm(p[1])

			if startvec[0,0] == finvec[0,0] and startvec[0,1] == finvec[0,1] and startvec[0,2] == finvec[0,2]:
				return(Rotation(np.eye(3)))
			else:
				#Versor perpendicular to both of thess
				perpvers=np.matrix(np.cross(startvec,finvec))/np.linalg.norm(np.matrix(np.cross(startvec,finvec)))
				angle=float(np.arccos(np.dot(startvec,finvec.T)))
				rot=Rotation(rodrigues=[perpvers[0,0],perpvers[0,1],perpvers[0,2],angle])
				return rot
		elif key == 'mat':
			mp = kwargs[key]
			return mp
		else:
			raise TypeError('cant understand orientation: ', kwargs)

'''
EXAMPLE:

def greet_me(**kwargs):
    if kwargs is not None:
        for key, value in kwargs.iteritems():
            print "%s == %s" %(key,value)
 
>>> greet_me(name="yasoob")
name == yasoob

'''

class Rotation(np.matrix):
	 """ A matrix containing an orientation vector for an atom

	 can be initialised either with euler=[theta, phi, psi] for a rotation matrix defined by three euler angles 
	 (defined as in http://mathworld.wolfram.com/EulerAngles.html)
	 or mat =[[a,b,c],[d,e,f],[g,h,i]] for a rotation  with elements abc, def, ghi

	 """

	 def __new__(cls, data=None, dtype=None, copy=True, **kwargs):
#		print "enter new:"
#		print "cls: ", cls
#		print "data: ", data
#		print "kwargs: ", kwargs
#		print "copy: ", copy
#		print "dtype: ", dtype
		if len(kwargs) > 0:
			 data2 = parse_input(**kwargs)
			 return np.matrix.__new__(cls, data2, dtype, copy)
		else:
			 return np.matrix.__new__(cls, data, dtype, copy)

	 def __array_finalize__(self, obj):
		pass


