import numpy as np
from numpy import cos, sin, matrix
from math import *
DEBUG='FALSE'

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
			#print 'Rot to align vector a with b using Rodrigues. NOTE: There is still a free axis of rotation about the vector from the origin to b.'
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
		elif key == 'doublealign':
			# Makes rotation First aligns p[0] with p[1], then rotates about vector from orgin to p[1] to align rotated p[2] with p[3]
			# Specify by rot=Rotation.Rotation(doublealign=[np.matrix([a0,a1,a2]),np.matrix([b0,b1,b2]),np.matrix([c0,c1,c2]),np.matrix([d0,d1,d2])])
			p = kwargs[key]
			rot1=Rotation(align=[p[0],p[1]])

			longaxisvers=p[1]/np.linalg.norm(p[1])

			startvecrot = (rot1*p[2].T).T
			finvec   = p[3]

			perpvect_start=startvecrot-(np.dot(startvecrot,longaxisvers.T)*longaxisvers)
			#vector in direction of rot1*p[2] from closest point on line in p[1] direction
			perpvers_start=perpvect_start/np.linalg.norm(perpvect_start)

			perpvect_fin=p[3]-(np.dot(p[3],longaxisvers.T)*longaxisvers)
			#vector in direction of p[3] from closest point on line in p[1] direction
			perpvers_fin=perpvect_fin/np.linalg.norm(perpvect_fin)


			#Rotation angle required about this vector

			#rotangle=float(np.arccos(np.dot(perpvers_start,perpvers_fin.T)))
			#rotvers=np.cross(perpvers_start,perpvers_fin)/np.linalg.norm(np.cross(perpvers_start,perpvers_fin))

			rot2=Rotation(align=[perpvers_start,perpvers_fin])

			rot=rot2*rot1

			if DEBUG == 'TRUE':
				print 'startvecrot-(np.dot(startvecrot,p[1].T)*p[1])'
				print 'startvecrot',startvecrot
				print '(np.dot(startvecrot,p[1].T)*p[1])',(np.dot(startvecrot,longaxisvers.T)*longaxisvers)
				print 'rot2', rot2
				print 'rot', rot
				print 'np.dot(perpvers_start,longaxisvers.T) (should be 0)', np.dot(perpvers_start,longaxisvers.T)
				print 'np.dot(perpvers_fin,longaxisvers.T) (should be 0)', np.dot(perpvers_fin,longaxisvers.T)
				print 'perpvers_start', perpvers_start
				print 'perpvers_fin', perpvers_fin

			if DEBUG == 'TRUE':
				perpvers=np.matrix(np.cross(perpvers_start,perpvers_fin))/np.linalg.norm(np.matrix(np.cross(perpvers_start,perpvers_fin)))
				#Versor perpendicular to both of these, should be in p[1] direction

				print 'p[0]', p[0]
				print 'p[1]', p[1]
				print 'p[2]', p[2]
				print 'p[3]', p[3]

				print 'perpvers_start', perpvers_start
				print 'perpvers_fin', perpvers_fin

				print 'rot1*p[0]', rot1*p[0].T
				print 'rot1*p[2]', rot1*p[2].T

				print 'rot*p[0]', rot*p[0].T
				print 'rot*p[2]', rot*p[2].T

				print 'perpvers, should align with p[1]', perpvers


			return rot


		elif key == 'mat':
			mp = kwargs[key]
			return mp
		else:
			raise TypeError('cant understand orientation: ', kwargs)


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


