import numpy as np

def parse_input( **kwargs):
	for key in kwargs:
		if key == 'iso' :
			p = kwargs[key]
			return  [ [p,0.,0.],[0.,p,0.],[0.,0.,p]] 
		elif key == 'noniso':
			mp = kwargs[key]
			return mp

		else:
			raise TypeError('cant understand polarizability', kwargs)

class Polarizability(np.matrix):
	 """ A matrix containing the polarization tensor for an atom

	 can be initialised either with iso=p for an isotropic tensor pI
	 or noniso =[[a,b,c],[d,e,f],[g,h,i]] for a non isotropic tensor with elements abc, def, ghi

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
	
	 def rotate(self, rot):
		rotated_array = rot.T * self * rot ; #SHERIDAN: check that this is the right way round!
							# I am too lazy to figure out if it should be rot^t * pol * rot
							# or rot *pol*rot^t or what - or whether it even makes a difference
		return Polarizability(noniso=rotated_array.tolist())

