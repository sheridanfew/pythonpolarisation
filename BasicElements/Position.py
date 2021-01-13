from numpy import matrix, sqrt
from Rotation import Rotation
import numpy as np

class Position(matrix):
	""" inherit matrix and call it a position
	
		
	"""
	pass
	
	def rotate(self, rot):
		if isinstance(rot, Rotation):
			print self
			print rot
			return (self*rot)
		raise TypeError('positions myust be rotated by rotations')
	def abs(self):
		return sqrt(np.dot(self, self.T)[0,0])
	
	def distance(self, p2):
		d = p2-self
		return d.abs()

