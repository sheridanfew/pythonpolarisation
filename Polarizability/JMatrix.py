import sys
sys.path.append('../')

from BasicElements import *
import numpy as np
from math import exp
from BasicElements.Register import GetRegister

def get_atoms_with_polarizability():
	n=0
	for i, atom in GetRegister('Atom'):
		if atom._haspol:	
			n+=1
	return n

class JMatrix(object):
	"""J Matrix creates the J matrix, as defined in Stern paper
	
	by looking up the register of atoms, it should automatically fill in the whole
	matrix when the object is initialised JMatrix._m will return the underlying matrix
	the JMatrix can be initialised with a certain cutoff 

	"""
	def __init__(self, **kwargs):
		natom   = get_atoms_with_polarizability()
		self._m = np.asmatrix(np.zeros((natom*3,natom*3)) ) # initialise an empty matrix
		self._jmtype = kwargs.get('jmtype', 'Stern') 
		self._fittype = kwargs.get('fittype', 'empirical') 
		tholeparameters = kwargs.get('tholeparameters','FALSE')

		if self._jmtype == 'Stern':   
			print 'STERN: NB. Should edit to allow atoms to have their own cutoff values. According to Stern, should sum the cutoff values for the two atoms involved, and compare to distance between those atoms. Atm this script just takes a single specified value for cutoff.'
			if kwargs.get('cutoff'):
				self._cutoff = kwargs.get('cutoff')       # this sets the cutoff. the default value I set is 0.
			if kwargs.get('screenradius'):
				self._cutoff = kwargs.get('screenradius')     
			else:
				self._cutoff = 0.0  
			self._screenradius=self._cutoff
		elif self._jmtype == 'TholeLin': #See Duijnen et al. Table 7 for these values
			#self._jmtype == 'TholeLinAniso' or self._jmtype == 'TholeLinIso' or
			if kwargs.get('screenradius'):
				self._screenradius=kwargs.get('screenradius') 
			elif self._fittype == 'mean':
				self._screenradius=1.6767 # Mean Fit
				if tholeparameters != 'FALSE':
					print 'Not set up mean Thole parameters yet! Only experimental.'
			elif self._fittype == 'components':
				self._screenradius=1.6623 # Component Fit
			elif self._fittype == 'empirical':
				self._screenradius=1.7278 # Empirical Fit
				#if tholeparameters != 'FALSE':
			elif not self._screenradius:
				print 'don\'t understand fittype and no screenradius defined'
				exit(2)
			#self._screenradius=1/self._screenradius
			print 'screenradius = ',self._screenradius				


		elif self._jmtype == 'TholeExp': #See Duijnen et al. Table 7 for these values
			print 'JM: TholeExp'
			self._screenradius=kwargs.get('screenradius', None) 
			if kwargs.get('screenradius'):
				self._screenradius=kwargs.get('screenradius') 
			elif self._fittype == 'mean':
				self._screenradius=2.4380 # Mean Fit
			elif self._fittype == 'components':
				self._screenradius=2.8990 # Component Fit
			elif self._fittype == 'empirical':
				print 'JM: TholeExp, empirical'
				self._screenradius=2.1304 # Empirical Fit
			elif self._fittype == 'extrascreen':
				self._screenradius=1.0 # Extra screen
			elif not self._screenradius:
				print 'don\'t understand fittype and no screenradius defined'
				exit(2)		
			#self._screenradius=1/self._screenradius
			print 'screenradius = ',self._screenradius

		else:
			print 'Didn\'t understand jmtype, should be Stern, TholeLin, or TholeExp'
			exit(2)

		print 'Generating JMatrix(', self._jmtype, ', screenradius = ', self._screenradius, ')'

		self.fill_matrix()

		#print 'J Matrix:'
		#print self._m



	def fill_matrix(self):
		""" this is automatically done when the JMatrix is initialised
		
		in principle you can call fill_matrix after the JMatrix has been initialised 
			
		"""
		row = 0
		for k, atom1 in GetRegister('Atom'):
			if np.linalg.det(atom1._pol)==0.0: # skip atoms that do not have polarizability
				print "atom does not have pol,(or pol is a singular matrix with no determinant), skipping  NB. This has caused problems so far! Python says makes Jmat a singular matrix..."
				continue	
			col = 0
			for k2, atom2 in GetRegister('Atom'):
				if np.linalg.det(atom2._pol)==0.0:
					print "atom does not have pol,skipping"
					continue
				self.fill_elements(row, col, atom1, atom2)
				col +=1
			row+=1
			
	
	def get_element(self,atom1, atom2):
		""" computes the submatrices that fit into the J_BB' matrix"""
		if atom1 is atom2 :
			""" the diagonal components J_BB are the inverse of the polarizability

			SHERIDAN: you can check this is true by looking up the definition on page 4731, second column just after equation 4

			"""
			return atom1._pol.I # Inverse matrix
		else: 
			""" the off-diagonal components J_BB' are described in  Stern equation 11. I think that there is an obvious typo in equaiton 11.
				Thole revisited equ 9
			"""
			dV = atom2._pos-atom1._pos
			d = dV.abs()
			I = np.mat(np.eye(3))          # this generates a3x3 identity matrix
			
			if self._jmtype == 'Stern':
				if d >= self._cutoff:
					return 1./d**3            * (I - 3. * (dV.T * dV )/d**2 )	#SHERIDAN: to understand how numpy does linear algebra
													# I suggest playing with it. the way i defined the positions
													# _pos.T * _pos gives you the "outer" product
													#_pos * _pos._T gives you the dot product
				else:
					return 1./self._cutoff**3 * (I - 3. * (dV.T * dV )/d**2 )

			elif self._jmtype == 'TholeExp':
				dV_vers=dV/d
				pol1_dV=np.linalg.norm(atom1._pol*dV_vers.T)
				# this is polarisability of atom 1 in direction connecting atoms 1 and 2. This will be just the iso value for atom 1 if polarisability is isotropic.
				pol2_dV=np.linalg.norm(atom2._pol*dV_vers.T)
				u=d/((pol1_dV*pol2_dV)**(1.0/6))
				#print 'u', u
				Tdiags= 1./d**3 * I * ( 1 - ( exp(-self._screenradius*u ) * (1 + self._screenradius * u + self._screenradius**2*u**2/2 ) ) ) 
				#Toffdiags= - 3. * (dV.T * dV )/d**5 * ( 1 - ( np.exp(-self._screenradius*d ) * (1 + self._screenradius * d + self._screenradius**2*d**2/2 + self._screenradius**3*d**3/6 ) ) )

				Toffdiags = - 3. * (dV.T * dV )*d**(-5)*(1 - (self._screenradius**3*u**3/6 + self._screenradius**2*u**2/2 + self._screenradius*u + 1)*np.exp(-self._screenradius*u) )

				#print 'Toffdiags'
				#print Toffdiags
				T = Tdiags + Toffdiags
				return T

			elif self._jmtype == 'TholeLin':
				#Uses polarisability in the bond direction for screening. This should usually be reasonable, as multiplied polarisabilities in perpendicular direction would need to be four times as large as those in bond directions to cause infinite polarisability perp. to bond direction. Usually bond direction is the one to worry about.
				dV_vers=dV/d
				pol1_dV=np.linalg.norm(atom1._pol*dV_vers.T)
				# this is polarisability of atom 1 in direction connecting atoms 1 and 2
				pol2_dV=np.linalg.norm(atom2._pol*dV_vers.T)
				# this is polarisability of atom 2 in direction connecting atoms 1 and 2
				s=self._screenradius*(pol1_dV*pol2_dV)**(1.0/6)

				v= d/s
				if d >= s:
						return 1./d**3            * (I - 3. * (dV.T * dV )/d**2 )
				else:
						return (4*v**3 - 3*v**4)/d**3 * I - 3. * (v**4/d**5) * (dV.T * dV )


			'''
			elif self._jmtype == 'TholeLinIso':
				s = self._screenradius*((atom1._pol[0,0]*atom2._pol[0,0])**(1.0/6))
				
				#print 'atom1._pol[0,0]', atom1._pol[0,0]
				#print 'atom2._pol[0,0]', atom2._pol[0,0]

				#print 'self._screenradius', self._screenradius

				#print '(atom1._pol[0,0]*atom2._pol[0,0])**(1.0/6)', (atom1._pol[0,0]*atom2._pol[0,0])**(1.0/6)

				v = d/s
				#print 'v = ',v,', d = ',d, ', s = ', s, '.'
				if d >= s:
						return 1./d**3            * (I - 3. * (dV.T * dV )/d**2 )
				else:
						#print 'MODIFIED INTERACTION'
						return (4*v**3 - 3*v**4)/d**3 * I - 3. * (v**4/d**5) * (dV.T * dV )
			'''

	'''
	OTHER IDEAS ABOUT GENERALISING THOLE TO ANISOTROPIC CASES:

				#Simplest option would be to use (Trace1/3 * Trace 2/3)^(1/6), but may not apture all subtleties
				# If dipoles sep. in cartesian axis (eg. x), and tensor is diagonal, would, I think make sense to take alpha_x for parallel part.
				# With anisotropic dipoles, it could also be the perpendicular part that blows up (requires perp polarisability ^ 1/6 over r to reach a higher value). Presumably 
				
		elif self._jmtype == 'TholeLin':
			polvec1=np.matrix([[atom1._pol[0,0],atom1._pol[1,1],atom1._pol[2,2]]])
			polvec2=np.matrix([[atom2._pol[0,0],atom2._pol[1,1],atom2._pol[2,2]]]) # these are the diagonal elements of the polarisability tensors of the two atoms.


			linpol1=
			linpol2= # As we have anisotropically polarisable dipoles, we are required to define our shape tensor in a different manner from Thole:
						# In fact, may make sense to modify polarisability in some directions and not others.
						# Most relevent would be to define an s for each part of the shape vector, where polarisabilities which are considered for each component of the tensor are those in which displacements are considered.
						# As there are two, eg. xzs in tensor, need to decide which x and which z we need to consider for each of these...
						# Makes sense to base this upon whether it is the component for influence of, eg., dip1_z on dip2_x or dip1_x on dip2_z. Can this be disentangled?
						# If aniso, certainly in general dip1_z.dip2_x doesn't equal dip1_x.dip2_z
						# My idea, from (2), and the above, if all below threshold:
			#T= np.multiply((dV.T * dV ),np.matrix([ [polvec1[0]*polvec2[0], polvec1[0]*polvec2[1], polvec1[0]*polvec2[2], 
			T= np.multiply((dV.T * dV ),(polvec1.T*polvec2))
						# However, to check element wise with this matrix whether below threshold:
			T=np.matrix(np.zeros((3,3)))
			for i in [0,1,2]:
				for j in [0,1,2]:
					s=self._screenradius*((polvec1.T[i]*polvec2.T[j]**(1/6)))
					v = d/s
					print 'Think about v definition for aniso. this may be fine, as absolute distance is taken in Thole, rather than distance in direction of component.'
					exit(2)
					if np.linalg.norm(dV[i])*np.linalg.norm(dV[j]) > s**2:
						if i == j:
							T[i,j]=1./d**3 - 3*dV[i]*dV[j]/d**5
						else:
							T[i,j]=-3*dV[i]*dV[j]/d**5
					else:
						if i == j:
							T[i,j]=(4*v**3 - 3*v**4)/d**3 - 3*(v**4)dV[i]*dV[j]/d**5
						else:
							T[i,j]=- 3*(v**4)dV[i]*dV[j]/d**5
			return T

		elif self._jmtype == 'TholeExp':
			T= 1./d**3            * ( I*(1 - (1+self._screenradius*d+(self._screenradius*d)**2/2)*exp(-self._screenradius*d))  - 3. * (dV.T * dV )/d**2 * (1+self._screenradius*d+(self._screenradius*d)**2/2+(self._screenradius*d)**3/6)*exp(-self._screenradius*d))
			print 'Think about terms, check this makes sense for aniso polarisability.'
			exit(2)
			return T
	'''

	
	def fill_elements(self,row, col, atom1, atom2):	
		"""this function computes the block element and adds it in at the correct row and col
			
		row and col refer to the position of the atom in the atom register. in order to convert
		this to an index for self._m, row and col must be multiplied by 3 - as each atom contributes
		three coordinates,x,y and z	
	
		"""
		mat = self.get_element(atom1, atom2) 
#		print mat
		for i in range(3):
			for j in range(3):
				self._m[row*3+i, col*3+j] = mat[i,j]

