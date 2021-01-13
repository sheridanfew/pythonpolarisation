from Molecule import Molecule
from MoleculeFactory import ReadMoleculeType
from MoleculeFactory import GetMolecule
from Rotation import Rotation as rot
from Position import Position
from ModifyPolarizability import *
#from RegisterMolecules import RegisterMolecules as _rm
from Register import RegisteredObject, GetRegister

import sys
import os

sys.path.append( os.path.dirname(__file__) + '/../Polarizability/')

from GetEnergyFromDips import *
from GetElectricField import *
from GetDipoles import *

import numpy as np

#eg. sexi_top=Crystal(mols_cen,mols_sur,TVs,cenpos,maxTVs,length[0],length[1],length[2])
#sexitopcrystal=Crystal(mols_cen=['sexi_cat_mola','sexi_neut_molb'],mols_sur=['sexi_neut_mola','sexi_neut_molb'],cenpos=[0,0,0],TVs,maxTVs,length[0],length[1],length[2],mols_outer=['sp_sexi_neut_mola','sp_sexi_neut_molb'],outer_maxTVs,outer_length[0],outer_length[1],outer_length[2],)


class Crystal ( RegisteredObject, list):
	"""A base crystal class
	"""

	def __init__(self, name,mols_cen,mols_sur,TVs,cenpos,length, **kwargs):
		list.__init__(self)
		self._type = "BaseClass"
		self._name = name
		self._cenpos =cenpos
		self._TVs =TVs

		self._TVangles = kwargs.get('TVangles')

		# Lists containing reorganisation energies from centrepos to different positions in the crystal. Different lists are cal. in different ways as in functions below.
		self._reorgs_staticpol=	[ [[[[] for i in range(length[2]+1)] for i in range(length[1]+1)] for i in range(length[0]+1)] for __ in xrange(len(mols_cen))]
		self._reorgs_shareq=	[ [[[[] for i in range(length[2]+1)] for i in range(length[1]+1)] for i in range(length[0]+1)] for __ in xrange(len(mols_cen))]
		self._reorgs_shareq_partial=	[ [[[{} for i in range(length[2]+1)] for i in range(length[1]+1)] for i in range(length[0]+1)] for __ in xrange(len(mols_cen))]
		#reorgs_shareq_partial[0][0][0][0][0.1] is the change in energy associated with moving 0.1 of the charge of the central mol to molecule at [0][0][0][0]

		#mols_outer only need be specified if want more coarse grained layer outside hi-res crystal (eg. single point polarisabilities for whole molecule rather than atomic res.)

		mols_outer = kwargs.get('mols_outer','none')

		self._maxTVs = kwargs.get('maxTVs','nomax')
		# maxTVs doesn't need to be specified, but allows the creation of a diamondlike shape. This will set a maximum number of TVs in any direction from cenpos. Otherwise, a cuboid will be created with lengths in a, b, c, as specified by max and min posns. (10000 assumed to be large enough that we will never reach this limit)
		outer_maxTVs = kwargs.get('outer_maxTVs',self._maxTVs)
		# if this is specified, this will allow upto maxTVs to be inner ring molecules, and above this outer molecules		


		outer_length = kwargs.get('outer_length',length)
		#Currently not implemented, use outer_maxTVs instead

		for i in np.arange(0,len(mols_cen),1):
				ReadMoleculeType(os.path.dirname(__file__) +  '/../Molecules/%s' % mols_cen[i])
				ReadMoleculeType(os.path.dirname(__file__) +  '/../Molecules/%s' % mols_sur[i])
				ReadMoleculeType(os.path.dirname(__file__) +  '/../Molecules/%s' % mols_outer[i])

		#print 'type(mols_cen)', type(mols_cen)
		self._mols=	[ [[[] for i in range(length[1]+1)] for i in range(length[0]+1)] for __ in xrange(len(mols_cen))]
		#print 'self._mols', self._mols
		#self._mols will contains all molecules.
		#mols[0] contains a list of all molecules in position a, mols[1] all mols in pos'n b, etc.
		#mols[0][x][y][z] contains molecule a in position x,y,z
		#mols may as such be iterated over in a number of ways to consider different molecules.
		#Need to make list the correct length in mols, a, b. We then append in c (see below)

		print 'Placing Molecules.'

		for a in np.arange(0,outer_length[0],1):
			for b in np.arange(0,outer_length[1],1):
				for c in np.arange(0,outer_length[2],1):
					if a == cenpos[0] and b == cenpos[1] and c == cenpos[2]:
						#central mol
						for i in np.arange(0,len(mols_cen),1):
							#print 'iabc', i, a, b, c
							#print 'Placing central Molecules:', mols_cen[i]
							self._mols[i][a][b].append(GetMolecule(os.path.dirname(__file__) +  '/../Molecules/%s' % mols_cen[i]))
							self._mols[i][a][b][c]().RedCoordMove(TVs=self._TVs, movevec=[a,b,c])
					#elif abs(a)-cenpos[0] <= (length[0]+1)/2 and (abs(b)-cenpos[1] <= (length[1]+1)/2 and abs(c)-cenpos[2] <= (length[2]+1)/2 and ( maxTVs == 'nomax' or (abs(a-cenpos[0]) + abs(b-cenpos[1]) + abs(c-cenpos[2])  <= maxTVs ):
					elif self._maxTVs == 'nomax' or (abs(a-cenpos[0]) + abs(b-cenpos[1]) + abs(c-cenpos[2]))  <= self._maxTVs:
						# Places molecule of this type if not in centres[0], and, if specified, less than MaxTVs translation vectors from centre unit cell, and still in fine resolution region.
						for i in np.arange(0,len(mols_sur),1):
							#print 'iabc', i, a, b, c
							#print 'Placing first ring Molecules:', mols_sur
							#print 'maxTVs', maxTVs
							self._mols[i][a][b].append(GetMolecule(os.path.dirname(__file__) +  '/../Molecules/%s' % mols_sur[i]))
							#print 'self._mols[i][a][b][c]', self._mols[i][a][b][c]
							#print 'type(self._mols[i][a][b][c])', type(self._mols[i][a][b][c])
							#print 'type(self._mols[i][a][b][c]())', type(self._mols[i][a][b][c]())
							#print 'self._TVs',self._TVs
							#print 'TVs', TVs
							self._mols[i][a][b][c]().RedCoordMove(TVs=self._TVs, movevec=[a,b,c])
					elif outer_maxTVs == 'nomax' or ( ( abs(a-cenpos[0]) + abs(b-cenpos[1]) + abs(c-cenpos[2]) ) <= outer_maxTVs):
						# Places molecule of this type if in coarse resolution region, and less than outer_TVs from centre
						for i in np.arange(0,len(mols_sur),1):
							#print 'iabc', i, a, b, c
							#print 'Placing second ring Molecules:', mols_sur
							#print 'outer_maxTVs', outer_maxTVs
							self._mols[i][a][b].append(GetMolecule(os.path.dirname(__file__) + '/../Molecules/%s' % mols_outer[i]))
							#print 'self._mols[i][a][b][c]', self._mols[i][a][b][c]
							#print 'type(self._mols[i][a][b][c])', type(self._mols[i][a][b][c])
							#print 'type(self._mols[i][a][b][c]())', type(self._mols[i][a][b][c]())
							self._mols[i][a][b][c]().RedCoordMove(TVs=self._TVs, movevec=[a,b,c])
					else:
						for i in np.arange(0,len(mols_sur),1):
							#print 'iabc', i, a, b, c
							#print 'empty'
							self._mols[i][a][b].append('empty')
		print 'Molecules placed.'
	def rotate(self, r):
		"""rotates the crystal by rotation r
		"""
		if not isinstance(r,rot):
			raise TypeError('must rotate molecules by a rotation')

		for a in range(0,len(self._mols[0])-1,1):
			for b in range(0,len(self._mols[0][0])-1,1):
				for c in range(0,len(self._mols[0][0][0])-1,1):
					for molincell in range(0,len(self._mols),1):
						if mols[molincell][a][b][c] != 'empty':
							self._mols[molincell][a][b][c].rotate(r)

		self._TVs=r*self._TVs
		print 'NB. Check rotation of TVs defined correctly)'

	def move(self, d):
		"""moves crystal by d


		"""
		for a in range(0,len(self._mols[0]-1),1):
			for b in range(0,len(self._mols[0][0]-1),1):
				for c in range(0,len(self._mols[0][0][0]-1),1):
					for molincell in range(0,len(self._mols),1):
						if mols[molincell][a][b][c] != 'empty':
							self._mols[molincell][a][b][c].move(d)		

	def print_posns(self):
		# new.dats for gnuplot with individual molecules
		if not os.path.exists('Dips_Posns_TVs'): os.makedirs('Dips_Posns_TVs')
		for a in range(0,len(self._mols[0])-1,1):
			for b in range(0,len(self._mols[0][0])-1,1):
				for c in range(0,len(self._mols[0][0][0])-1,1):
					for molincell in range(0,len(self._mols),1):
						if self._mols[molincell][a][b][c] != 'empty':
							f = open('Dips_Posns_TVs/%s_pos_%s%s%s_mol%s.dat' % (self._name,a,b,c,molincell), 'w')
							f.write('Molecule Starts:\n')
							for molecule in self._mols[molincell][a][b][c]():
								#print 'making posns file'
								molstr=str(molecule())
								for atom in GetRegister('Atom'):
									astr=str(atom)
									if molstr in astr:
										f.write(astr)
										f.write('\n')
							f.write('Molecule Ends.')
							f.flush()
							f.close()

	def ModifyPolarizabilityCry(self,**kwargs):
		print 'Modifying Polarisability...'
		jmtype=kwargs.get('jmtype', '')
		fittype=kwargs.get('fittype', '')
		for a in range(0,len(self._mols[0]),1):
			print 'a', a
			for b in range(0,len(self._mols[0][0]),1):
				print 'b', b
				for c in range(0,len(self._mols[0][0][0]),1):
					print 'c', c
					for molincell in range(0,len(self._mols),1):
						if self._mols[molincell][a][b][c] != 'empty':
							ModifyPolarizability.ModifyPolarizability(molecule=self._mols[molincell][a][b][c](),jmtype=jmtype,fittype=fittype)

	def calc_reorg_staticpol(self,a1,b1,c1,molincell1,a2,b2,c2,molincell2,dips,oldUqd):
		"""calc_reorg_staticpol(self,a1,b1,c1,molincell1,a2,b2,c2,molincell2,Efield,dips=dips._m,oldUqd):
		Calculates change in q-d interation energy for given set of E-fields and j matrix associated with moving charge of mol from a1,b1,c1,molincell1 to pos 2 whilst keeping polarisation of surroundings constant. Eqd, charge-dip interaction in initial configuration, must be specified.
		Recorded in list for mol2, mol1 is assumed to be central mol.
		"""
		print 'Attempting to swapping charges of molecule ', molincell1, 'at', a1,',', b1,',', c1,' and molecule ', molincell2, 'at', a2,',', b2,',', c2,'.'
		print 'cenpos', self._cenpos
		print 'maxTVs: ', self._maxTVs
		print 'nVs: ', (abs(a2-self._cenpos[0]) + abs(b2-self._cenpos[1]) + abs(c2-self._cenpos[2]))

		if a2 == a1 and b2 == b1 and c2 == c1 and molincell2 == molincell1:
			print 'Same Molecule'
			self._reorgs_staticpol[molincell2][a2][b2][c2]=0.0
			return self._reorgs_staticpol[molincell2][a2][b2][c2]
		elif self._reorgs_staticpol[molincell2][a2][b2][c2]:
			print 'Already Calculated'
			return self._reorgs_staticpol[molincell2][a2][b2][c2]
		elif self._mols[molincell2][a2][b2][c2] == 'empty':
			print 'No molecule in 2nd position.'
			return 'N/A'
		elif (abs(a2-self._cenpos[0]) + abs(b2-self._cenpos[1]) + abs(c2-self._cenpos[2]))  > self._maxTVs:
			print 'Outside of maxTVs for central moltype. Can only currently calculate roerg energies with molecules of the same type.'
		else:
			for cen_atom, x_atom in zip(self._mols[molincell1][a1][b1][c1](), self._mols[molincell2][a2][b2][c2]()):
				cen_initial=cen_atom()._crg
				x_initial=x_atom()._crg
				cen_atom()._crg=x_initial
				x_atom()._crg=cen_initial

			newUqd=np.multiply(get_U_qdip(dips=dips,Efield=get_electric_field(E0=np.matrix([0.,0.,0.]))),27.211)

			print 'Olduqd ',oldUqd,', Newuqd', newUqd,'.'

			self._reorgs_staticpol[molincell2][a2][b2][c2]=float(newUqd)-float(oldUqd)

			#Put charge back
			for cen_atom, x_atom in zip(self._mols[molincell1][a1][b1][c1](), self._mols[molincell2][a2][b2][c2]()):
				cen_initial=cen_atom()._crg
				x_initial=x_atom()._crg
				cen_atom()._crg=x_initial
				x_atom()._crg=cen_initial
			print 'Success.'
			return self._reorgs_staticpol[molincell2][a2][b2][c2]

	def print_reorgs_staticpol(self):
		lengths=[len(self._mols[0]),len(self._mols[0][0]),len(self._mols[0][0][0])+1]
		print 'lengths', lengths
		names=['a','b','c']

		if not os.path.exists('ReorgPlanesStaticPol'): os.makedirs('ReorgPlanesStaticPol')

		for view, col, row in zip(names,[1,2,0],[2,0,1]):
		# col, row here refer to position index in columns or rows
		# b\a
			thirdindex=[v for v in [0,1,2] if not ( v == row or v == col)][0]
			f = open('ReorgPlanesStaticPol/%s_reorgs_staticpol_planes_CENTRE_%s.dat' % (self._name,view), 'w')
			f.write(names[row] + '\\' + names[col] + '\tmol')
			for colpos in range(0,lengths[col]-1,1):
				f.write('\t' + str(colpos))
			for rowpos in range(0,lengths[row]-1,1):
					for molincell in range(0,len(self._mols),1):
						f.write('\n' + str(rowpos) + '\t' + str(molincell))
						for colpos in range(0,lengths[col]-1,1):
							pos=[0,0,0]
							pos[row]=rowpos
							pos[col]=colpos
							pos[thirdindex]=self._cenpos[thirdindex]
							f.write('\t' + str(self._reorgs_staticpol[molincell][pos[0]][pos[1]][pos[2]]))
			f.flush()
			f.close()

			# Same for non-central planes
			for planepos in range(0,lengths[thirdindex]-1,1):
				f = open('ReorgPlanesStaticPol/%s_reorgs_staticpol_planes_%s_%s.dat' % (self._name,planepos,view), 'w')
				f.write(names[row] + '\\' + names[col] + '\tmol')
				for colpos in range(0,lengths[col]-1,1):
					f.write('\t' + str(colpos))
				for rowpos in range(0,lengths[row]-1,1):
						for molincell in range(0,len(self._mols),1):
							f.write('\n' + str(rowpos) + '\t' + str(molincell))
							for colpos in range(0,lengths[col]-1,1):
								pos=[0,0,0]
								pos[row]=rowpos
								pos[col]=colpos
								pos[thirdindex]=planepos
								f.write('\t' + str(self._reorgs_staticpol[molincell][pos[0]][pos[1]][pos[2]]))

	def calc_reorg_shareq(self,a1,b1,c1,molincell1,a2,b2,c2,molincell2,jm,oldUqd):
		"""calc_reorg_shareq(self,a1,b1,c1,molincell1,a2,b2,c2,molincell2,Efield,dips=dips._m,oldUqd):
		Calculates change in q-d interation energy for given set of E-fields and j matrix associated with splitting charge of mol at a1,b1,c1,molincell1 equally between this site and pos 2, allowing dipoles of surroundings to relax. oldUqd, charge-dip interaction in initial configuration, must be specified.
		Recorded in list for mol2, mol1 is assumed to be central mol.
		"""
		print 'Attempting to swapping charges of molecule ', molincell1, 'at', a1,',', b1,',', c1,' and molecule ', molincell2, 'at', a2,',', b2,',', c2,'.'
		print 'cenpos', self._cenpos
		print 'maxTVs: ', self._maxTVs
		print 'nVs: ', (abs(a2-self._cenpos[0]) + abs(b2-self._cenpos[1]) + abs(c2-self._cenpos[2]))

		if a2 == a1 and b2 == b1 and c2 == c1 and molincell2 == molincell1:
			print 'Same Molecule'
			self._reorgs_shareq[molincell2][a2][b2][c2]=0.0
			return self._reorgs_shareq[molincell2][a2][b2][c2]
		elif self._reorgs_shareq[molincell2][a2][b2][c2]:
			print 'Already Calculated'
			return self._reorgs_shareq[molincell2][a2][b2][c2]
		elif self._mols[molincell2][a2][b2][c2] == 'empty':
			print 'No molecule in 2nd position.'
			return 'N/A'
		elif (abs(a2-self._cenpos[0]) + abs(b2-self._cenpos[1]) + abs(c2-self._cenpos[2]))  > self._maxTVs:
			print 'Outside of maxTVs for central moltype. Can only currently calculate roerg energies with molecules of the same type.'
		else:
			for cen_atom, x_atom in zip(self._mols[molincell1][a1][b1][c1](), self._mols[molincell2][a2][b2][c2]()):
				cen_initial=cen_atom()._crg
				x_initial=x_atom()._crg
				#print 'cen_initial', cen_initial
				#print 'x_initial', x_initial
				cen_atom()._crg=(x_initial+cen_initial)/2.0
				x_atom()._crg=(x_initial+cen_initial)/2.0
				#print 'cen_swapped: ',cen_atom()._crg
 				#print 'x_swapped: ',x_atom()._crg
				#print 'charges split:', cen_atom()._crg, cen_atom()._crg
			#print 'charges still split:', cen_atom()._crg, cen_atom()._crg
			Efield = get_electric_field(E0=np.matrix([0.,0.,0.]))

			d = get_dipoles(Efield=Efield,jm=jm)
			split_d = split_dipoles_onto_atoms(d)

			f = open('Dips_Posns_TVs/%s_dipoles_q%s%s%s_mol%s.dat' % (self._name,a2,b2,c2,molincell2), 'w')
			for dd in split_d:
				dstr=str(dd)
				f.write(dstr)
				f.write('\n')
			f.flush()
			f.close()

			newUqd=np.multiply(get_U_qdip(dips=d,Efield=Efield),27.211)
			#newUdd=np.multiply(get_U_dipdip(jm=jm,dips=d.T),27.211)
			print 'Olduqd ',oldUqd,', Newuqd', newUqd,'.'
			#print 'newUdd,', newUdd

			self._reorgs_shareq[molincell2][a2][b2][c2]=(float(newUqd)-float(oldUqd))/2.0
			# Making use of Udd = -0.5 Uqd -> delta U = 0.5 * delta Uqd

			#Put charge back
			for cen_atom, x_atom in zip(self._mols[molincell1][a1][b1][c1](), self._mols[molincell2][a2][b2][c2]()):
				cen_atom()._crg=cen_initial
				x_atom()._crg=x_initial
			#print 'Success.'
			#print 'charges back:', cen_atom()._crg, x_atom()._crg
			return self._reorgs_shareq[molincell2][a2][b2][c2]

	def calc_reorg_shareq_partial(self,a1,b1,c1,molincell1,a2,b2,c2,molincell2,jm,oldUqd,q):
		"""calc_reorg_shareq(self,a1,b1,c1,molincell1,a2,b2,c2,molincell2,Efield,dips=dips._m,oldUqd):
		Calculates change in q-d interation energy for given set of E-fields and j matrix associated with transferring a factor, q, of charge of mol at a1,b1,c1,molincell1 to pos 2, allowing dipoles of surroundings to relax. oldUqd, charge-dip interaction in initial configuration, must be specified.
		Recorded in list for mol2, mol1 is assumed to be central mol.
		"""
		print 'Attempting to put factor', q, ' of charge of molecule ', molincell1, 'at', a1,',', b1,',', c1,' on molecule ', molincell2, 'at', a2,',', b2,',', c2,'.'
		print 'cenpos', self._cenpos
		print 'maxTVs: ', self._maxTVs
		print 'nVs: ', (abs(a2-self._cenpos[0]) + abs(b2-self._cenpos[1]) + abs(c2-self._cenpos[2]))

		if a2 == a1 and b2 == b1 and c2 == c1 and molincell2 == molincell1:
			print 'Same Molecule'
			self._reorgs_shareq_partial[molincell2][a2][b2][c2][q]=0.0
			return self._reorgs_shareq_partial[molincell2][a2][b2][c2][q]
		elif q in self._reorgs_shareq_partial[molincell2][a2][b2][c2]:
			print 'Already Calculated'
			return self._reorgs_shareq_partial[molincell2][a2][b2][c2][q]
		elif self._mols[molincell2][a2][b2][c2] == 'empty':
			print 'No molecule in 2nd position.'
			self._reorgs_shareq_partial[molincell2][a2][b2][c2][q]='N/A'
			return 'N/A'
		elif (abs(a2-self._cenpos[0]) + abs(b2-self._cenpos[1]) + abs(c2-self._cenpos[2]))  > self._maxTVs:
			print 'Outside of maxTVs for central moltype. Can only currently calculate roerg energies with molecules of the same type.'
		else:
			for cen_atom, x_atom in zip(self._mols[molincell1][a1][b1][c1](), self._mols[molincell2][a2][b2][c2]()):
				cen_initial=cen_atom()._crg
				x_initial=x_atom()._crg
				print 'cen_initial', cen_initial
				print 'x_initial', x_initial
				cen_atom()._crg=(1.0-q)*(x_initial+cen_initial)
				x_atom()._crg=q*(x_initial+cen_initial)
				print 'cen_swapped: ',cen_atom()._crg
 				print 'x_swapped: ',x_atom()._crg
			print 'charges still split:', cen_atom()._crg, cen_atom()._crg
			Efield = get_electric_field(E0=np.matrix([0.,0.,0.]))

			d = get_dipoles(Efield=Efield,jm=jm)
			split_d = split_dipoles_onto_atoms(d)

			f = open('Dips_Posns_TVs/%s_dipoles_q%s%s%s_mol%s.dat' % (self._name,a2,b2,c2,molincell2), 'w')
			for dd in split_d:
				dstr=str(dd)
				f.write(dstr)
				f.write('\n')
			f.flush()
			f.close()

			newUqd=np.multiply(get_U_qdip(dips=d,Efield=Efield),27.211)
			#newUdd=np.multiply(get_U_dipdip(jm=jm,dips=d.T),27.211)
			print 'Olduqd ',oldUqd,', Newuqd', newUqd,'.'
			#print 'newUdd,', newUdd

			self._reorgs_shareq_partial[molincell2][a2][b2][c2][q]=(float(newUqd)-float(oldUqd))/2.0
			# Making use of Udd = -0.5 Uqd -> delta U = 0.5 * delta Uqd

			#Put charge back
			for cen_atom, x_atom in zip(self._mols[molincell1][a1][b1][c1](), self._mols[molincell2][a2][b2][c2]()):
				cen_atom()._crg=cen_initial
				x_atom()._crg=x_initial
			print 'Success.'
			print 'charges back:', cen_atom()._crg, x_atom()._crg
			return self._reorgs_shareq_partial[molincell2][a2][b2][c2][q]


	def print_reorgs_shareq(self):
		lengths=[len(self._mols[0]),len(self._mols[0][0]),len(self._mols[0][0][0])+1]
		print 'lengths', lengths
		names=['a','b','c']

		if not os.path.exists('ReorgPlanes_shareq'): os.makedirs('ReorgPlanes_shareq')

		for view, col, row in zip(names,[1,2,0],[2,0,1]):
		# col, row here refer to position index in columns or rows
		# b\a
			thirdindex=[v for v in [0,1,2] if not ( v == row or v == col)][0]
			f = open('ReorgPlanes_shareq/%s_reorgs_shareq_planes_CENTRE_%s.dat' % (self._name,view), 'w')
			f.write(names[row] + '\\' + names[col] + '\tmol')
			for colpos in range(0,lengths[col]-1,1):
				f.write('\t' + str(colpos))
			for rowpos in range(0,lengths[row]-1,1):
					for molincell in range(0,len(self._mols),1):
						f.write('\n' + str(rowpos) + '\t' + str(molincell))
						for colpos in range(0,lengths[col]-1,1):
							pos=[0,0,0]
							pos[row]=rowpos
							pos[col]=colpos
							pos[thirdindex]=self._cenpos[thirdindex]
							f.write('\t' + str(self._reorgs_shareq[molincell][pos[0]][pos[1]][pos[2]]))
			f.flush()
			f.close()

			# Same for non-central planes
			for planepos in range(0,lengths[thirdindex]-1,1):
				f = open('ReorgPlanes_shareq/%s_reorgs_shareq_planes_%s_%s.dat' % (self._name,planepos,view), 'w')
				f.write(names[row] + '\\' + names[col] + '\tmol')
				for colpos in range(0,lengths[col]-1,1):
					f.write('\t' + str(colpos))
				for rowpos in range(0,lengths[row]-1,1):
						for molincell in range(0,len(self._mols),1):
							f.write('\n' + str(rowpos) + '\t' + str(molincell))
							for colpos in range(0,lengths[col]-1,1):
								pos=[0,0,0]
								pos[row]=rowpos
								pos[col]=colpos
								pos[thirdindex]=planepos
								f.write('\t' + str(self._reorgs_shareq[molincell][pos[0]][pos[1]][pos[2]]))		


