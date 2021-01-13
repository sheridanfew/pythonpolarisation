import sys
sys.path.append('../')

from BasicElements import *
from BasicElements.MoleculeFactory import *
#from BasicElements.ModifyPolarizability import *
#Seems not to work when imported, don't really understand why.
from Polarizability import *
from Polarizability.GetDipoles import *

import numpy as np
from BasicElements.Register import GetRegister


def eta( namefile, calpol, **kwargs, cut):
    
	ReadMoleculeType(namefile)
	mol = GetMolecule(namefile)
    
	def ModifyPolarizability(molecule, C, N, H, S):
		""" takes a mol and changes the isotropic polarizability
		for each atom
		"""
		C = Polarizability(iso=C)
		H = Polarizability(iso=H)
		N = Polarizability(iso=N)
		S = Polarizability(iso=S)
		for atom in molecule:
			if (atom()._elname=="C"):
				atom()._pol= C
			elif (atom()._elname=="H"):
				atom()._pol= H
			elif (atom()._elname=="N"):
				atom()._pol= N
			elif (atom()._elname=="S"):
				atom()._pol= S

	pypol=np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
	etamat=np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])    

	for i in np.arange(0. ,2.1 ,1. ):
		E0 = np.matrix([0.,0.,0.])
		E0[0,i]=1.
		ModifyPolarizability(mol(), C=c,S=s,H=h,N=n )
		d = get_dipoles(E0=E0,cutoff=cut)
		split_d = split_dipoles_onto_atoms(d)
		tot = np.matrix([0.,0.,0.])
       		for dd in split_d:
        		tot += dd
		print 'tot'
		print tot
       		pypol.T[i] = tot
    
	etamat=np.multiply((pypol-calpol),(pypol-calpol))/np.multiply(calpol,calpol)
    
	for i in np.arange(0,3,1):
		for j in np.arange(0,3,1):
                	if calpol[i,j] == 0:
                    		etamat[i,j] = '0'

	eta=0
	eta = etamat[0,0]+etamat[1,0]+etamat[1,1]+etamat[2,0]+etamat[2,1]+etamat[2,2]
	print 'c, h, s, cut'
	print c, h, s, cut
	print 'calpol:'
	print calpol
	print 'pypol:'
	print pypol
	print 'etamat:'
	print etamat
	print 'eta:'
	print eta
	mol().__del__
	return eta

def abs_square_dif_fit( namefile, calpol, c, h, n, s, cut):
    
	ReadMoleculeType(namefile)
	mol = GetMolecule(namefile)
    
	def ModifyPolarizability(molecule, C, N, H, S):
		""" takes a mol and changes the isotropic polarizability
		for each atom
		"""
		C = Polarizability(iso=C)
		H = Polarizability(iso=H)
		N = Polarizability(iso=N)
		S = Polarizability(iso=S)
		for atom in molecule:
			if (atom()._elname=="C"):
				atom()._pol= C
			elif (atom()._elname=="H"):
				atom()._pol= H
			elif (atom()._elname=="N"):
				atom()._pol= N
			elif (atom()._elname=="S"):
				atom()._pol= S

	pypol=np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
	etamat=np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])    

	for i in np.arange(0. ,2.1 ,1. ):
		E0 = np.matrix([0.,0.,0.])
		E0[0,i]=1.
		ModifyPolarizability(mol(), C=c,S=s,H=h,N=n )
		d = get_dipoles(E0=E0,cutoff=cut)
		split_d = split_dipoles_onto_atoms(d)
		tot = np.matrix([0.,0.,0.])
       		for dd in split_d:
        		tot += dd
		print 'tot'
		print tot
       		pypol.T[i] = tot
    
	etamat=np.multiply((pypol-calpol),(pypol-calpol))
    
	for i in np.arange(0,3,1):
		for j in np.arange(0,3,1):
                	if calpol[i,j] == 0:
                    		etamat[i,j] = '0'

	eta=0
	eta = etamat[0,0]+etamat[1,0]+etamat[1,1]+etamat[2,0]+etamat[2,1]+etamat[2,2]
	print 'c, h, s, cut'
	print c, h, s, cut
	print 'calpol:'
	print calpol
	print 'pypol:'
	print pypol
	print 'etamat:'
	print etamat
	print 'absdiffit:'
	print eta
	mol().__del__
	return eta
