import sys
sys.path.append('../')

from BasicElements import *
from BasicElements.MoleculeFactory import *
#from BasicElements.ModifyPolarizability import *
#Seems not to work when imported, don't really understand why.
from Polarizability import *
from Polarizability.GetDipoles import *
from Polarizability.JMatrix import JMatrix

import numpy as np
from BasicElements.Register import GetRegister


def eta( namefile, calpol, cut, **kwargs):
	""" namefile is path to file containing molecule specifications, calpol is molecular polarisability we are calibrating to, cut is cutoff for interation of polarisable dipoles (might be good to make dif for dif atom types at some point) """
	""" kwargs specify desired atomic polarisabilities, will change polarisability of each atom of this type, eg. C=3.0 """
	""" kwargs also specify costslope (universal, set to 0 if no cost desired) and cost centre (for each atom type), specified as eg. C_costcentre=3.5 """
    	'''Fit parameter here is difference in squared differences of each component of polarisability tensor normalised by square of calibrated value'''
	ReadMoleculeType(namefile)
	mol = GetMolecule(namefile)
    

	print "keywords understood:"
	for key in kwargs:
       		print "%s: %s" % (key, kwargs[key])

	""" For cost function, will scale fitting parameter by costslope*(pol-costcentre)"""
	costslope = kwargs.get('costslope',0.0)

	""" ModifyPolarizability """

	print namefile

#New:
	for atom in mol():
      		for key in kwargs:
			if (atom()._elname.upper()==key.upper()):
				b=np.matrix(atom()._pol)
				mag=np.sqrt(b[0,0]*b[0,0] + b[1,1]*b[1,1] + b[2,2]*b[2,2])
				newpol=np.divide(kwargs[key]*b,mag)
#				print 'newpol', newpol
				atom()._pol=newpol

	for atom in mol():
      		for key in kwargs:
			if (atom()._elname.upper()==key.upper()):
				b=np.matrix(atom()._pol)
				mag=np.sqrt(b[0,0]*b[0,0] + b[1,1]*b[1,1] + b[2,2]*b[2,2])
				newpol=np.divide(kwargs[key]*b,mag)
				atom()._pol=newpol


	costfactor=1.0
      	for costkey in kwargs:
		for atomkey in kwargs:
			if (costkey==atomkey,'_costcentre'):
#				print 'type: costslope, costfactor, kwargs[costkey], kwargs[atomkey]'
#				t_cs=type(costslope)
#				t_cf=type(costfactor)
#				t_kck=type(kwargs[costkey])
#				t_kak=type(kwargs[atomkey])
#				print t_cs, t_cf, t_kck, t_kak

				costfactor = costfactor * (1.0 + (costslope*((kwargs[costkey]-kwargs[atomkey]) ** 2.0)))

	print 'costfactor', costfactor

	



#	OLD for iso
#	
#	I = np.mat(np.eye(3))
#
#	for atom in mol():
#      		for key in kwargs:
#			if (atom()._elname.upper()==key.upper()):
#				atom()._pol=kwargs[key]*I



#	print mol()
#	for atom in mol():
#		print atom()
#		print atom()._pol

	
	'''pypol will give polarisability in python model, etamat will give goodness of fit value for each component of the polarisability tensor '''
	pypol=np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
	etamat=np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])    

	jm=JMatrix(jmtype='Stern',cutoff=cut)

	for i in np.arange(0. ,2.1 ,1. ):
		E0 = np.matrix([0.,0.,0.])
		E0[0,i]=1.
#		for atom in mol():
#			print atom()
#			print atom()._pol
		d = get_dipoles(E0=E0,jm=jm._m)
		split_d = split_dipoles_onto_atoms(d)
		tot = np.matrix([0.,0.,0.])
       		for dd in split_d:
        		tot += dd
		print 'tot'
		print tot
       		pypol.T[i] = tot
    
	etamat=np.multiply((pypol-calpol),(pypol-calpol))
    
#	If calibration polarisability is 0, eta not appropriate fit here (and previous funct will have made undefined), set value in etamat to 0
	for i in np.arange(0,3,1):
		for j in np.arange(0,3,1):
                	if calpol[i,j] == 0:
                    		etamat[i,j] = '0'

	''' eta sums components of etamat to give overall goodness of fit '''
	eta=0
	eta = etamat[0,0]+etamat[1,0]+etamat[1,1]+etamat[2,0]+etamat[2,1]+etamat[2,2]
	eta_w_cost = costfactor*eta

	print 'eta_w_cost', eta_w_cost

#	print mol()
#	print 'c, h, s, cut'
#	print c, h, s, cut
	print 'calpol:'
	print calpol
	print 'pypol:'
	print pypol
#	print 'etamat:'
#	print etamat
	py_over_cal=pypol/calpol
	print 'py_over_cal'
	print py_over_cal
	print 'py over cal three linear compnents:'
	print namefile, py_over_cal[0,0], py_over_cal[1,1], py_over_cal[2,2], 
	print 'eta:'
	print eta
#	t_e=type(eta)
#	print 'type eta:', t_e
#	t_ewc=type(eta_w_cost)
#	print 'type eta_w_cost:', t_ewc
	print 'eta_w_cost:'
	print eta_w_cost
	mol().__del__
	return eta_w_cost

def abs_square_dif_fit( namefile, calpol, cut, **kwargs):
    	'''similar to eta function, for explanation of terms see eta. Fit parameter here is difference square of difference between each component of python and calibration polarisability tensors, NOT normalised by square of calibrated value => favours getting polarisability right in highly polarisable directions, possibly at the cost of less polarisable directions.'''
	ReadMoleculeType(namefile)
	mol = GetMolecule(namefile)
    
	""" ModifyPolarizability """

	print namefile

	print "keywords understood:"
	for key in kwargs:
       		print "%s: %s" % (key, kwargs[key])
	

	for atom in mol():
      		for key in kwargs:
			if (atom()._elname.upper()==key.upper()):
				b=np.matrix(atom()._pol)
				mag=np.sqrt(b[0,0]*b[0,0] + b[1,1]*b[1,1] + b[2,2]*b[2,2])
				newpol=np.divide(kwargs[key]*b,mag)
				atom()._pol=newpol

	pypol=np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
	etamat=np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])    

	jm=JMatrix(jmtype='Stern',cutoff=cut)

	for i in np.arange(0. ,2.1 ,1. ):
		E0 = np.matrix([0.,0.,0.])
		E0[0,i]=1.
#		for atom in mol():
#			print atom()
#			print atom()._pol
		d = get_dipoles(E0=E0,jm=jm._m)
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
	t_e=type(eta)
	print 'type eta:', t_e
#	print 'c, h, s, cut'
#	print c, h, s, cut
	print 'calpol:'
	print calpol
	print 'pypol:'
	print pypol
#	print 'etamat:'
#	print etamat
	py_over_cal=pypol/calpol
	print 'py_over_cal'
	print py_over_cal
	print 'absdiffit:'
	print eta
	mol().__del__
	return eta
