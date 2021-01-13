import sys
sys.path.append('../')
from BasicElements import *
from BasicElements.Register import GetRegister
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from BasicElements.Crystal import *
from BasicElements.ModifyPolarizability import *
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from Polarizability import *
from Polarizability.GetEnergyFromDips import *
from Polarizability.JMatrix import JMatrix
import numpy as np
from math import *
from time import gmtime, strftime
import os


basename='nTs_TholeTests'
f = open('Polarisabilities_nT_long.csv', 'w')
f.write(basename)
f.write('\nLength\tLong_calpol\tTholeExp_empirical\tTholeExp_mean\tTholeExp_components\tTholeLin_empirical\tTholeLin_mean\tTholeLin_components')

g = open('Polarisabilities_nT_trace.csv', 'w')
g.write(basename)
g.write('\nLength\tTrace_calpol\tTholeExp_empirical\tTholeExp_mean\tTholeExp_components\tTholeLin_empirical\tTholeLin_mean\tTholeLin_components')

h = open('Polarisabilities_nT_long_ratios.csv', 'w')
h.write(basename)
h.write('\nLength\tLong_calpol\tTholeExp_empirical_ratio\tTholeExp_mean_ratio\tTholeExp_components_ratio\tTholeLin_empirical_ratio\tTholeLin_mean_ratio\tTholeLin_components_ratio')

j = open('Polarisabilities_nT_trace_ratios.csv', 'w')
j.write(basename)
j.write('\nLength\tTrace_calpol\tTholeExp_empirical_ratio\tTholeExp_mean_ratio\tTholeExp_components_ratio\tTholeLin_empirical_ratio\tTholeLin_mean_ratio\tTholeLin_components_ratio')


print strftime("%a, %d %b %Y %X +0000", gmtime())
for n in range(1,9,1):
	name= str( 'thio_' + str(n) + 'T_E_100' )
	print 'N', n, ' Namefile: ', name
	namefile= str( 'thio_' + str(n) + 'T_neut_fit_aniso_no_charge.xyz' )
	ReadMoleculeType('../Molecules/' + namefile)
	mol = GetMolecule('../Molecules/' + namefile)
	exec( 'from Molecules.pol_thio_'  + str(n) + 'T_neut import pol as calpol' )
	calpoldiag=list(np.diag(calpol))

	calpol_max = max(calpoldiag)
	calpol_max_index = calpoldiag.index(calpol_max)
	calpol_tr = np.trace(calpol)

	f.write('\n%s\t%s' % (n,str(calpol_max))) 
	h.write('\n%s\t%s' % (n,str(calpol_max)))

	g.write('\n%s\t%s' % (n,str(calpol_tr)))
	j.write('\n%s\t%s' % (n,str(calpol_tr)))


	for jmtype in [ 'TholeExp', 'TholeLin' ]:
		for fittype in ['empirical', 'mean', 'components']:
			ModifyPolarizability(molecule=mol(),jmtype=jmtype,fittype=fittype)
			print 'mol()'			
			print mol()
			jm=JMatrix(jmtype=jmtype,fittype=fittype)

			pypol=np.matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])

			for i in np.arange(0. ,2.1 ,1. ):
					 E0 = np.matrix([0.,0.,0.])
					 E0[0,i]=1.
					 d = get_dipoles(E0=E0,jm=jm._m)
					 split_d = split_dipoles_onto_atoms(d)
					 tot = np.matrix([0.,0.,0.])
					 for dd in split_d:
						 tot += dd
					 print 'tot'
					 print tot
					 pypol.T[i] = tot

			pypoldiag=list(np.diag(pypol))
			pypol_tr=np.trace(pypol)
			ratios=np.divide(pypol,calpol)

			f.write('\t%s' % (str(pypoldiag[calpol_max_index]))) 
			h.write('\t%s' % (str(ratios[calpol_max_index,calpol_max_index])))

			g.write('\t%s' % (str(np.trace(pypol))))
			j.write('\t%s' % (str(pypol_tr/calpol_tr)))

f.flush()
f.close()
g.flush()
g.close()
h.flush()
h.close()
j.flush()
j.close()


'''
			# print dipoles
			if not os.path.exists('Dips_Posns_TVs'): os.makedirs('Dips_Posns_TVs')
			f = open('Dips_Posns_TVs/%s_dipoles.dat' % name, 'w')
			for dd in split_d:
				dstr=str(dd)
				f.write(dstr)
				f.write('\n')
			f.flush()
			f.close()

			f = open('Dips_Posns_TVs/%s_posns.dat' % name, 'w')
			f.write('Molecule Starts:\n')
			for atom in GetRegister('Atom'):
				astr=str(atom)
				f.write(astr)
				f.write('\n')
			f.write('Molecule Ends.')
			f.flush()
			f.close()

print 'Job Completed Successfully.
'''
