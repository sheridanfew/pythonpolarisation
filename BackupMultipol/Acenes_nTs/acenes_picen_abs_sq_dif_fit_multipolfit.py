import sys
sys.path.append('../')
from BasicElements import Molecule,Position, Polarizability
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from BasicElements import *
from Polarizability.etacalc import *
import numpy as np
#from scipy.optimize import minimize
import minuit
import time

from Molecules.pol_BENZENE_picen_fit_aniso import pol
HFpol_BENZENE_picen_fit_aniso=pol
from Molecules.pol_NAPHTA_picen_fit_aniso import pol
HFpol_NAPHTA_picen_fit_aniso=pol
from Molecules.pol_TETCEN_picen_fit_aniso import pol
HFpol_TETCEN_picen_fit_aniso=pol
from Molecules.pol_ANTCEN_picen_fit_aniso import pol
HFpol_ANTCEN_picen_fit_aniso=pol
from Molecules.pol_Pc_picen_fit_aniso import pol
HFpol_Pc_picen_fit_aniso=pol
def etatot( cut,C,H,single,pi):
	etalist=[abs_square_dif_fit('../Molecules/BENZENE_picen_fit_aniso_no_charge.xyz',HFpol_BENZENE_picen_fit_aniso,cut,C=C,H=H,single=single,pi=pi), abs_square_dif_fit('../Molecules/NAPHTA_picen_fit_aniso_no_charge.xyz',HFpol_NAPHTA_picen_fit_aniso,cut,C=C,H=H,single=single,pi=pi), abs_square_dif_fit('../Molecules/TETCEN_picen_fit_aniso_no_charge.xyz',HFpol_TETCEN_picen_fit_aniso,cut,C=C,H=H,single=single,pi=pi), abs_square_dif_fit('../Molecules/ANTCEN_picen_fit_aniso_no_charge.xyz',HFpol_ANTCEN_picen_fit_aniso,cut,C=C,H=H,single=single,pi=pi), abs_square_dif_fit('../Molecules/Pc_picen_fit_aniso_no_charge.xyz',HFpol_Pc_picen_fit_aniso,cut,C=C,H=H,single=single,pi=pi)]
	etatot=0
	for i in np.arange(0,len(etalist),1):
		etatot=etatot+etalist[i]
	print 'etatot'
	print etatot
	return etatot
m=minuit.Minuit(etatot, cut=8.0, fix_cut=True ,C=10.0,limit_C=(5.0100,15.0),H=0.001,fix_H=True,single=5.0,limit_single=(0.0100,10.0),pi=15.0,limit_pi=(5.0100,25.0))
m.printMode = 1
m.migrad()

