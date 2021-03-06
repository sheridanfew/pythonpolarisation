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

from Molecules.pol_thio_1T_neut_fit_aniso import pol
HFpol_thio_1T_neut_fit_aniso=pol
from Molecules.pol_thio_2T_neut_fit_aniso import pol
HFpol_thio_2T_neut_fit_aniso=pol
from Molecules.pol_thio_3T_neut_fit_aniso import pol
HFpol_thio_3T_neut_fit_aniso=pol
from Molecules.pol_thio_4T_neut_fit_aniso import pol
HFpol_thio_4T_neut_fit_aniso=pol
from Molecules.pol_thio_5T_neut_fit_aniso import pol
HFpol_thio_5T_neut_fit_aniso=pol
from Molecules.pol_thio_6T_neut_fit_aniso import pol
HFpol_thio_6T_neut_fit_aniso=pol
from Molecules.pol_thio_7T_neut_fit_aniso import pol
HFpol_thio_7T_neut_fit_aniso=pol
from Molecules.pol_thio_8T_neut_fit_aniso import pol
HFpol_thio_8T_neut_fit_aniso=pol
def etatot( cut,C,H,S,single,pi):
	etalist=[abs_square_dif_fit('../Molecules/thio_1T_neut_fit_aniso_no_charge.xyz',HFpol_thio_1T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_2T_neut_fit_aniso_no_charge.xyz',HFpol_thio_2T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_3T_neut_fit_aniso_no_charge.xyz',HFpol_thio_3T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_4T_neut_fit_aniso_no_charge.xyz',HFpol_thio_4T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_5T_neut_fit_aniso_no_charge.xyz',HFpol_thio_5T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_6T_neut_fit_aniso_no_charge.xyz',HFpol_thio_6T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_7T_neut_fit_aniso_no_charge.xyz',HFpol_thio_7T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_8T_neut_fit_aniso_no_charge.xyz',HFpol_thio_8T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi)]
	etatot=0
	for i in np.arange(0,len(etalist),1):
		etatot=etatot+etalist[i]
	print 'etatot'
	print etatot
	return etatot
m=minuit.Minuit(etatot, cut=8.0, fix_cut=True ,C=10.0,limit_C=(5.0100,15.0),H=0.001,fix_H=True,S=15.0,limit_S=(10.0100,20.0),single=5.0,limit_single=(0.0100,10.0),pi=15.0,limit_pi=(5.0100,25.0))
m.printMode = 1
m.migrad()

