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

from Molecules.pol_thio_1T_neut import pol
HFpol_thio_1T_neut_aniso=pol
from Molecules.pol_thio_2T_neut import pol
HFpol_thio_2T_neut_aniso=pol
from Molecules.pol_thio_3T_neut import pol
HFpol_thio_3T_neut_aniso=pol
from Molecules.pol_thio_4T_neut import pol
HFpol_thio_4T_neut_aniso=pol
from Molecules.pol_thio_5T_neut import pol
HFpol_thio_5T_neut_aniso=pol
from Molecules.pol_thio_6T_neut import pol
HFpol_thio_6T_neut_aniso=pol
from Molecules.pol_thio_7T_neut import pol
HFpol_thio_7T_neut_aniso=pol
from Molecules.pol_thio_8T_neut import pol
HFpol_thio_8T_neut_aniso=pol
def abs_sq_dif_fittot( cut,C,H,S,single,pi):
	abs_sq_dif_fitlist=[abs_square_dif_fit('../Molecules/thio_1T_neut_aniso_no_charge.xyz',HFpol_thio_1T_neut_aniso,cut=cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_2T_neut_aniso_no_charge.xyz',HFpol_thio_2T_neut_aniso,cut=cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_3T_neut_aniso_no_charge.xyz',HFpol_thio_3T_neut_aniso,cut=cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_4T_neut_aniso_no_charge.xyz',HFpol_thio_4T_neut_aniso,cut=cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_5T_neut_aniso_no_charge.xyz',HFpol_thio_5T_neut_aniso,cut=cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_6T_neut_aniso_no_charge.xyz',HFpol_thio_6T_neut_aniso,cut=cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_7T_neut_aniso_no_charge.xyz',HFpol_thio_7T_neut_aniso,cut=cut,C=C,H=H,S=S,single=single,pi=pi), abs_square_dif_fit('../Molecules/thio_8T_neut_aniso_no_charge.xyz',HFpol_thio_8T_neut_aniso,cut=cut,C=C,H=H,S=S,single=single,pi=pi)]
	abs_sq_dif_fittot=0
	for i in np.arange(0,len(abs_sq_dif_fitlist),1):
		abs_sq_dif_fittot=abs_sq_dif_fittot+abs_sq_dif_fitlist[i]
	print 'abs_sq_dif_fittot'
	print abs_sq_dif_fittot
	return abs_sq_dif_fittot
m=minuit.Minuit(abs_sq_dif_fittot, cut=8.0, fix_cut=True ,C=12.1626349101,fix_C=True,H=0.001,fix_H=True,S=40.0,limit_S=(25.0100,80.0),single=0.20929801812,fix_single=True,pi=8.02055851811,fix_pi=True)
m.printMode = 1
m.migrad()

