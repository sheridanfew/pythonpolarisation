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
def etatot( cut,C,H,S,single,pi):
	etalist=[eta('../Molecules/thio_1T_neut_fit_aniso_no_charge.xyz',HFpol_thio_1T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/thio_2T_neut_fit_aniso_no_charge.xyz',HFpol_thio_2T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/thio_3T_neut_fit_aniso_no_charge.xyz',HFpol_thio_3T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/thio_4T_neut_fit_aniso_no_charge.xyz',HFpol_thio_4T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/thio_5T_neut_fit_aniso_no_charge.xyz',HFpol_thio_5T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/thio_6T_neut_fit_aniso_no_charge.xyz',HFpol_thio_6T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/thio_7T_neut_fit_aniso_no_charge.xyz',HFpol_thio_7T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/thio_8T_neut_fit_aniso_no_charge.xyz',HFpol_thio_8T_neut_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/BENZENE_picen_fit_aniso_no_charge.xyz',HFpol_BENZENE_picen_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/NAPHTA_picen_fit_aniso_no_charge.xyz',HFpol_NAPHTA_picen_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/TETCEN_picen_fit_aniso_no_charge.xyz',HFpol_TETCEN_picen_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/ANTCEN_picen_fit_aniso_no_charge.xyz',HFpol_ANTCEN_picen_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi), eta('../Molecules/Pc_picen_fit_aniso_no_charge.xyz',HFpol_Pc_picen_fit_aniso,cut,C=C,H=H,S=S,single=single,pi=pi)]
	etatot=0
	for i in np.arange(0,len(etalist),1):
		etatot=etatot+etalist[i]
	print 'etatot'
	print etatot
	return etatot
