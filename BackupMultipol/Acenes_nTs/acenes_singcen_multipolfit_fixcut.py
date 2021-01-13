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

from Molecules.pol_BENZENE_singcen_fit_aniso import pol
HFpol_BENZENE_singcen_fit_aniso=pol
from Molecules.pol_NAPHTA_singcen_fit_aniso import pol
HFpol_NAPHTA_singcen_fit_aniso=pol
from Molecules.pol_TETCEN_singcen_fit_aniso import pol
HFpol_TETCEN_singcen_fit_aniso=pol
from Molecules.pol_ANTCEN_singcen_fit_aniso import pol
HFpol_ANTCEN_singcen_fit_aniso=pol
from Molecules.pol_Pc_singcen_fit_aniso import pol
HFpol_Pc_singcen_fit_aniso=pol
def etatot( cut,C,H,single,pi):
	etalist=[eta('../Molecules/BENZENE_singcen_fit_aniso_no_charge.xyz',HFpol_BENZENE_singcen_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/NAPHTA_singcen_fit_aniso_no_charge.xyz',HFpol_NAPHTA_singcen_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/TETCEN_singcen_fit_aniso_no_charge.xyz',HFpol_TETCEN_singcen_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/ANTCEN_singcen_fit_aniso_no_charge.xyz',HFpol_ANTCEN_singcen_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/Pc_singcen_fit_aniso_no_charge.xyz',HFpol_Pc_singcen_fit_aniso,cut,C=C,H=H,single=single,pi=pi)]
	etatot=0
	for i in np.arange(0,len(etalist),1):
		etatot=etatot+etalist[i]
	print 'etatot'
	print etatot
	return etatot
