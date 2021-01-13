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

from Molecules.pol_cyclopropane import pol
HFpol_cyclopropane=pol
from Molecules.pol_cyclohexane import pol
HFpol_cyclohexane=pol
from Molecules.pol_cyclopentane import pol
HFpol_cyclopentane=pol
from Molecules.pol_chain_dodecane import pol
HFpol_chain_dodecane=pol
from Molecules.pol_chain_ethane import pol
HFpol_chain_ethane=pol
from Molecules.pol_chain_hexane import pol
HFpol_chain_hexane=pol
from Molecules.pol_chain_methane import pol
HFpol_chain_methane=pol
from Molecules.pol_chain_propane import pol
HFpol_chain_propane=pol
from Molecules.pol_BENZENE_picen import pol
HFpol_BENZENE_picen=pol
from Molecules.pol_NAPHTA_picen import pol
HFpol_NAPHTA_picen=pol
from Molecules.pol_TETCEN_picen import pol
HFpol_TETCEN_picen=pol
from Molecules.pol_ANTCEN_picen import pol
HFpol_ANTCEN_picen=pol
from Molecules.pol_Pc_picen import pol
HFpol_Pc_picen=pol
from Molecules.pol_Benzene_MOPAC import pol
HFpol_Benzene_MOPAC=pol
from Molecules.pol_methane_MOPAC import pol
HFpol_methane_MOPAC=pol

def etatot( screen,C,H,single,pi,dpi,double):
	etalist=[tholefit('../Molecules/py_cyclopropane.xyz',HFpol_cyclopropane,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen),tholefit('../Molecules/py_cyclopentane.xyz',HFpol_cyclopentane,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen),tholefit('../Molecules//py_cyclohexane.xyz',HFpol_cyclohexane,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen),tholefit('../Molecules/chain_dodecane_fit_aniso_no_charge.xyz',HFpol_chain_dodecane,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen), tholefit('../Molecules/chain_ethane_fit_aniso_no_charge.xyz',HFpol_chain_ethane,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen),  tholefit('../Molecules/chain_hexane_fit_aniso_no_charge.xyz',HFpol_chain_hexane,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen), tholefit('../Molecules/chain_methane_fit_aniso_no_charge.xyz',HFpol_chain_methane,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen),  tholefit('../Molecules/chain_propane_fit_aniso_no_charge.xyz',HFpol_chain_propane,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen), tholefit('../Molecules/BENZENE_picen_fit_aniso_no_charge.xyz',HFpol_BENZENE_picen,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen), tholefit('../Molecules/py_Benzene_MOPAC.xyz',HFpol_Benzene_MOPAC,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen), tholefit('../Molecules/py_methane_MOPAC.xyz',HFpol_methane_MOPAC,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen)]

#, tholefit('../Molecules/NAPHTA_picen_fit_aniso_no_charge.xyz',HFpol_NAPHTA_picen,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen), tholefit('../Molecules/TETCEN_picen_fit_aniso_no_charge.xyz',HFpol_TETCEN_picen,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen), tholefit('../Molecules/ANTCEN_picen_fit_aniso_no_charge.xyz',HFpol_ANTCEN_picen,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen), tholefit('../Molecules/Pc_picen_fit_aniso_no_charge.xyz',HFpol_Pc_picen,C=C,H=H,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen)

	etatot=0
	for i in np.arange(0,len(etalist),1):
		etatot=etatot+etalist[i]
	print 'etatot'
	print etatot
	return etatot


#etatot(screen=1.68384957924, H=1.3849, C=6.5780, single=0.000001, double=0.000001,pi=0.000001,dpi=0.0001) # Thole Exp Components

#etatot(screen=2.8990, H=1.3849, C=6.5780, single=0.000001, double=0.000001,pi=0.000001,dpi=0.0001) # Thole Exp Components

etatot(screen=2.1304, H=2.7927, C=8.6959, single=0.000001, double=0.000001,pi=0.000001,dpi=0.000001) # Thole Exp Empirical

#etatot(screen=1.7278, H=3.5020, C=10.1756, single=0.0001, double=0.0001,pi=0.0001,dpi=0.0001) # Thole Lin Empirical

#etatot(screen=1.6623, H=0.5974, C=7.7764, single=0.0001, double=0.0001,pi=0.0001,dpi=0.0001) # Thole Lin Components


#'H':2.7927, 'C':8.6959, 'N': 6.5565, 'O': 5.7494, 'F': 3.0013, 'S': 16.6984, 'Cl':16.1979, 'Br': 23.5714, 'I':36.9880, 'screen':2.1304}\n"

#m=minuit.Minuit(etatot, H=0.001,fix_H=True,limit_single=(0.0100,20.0),limit_pi=(0.0100,20.0),limit_dpi=(0.0100,20.0),limit_S=(0.0100,20.0), screen=2.8990, fix_screen=True)
#m.scan(("C",5,0.0100,20.0),("single",5,0.0100,20.0),("pi",5,0.0100,20.0),("dpi",5,0.0100,20.0),("S",5,0.0100,20.0),)
#m.printMode = 1
#m.migrad()
