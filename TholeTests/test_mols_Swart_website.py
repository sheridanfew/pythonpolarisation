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

from Molecules.pol_Water_SWART import pol
HFpol_Water_SWART=pol
from Molecules.pol_CO_SWART import pol
HFpol_CO_SWART=pol

def etatot( screen,C,H,O,single,pi,dpi,double):
	etalist=[tholefit('../Molecules/py_Water_SWART.xyz',HFpol_Water_SWART,C=C,H=H,O=O,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen),tholefit('../Molecules/py_CO_SWART.xyz',HFpol_CO_SWART,C=C,H=H,O=O,single=single,pi=pi,dpi=dpi,jmtype='TholeExp',screenradius=screen)]

#, tholefit('../Molecules/NAPHTA_picen_fit_aniso_no_charge.xyz',HFpol_NAPHTA_picen,C=C,H=H,O=O,single=single,pi=pi,dpi=dpi,jmtype='TholeLinIso',screenradius=screen), tholefit('../Molecules/TETCEN_picen_fit_aniso_no_charge.xyz',HFpol_TETCEN_picen,C=C,H=H,O=O,single=single,pi=pi,dpi=dpi,jmtype='TholeLinIso',screenradius=screen), tholefit('../Molecules/ANTCEN_picen_fit_aniso_no_charge.xyz',HFpol_ANTCEN_picen,C=C,H=H,O=O,single=single,pi=pi,dpi=dpi,jmtype='TholeLinIso',screenradius=screen), tholefit('../Molecules/Pc_picen_fit_aniso_no_charge.xyz',HFpol_Pc_picen,C=C,H=H,O=O,single=single,pi=pi,dpi=dpi,jmtype='TholeLinIso',screenradius=screen)

	etatot=0
	for i in np.arange(0,len(etalist),1):
		etatot=etatot+etalist[i]
	print 'etatot'
	print etatot
	return etatot

#etatot(screen=2.1304, H=3.5020, C=10.1756, O=6.39400 ,single=0.000001, double=0.000001,pi=0.000001,dpi=0.000001) # Thole Linear as on Swart Website (a differes from paper)

etatot(screen=2.1304, H=2.7927, C=8.6959, O=5.74940 ,single=0.000001, double=0.000001,pi=0.000001,dpi=0.000001) # Thole Exp Empirical

#etatot(screen=1.68384957924, H=1.3849, C=6.5780, single=0.000001, double=0.000001,pi=0.000001,dpi=0.0001) # Thole Exp Components

#etatot(screen=2.8990, H=1.3849, C=6.5780, single=0.000001, double=0.000001,pi=0.000001,dpi=0.0001) # Thole Exp Components

#etatot(screen=2.1304, H=2.7927, C=8.6959, single=0.000001, double=0.000001,pi=0.000001,dpi=0.000001) # Thole Exp Empirical

#etatot(screen=1.7278, H=3.5020, C=10.1756, single=0.0001, double=0.0001,pi=0.0001,dpi=0.0001) # Thole Lin Empirical

#etatot(screen=1.6623, H=0.5974, C=7.7764, single=0.0001, double=0.0001,pi=0.0001,dpi=0.0001) # Thole Lin Components


#'H':2.7927, 'C':8.6959, 'N': 6.5565, 'O': 5.7494, 'F': 3.0013, 'S': 16.6984, 'Cl':16.1979, 'Br': 23.5714, 'I':36.9880, 'screen':2.1304}\n"

#m=minuit.Minuit(etatot, H=0.001,fix_H=True,limit_single=(0.0100,20.0),limit_pi=(0.0100,20.0),limit_dpi=(0.0100,20.0),limit_S=(0.0100,20.0), screen=2.8990, fix_screen=True)
#m.scan(("C",5,0.0100,20.0),("single",5,0.0100,20.0),("pi",5,0.0100,20.0),("dpi",5,0.0100,20.0),("S",5,0.0100,20.0),)
#m.printMode = 1
#m.migrad()
