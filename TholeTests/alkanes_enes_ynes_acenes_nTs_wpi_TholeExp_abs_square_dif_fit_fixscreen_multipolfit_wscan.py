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

from Molecules.pol_chain_1_3_5_7_9_undecapentaene_pi import pol
HFpol_chain_1_3_5_7_9_undecapentaene_pi=pol
from Molecules.pol_chain_1_3_5_7_9_undecapentayne_doublepi import pol
HFpol_chain_1_3_5_7_9_undecapentayne_doublepi=pol
from Molecules.pol_chain_1_3_5_7_nonatetraene_pi import pol
HFpol_chain_1_3_5_7_nonatetraene_pi=pol
from Molecules.pol_chain_1_3_5_7_nonatetrayne_doublepi import pol
HFpol_chain_1_3_5_7_nonatetrayne_doublepi=pol
from Molecules.pol_chain_1_3_5_heptatriene_pi import pol
HFpol_chain_1_3_5_heptatriene_pi=pol
from Molecules.pol_chain_1_3_5_heptatriyne_doublepi import pol
HFpol_chain_1_3_5_heptatriyne_doublepi=pol
from Molecules.pol_chain_1_3_pentadiene_pi import pol
HFpol_chain_1_3_pentadiene_pi=pol
from Molecules.pol_chain_1_3_pentadiyne_doublepi import pol
HFpol_chain_1_3_pentadiyne_doublepi=pol
from Molecules.pol_chain_2_4_6_8_10_dodecapentaene_pi import pol
HFpol_chain_2_4_6_8_10_dodecapentaene_pi=pol
from Molecules.pol_chain_2_4_6_8_10_dodecapentayne_doublepi import pol
HFpol_chain_2_4_6_8_10_dodecapentayne_doublepi=pol
from Molecules.pol_chain_2_4_6_8_decatetraene_pi import pol
HFpol_chain_2_4_6_8_decatetraene_pi=pol
from Molecules.pol_chain_2_4_6_8_decatetrayne_doublepi import pol
HFpol_chain_2_4_6_8_decatetrayne_doublepi=pol
from Molecules.pol_chain_2_4_6_octatriene_pi import pol
HFpol_chain_2_4_6_octatriene_pi=pol
from Molecules.pol_chain_2_4_6_octatriyne_doublepi import pol
HFpol_chain_2_4_6_octatriyne_doublepi=pol
from Molecules.pol_chain_2_4_hexadiene_pi import pol
HFpol_chain_2_4_hexadiene_pi=pol
from Molecules.pol_chain_2_4_hexadiyne_doublepi import pol
HFpol_chain_2_4_hexadiyne_doublepi=pol
from Molecules.pol_chain_butane import pol
HFpol_chain_butane=pol
from Molecules.pol_chain_butene_pi import pol
HFpol_chain_butene_pi=pol
from Molecules.pol_chain_butyne_doublepi import pol
HFpol_chain_butyne_doublepi=pol
from Molecules.pol_chain_decane import pol
HFpol_chain_decane=pol
from Molecules.pol_chain_dodecane import pol
HFpol_chain_dodecane=pol
from Molecules.pol_chain_ethane import pol
HFpol_chain_ethane=pol
from Molecules.pol_chain_heptane import pol
HFpol_chain_heptane=pol
from Molecules.pol_chain_hexane import pol
HFpol_chain_hexane=pol
from Molecules.pol_chain_methane import pol
HFpol_chain_methane=pol
from Molecules.pol_chain_nonane import pol
HFpol_chain_nonane=pol
from Molecules.pol_chain_octane import pol
HFpol_chain_octane=pol
from Molecules.pol_chain_pentane import pol
HFpol_chain_pentane=pol
from Molecules.pol_chain_propane import pol
HFpol_chain_propane=pol
from Molecules.pol_chain_undecane import pol
HFpol_chain_undecane=pol
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
from Molecules.pol_thio_1T_neut import pol
HFpol_thio_1T_neut=pol
from Molecules.pol_thio_2T_neut import pol
HFpol_thio_2T_neut=pol
from Molecules.pol_thio_3T_neut import pol
HFpol_thio_3T_neut=pol
from Molecules.pol_thio_4T_neut import pol
HFpol_thio_4T_neut=pol
from Molecules.pol_thio_5T_neut import pol
HFpol_thio_5T_neut=pol
from Molecules.pol_thio_6T_neut import pol
HFpol_thio_6T_neut=pol
from Molecules.pol_thio_7T_neut import pol
HFpol_thio_7T_neut=pol
from Molecules.pol_thio_8T_neut import pol
HFpol_thio_8T_neut=pol
def etatot( screen,C,H,single,pi,dpi,S):
	etalist=[abs_square_dif_fit('../Molecules/chain_1_3_5_7_9_undecapentaene_pi_fit_aniso_no_charge.xyz',HFpol_chain_1_3_5_7_9_undecapentaene_pi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_1_3_5_7_9_undecapentayne_doublepi_fit_aniso_no_charge.xyz',HFpol_chain_1_3_5_7_9_undecapentayne_doublepi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_1_3_5_7_nonatetraene_pi_fit_aniso_no_charge.xyz',HFpol_chain_1_3_5_7_nonatetraene_pi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_1_3_5_7_nonatetrayne_doublepi_fit_aniso_no_charge.xyz',HFpol_chain_1_3_5_7_nonatetrayne_doublepi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_1_3_5_heptatriene_pi_fit_aniso_no_charge.xyz',HFpol_chain_1_3_5_heptatriene_pi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_1_3_5_heptatriyne_doublepi_fit_aniso_no_charge.xyz',HFpol_chain_1_3_5_heptatriyne_doublepi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_1_3_pentadiene_pi_fit_aniso_no_charge.xyz',HFpol_chain_1_3_pentadiene_pi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_1_3_pentadiyne_doublepi_fit_aniso_no_charge.xyz',HFpol_chain_1_3_pentadiyne_doublepi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_2_4_6_8_10_dodecapentaene_pi_fit_aniso_no_charge.xyz',HFpol_chain_2_4_6_8_10_dodecapentaene_pi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_2_4_6_8_10_dodecapentayne_doublepi_fit_aniso_no_charge.xyz',HFpol_chain_2_4_6_8_10_dodecapentayne_doublepi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_2_4_6_8_decatetraene_pi_fit_aniso_no_charge.xyz',HFpol_chain_2_4_6_8_decatetraene_pi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_2_4_6_8_decatetrayne_doublepi_fit_aniso_no_charge.xyz',HFpol_chain_2_4_6_8_decatetrayne_doublepi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_2_4_6_octatriene_pi_fit_aniso_no_charge.xyz',HFpol_chain_2_4_6_octatriene_pi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_2_4_6_octatriyne_doublepi_fit_aniso_no_charge.xyz',HFpol_chain_2_4_6_octatriyne_doublepi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_2_4_hexadiene_pi_fit_aniso_no_charge.xyz',HFpol_chain_2_4_hexadiene_pi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_2_4_hexadiyne_doublepi_fit_aniso_no_charge.xyz',HFpol_chain_2_4_hexadiyne_doublepi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_butane_fit_aniso_no_charge.xyz',HFpol_chain_butane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_butene_pi_fit_aniso_no_charge.xyz',HFpol_chain_butene_pi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_butyne_doublepi_fit_aniso_no_charge.xyz',HFpol_chain_butyne_doublepi,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_decane_fit_aniso_no_charge.xyz',HFpol_chain_decane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_dodecane_fit_aniso_no_charge.xyz',HFpol_chain_dodecane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_ethane_fit_aniso_no_charge.xyz',HFpol_chain_ethane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_heptane_fit_aniso_no_charge.xyz',HFpol_chain_heptane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_hexane_fit_aniso_no_charge.xyz',HFpol_chain_hexane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_methane_fit_aniso_no_charge.xyz',HFpol_chain_methane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_nonane_fit_aniso_no_charge.xyz',HFpol_chain_nonane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_octane_fit_aniso_no_charge.xyz',HFpol_chain_octane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_pentane_fit_aniso_no_charge.xyz',HFpol_chain_pentane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_propane_fit_aniso_no_charge.xyz',HFpol_chain_propane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/chain_undecane_fit_aniso_no_charge.xyz',HFpol_chain_undecane,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/BENZENE_picen_fit_aniso_no_charge.xyz',HFpol_BENZENE_picen,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/NAPHTA_picen_fit_aniso_no_charge.xyz',HFpol_NAPHTA_picen,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/TETCEN_picen_fit_aniso_no_charge.xyz',HFpol_TETCEN_picen,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/ANTCEN_picen_fit_aniso_no_charge.xyz',HFpol_ANTCEN_picen,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/Pc_picen_fit_aniso_no_charge.xyz',HFpol_Pc_picen,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/thio_1T_neut_fit_aniso_no_charge.xyz',HFpol_thio_1T_neut,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/thio_2T_neut_fit_aniso_no_charge.xyz',HFpol_thio_2T_neut,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/thio_3T_neut_fit_aniso_no_charge.xyz',HFpol_thio_3T_neut,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/thio_4T_neut_fit_aniso_no_charge.xyz',HFpol_thio_4T_neut,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/thio_5T_neut_fit_aniso_no_charge.xyz',HFpol_thio_5T_neut,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/thio_6T_neut_fit_aniso_no_charge.xyz',HFpol_thio_6T_neut,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/thio_7T_neut_fit_aniso_no_charge.xyz',HFpol_thio_7T_neut,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen), abs_square_dif_fit('../Molecules/thio_8T_neut_fit_aniso_no_charge.xyz',HFpol_thio_8T_neut,C=C,H=H,single=single,pi=pi,dpi=dpi,S=S,jmtype='TholeExp',screenradius=screen)]
	etatot=0
	for i in np.arange(0,len(etalist),1):
		etatot=etatot+etalist[i]
	print 'etatot'
	print etatot
	return etatot
m=minuit.Minuit(etatot, limit_C=(0.0100,20.0),H=0.001,fix_H=True,limit_single=(0.0100,20.0),limit_pi=(0.0100,20.0),limit_dpi=(0.0100,20.0),limit_S=(0.0100,20.0), screen=2.8990, fix_screen=True)
m.scan(("C",5,0.0100,20.0),("single",5,0.0100,20.0),("pi",5,0.0100,20.0),("dpi",5,0.0100,20.0),("S",5,0.0100,20.0),)
m.printMode = 1
m.migrad()
