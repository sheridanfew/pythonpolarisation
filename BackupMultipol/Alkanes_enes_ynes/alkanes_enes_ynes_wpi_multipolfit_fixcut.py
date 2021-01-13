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

from Molecules.pol_chain_1_3_5_7_9_undecapentaene_pi_fit_aniso import pol
pol_chain_1_3_5_7_9_undecapentaene_pi_fit_aniso=pol
from Molecules.pol_chain_1_3_5_7_9_undecapentayne_doublepi_fit_aniso import pol
pol_chain_1_3_5_7_9_undecapentayne_doublepi_fit_aniso=pol
from Molecules.pol_chain_1_3_5_7_nonatetraene_pi_fit_aniso import pol
pol_chain_1_3_5_7_nonatetraene_pi_fit_aniso=pol
from Molecules.pol_chain_1_3_5_7_nonatetrayne_doublepi_fit_aniso import pol
pol_chain_1_3_5_7_nonatetrayne_doublepi_fit_aniso=pol
from Molecules.pol_chain_1_3_5_heptatriene_pi_fit_aniso import pol
pol_chain_1_3_5_heptatriene_pi_fit_aniso=pol
from Molecules.pol_chain_1_3_5_heptatriyne_doublepi_fit_aniso import pol
pol_chain_1_3_5_heptatriyne_doublepi_fit_aniso=pol
from Molecules.pol_chain_1_3_pentadiene_pi_fit_aniso import pol
pol_chain_1_3_pentadiene_pi_fit_aniso=pol
from Molecules.pol_chain_1_3_pentadiyne_doublepi_fit_aniso import pol
pol_chain_1_3_pentadiyne_doublepi_fit_aniso=pol
from Molecules.pol_chain_2_4_6_8_10_dodecapentaene_pi_fit_aniso import pol
pol_chain_2_4_6_8_10_dodecapentaene_pi_fit_aniso=pol
from Molecules.pol_chain_2_4_6_8_10_dodecapentayne_doublepi_fit_aniso import pol
pol_chain_2_4_6_8_10_dodecapentayne_doublepi_fit_aniso=pol
from Molecules.pol_chain_2_4_6_8_decatetraene_pi_fit_aniso import pol
pol_chain_2_4_6_8_decatetraene_pi_fit_aniso=pol
from Molecules.pol_chain_2_4_6_8_decatetrayne_doublepi_fit_aniso import pol
pol_chain_2_4_6_8_decatetrayne_doublepi_fit_aniso=pol
from Molecules.pol_chain_2_4_6_octatriene_pi_fit_aniso import pol
pol_chain_2_4_6_octatriene_pi_fit_aniso=pol
from Molecules.pol_chain_2_4_6_octatriyne_doublepi_fit_aniso import pol
pol_chain_2_4_6_octatriyne_doublepi_fit_aniso=pol
from Molecules.pol_chain_2_4_hexadiene_pi_fit_aniso import pol
pol_chain_2_4_hexadiene_pi_fit_aniso=pol
from Molecules.pol_chain_2_4_hexadiyne_doublepi_fit_aniso import pol
pol_chain_2_4_hexadiyne_doublepi_fit_aniso=pol
from Molecules.pol_chain_butane_fit_aniso import pol
pol_chain_butane_fit_aniso=pol
from Molecules.pol_chain_butene_pi_fit_aniso import pol
pol_chain_butene_pi_fit_aniso=pol
from Molecules.pol_chain_butyne_doublepi_fit_aniso import pol
pol_chain_butyne_doublepi_fit_aniso=pol
from Molecules.pol_chain_decane_fit_aniso import pol
pol_chain_decane_fit_aniso=pol
from Molecules.pol_chain_dodecane_fit_aniso import pol
pol_chain_dodecane_fit_aniso=pol
from Molecules.pol_chain_ethane_fit_aniso import pol
pol_chain_ethane_fit_aniso=pol
from Molecules.pol_chain_heptane_fit_aniso import pol
pol_chain_heptane_fit_aniso=pol
from Molecules.pol_chain_hexane_fit_aniso import pol
pol_chain_hexane_fit_aniso=pol
from Molecules.pol_chain_methane_fit_aniso import pol
pol_chain_methane_fit_aniso=pol
from Molecules.pol_chain_nonane_fit_aniso import pol
pol_chain_nonane_fit_aniso=pol
from Molecules.pol_chain_octane_fit_aniso import pol
pol_chain_octane_fit_aniso=pol
from Molecules.pol_chain_pentane_fit_aniso import pol
pol_chain_pentane_fit_aniso=pol
from Molecules.pol_chain_propane_fit_aniso import pol
pol_chain_propane_fit_aniso=pol
from Molecules.pol_chain_undecane_fit_aniso import pol
pol_chain_undecane_fit_aniso=pol
def etatot( cut,C,H,single,pi):
	etalist=[eta('../Molecules/chain_1_3_5_7_9_undecapentaene_pi_fit_aniso_no_charge.xyz',pol_chain_1_3_5_7_9_undecapentaene_pi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_1_3_5_7_9_undecapentayne_doublepi_fit_aniso_no_charge.xyz',pol_chain_1_3_5_7_9_undecapentayne_doublepi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_1_3_5_7_nonatetraene_pi_fit_aniso_no_charge.xyz',pol_chain_1_3_5_7_nonatetraene_pi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_1_3_5_7_nonatetrayne_doublepi_fit_aniso_no_charge.xyz',pol_chain_1_3_5_7_nonatetrayne_doublepi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_1_3_5_heptatriene_pi_fit_aniso_no_charge.xyz',pol_chain_1_3_5_heptatriene_pi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_1_3_5_heptatriyne_doublepi_fit_aniso_no_charge.xyz',pol_chain_1_3_5_heptatriyne_doublepi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_1_3_pentadiene_pi_fit_aniso_no_charge.xyz',pol_chain_1_3_pentadiene_pi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_1_3_pentadiyne_doublepi_fit_aniso_no_charge.xyz',pol_chain_1_3_pentadiyne_doublepi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_2_4_6_8_10_dodecapentaene_pi_fit_aniso_no_charge.xyz',pol_chain_2_4_6_8_10_dodecapentaene_pi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_2_4_6_8_10_dodecapentayne_doublepi_fit_aniso_no_charge.xyz',pol_chain_2_4_6_8_10_dodecapentayne_doublepi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_2_4_6_8_decatetraene_pi_fit_aniso_no_charge.xyz',pol_chain_2_4_6_8_decatetraene_pi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_2_4_6_8_decatetrayne_doublepi_fit_aniso_no_charge.xyz',pol_chain_2_4_6_8_decatetrayne_doublepi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_2_4_6_octatriene_pi_fit_aniso_no_charge.xyz',pol_chain_2_4_6_octatriene_pi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_2_4_6_octatriyne_doublepi_fit_aniso_no_charge.xyz',pol_chain_2_4_6_octatriyne_doublepi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_2_4_hexadiene_pi_fit_aniso_no_charge.xyz',pol_chain_2_4_hexadiene_pi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_2_4_hexadiyne_doublepi_fit_aniso_no_charge.xyz',pol_chain_2_4_hexadiyne_doublepi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_butane_fit_aniso_no_charge.xyz',pol_chain_butane_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_butene_pi_fit_aniso_no_charge.xyz',pol_chain_butene_pi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_butyne_doublepi_fit_aniso_no_charge.xyz',pol_chain_butyne_doublepi_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_decane_fit_aniso_no_charge.xyz',pol_chain_decane_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_dodecane_fit_aniso_no_charge.xyz',pol_chain_dodecane_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_ethane_fit_aniso_no_charge.xyz',pol_chain_ethane_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_heptane_fit_aniso_no_charge.xyz',pol_chain_heptane_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_hexane_fit_aniso_no_charge.xyz',pol_chain_hexane_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_methane_fit_aniso_no_charge.xyz',pol_chain_methane_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_nonane_fit_aniso_no_charge.xyz',pol_chain_nonane_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_octane_fit_aniso_no_charge.xyz',pol_chain_octane_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_pentane_fit_aniso_no_charge.xyz',pol_chain_pentane_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_propane_fit_aniso_no_charge.xyz',pol_chain_propane_fit_aniso,cut,C=C,H=H,single=single,pi=pi), eta('../Molecules/chain_undecane_fit_aniso_no_charge.xyz',pol_chain_undecane_fit_aniso,cut,C=C,H=H,single=single,pi=pi)]
	etatot=0
	for i in np.arange(0,len(etalist),1):
		etatot=etatot+etalist[i]
	print 'etatot'
	print etatot
	return etatot
