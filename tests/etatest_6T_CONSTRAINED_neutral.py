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
#import minuit
import time

from Molecules.mol_6T_CONSTRAINED_neutral_HFpol import HFpol
HFpol_6T_CONSTRAINED_neutral=HFpol

def etacalc( c, h, s, n, cut):
	etacalc=eta('../Molecules/6T_CONSTRAINED_neutral_no_charge.xyz',HFpol_6T_CONSTRAINED_neutral, c=c, h=h, s=s, n=n, cut=cut)
	return etacalc

etacalc(c=10.0, h=0., s=10.0, n=5.155, cut=5.0)

