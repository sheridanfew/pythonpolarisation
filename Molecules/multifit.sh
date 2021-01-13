FitScreen="FALSE"
FITSCREENONLY="FALSE"
FitDist=0.5
PolHydrogens="FALSE"
STARTPARS="empirical"
NOBONDS="FALSE"
POLFILE=""

function USAGE()
{
 cat << EOF
Makes opts of atomic pol parameters for set of mols calibrated to HF results. Should be run from Multimol_Pol_Determination folder

USAGE:  multipolfit.sh -n nTs 1T 2T 3T 4T...

OPTIONS:
	-a			Fit screening parameter
	-b			No bonds, default $NOBONDS
	-c			Fit only screening parameter.
	-d			Distance from start value to scan (in proportion of start value, default 0.5)
	-f			Fit parameter to use (eta, abs_square_dif_fit, tholefit)
	-h			Make hydrogens polarisable and fit, default $PolHydrogens
	-n			name of series
	-p			name of polfile (assumes pol_$( echo \${file} | sed 's/_aniso_w_connectivity//') if not stated)
	-j 			JM Type (Stern, TholeLinAniso, TholeExp)
	-s			Choose start parameters to use for fit (either van Duijnen's fit to components of calculated pol. tensors, 'components', or to an empirical fit, 'empirical', default is '$STARTPARS'.)

EOF
}



while getopts "abcd:f:hj:n:p:s:?" Option
do
    case $Option in
	a	  ) FitScreen="TRUE";;
	b	  ) NOBONDS="TRUE";;
	c	  ) FITSCREENONLY="TRUE";;
	d	  ) FitDist="0.5";;
	f    ) FitType=$OPTARG;;
	h    ) PolHydrogens="TRUE";;
	j    ) JMType=$OPTARG;;
	n    ) name="$OPTARG";;
	p	  ) POLFILE="$OPTARG";;
	s    ) STARTPARS=$OPTARG;;
        ?    ) USAGE
               exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               USAGE   # DEFAULT
    esac
done

shift $(($OPTIND - 1))

cat << EOF > outputfile.py
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

fitdist=${FitDist}

EOF

if [[ $STARTPARS == 'components' ]]; then
	if [[ $JMType == 'TholeLinAniso' ]]; then echo -en "cenpos={'H':0.5974, 'C':7.7764, 'N': 4.7575, 'O': 2.6442, 'F': 0.8149, 'S': 11.0208, 'Se': 11.0208, 'Cl':7.5682, 'Br': 12.7339, 'I':20.2867, 'screen':1.6623}\n" >> outputfile.py
	elif [[ $JMType == 'TholeExp' ]]; then echo -en "cenpos={'H':1.3849, 'C':6.5780, 'N': 2.9855, 'O': 2.3184, 'F': 1.4906, 'S': 10.8651, 'Se': 10.8651, 'Cl':6.7129, 'Br': 12.3384, 'I':20.0081, 'screen':2.8990}\n" >> outputfile.py
	fi
fi

if [[ $STARTPARS == 'empirical' ]]; then
	if [[ $JMType == 'TholeLinAniso' ]]; then echo -en "cenpos={'H':3.5020, 'C':10.01756, 'N': 7.6048, 'O': 6.3940, 'F': 2.9413, 'S': 19.7422, 'Se': 19.7422, 'Cl':16.1145, 'Br': 22.6671, 'I':34.2434, 'screen':1.7278}\n" >> outputfile.py
	elif [[ $JMType == 'TholeExp' ]]; then echo -en "cenpos={'H':2.7927, 'C':8.6959, 'N': 6.5565, 'O': 5.7494, 'F': 3.0013, 'S': 16.6984, 'Se': 19.7422,  'Cl':16.1979, 'Br': 23.5714, 'I':36.9880, 'screen':2.1304}\n" >> outputfile.py
	fi
fi

if [[ $JMType == 'Stern' ]]; then echo -en "cenpos={'H':2.0, 'C':8.0, 'N': (6.74833515*0.764), 'O': (6.74833515*0.405), 'F': 3.58, 'S': 23.2, 'Se': 19.7422, 'Cl':17.6, 'Br': 25.6, 'I':42.6, 'screen':2.8990}\n" >> outputfile.py
fi

 echo -en "cenpos=dict(cenpos, **{'single':5.0, 'pi':10.0, 'double': 10.0, 'triple': 10.0, 'dpi':10.0})\n" >> outputfile.py


for file in $@ 
do
 for atomtype in C H N O S Se Si
 do
  #Sets fitting criteria for individual atom types. Only does this if atom is present in molecule being considered, and not already put in this atom type for another molecule.
  if [[ $(grep -c "elname'\:'$atomtype'" ${file}_aniso_no_charge.xyz) -ge 1 && $(echo "$etaopts" | grep -c "${atomtype}=${atomtype}") == "0" ]]
  then

	etacalcopts="${etacalcopts},${atomtype}"
	etaopts="${etaopts},${atomtype}=${atomtype}"
	mvaluesout="${mvaluesout}m.values['${atomtype}'],"


	if [[ $atomtype == 'H' && $PolHydrogens == "FALSE" ]]
	then
	  minuitopts="${minuitopts},H=0.0001,fix_H=True"
	  minuitopts_wscan="${minuitopts_wscan}H=0.0001,fix_H=True,"	
	elif [[ $FITSCREENONLY == 'TRUE' ]]
	then
	  minuitopts="${minuitopts},${atomtype}=cenpos['${atomtype}'],fix_${atomtype}=True"
	  minuitopts_wscan="${minuitopts_wscan}${atomtype}=cenpos['${atomtype}'],fix_${atomtype}=True,"
	else
	  minuitopts="${minuitopts},${atomtype}=cenpos['${atomtype}'],limit_${atomtype}=((cenpos['${atomtype}']*(1.0-fitdist)),(cenpos['${atomtype}']*(1.0+fitdist)))"
	  minuitopts_wscan="${minuitopts_wscan}limit_${atomtype}=((cenpos['${atomtype}']*(1.0-fitdist)),(cenpos['${atomtype}']*(1.0+fitdist))),"
	  minuitscanopts="${minuitscanopts}(\"${atomtype}\",5,(cenpos['${atomtype}']*(1.0-fitdist)),(cenpos['${atomtype}']*(1.0+fitdist))),"	
	fi
  else
	mvaluesout="${mvaluesout}'N/A ',"
  fi
 done

 for atomtype in single pi double triple dpi
 do
  #Sets fitting criteria for individual bond types. Only does this if bond is present in molecule being considered, and not already put in this bond type for another molecule.
  if [[ $(grep -c "elname'\:'$atomtype'" ${file}_aniso_no_charge.xyz) -ge 1 && $(echo "$etaopts" | grep -c "${atomtype}=${atomtype}") == "0" ]]
  then

	etacalcopts="${etacalcopts},${atomtype}"
	etaopts="${etaopts},${atomtype}=${atomtype}"
	mvaluesout="${mvaluesout}m.values['${atomtype}'],"

	if [[ $NOBONDS == "TRUE" ]]
	then
	  minuitopts="${minuitopts},${atomtype}=0.0001,fix_${atomtype}=True"
	  minuitopts_wscan="${minuitopts_wscan}${atomtype}=0.0001,fix_${atomtype}=True,"	
	elif [[ $FITSCREENONLY == 'TRUE' ]]
	then
	  minuitopts="${minuitopts},${atomtype}=cenpos['${atomtype}'],fix_${atomtype}=True"
	  minuitopts_wscan="${minuitopts_wscan}${atomtype}=cenpos['${atomtype}'],fix_${atomtype}=True,"
	else
	  minuitopts="${minuitopts},${atomtype}=cenpos['${atomtype}'],limit_${atomtype}=(0.0001,(cenpos['${atomtype}']*2))"
	  minuitopts_wscan="${minuitopts_wscan}limit_${atomtype}=(0.0001,(cenpos['${atomtype}']*2)),"
	  minuitscanopts="${minuitscanopts}(\"${atomtype}\",5,0.0001,(cenpos['${atomtype}']*2)),"	
	fi
  else
	mvaluesout="${mvaluesout}'N/A ',"
  fi
 done
done


#NB Empirical exponential fit seems to be the one they are happiest with

for file in $@ 
do
echo -en "from Molecules.pol_$( echo ${file} | sed 's/_aniso_w_connectivity//') import pol\nHFpol_${file}=pol\n" >> outputfile.py
done

echo -en "def etatot( screen$etacalcopts):\n\tetalist=[" >> outputfile.py

# ./multifitpol.sh -n nTnew thio_1T_neut thio_2T_neut_fit_aniso thio_3T_neut_fit_aniso thio_4T_neut_fit_aniso thio_5T_neut_fit_aniso thio_6T_neut_fit_aniso thio_7T_neut_fit_aniso thio_8T_neut_fit_aniso

for file in $@ 
do
  if [[ $file != $1 ]]; then echo -en ", " >> outputfile.py; fi
  echo -en "${FitType}('../Molecules/${file}_aniso_no_charge.xyz',HFpol_${file}$etaopts,jmtype='$JMType',screenradius=screen)" >> outputfile.py
done

cat << EOF >> outputfile.py
]
	etatot=0
	for i in np.arange(0,len(etalist),1):
		etatot=etatot+etalist[i]
	print 'etatot'
	print etatot
	return etatot
EOF


basename="${name}_${JMType}_${FitType}_"

if [[ $NOBONDS == "TRUE" ]]; then basename="${basename}nobonds_"; fi
if [[ $PolHydrogens == "FALSE" ]]; then basename="${basename}noH_"; 
else basename="${basename}wH_";  fi
if [[ $FITSCREENONLY == "TRUE" ]]; then basename="${basename}screenonly_"; fi

basename="${basename}${STARTPARS}start_"


if [[ $FitScreen = "FALSE" ]]
then
	#make version without fitscreen, without scan
	cp outputfile.py ../Multimol_Pol_Determination/${basename}fixscreen_multipolfit.py
	cat << EOF >> ../Multimol_Pol_Determination/${basename}fixscreen_multipolfit.py
m=minuit.Minuit(etatot, screen=cenpos['screen'], fix_screen=True $minuitopts)
#m.printMode = 1
m.migrad()
EOF

	#make version without fitscreen, with scan
	cp outputfile.py ../Multimol_Pol_Determination/${basename}fixscreen_multipolfit_wscan.py
	cat << EOF >> ../Multimol_Pol_Determination/${basename}fixscreen_multipolfit_wscan.py
m=minuit.Minuit(etatot, ${minuitopts_wscan} screen=cenpos['screen'], fix_screen=True)
m.scan($minuitscanopts)
#m.printMode = 1
m.migrad()
EOF

else
	#make version with fitscreen, without scan
	cp outputfile.py ../Multimol_Pol_Determination/${basename}fitscreen_multipolfit.py
	cat << EOF >> ../Multimol_Pol_Determination/${basename}fitscreen_multipolfit.py
m=minuit.Minuit(etatot, screen=cenpos['screen'], limit_screen=((cenpos['screen']*(1.0-fitdist)),(cenpos['screen']*(1.0+fitdist)))$minuitopts)
#m.printMode = 1
m.migrad()
EOF

	#make version without fitscreen, with scan
	cp outputfile.py ../Multimol_Pol_Determination/${basename}fitscreen_multipolfit_wscan.py
	cat << EOF >> ../Multimol_Pol_Determination/${basename}fitscreen_multipolfit_wscan.py
m=minuit.Minuit(etatot, ${minuitopts_wscan} limit_screen=((cenpos['screen']*(1.0-fitdist)),(cenpos['screen']*(1.0+fitdist))))
m.scan($minuitscanopts ("screen",5,(cenpos['screen']*(1.0-fitdist)),(cenpos['screen']*(1.0+fitdist))))
#m.printMode = 1
m.migrad()
EOF

fi

rm outputfile.py
