#Script to convert gaussian output files to molecules in python polarisation model
#
#Sheridan Few

file=""
log=""
cif_nos_log_orient_file=""
connectivity_file=""
iso="FALSE"
specifiedname=""

function USAGE()
{
 cat << EOF
Converts log files from gaussian to xyzs for molecule factory in James Python Polarisation script.

Takes as input a gaussian log file for charges and polarisability, and an xyz with molecular orientation as in the cif file. Exports two xyz files, one with the orientation of the gaussian log (for calibrating atomic polarisabilities to molecular polarisability from gaussian) and one with the orientation of the xyz_with_cif_structure.xyz, but atomic charges from the gaussian file.

Can also specify a third file if atom numbers of these two do not match up (see option -a)

Converts to Bohr automatically, but can leave in A, and add "'convert_to_bohr':True," to each line like:

{'convert_to_bohr':True,'elname':'C', 'pos':Position([ -3.2281290, 0.1935640, -0.0000590]), 'crg':-0.1441700, 'pol':Polarizability(iso=1.0)}

USAGE: ./gauss_log_to_py_xyz.sh -b ${struct}.xyz -n ${struct} -y ${struct}_connectivity.dat ${struct}_neut_HF631PLUSGSTAR.log
	(exports blah_mul.xyz, blah_chelpg.xyz, blah_no_charge.py pol_blah.py)

OPTIONS:
	-a			(Only required if extracting atomic charges from a log with different atom numbering from cif structure. If not specified, then script assumes numbering of xyz from cif file for structure, and gaussian .log are ther same) xyz w. atom nos corresponding to cif file, but orientation as in log file.
	-b 			xyz file with final coords (as in original cif file, for making crystals)
	-i			Isotropic (Will assign atoms by no. of bonds, aniso will treat all atoms of same element same iso polarisability, but adds anisotropically polarisable bonds)
	-k			Alkene: Make all bonds between pairs of carbons pi bonds.
	-l			(redundant here) .log file to use for charges (if different from main file, used for cryst structure)
	-n			Specify name
	-y			Connectivity data file (can get by getting gjf from cif, then opening&resaving as neutral (0 1) in gaussview, extract connectivity lines at end. can use extract_connectivity.sh to get this. NB. Should check gaussview's assignments and edit in connectivity file.) All atoms of same element will be treated as identical if this option is not specified.)

EOF
}



while getopts "a:b:ikl:n:y:?" Option
do
    case $Option in
		a    ) cif_nos_log_orient_file=$OPTARG;;
		b    ) cif_xyz_file=$OPTARG;;
		i    ) iso="TRUE";;
		k    ) ALKENE="TRUE";;
		l    ) log=$OPTARG;;
		n    ) specifiedname=$OPTARG;;
		y    ) connectivity_file=$OPTARG;;
		?    ) USAGE
		       exit 0;;
		*    ) echo ""
		       echo "Unimplemented option chosen."
		       USAGE   # DEFAULT
    esac
done

shift $(($OPTIND - 1))



for file in "$@"
do
#Check files exist
	if [ ! -f ${file} ]
	then
               	echo "file $file does not exist."
		exit 2 
	elif [[ ! -f ${log} && $log != "" ]]
	then
               	echo "file $log does not exist."
		exit 2 
	elif [[ ! -f ${cif_nos_log_orient_file} && $cif_nos_log_orient_file != "" ]]
	then
               	echo "file $cif_nos_log_orient_file does not exist."
		exit 2 
	fi

	nstarts=$(grep -c 'Will use up to' $file )
        nfins=$(grep -c 'Normal termination' $file )

	if [[ $nstarts != $nfins ]]
	then
	  echo "Job $file did not terminate correctly. Will not use coords. Exiting"
	  exit 2
	fi

#	rename 's/*/STAR/g' $file
#	rename 's/+/PLUS/g' $file

done


for file in "$@"
do
#name based on whether orienting differently
  echo "Processing $file, iso = $iso"

  basename="$(echo "${cif_xyz_file}" | sed 's/.log//g' | sed 's/.xyz/_/g' | sed 's/-/_/g' )"
  if [[ $specifiedname != "" ]]
  then
    name="${specifiedname}"
  else
    if [[ $cif_nos_log_orient_file == "" ]]
    then
      name="${basename}"
    else
      name="$(echo "${cif_xyz_file}_${cif_nos_log_orient_file}_${basename}" | sed 's/.log//g' | sed 's/.xyz/_/g' | sed 's/-/_/g' )"
    fi
  fi

  if [[ $iso == "TRUE" ]]
  then
    name="${name}_iso"
  else
    name="${name}_aniso"
  fi

#  if [[ $connectivity_file != "" ]]
#  then
#    name=$(echo "${name}_w_connectivity")
#  fi

name_for_file="${name}_cut_$(echo "$cutset" | sed "s/\./_/")_$(date +"%Y_%m_%d")"
if [[ -f ${name}_mul.xyz || -f ${name}_chelpg.xyz || -f ${name}_no_charge.xyz || -f ${name}_cifstuct_mul.xyz || -f ${name}_cifstuct_chelpg.xyz || -f ${name}_cifstuct_no_charge.xyz ]]
then
    echo "xyz of this name, ${name}_(mul/chelpg/no_charge).xyz already exists. Exiting."
    exit 2
fi

### CONNECTIVITY 

echo "Recording connectivity"

if [[ $connectivity_file != "" ]]
then
	#PART FOR ISOTROPIC MODEL
  	#Results in array NBonds where NBonds[i] is no. of bonds on atom i.
	connectivity_count=1
	#Connectivity count keeps track of current line (&corresponding atom no) analysing in connectivity file 
	if [[ $iso == "TRUE" ]]
	then
	while read line
	do
	  #Number of bonds to this atom calculated from no. of entried on line
	  Nbondsnew=$(( ( $(echo $line | wc -w) - 1 ) / 2))
#	  echo "Nbondsnew($connectivity_count)=$Nbondsnew"
	  NBonds[$connectivity_count]=$(( ${NBonds[${connectivity_count}]} + $Nbondsnew ))
	  #One bond added to bond count for each atom which this atom is bonded to
	  NBonds[$(echo $line | awk '{printf "%s",$2}')]=$(( ${NBonds[$(echo $line | awk '{printf "%s",$2}')]} + 1 ))
	  NBonds[$(echo $line | awk '{printf "%s",$4}')]=$(( ${NBonds[$(echo $line | awk '{printf "%s",$4}')]} + 1 ))
	  NBonds[$(echo $line | awk '{printf "%s",$6}')]=$(( ${NBonds[$(echo $line | awk '{printf "%s",$6}')]} + 1 ))
	  NBonds[$(echo $line | awk '{printf "%s",$8}')]=$(( ${NBonds[$(echo $line | awk '{printf "%s",$8}')]} + 1 ))
	  #If there aren't entries for $4,6,8 contributions will all go to array member NBonds[0], not called in this script so no problem)
	  connectivity_count=$(( $connectivity_count + 1 ))
	done < $connectivity_file
	fi

	if [[ $iso == "FALSE" ]]
	then
	#PART FOR ANISOTROPIC MODEL
  	#Results in arrays:
	# bond_partners, tracking bond partners	    # eg. if bond n is between atoms 1 and 2, bond_partners[n]="1 2"
	# bond_type, tracking bond type		    # eg. if bond n is single, bond_type[n]="single"
	#
	# string bond_count_total, total no. of bonds.

	connectivity_count=1
	#Connectivity count keeps track of current line (&corresponding atom no) analysing in connectivity file 
	bond_count=0
	#Keeps track and gives index to all bonds. Set to 0 initially so that final number is total # bonds.
	while read line
	do
	  nentries=$(( ($(echo $line | awk '{printf "%d",NF}') - 1) /2 ))
	  #number of bonds on line
#	  echo "line number: $connectivity_count, line: $line, nentries: $nentries"
	  entry=1
	  #Keeps track of bond number to this atom (line)
	  if [[ $entry -le $nentries ]]
	  then
	   while [[ $entry -le $nentries ]]
	   do
#	    echo "entry: $entry"
	    bond_count=$(( $bond_count + 1 ))

	    #DEFINE BOND PARTNERS
	    bondpartner_argno=$(( 2 * $entry ))
	    bond_partners[$bond_count]="$connectivity_count $(echo $line | cut -d " " -f $bondpartner_argno )"
	    # eg. if bond n is between atoms 1 and 2, bond_partners[n]="1 2"
#	    echo "bond_partners[$bond_count]=${bond_partners[$bond_count]}"

	    #DEFINE BOND TYPE (NB. Connectivity file from Gaussian can be wrong, likes to make 1.5s 2.0/1.0s, should be checked manually )

	    bondtype_argno=$(( $bondpartner_argno + 1 ))
	    bond_type_gauss[$bond_count]="$(echo $line | awk -v N=$bondtype_argno '{print $N}')"

		if [ ${bond_type_gauss[$bond_count]} == '1.0' ]
		then
		  bond_type[$bond_count]="single"
		elif [ ${bond_type_gauss[${bond_count}]} == '1.5' ]
		then
		  bond_type[$bond_count]="pi"
		elif [ ${bond_type_gauss[${bond_count}]} == '2.0' ]
		then
		  bond_type[$bond_count]="double"
		elif [ ${bond_type_gauss[${bond_count}]} == '3.0' ]
		then
		  bond_type[$bond_count]="triple"
		else
		  echo "Can't understand bond type for bond $entry to atom $connectivity_count. Not single, pi, double, or triple."
		  echo "bond_type_gauss[${bond_count}] = ${bond_type_gauss[${bond_count}]}"
		  echo "line: $line"
		  echo -en "\n"
		  echo $line | cut -d " " -f $bondtype_argno
		  echo -en "\n"
		  echo "bond_type[${bond_count}] = ${bond_type[${bond_count}]}"
		  exit 2
		fi
		  echo "bond_type[${bond_count}] = ${bond_type[${bond_count}]}"
	    	entry=$(( $entry + 1 ))
	   done
	  fi
	  connectivity_count=$(( $connectivity_count + 1 ))
	done < $connectivity_file
	bond_count_total=$bond_count
	fi
fi

#Polarisability

if [[ $file == *log ]]
then
	echo "Recording desired Polarisability from $file"	
	linestartpol=$(( $(grep -n 'SCF Polarizability' "${file}" | tail -n 1 | awk -F":" 'NR==1{print $1}') + 2 ))
	echo $linestartpol
	pol[00]=$(cat ${file} | sed -n ${linestartpol}p | awk '{printf "%s",$2}' | sed 's/D/E/g' )
	pol[01]=$(cat ${file} | sed -n $(( ${linestartpol} + 1 ))p | awk '{printf "%s",$2}' | sed 's/D/E/g' )
	pol[02]=$(cat ${file} | sed -n $(( ${linestartpol} + 2 ))p | awk '{printf "%s",$2}' | sed 's/D/E/g' )
	pol[11]=$(cat ${file} | sed -n $(( ${linestartpol} + 1 ))p | awk '{printf "%s",$3}' | sed 's/D/E/g' )
	pol[12]=$(cat ${file} | sed -n $(( ${linestartpol} + 2 ))p | awk '{printf "%s",$3}' | sed 's/D/E/g' )
	pol[22]=$(cat ${file} | sed -n $(( ${linestartpol} + 2 ))p | awk '{printf "%s",$4}' | sed 's/D/E/g' )

	cat ${file} | sed -n $(( ${linestartpol} + 2 ))p | awk '{printf "%s",$2}' | sed 's/D/E/g'

	for i in 00 01 02 11 12 22
	do
		echo $i
		echo ${pol[${i}]}
	done

	cat << EOF > pol_${name}.py
import numpy as np
pol=np.matrix([[${pol[00]},0.0E+00,0.000000E+00],[${pol[01]},${pol[11]},0.000000E+00],[${pol[02]},${pol[12]},${pol[22]}]])
EOF

fi

#COORDINATES


####################################################################################################################################
#GET COORDS
#
# As long as gaussian opt was carried out with input atom numbers as in .cif file (can be converted to a gjf using openbabel), the only relevent coords here are the (x,y,z,atomtype)_cif 

#This section is messy as it was made with intention of dealing with different atom numbers in .cif and .log. This hasn't been a continuing issue, as I've used openbabel to convert .cif files to .gjf, which I have then used for Gaussian input files (maintaining correspondence between atom numbers). In this case, a gjf should be made with the cif, then reoriented using gaussview (symmetrize function) to roughly align with that in the log file. An xyz should then be made using these coords. Gaussview will probs put same way as log automatically, but if not, I have tools for rotating the xyz file. This script then compares this file in the log orientation with cif atom numbers to the output log file (by finding atoms at shortest physical distance), establishes correspondence, and makes an xyz readable by python pol. model with the cif coords, and charges from the correctly correspnding atoms in the log file.

#A correspondence is required to assign Mulliken/Chelpg charges extracted from Gaussian to atoms in the polarisation model.

########################################################################################################################################

#Atom pos'ns from coords.tmp, cif_xyz_file, cif_nos_log_orient_file to arrays

#get coords

echo "Making coords file"

if [[ $file == *log ]]
then
	get_fin_pos.awk $file > coords.tmp
	Natoms=$( grep 'NAtoms=' "${file}" | sed -n 1p | awk '{printf "%1.f",$2}')
elif [[ $file == *gjf ]]
then
	hashline=`expr $(grep -n '#' "${file}" | awk -F":" 'NR==1{print $1}') + 2`
	sed -n $(( $hashline + 5)),$(( $(cat $file | wc -l) - 1 ))p $file > coords.tmp
	trimnewlines.sh coords.tmp
	Natoms=$(cat coords.tmp | wc -l)
else
	echo "Can't read file type! Should be .log or .gjf"
	exit 2
fi

echo "Importing coords into arrays"

x=1
while [ $x -le $Natoms ]
do
# echo "Recording position of atom ${x} in $file"	
  posnline="${x}"
  x_log[$x]=$(cat coords.tmp | sed -n ${posnline}p | awk '{printf "%.7f",$2}')
  x_logbohr[$x]=$(echo "${x_log[${x}]} / 0.5291772" | bc -l )
  y_log[$x]=$(cat coords.tmp | sed -n ${posnline}p | awk '{printf "%.7f",$3}')
  y_logbohr[$x]=$(echo "${y_log[${x}]} / 0.5291772" | bc -l )	
  z_log[$x]=$(cat coords.tmp | sed -n ${posnline}p | awk '{printf "%.7f",$4}')
  z_logbohr[$x]=$(echo "${z_log[${x}]} / 0.5291772" | bc -l )
  atomtype_log[$x]=$(cat coords.tmp | sed -n ${posnline}p | awk '{printf "%s",$1}')
# echo "position of atom ${x} = ${x_log[${x}]},${y_log[${x}]},${z_log[${x}]}"

  if [[ $cif_xyz_file != "" ]]
  then
#   echo "Recording position of atom ${x} in $cif_xyz_file"
    x_cif[$x]=$(cat $cif_xyz_file | sed -n ${x}p | awk '{printf "%.7f",$2}')
    x_cifbohr[$x]=$(echo "${x_cif[${x}]} / 0.5291772" | bc -l )
    y_cif[$x]=$(cat $cif_xyz_file | sed -n ${x}p | awk '{printf "%.7f",$3}')
    y_cifbohr[$x]=$(echo "${y_cif[${x}]} / 0.5291772" | bc -l )	
    z_cif[$x]=$(cat $cif_xyz_file | sed -n ${x}p | awk '{printf "%.7f",$4}')
    z_cifbohr[$x]=$(echo "${z_cif[${x}]} / 0.5291772" | bc -l )
    atomtype_cif[$x]=$(cat $cif_xyz_file | sed -n ${posnline}p | awk '{printf "%s",$1}')
#   echo "position of atom ${x} = ${x_cif[${x}]},${y_cif[${x}]},${z_cif[${x}]}"
  fi
    

  if [[ $cif_nos_log_orient_file != "" ]]
  then
#   echo "Recording position of atom ${x} in $cif_nos_log_orient_file"
    x_cif_nos_log_orient[$x]=$(cat $cif_nos_log_orient_file | sed -n ${x}p | awk '{printf "%.7f",$2}')
#   x_cif_nos_log_orientbohr[$x]=$(echo "${x_cif_nos_log_orient[${x}]} / 0.5291772" | bc -l )
    y_cif_nos_log_orient[$x]=$(cat $cif_nos_log_orient_file | sed -n ${x}p | awk '{printf "%.7f",$3}')
#   y_cif_nos_log_orientbohr[$x]=$(echo "${y_cif_nos_log_orient[${x}]} / 0.5291772" | bc -l )	
    z_cif_nos_log_orient[$x]=$(cat $cif_nos_log_orient_file | sed -n ${x}p | awk '{printf "%.7f",$4}')
#   z_cif_nos_log_orientbohr[$x]=$(echo "${z_cif_nos_log_orient[${x}]} / 0.5291772" | bc -l )
    atomtype_cif_nos_log_orient[$x]=$(cat $cif_nos_log_orient_file | sed -n ${posnline}p | awk '{printf "%s",$1}')
#   echo "position of atom ${x} = ${x_cif_nos_log_orient[${x}]},${y_cif_nos_log_orient[${x}]},${z_cif_nos_log_orient[${x}]}"
  fi

  x=$(( $x + 1 ))
done

#Only need to worry about if using cif_nos_log_orient_file

if [[ $cif_nos_log_orient_file == "" ]]
then
  echo "no cif_nos_log_orient_file, => correspondence is array such that corrsspndence[n]='n'"
  i=1
  while [ $i -le $Natoms ]
  do
    correspondence[$i]="$i"
    i=$(( $i + 1 ))
  done
else
  # Work out correspondence between atom nos in cif and gausslog files (necessary to assign pointcharges to cif structure if gaussian log atom nos differ from those in cif structure file)
  echo "found cif_nos_log_orient_file, => working out correspondence between atom nos in cif and gausslog files"
  i=1
  while [ $i -le $Natoms ]
  do
    mindistabs="100000"
    j=1
    while [ $j -le $Natoms ]
    do
      distabsij=$(echo "( ((${x_log[$i]})-(${x_cif_nos_log_orient[$j]}))*((${x_log[$i]})-(${x_cif_nos_log_orient[$j]})) ) + ( ((${y_log[$i]})-(${y_cif_nos_log_orient[$j]}))*((${y_log[$i]})-(${y_cif_nos_log_orient[$j]})) ) + ( ((${z_log[$i]})-(${z_cif_nos_log_orient[$j]}))*((${z_log[$i]})-(${z_cif_nos_log_orient[$j]})) )" | bc -l)

      smaller=$(echo "$distabsij $mindistabs" | awk '{if ($1 < $2) print "TRUE"; else print "FALSE"}')

      if [[ $smaller == "TRUE" ]]
      then
	  mindistabs="$distabsij"
	  correspondence[$j]="$i"
      fi
      j=$(( $j + 1 ))
    done
    i=$(( $i + 1 ))
  done

  # Check correspondence is correct

  i=1
  while [ $i -le $Natoms ]
  do
    echo "${atomtype_cif_nos_log_orient[$i]} ${atomtype_log[${correspondence[$i]}]}"
    if [[ "${atomtype_cif_nos_log_orient[$i]}" != "${atomtype_log[${correspondence[$i]}]}" ]]
    then
      echo "correspondence wrong! Make sure orientation of cif_nos_log_orient_file is same as gaussian log (and maybe that this part of script works...!)"
      exit 2
    fi
    i=$(( $i + 1 ))
  done
fi

################################################################################################################################################################################################


#MULLIKEN CHARGES
#Results in array qmul giving Mulliken charges	# eg. if mul charge on atom n 0.01, qmul[n]="0.01"

echo "Recording Mulliken Charges"

if [[ $log = "" ]]
then
  chargesfile=$file
else
  chargesfile=$log
fi

linestartq=`expr $(grep -n 'Mulliken atomic charges:' "${chargesfile}" | tail -n 1 | awk -F":" 'NR==1{print $1}') + 2`

x=1
while [ $x -le $Natoms ]
do
#  echo "Recording Mulliken charge of atom ${x}"  
  chargeline=$(( ${linestartq} + ${x} - 1 ))
  qmul[$x]=$(cat ${chargesfile} | sed -n ${chargeline}p | awk '{printf "%.7f",$3}')
#  atomtype_log[$x]=$(cat ${chargesfile} | sed -n ${chargeline}p | awk '{printf "%s",$2}')
#  echo "Mulliken charge on atom ${x} = ${qmul[${x}]}"
  x=$(( $x + 1 ))
done

# CHELPG CHARGES	
#Results in array qchelpg giving chelpg charges	# eg. if chelpg charge on atom n 0.01, qchelpg[n]="0.01"


chelpg=$(grep -c 'Fitting point charges to electrostatic potential' $file )
if [[ $chelpg -ge 1 ]]
then
  echo "Recording Chelpg Charges"
  linestartq=`expr $(grep -n 'Fitting point charges to electrostatic potential' "${chargesfile}" | awk -F":" 'NR==1{print $1}') + 4`
  x=1
  while [ $x -le $Natoms ]
  do
  	if [[ $chelpg -ge 1 ]]
  		then
#   		  echo "Recording chelpg charge of atom ${x}"	
  		  chargeline=$(( ${linestartq} + ${x} - 1 ))
  		  qchelpg[$x]=$(cat ${chargesfile} | sed -n ${chargeline}p | awk '{printf "%.7f",$3}')
#  		  echo "Chelpg charge on atom ${x} = ${qchelpg1[${x}]}"
  	fi
    	x=$(( $x + 1 ))
  done
else
  echo "No Chelpg Charges found in file $file" 

fi

#EXPORT to xyz files readable by python pol model

#Results in files:
#[molname]_[mul,chelpg,no_charge].xyz	with orientations as in log file, to use for calibrating atomic polarisabilities to polarisability extracted from Gaussian
#[molname]_cifstruct_[mul,chelpg,no_charge].xyz	with orientations as in cif file, to use to build a crystal once atomic polarisabilities have been appropriately determined

#ADD ATOMS TO FILES

echo "Adding atoms to python .xyz files"

#EXPORT File with orientations from gaussian log:
x=1
while [ $x -le $Natoms ]
do

#correspondence is array such that corrsspndence[n]="n" if there is no cif_nos_log_orient_file

  atomtype_log_ID=$(echo "${atomtype_log[$x]}" | tr -d "\n" | od -An -t dC | sed 's/ //g')

  cat << EOF >> ${name}_mul.xyz
{'elname':'${atomtype_log[$x]}${NBonds[$x]}', 'pos':Position([ ${x_logbohr[${x}]}, ${y_logbohr[${x}]}, ${z_logbohr[${x}]}]), 'crg':${qmul[${correspondence[$x]}]}, 'pol':Polarizability(iso=1.0)}
EOF

	  cat << EOF >> ${name}_no_charge.xyz
{'elname':'${atomtype_log[$x]}${NBonds[$x]}', 'pos':Position([ ${x_logbohr[${x}]}, ${y_logbohr[${x}]}, ${z_logbohr[${x}]}]), 'crg':0., 'pol':Polarizability(iso=1.0)}
EOF

	if [[ $chelpg -ge 1 ]]
		then
		  cat << EOF >> ${name}_chelpg.xyz
{'elname':'${atomtype_log[$x]}${NBonds[$x]}', 'pos':Position([ ${x_logbohr[${x}]}, ${y_logbohr[${x}]}, ${z_logbohr[${x}]}]), 'crg':${qchelpg[${correspondence[$x]}]}, 'pol':Polarizability(iso=1.0)}
EOF

	fi


#EXPORT File, with orientations from cif:

	if [[ $cif_xyz_file != "" ]]
	then
	
	  atomtype_cif_ID=$(echo "${atomtype_cif[$x]}" | tr -d "\n" | od -An -t dC | sed 's/ //g')

	  cat << EOF >> ${name}_cifstruct_mul.xyz
{'elname':'${atomtype_cif[$x]}${NBonds[$x]}', 'pos':Position([ ${x_cifbohr[${x}]}, ${y_cifbohr[${x}]}, ${z_cifbohr[${x}]}]), 'crg':${qmul[${correspondence[$x]}]}, 'pol':Polarizability(iso=1.0)}
EOF

	  cat << EOF >> ${name}_cifstruct_no_charge.xyz
{'elname':'${atomtype_cif[$x]}${NBonds[$x]}', 'pos':Position([ ${x_cifbohr[${x}]}, ${y_cifbohr[${x}]}, ${z_cifbohr[${x}]}]), 'crg':0., 'pol':Polarizability(iso=1.0)}
EOF

	  if [[ $chelpg -ge 1 ]]
		then
		  cat << EOF >> ${name}_cifstruct_chelpg.xyz
{'elname':'${atomtype_cif[$x]}${NBonds[$x]}', 'pos':Position([ ${x_cifbohr[${x}]}, ${y_cifbohr[${x}]}, ${z_cifbohr[${x}]}]), 'crg':${qchelpg[${correspondence[$x]}]}, 'pol':Polarizability(iso=1.0)}
EOF
	  fi
	fi

  x=$(( $x + 1 ))
done

#ADD BONDS TO FILES

if [[ $iso == "FALSE" ]]
then

echo "Determining bond properties"

bond_count=1
while [[ $bond_count -le $bond_count_total ]]
do
	#index of bond_partners for bond # bond_count
	bondpartnera=$(echo "${bond_partners[$bond_count]}" | awk '{printf "%s",$1}' )
	bondpartnerb=$(echo "${bond_partners[$bond_count]}" | awk '{printf "%s",$2}' )

	if [[ $ALKENE == "TRUE" ]]; then
		if [ ${atomtype_log[${bondpartnera}]} == 'C' ] &&  [ ${atomtype_log[${bondpartnerb}]} == 'C' ]; then
			echo "Connectivity from connectivity.dat file overriden as told that this mol is an alkene. Each C-C bond will be assumed to be 1.5, pi"
			bond_type[${bond_count}]='pi'
		fi
	fi

	if [[ $ALKYNE == "TRUE" ]]; then
		if [ ${atomtype_log[${bondpartnera}]} == 'C' ] &&  [ ${atomtype_log[${bondpartnerb}]} == 'C' ]; then
			echo "Connectivity from connectivity.dat file overriden as told that this mol is an alkene. Each C-C bond will be assumed to be 1.5, pi"
			bond_type[${bond_count}]='dpi'
		fi
	fi
	#bond_type ID (for polguess array)
	bond_type_ID=$(echo "${bond_type[$bond_count]}" | tr -d "\n" | od -An -t dC | sed 's/ //g')


#BOND PROPERTIES FOR LOG ORIENT FILE

	#Bond x y z pos'n halfway between atoms
	x_bondpos_log_bohr[$bond_count]=$(echo "( ${x_logbohr[$bondpartnera]} + ${x_logbohr[$bondpartnerb]} ) * 0.5" | bc -l)
	y_bondpos_log_bohr[$bond_count]=$(echo "( ${y_logbohr[$bondpartnera]} + ${y_logbohr[$bondpartnerb]} ) * 0.5" | bc -l)
	z_bondpos_log_bohr[$bond_count]=$(echo "( ${z_logbohr[$bondpartnera]} + ${z_logbohr[$bondpartnerb]} ) * 0.5" | bc -l)

	#Bondlength in x, y, z
	x_bondlength_log_bohr[$bond_count]=$(echo "${x_logbohr[$bondpartnera]} - ${x_logbohr[$bondpartnerb]}" | bc -l)
	y_bondlength_log_bohr[$bond_count]=$(echo "${y_logbohr[$bondpartnera]} - ${y_logbohr[$bondpartnerb]}" | bc -l)
	z_bondlength_log_bohr[$bond_count]=$(echo "${z_logbohr[$bondpartnera]} - ${z_logbohr[$bondpartnerb]}" | bc -l)


	#Magnitude of bondlength
	mag_bondlength_log_bohr[$bond_count]=$( echo "sqrt( ${x_bondlength_log_bohr[$bond_count]}^2 + ${y_bondlength_log_bohr[$bond_count]}^2 + ${z_bondlength_log_bohr[$bond_count]}^2 )" | bc -l )

	#Versor along bond direction
	x_bondversor_log[$bond_count]=$(echo "${x_bondlength_log_bohr[$bond_count]} / ${mag_bondlength_log_bohr[$bond_count]}" | bc -l)
	y_bondversor_log[$bond_count]=$(echo "${y_bondlength_log_bohr[$bond_count]} / ${mag_bondlength_log_bohr[$bond_count]}" | bc -l)
	z_bondversor_log[$bond_count]=$(echo "${z_bondlength_log_bohr[$bond_count]} / ${mag_bondlength_log_bohr[$bond_count]}" | bc -l)

	cat > getbondpol.py << EOF 
#Polarisability directed along x, rotate to direction of bond.
import sys
sys.path.append('../')
import numpy as np
from BasicElements import Rotation

startvec=np.matrix([[1.0, 0.0, 0.0]])
finvec_unnorm=np.matrix([[${x_bondversor_log[${bond_count}]}, ${y_bondversor_log[${bond_count}]}, ${z_bondversor_log[${bond_count}]}]])
finvec=finvec_unnorm/np.linalg.norm(finvec_unnorm)
polstart=1.0*np.matrix([[1.0, 0.0, 0.0],[0.0, 0.001, 0.0],[0.0, 0.0, 0.001]])

rot=Rotation(align=[startvec,finvec])

polrot=rotated_array = rot.T * polstart * rot

for i in [0,1,2]:
	for j in [0,1,2]:
		print polrot[i,j]
EOF

	python getbondpol.py > bondpol.out
	pol_component_count=1
	while read pol_component
	do
		bondpol[${pol_component_count}]=$pol_component
		pol_component_count=$(( $pol_component_count + 1 ))
	done < bondpol.out

	rm getbondpol.py bondpol.out

	# If any diagonals are zero, the determinant is zero, so the matrix is non invertable, so the procedure to calc J matrix fails.
#	if [[ ${bondpol[1]} == 0 ]]; then bondpol[1]="0.0001"; fi
#	if [[ ${bondpol[5]} == 0 ]]; then bondpol[5]="0.0001"; fi
#	if [[ ${bondpol[9]} == 0 ]]; then bondpol[9]="0.0001"; fi


	cat << EOF >> bonds_log.xyz
{'elname':'${bond_type[$bond_count]}', 'pos':Position([ ${x_bondpos_log_bohr[${bond_count}]}, ${y_bondpos_log_bohr[${bond_count}]}, ${z_bondpos_log_bohr[${bond_count}]}]), 'crg':'0.0', 'pol':Polarizability(noniso =[ [${bondpol[1]},${bondpol[2]},${bondpol[3]}], [${bondpol[4]},${bondpol[5]},${bondpol[6]}], [${bondpol[7]}, ${bondpol[8]}, ${bondpol[9]} ]])}
EOF

#IDENTICAL FOR CIF ORIENT FILE
 	if [[ $cif_xyz_file != "" ]]
  	then
	#Bond x y z pos'n halfway between atoms
	x_bondpos_cif_bohr[$bond_count]=$(echo "( ${x_cifbohr[$bondpartnera]} + ${x_cifbohr[$bondpartnerb]} ) * 0.5" | bc -l)
	y_bondpos_cif_bohr[$bond_count]=$(echo "( ${y_cifbohr[$bondpartnera]} + ${y_cifbohr[$bondpartnerb]} ) * 0.5" | bc -l)
	z_bondpos_cif_bohr[$bond_count]=$(echo "( ${z_cifbohr[$bondpartnera]} + ${z_cifbohr[$bondpartnerb]} ) * 0.5" | bc -l)

	#Bondlength in x, y, z
	x_bondlength_cif_bohr[$bond_count]=$(echo "${x_cifbohr[$bondpartnera]} - ${x_cifbohr[$bondpartnerb]}" | bc -l)
	y_bondlength_cif_bohr[$bond_count]=$(echo "${y_cifbohr[$bondpartnera]} - ${y_cifbohr[$bondpartnerb]}" | bc -l)
	z_bondlength_cif_bohr[$bond_count]=$(echo "${z_cifbohr[$bondpartnera]} - ${z_cifbohr[$bondpartnerb]}" | bc -l)


	#Magnitude of bondlength
	mag_bondlength_cif_bohr[$bond_count]=$( echo "sqrt( ${x_bondlength_cif_bohr[$bond_count]}^2 + ${y_bondlength_cif_bohr[$bond_count]}^2 + ${z_bondlength_cif_bohr[$bond_count]}^2 )" | bc -l )

	#Versor along bond direction
	x_bondversor_cif[$bond_count]=$(echo "${x_bondlength_cif_bohr[$bond_count]} / ${mag_bondlength_cif_bohr[$bond_count]}" | bc -l)
	y_bondversor_cif[$bond_count]=$(echo "${y_bondlength_cif_bohr[$bond_count]} / ${mag_bondlength_cif_bohr[$bond_count]}" | bc -l)
	z_bondversor_cif[$bond_count]=$(echo "${z_bondlength_cif_bohr[$bond_count]} / ${mag_bondlength_cif_bohr[$bond_count]}" | bc -l)

	cat > getbondpol.py << EOF 
#Polarisability directed along x, rotate to direction of bond.
import sys
sys.path.append('../')
import numpy as np
from BasicElements import Rotation

startvec=np.matrix([[1.0, 0.0, 0.0]])
finvec_unnorm=np.matrix([[${x_bondversor_cif[${bond_count}]}, ${y_bondversor_cif[${bond_count}]}, ${z_bondversor_cif[${bond_count}]}]])
finvec=finvec_unnorm/np.linalg.norm(finvec_unnorm)
polstart=1.0*np.matrix([[1.0, 0.0, 0.0],[0.0, 0.001, 0.0],[0.0, 0.0, 0.001]])

rot=Rotation(align=[startvec,finvec])

polrot=rotated_array = rot.T * polstart * rot

for i in [0,1,2]:
	for j in [0,1,2]:
		print polrot[i,j]
EOF

	python getbondpol.py > bondpol.out
	pol_component_count=1
	while read pol_component
	do
		bondpol[${pol_component_count}]=$pol_component
		pol_component_count=$(( $pol_component_count + 1 ))
	done < bondpol.out

	rm getbondpol.py bondpol.out

	# If any diagonals are zero, the determinant is zero, so the matrix is non invertable, so the procedure to calc J matrix fails.
#	if [[ ${bondpol[1]} == 0 ]]; then bondpol[1]="0.0001"; fi
#	if [[ ${bondpol[5]} == 0 ]]; then bondpol[5]="0.0001"; fi
#	if [[ ${bondpol[9]} == 0 ]]; then bondpol[9]="0.0001"; fi

	cat << EOF >> bonds_cif.xyz
{'elname':'${bond_type[$bond_count]}', 'pos':Position([ ${x_bondpos_cif_bohr[${bond_count}]}, ${y_bondpos_cif_bohr[${bond_count}]}, ${z_bondpos_cif_bohr[${bond_count}]}]), 'crg':'0.0', 'pol':Polarizability(noniso =[ [${bondpol[1]},${bondpol[2]},${bondpol[3]}], [${bondpol[4]},${bondpol[5]},${bondpol[6]}], [${bondpol[7]}, ${bondpol[8]}, ${bondpol[9]} ]])}
EOF
	fi

bond_count=$(( $bond_count + 1 ))
done

#Put bonds in files

echo "Adding bonds to python .xyz files"

for charge in mul chelpg no_charge
do
  cat bonds_log.xyz >> ${name}_${charge}.xyz
  if [[ $cif_xyz_file != "" ]]
  then
    cat bonds_cif.xyz >> ${name}_cifstruct_${charge}.xyz
  fi
done

rm bonds_log.xyz
rm bonds_cif.xyz

rm coords.tmp

fi

echo "Done!"
done

