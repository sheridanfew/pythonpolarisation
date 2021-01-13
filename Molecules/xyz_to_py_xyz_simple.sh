#Script to convert gaussian output files to molecules in python polarisation model
#
#Sheridan Few

specifiedname="NONE"
BOHR="FALSE"

function USAGE()
{
 cat << EOF
Makes python readable xyz from stndard molecular xyz file. Converts to Bohr automatically.xyz

USAGE:  xyz_to_py_xyz_simple.sh mol.xyz

generates

py_mol.xyz

OPTIONS:
	-b	xyz already in Bohr rather than Angstroms (default $BOHR)		
	-n	Specify name

EOF
}



while getopts "bn:?" Option
do
    case $Option in
	b    ) BOHR="TRUE";;		
        n    ) specifiedname=$OPTARG;;
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
	fi
done

for file in "$@"
do

	if [[ $specifiedname == "NONE" ]]; then basename=$( basename $file | cut -d'.' -f1)
	else basename="$specifiedname"
	fi

	if [ -f py_${basename}.xyz ]
	then
		echo "py_${basename}.xyz already exists. Remove? (y/n)"
		read REMOVE
	
		if [[ $REMOVE = "y" ]]; then rm py_${basename}.xyz
		else exit 0 
		fi
	fi

	#COORDINATES

	echo "Importing coords into arrays"

	x=1
	while read h
	do
	# echo "Recording position of atom ${x} in $file"	
	  posnline="${x}"
	  x[$x]=$(echo $h | awk '{printf "%.7f",$2}')
	  x_bohr[$x]=$(echo "${x[${x}]} / 0.5291772" | bc -l )
	  y[$x]=$(echo $h | awk '{printf "%.7f",$3}')
	  y_bohr[$x]=$(echo "${y[${x}]} / 0.5291772" | bc -l )	
	  z[$x]=$(echo $h | awk '{printf "%.7f",$4}')
	  z_bohr[$x]=$(echo "${z[${x}]} / 0.5291772" | bc -l )
	  atomtype[$x]=$(echo $h | awk '{printf "%s",$1}')
	  x=$(( $x + 1 ))
	done < $file

	NAtoms=$(( $x - 1 ))

	#EXPORT to xyz files readable by python pol model

	#ADD ATOMS TO FILES

	echo "Adding atoms to python .xyz file"
	x=1
	while [[ $x -le $NAtoms ]]
	do
	  if [[ $BOHR == "FALSE" ]]; then
	  cat << EOF >> py_${basename}.xyz
{'elname':'${atomtype[$x]}', 'pos':Position([ ${x_bohr[${x}]}, ${y_bohr[${x}]}, ${z_bohr[${x}]}]), 'crg':0.0, 'pol':Polarizability(iso=1.0)}
EOF
	  else
          cat << EOF >> py_${basename}.xyz
{'elname':'${atomtype[$x]}', 'pos':Position([ ${x[${x}]}, ${y[${x}]}, ${z[${x}]}]), 'crg':0.0, 'pol':Polarizability(iso=1.0)}
EOF

	  fi
	  x=$(( $x + 1 ))
	done
done
