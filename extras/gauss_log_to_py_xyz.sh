function USAGE()
{
 cat << EOF
Converts log files from gaussian to xyzs for molecule factory in James Python Polarisation script.

Converts to Bohr automatically, but can leave in A, and add "'convert_to_bohr':True," to each line like:

{'convert_to_bohr':True,'elname':'C', 'pos':Position([ -3.2281290, 0.1935640, -0.0000590]), 'crg':-0.1441700, 'pol':Polarizability(iso=5)}

USAGE:  gauss_log-to_py_yyz.sh FOO.log BAR
	(exports bar_mul.xyz, bar_chelpg.xyz)

OPTIONS:
	-p	polarisability of atoms (1 by default)


EOF
}

polarisability="1"

while getopts "?" Option
do
    case $Option in
	p    ) polarisability=$OPTARG;;
        ?    ) USAGE
               exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               USAGE   # DEFAULT
    esac
done

shift $(($OPTIND - 1))


	for file in "$1"
		do

		if [ ! -f ${file} ]
			then
                	echo "file $@ does not exist."
			exit 2 

		fi
	done

file=$1
name=$2


                        mol1Natoms=$( grep 'NAtoms=' "${file}" | sed -n 1p | awk '{printf "%1.f",$2}')

#Mulliken charges


	echo "Recording Mulliken Charges & atom types"

                        linestartq=`expr $(grep -n 'Mulliken atomic charges' "${file}" | awk -F":" 'NR==1{print $1}') + 2`
                        lineendq=`expr ${linestartq} + ${mol1Natoms} - 1`

                        x=1
                        while [ $x -le $mol1Natoms ]
                                do
                                  echo "Recording Mulliken charge of atom ${x}"  
                                  chargeline=$(( ${linestartq} + ${x} - 1 ))
                                  qmul1[$x]=$(cat ${file} | sed -n ${chargeline}p | awk '{printf "%.7f",$3}')
				  atomtype[$x]=$(cat ${file} | sed -n ${chargeline}p | awk '{printf "%s",$2}')
                                  echo "Mulliken charge on atom ${x} = ${qmul1[${x}]}"
                                  x=$(( $x + 1 ))
                                done

# Chelpg Charges & atom posns

			echo "Recording Chelpg Charges & atom posns"

			linestartq=`expr $(grep -n 'Fitting point charges to electrostatic potential' "${file}" | awk -F":" 'NR==1{print $1}') + 4`

			lineendq=`expr ${linestart} + ${mol1Natoms} - 1`

			linestartposn=`expr $(grep -n 'Standard orientation:' "${file}" | awk -F":" 'NR==1{print $1}') + 5`

# NB. Used to be linestartposn=`expr $(grep -n 'Z-Matrix orientation' "${file}" | awk -F":" 'NR==1{print $1}') + 5`

		#			echo "mol1Natoms = $mol1Natoms"			

			x=1
				while [ $x -le $mol1Natoms ]
					do
			 		  echo "Recording chelpg charge of atom ${x}"	
					  chargeline=$(( ${linestartq} + ${x} - 1 ))
					  qchelpg1[$x]=$(cat ${file} | sed -n ${chargeline}p | awk '{printf "%.7f",$3}')
					  echo "Chelpg charge on atom ${x} = ${qchelpg1[${x}]}"

					  echo "Recording position of atom ${x}"	
					  posnline=$(( ${linestartposn} + ${x} - 1 ))
					  x1[$x]=$(cat ${file} | sed -n ${posnline}p | awk '{printf "%.7f",$4}')
					  x1bohr[$x]=$(echo "${x1[${x}]} * 0.5291772" | bc -l )

					  y1[$x]=$(cat ${file} | sed -n ${posnline}p | awk '{printf "%.7f",$5}')
					  y1bohr[$x]=$(echo "${y1[${x}]} * 0.5291772" | bc -l )
	
					  z1[$x]=$(cat ${file} | sed -n ${posnline}p | awk '{printf "%.7f",$6}')
					  z1bohr[$x]=$(echo "${z1[${x}]} * 0.5291772" | bc -l )
					  echo "position of atom ${x} = ${x1[${x}]},${y1[${x}]},${z1[${x}]}"
#ROTATION OPTION
#					  x1[$x]=${x1[${x}]}
#					  y1[$x]=${y1[${x}]}
#					  z1[$x]=${z1[${x}]}
						 

		#			  echo -en "${charge_ring[${x}]}\t" >> coulomb.csv
				  	  x=$(( $x + 1 ))
					done
#Export:

			x=1
				while [ $x -le $mol1Natoms ]
					do
					  cat << EOF >> ${name}_mul.xyz
{'elname':'${atomtype[$x]}', 'pos':Position([ ${x1bohr[${x}]}, ${y1bohr[${x}]}, ${z1bohr[${x}]}]), 'crg':${qmul1[${x}]}, 'pol':Polarizability(iso=${polarisability})}
EOF
					  cat << EOF >> ${name}_chelpg.xyz
{'elname':'${atomtype[$x]}', 'pos':Position([ ${x1bohr[${x}]}, ${y1bohr[${x}]}, ${z1bohr[${x}]}]), 'crg':${qmul1[${x}]}, 'pol':Polarizability(iso=${polarisability})}
EOF
					  cat << EOF >> ${name}_no_charge.xyz
{'elname':'${atomtype[$x]}', 'pos':Position([ ${x1bohr[${x}]}, ${y1bohr[${x}]}, ${z1bohr[${x}]}]), 'crg':0., 'pol':Polarizability(iso=${polarisability})}
EOF
				  	  x=$(( $x + 1 ))
					done


