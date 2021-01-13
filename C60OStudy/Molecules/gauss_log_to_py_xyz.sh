file=""
logfile=""
cif_nos_log_orient_file=""
connectivity_file=""

polguess[$(echo "C" | tr -d "\n" | od -An -t dC | sed 's/ //g')]="9.0799"
polguess[$(echo "H" | tr -d "\n" | od -An -t dC | sed 's/ //g')]="0.001"
polguess[$(echo "N" | tr -d "\n" | od -An -t dC | sed 's/ //g')]="9.855"
polguess[$(echo "O" | tr -d "\n" | od -An -t dC | sed 's/ //g')]="5.160"
polguess[$(echo "S" | tr -d "\n" | od -An -t dC | sed 's/ //g')]="15.4061"
polguess[$(echo "Si" | tr -d "\n" | od -An -t dC | sed 's/ //g')]="10.0"
cutset="5.0"



function USAGE()
{
 cat << EOF
Converts log files from gaussian to xyzs for molecule factory in James Python Polarisation script.

Takes as input a gaussian log file for charges and polarisability, and an xyz with molecular orientation as in the cif file. Exports two xyz files, one with the orientation of the gaussian log (for calibrating atomic polarisabilities to molecular polarisability from gaussian) and one with the orientation of the xyz_with_cif_structure.xyz, but atomic charges from the gaussian file.

Can also specify a third file if atom numbers of these two do not match up (see option -a)

Converts to Bohr automatically, but can leave in A, and add "'convert_to_bohr':True," to each line like:

{'convert_to_bohr':True,'elname':'C', 'pos':Position([ -3.2281290, 0.1935640, -0.0000590]), 'crg':-0.1441700, 'pol':Polarizability(iso=5)}

USAGE:  gauss_log-to_py_xyz.sh -b xyz_with_cif_structure.xyz -c 4.0 -h 1.0 -y connectivity.dat gaussfile.log
	(exports blah_mul.xyz, blah_${poptype}.xyz, blah_no_charge.py pol_blah.py)

OPTIONS:
	-a			(Only required if extracting atomic charges from a log with different atom numbering from cif structure. If not specified, then script assumes numbering of xyz from cif file for structure, and gaussian .log are ther same) xyz w. atom nos corresponding to cif file, but orientation as in log file.
	-b 			xyz file with final coords (as in original cif file, for making crystals)
	-c,h,i(Si),n,s,o	Guess of polarisability of each of these elements, goes into minuit fit jobs (in Singlemol_pol_determination) and polarisability in xyz file.		
	-l			(redundant here) .log file to use for charges (if different from main file, used for cryst structure)
	-y			Connectivity data file (can get by getting gjf from cif, then opening&resaving as neutral (0 1) in gaussview, extract connectivity lines at end.) All atoms of same element will be treated as identical if this option is not specified.)
	-z			Set cutoff value

DEFAULTS:
cpolguess=${polguess[$(echo "C" | tr -d "\n" | od -An -t dC | sed 's/ //g')]} au
Hpolguess=${polguess[$(echo "H" | tr -d "\n" | od -An -t dC | sed 's/ //g')]} au
Npolguess=${polguess[$(echo "N" | tr -d "\n" | od -An -t dC | sed 's/ //g')]} au
Opolguess=${polguess[$(echo "O" | tr -d "\n" | od -An -t dC | sed 's/ //g')]} au
Spolguess=${polguess[$(echo "S" | tr -d "\n" | od -An -t dC | sed 's/ //g')]} au
cutset=$cutset Bohr

EOF
}



while getopts "a:b:c:h:i:l:n:o:s:t:y:z:?" Option
do
    case $Option in
	a    ) cif_nos_log_orient_file=$OPTARG;;
	b    ) cif_xyz_file=$OPTARG;;
	c    ) polguess[$(echo "C" | tr -d "\n" | od -An -t dC | sed 's/ //g')]=$OPTARG;;
	h    ) polguess[$(echo "H" | tr -d "\n" | od -An -t dC | sed 's/ //g')]=$OPTARG;;
	i    ) polguess[$(echo "Si" | tr -d "\n" | od -An -t dC | sed 's/ //g')]=$OPTARG;;
	l    ) logfile=$OPTARG;;
	n    ) polguess[$(echo "N" | tr -d "\n" | od -An -t dC | sed 's/ //g')]=$OPTARG;;
	o    ) polguess[$(echo "O" | tr -d "\n" | od -An -t dC | sed 's/ //g')]=$OPTARG;;
	s    ) polguess[$(echo "S" | tr -d "\n" | od -An -t dC | sed 's/ //g')]=$OPTARG;;
	y    ) connectivity_file=$OPTARG;;
	z    ) cutset=$OPTARG;;
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
	elif [[ ! -f ${logfile} && $logfile != "" ]]
	then
               	echo "file $logfile does not exist."
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
  if [[ $cif_nos_log_orient_file == "" ]]
  then
    name=$(echo "${file}" | sed 's/.log//g' | sed 's/.xyz/_/g' | sed 's/-/_/g' )
  else
    name=$(echo "${cif_xyz_file}_${cif_nos_log_orient_file}_${file}" | sed 's/.log//g' | sed 's/.xyz/_/g' | sed 's/-/_/g' )
  fi

#  if [[ $connectivity_file != "" ]]
#  then
#    name=$(echo "${name}_w_connectivity")
#  fi

poptype=""
if [[ $(grep -c -i 'pop=chelpg' $file) -ge 1 ]]
then
	poptype="chelpg"
elif [[ $(grep -c -i 'pop=ESP' $file) -ge 1 ]]
then
	poptype="ESP"
fi


name_for_file="$name_$(date +"%Y_%m_%d")"
if [[ -f ${name}_mul.xyz || -f ${name}_${poptype}.xyz || -f ${name}_no_charge.xyz || -f ${name}_cifstuct_mul.xyz || -f ${name}_cifstuct_${poptype}.xyz || -f ${name}_cifstuct_no_charge.xyz ]]
then
    echo "xyz of this name, ${name}_(mul/${poptype}/no_charge).xyz already exists. Exiting."
    exit 2
fi

#get coords
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

#connectivity (/no. of bonds)
if [[ $connectivity_file != "" ]]
then
	connectivity_count=1 
	while read line
	do
	  Nbondsnew=$(( ( $(echo $line | wc -w) - 1 ) / 2))
	  echo "Nbondsnew($connectivity_count)=$Nbondsnew"
	  NBonds[$connectivity_count]=$(( ${NBonds[${connectivity_count}]} + $Nbondsnew ))
	  NBonds[$(echo $line | awk '{printf "%s",$2}')]=$(( ${NBonds[$(echo $line | awk '{printf "%s",$2}')]} + 1 ))
	  NBonds[$(echo $line | awk '{printf "%s",$4}')]=$(( ${NBonds[$(echo $line | awk '{printf "%s",$4}')]} + 1 ))
	  NBonds[$(echo $line | awk '{printf "%s",$6}')]=$(( ${NBonds[$(echo $line | awk '{printf "%s",$6}')]} + 1 ))
	  NBonds[$(echo $line | awk '{printf "%s",$8}')]=$(( ${NBonds[$(echo $line | awk '{printf "%s",$8}')]} + 1 ))
	  #If there aren't entries for $4,6,8 contributions will all go to array member NBonds[0], not called in this script so no problem)
	  connectivity_count=$(( $connectivity_count + 1 ))
	done < $connectivity_file
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

#Atom pos'ns from coords.tmp, cif_xyz_file, cif_nos_log_orient_file to arrays

x=1
while [ $x -le $Natoms ]
do
  echo "Recording position of atom ${x} in $file"	
  posnline="${x}"
  x1[$x]=$(cat coords.tmp | sed -n ${posnline}p | awk '{printf "%.7f",$2}')
  x1bohr[$x]=$(echo "${x1[${x}]} / 0.5291772" | bc -l )
  y1[$x]=$(cat coords.tmp | sed -n ${posnline}p | awk '{printf "%.7f",$3}')
  y1bohr[$x]=$(echo "${y1[${x}]} / 0.5291772" | bc -l )	
  z1[$x]=$(cat coords.tmp | sed -n ${posnline}p | awk '{printf "%.7f",$4}')
  z1bohr[$x]=$(echo "${z1[${x}]} / 0.5291772" | bc -l )
  atomtype1[$x]=$(cat coords.tmp | sed -n ${posnline}p | awk '{printf "%s",$1}')
  echo "position of atom ${x} = ${x1[${x}]},${y1[${x}]},${z1[${x}]}"

  if [[ $cif_xyz_file != "" ]]
  then
    echo "Recording position of atom ${x} in $cif_xyz_file"
    x3[$x]=$(cat $cif_xyz_file | sed -n ${x}p | awk '{printf "%.7f",$2}')
    x3bohr[$x]=$(echo "${x3[${x}]} / 0.5291772" | bc -l )
    y3[$x]=$(cat $cif_xyz_file | sed -n ${x}p | awk '{printf "%.7f",$3}')
    y3bohr[$x]=$(echo "${y3[${x}]} / 0.5291772" | bc -l )	
    z3[$x]=$(cat $cif_xyz_file | sed -n ${x}p | awk '{printf "%.7f",$4}')
    z3bohr[$x]=$(echo "${z3[${x}]} / 0.5291772" | bc -l )
    atomtype3[$x]=$(cat $cif_xyz_file | sed -n ${posnline}p | awk '{printf "%s",$1}')
    echo "position of atom ${x} = ${x3[${x}]},${y3[${x}]},${z3[${x}]}"
  fi
    

  if [[ $cif_nos_log_orient_file != "" ]]
  then
    echo "Recording position of atom ${x} in $cif_nos_log_orient_file"
    x2[$x]=$(cat $cif_nos_log_orient_file | sed -n ${x}p | awk '{printf "%.7f",$2}')
#   x2bohr[$x]=$(echo "${x2[${x}]} / 0.5291772" | bc -l )
    y2[$x]=$(cat $cif_nos_log_orient_file | sed -n ${x}p | awk '{printf "%.7f",$3}')
#   y2bohr[$x]=$(echo "${y2[${x}]} / 0.5291772" | bc -l )	
    z2[$x]=$(cat $cif_nos_log_orient_file | sed -n ${x}p | awk '{printf "%.7f",$4}')
#   z2bohr[$x]=$(echo "${z2[${x}]} / 0.5291772" | bc -l )
    atomtype2[$x]=$(cat $cif_nos_log_orient_file | sed -n ${posnline}p | awk '{printf "%s",$1}')
    echo "position of atom ${x} = ${x2[${x}]},${y2[${x}]},${z2[${x}]}"
  fi

  x=$(( $x + 1 ))
done

#Only need to worry about if using cif_nos_log_orient_file

if [[ $cif_nos_log_orient_file == "" ]]
then
  i=1
  while [ $i -le $Natoms ]
  do
    correspondence[$i]="$i"
    i=$(( $i + 1 ))
  done
else
  # Work out correspondence between atom nos in cif and gausslog files (necessary to assign pointcharges to cif structure if gaussian log atom nos differ from those in cif structure file)
  i=1
  while [ $i -le $Natoms ]
  do
    mindistabs="100000"
    j=1
    while [ $j -le $Natoms ]
    do
      distabsij=$(echo "( ((${x1[$i]})-(${x2[$j]}))*((${x1[$i]})-(${x2[$j]})) ) + ( ((${y1[$i]})-(${y2[$j]}))*((${y1[$i]})-(${y2[$j]})) ) + ( ((${z1[$i]})-(${z2[$j]}))*((${z1[$i]})-(${z2[$j]})) )" | bc -l)

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
    echo "${atomtype2[$i]} ${atomtype1[${correspondence[$i]}]}"
    if [[ "${atomtype2[$i]}" != "${atomtype1[${correspondence[$i]}]}" ]]
    then
      echo "correspondence wrong! Make sure orientation of cif_nos_log_orient_file is same as gaussian log (and maybe that this part of script works...!)"
      exit 2
    fi
    i=$(( $i + 1 ))
  done
fi




#Mulliken charges

echo "Recording Mulliken Charges & atom types"

if [[ $logfile = "" ]]
then
  chargesfile=$file
else
  chargesfile=$logfile
fi

linestartq=`expr $(grep -n 'Mulliken atomic charges:' "${chargesfile}" | tail -n 1 | awk -F":" 'NR==1{print $1}') + 2`
lineendq=`expr ${linestartq} + ${Natoms} - 1`

x=1
while [ $x -le $Natoms ]
do
  echo "Recording Mulliken charge of atom ${x}"  
  chargeline=$(( ${linestartq} + ${x} - 1 ))
  qmullog[$x]=$(cat ${chargesfile} | sed -n ${chargeline}p | awk '{printf "%.7f",$3}')
#  atomtype1[$x]=$(cat ${chargesfile} | sed -n ${chargeline}p | awk '{printf "%s",$2}')
  echo "Mulliken charge on atom ${x} = ${qmullog[${x}]}"
  x=$(( $x + 1 ))
done

# ${poptype} Charges

echo "Recording ${poptype} Charges & atom posns"

if [[ $poptype != "" ]]
then
  if [[ $poptype == "chelpg" ]]
  then
    linestartq=`expr $(grep -n 'Fitting point charges to electrostatic potential' "${chargesfile}" | awk -F":" 'NR==1{print $1}') + 4`
  elif [[ $poptype == "ESP" ]]
  then
    linestartq=`expr $(grep -n 'Fitting point charges to electrostatic potential' "${chargesfile}" | awk -F":" 'NR==1{print $1}') + 4`
  fi

  lineendq=`expr ${linestart} + ${Natoms} - 1`

  x=1
  while [ $x -le $Natoms ]
  do
   		  echo "Recording ${poptype} charge of atom ${x}"	
  		  chargeline=$(( ${linestartq} + ${x} - 1 ))
  		  qpoplog[$x]=$(cat ${chargesfile} | sed -n ${chargeline}p | awk '{printf "%.7f",$3}')
  		  echo "${poptype} charge on atom ${x} = ${qpoplog[${x}]}"
    	x=$(( $x + 1 ))
  done
fi

#Export to xyz file readable by python, with orientations from gaussian log:

echo "correspondence[1]=${correspondence[1]}"
echo "qmullog[${correspondence[1]}]=${qmullog[${correspondence[1]}]}"

			x=1
				while [ $x -le $Natoms ]
				  do
					  cat << EOF >> ${name}_mul.xyz
{'elname':'${atomtype1[$x]}${NBonds[$x]}', 'pos':Position([ ${x1bohr[${x}]}, ${y1bohr[${x}]}, ${z1bohr[${x}]}]), 'crg':${qmullog[${correspondence[$x]}]}, 'pol':Polarizability(iso=${polguess[$(echo "${atomtype1[$x]}" | tr -d "\n" | od -An -t dC | sed 's/ //g')]})}
EOF
					if [[ $poptype != "" ]]
						then
						  cat << EOF >> ${name}_${poptype}.xyz
{'elname':'${atomtype1[$x]}${NBonds[$x]}', 'pos':Position([ ${x1bohr[${x}]}, ${y1bohr[${x}]}, ${z1bohr[${x}]}]), 'crg':${qpoplog[${correspondence[$x]}]}, 'pol':Polarizability(iso=${polguess[$(echo "${atomtype1[$x]}" | tr -d "\n" | od -An -t dC | sed 's/ //g')]})}
EOF
					fi
					  cat << EOF >> ${name}_no_charge.xyz
{'elname':'${atomtype1[$x]}${NBonds[$x]}', 'pos':Position([ ${x1bohr[${x}]}, ${y1bohr[${x}]}, ${z1bohr[${x}]}]), 'crg':0., 'pol':Polarizability(iso=${polguess[$(echo "${atomtype1[$x]}" | tr -d "\n" | od -An -t dC | sed 's/ //g')]})}
EOF

#Export to xyz file readable by python, with orientations from cif:

					if [[ $cif_xyz_file != "" ]]
					then
					  cat << EOF >> ${name}_cifstruct_mul.xyz
{'elname':'${atomtype3[$x]}${NBonds[$x]}', 'pos':Position([ ${x3bohr[${x}]}, ${y3bohr[${x}]}, ${z3bohr[${x}]}]), 'crg':${qmullog[${correspondence[$x]}]}, 'pol':Polarizability(iso=${polguess[$(echo "${atomtype3[$x]}" | tr -d "\n" | od -An -t dC | sed 's/ //g')]})}
EOF
					  if [[ $poptype != ""  ]]
						then
						  cat << EOF >> ${name}_cifstruct_${poptype}.xyz
{'elname':'${atomtype3[$x]}${NBonds[$x]}', 'pos':Position([ ${x3bohr[${x}]}, ${y3bohr[${x}]}, ${z3bohr[${x}]}]), 'crg':${qpoplog[${correspondence[$x]}]}, 'pol':Polarizability(iso=${polguess[$(echo "${atomtype3[$x]}" | tr -d "\n" | od -An -t dC | sed 's/ //g')]})}
EOF
					  fi
					  cat << EOF >> ${name}_cifstruct_no_charge.xyz
{'elname':'${atomtype3[$x]}${NBonds[$x]}', 'pos':Position([ ${x3bohr[${x}]}, ${y3bohr[${x}]}, ${z3bohr[${x}]}]), 'crg':0., 'pol':Polarizability(iso=${polguess[$(echo "${atomtype3[$x]}" | tr -d "\n" | od -An -t dC | sed 's/ //g')]})}
EOF
					fi

				  	x=$(( $x + 1 ))
				  done

#Options for minuit fit to (determine appropriate values for atomic polarisabilities to reproduce molecular polarisabilities from gaussian log)

minuitopts=""
minuitopts_wscan=""
minuitscanopts=""
etacalcopts=""
etaopts=""
mvaluesout="m.fval, m.ncalls,m.values['cut'],Uqq,Uqd,Udd,energyev,"

: <<'EXAMPLE'

https://code.google.com/p/pyminuit/wiki/GettingStartedGuide

>>> m = minuit.Minuit(lambda x, y: (x-1)**2 + (y-2)**2, x=3, y=4)
>>> m.scan(("x", 30, -3, 7), ("y", 30, -3, 7), output=False)
>>> m.values
{'y': 2.1666666666666661, 'x': 1.1666666666666663}
>>> m.migrad()
>>> m.values
{'y': 2.0000000000041425, 'x': 1.0000000000042153}

EXAMPLE


for atomtype in C H N O S Se Si C1 C2 C3 C4 H1 N1 N2 N3 N4 O1 O2 O3 O4 S1 S2 S3 S4 Se1 Se2 Se3 Se4 Si1 Si2 Si3 Si4
do
#Sets fitting criteria for individual atom types (could add more to this list easily)
  if [[ $(grep -c "elname'\:'$atomtype'" ${name}_no_charge.xyz) -ge 1 ]]
  then
	atomtype[${atomtypecount}]="$atomtype"
	
	etacalcopts="${etacalcopts},${atomtype}"
	etaopts="${etaopts},${atomtype}=${atomtype}"
	mvaluesout="${mvaluesout}m.values['${atomtype}'],"
	if [[ $atomtype == C[1-4] || $atomtype == "C" ]]
	then
	  echo "${atomtype}pol: ${polguess[$(echo "${atomtype}" | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]}, C atom no: $(echo "${atomtype}" | tr -d "\n" | od -An -t dC | sed 's/ //g')"
	  minuitopts="${minuitopts}${atomtype}=${polguess[$(echo "${atomtype}" | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]},limit_${atomtype}=(1,14),"
	  minuitopts_wscan="${minuitopts_wscan}limit_${atomtype}=(1,14),"
	  minuitscanopts="${minuitscanopts}(\"${atomtype}\",5,1,14),"
#	  minuitopts="${minuitopts}${atomtype}=${polguess[$(echo "${atomtype}" | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]},fix_${atomtype}=True,"
#	  minuitopts_wscan="${minuitopts_wscan}${atomtype}=${polguess[$(echo "${atomtype}" | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]},fix_${atomtype}=True,"	
	elif [[ $atomtype == H[1-4] || $atomtype == "H" ]]
	then
	  minuitopts="${minuitopts}${atomtype}=${polguess[$(echo "${atomtype}" | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]},fix_${atomtype}=True,"
	  minuitopts_wscan="${minuitopts_wscan}${atomtype}=${polguess[$(echo "${atomtype}"  | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]},fix_${atomtype}=True,"	
	elif [[ $atomtype == N[1-4] || $atomtype == "N" ]]
	then
	  minuitopts="${minuitopts}${atomtype}=${polguess[$(echo "${atomtype}" | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]},limit_${atomtype}=(1,20),"
	  minuitopts_wscan="${minuitopts_wscan}limit_${atomtype}=(1,20),"
	  minuitscanopts="${minuitscanopts}(\"${atomtype}\",5,4,20),"	
	elif [[ $atomtype == O[1-4] || $atomtype == "O" ]]
	then
	  minuitopts="${minuitopts}${atomtype}=${polguess[$(echo "${atomtype}" | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]},limit_${atomtype}=(1,30),"
	  minuitopts_wscan="${minuitopts_wscan}limit_${atomtype}=(1,30),"
	  minuitscanopts="${minuitscanopts}(\"${atomtype}\",7,2,24),"	
	elif [[ $atomtype == S[1-4] || $atomtype == "S" ]]
	then
	  minuitopts="${minuitopts}${atomtype}=${polguess[$(echo "${atomtype}" | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]},limit_${atomtype}=(1,30),"
	  minuitopts_wscan="${minuitopts_wscan}limit_${atomtype}=(1,30),"
	  minuitscanopts="${minuitscanopts}(\"${atomtype}\",7,4,24),"	
	elif [[ $atomtype == Si[1-4] || $atomtype == "Si" ]]
	then
	  minuitopts="${minuitopts}${atomtype}=${polguess[$(echo "${atomtype}" | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]},limit_${atomtype}=(1,30),"
	  minuitopts_wscan="${minuitopts_wscan}limit_${atomtype}=(1,30),"
	  minuitscanopts="${minuitscanopts}(\"${atomtype}\",7,4,24),"	
	else
	  echo "WARNING: Unrecognised atom type $atomtype, assuming limits for pol. calibration as for carbon. Can edit script to add faculties for dealing with ${atomtype}."
	  echo "${atomtype}pol: ${polguess[$(echo "${atomtype}" | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]}, atom no: $(echo "${atomtype}" | tr -d "\n" | od -An -t dC | sed 's/ //g')"
	  minuitopts="${minuitopts}${atomtype}=${polguess[$(echo "${atomtype}" | sed 's/[0-9]*//g' | tr -d "\n" | od -An -t dC | sed 's/ //g')]},limit_${atomtype}=(1,12),"
	  minuitopts_wscan="${minuitopts_wscan}limit_${atomtype}=(1,12),"
	  minuitscanopts="${minuitscanopts}(\"${atomtype}\",5,1,12),"
	fi
	atomtypecount=$(( $atomtypecount + 1 ))
  else
	mvaluesout="${mvaluesout}'N/A ',"
  fi
done

natomtypes=$atomtypecount-1

#Make etafit python file

cat << EOF > ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py
import sys
sys.path.append('../')
from BasicElements import Molecule,Position, Polarizability
from BasicElements.MoleculeFactory import ReadMoleculeType
from BasicElements.MoleculeFactory import GetMolecule
from Polarizability.GetDipoles import get_dipoles,split_dipoles_onto_atoms
from BasicElements import *
from Polarizability.etacalc import *
from Polarizability.GetEnergy import *
import numpy as np
#from scipy.optimize import minimize
import minuit

from Molecules.pol_${name} import pol
pol_${name}=pol

def etacalc( cut$etacalcopts):
	etacalc=eta('../Molecules/${name}_no_charge.xyz',pol_${name}, cut$etaopts)
	return etacalc

m=minuit.Minuit(etacalc, ${minuitopts} cut=${cutset}, fix_cut=True)
m.printMode = 1
m.migrad()

EOF

#Make modify polarisability function in order to output energy of mol resulting from chosen atomic polarisability values (make sure nothing crazy like +ve overall energies which can result from too small a cutoff/too large polarisabilities too close. See Appendix to Pc polarisability Castet paper/Cleaven's work)

echo -en 'def ModifyPolarizability(molecule, ' >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py

atomtypecount=1
while [[ $atomtypecount -lt $natomtypes ]]; do echo -en "${atomtype[${atomtypecount}]}," >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py; atomtypecount=$(( $atomtypecount + 1 )); done

echo -en "${atomtype[${natomtypes}]}):\n" >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py
echo -en  "\t\"\"\" takes a mol and changes the isotropic polarizability\n\tfor each atom \"\"\"\n" >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py

atomtypecount=1
while [[ $atomtypecount -le $natomtypes ]]; do echo -en "\t${atomtype[${atomtypecount}]} = Polarizability(iso=${atomtype[${atomtypecount}]})\n" >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py; atomtypecount=$(( $atomtypecount + 1 )); done

echo -en "\tfor atom in molecule:\n" >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py

atomtypecount=1
while [[ $atomtypecount -le $natomtypes ]]; do echo -en "\t\tif (atom()._elname==\"${atomtype[${atomtypecount}]}\"):\n\t\t\tatom()._pol= ${atomtype[${atomtypecount}]}\n" >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py; atomtypecount=$(( $atomtypecount + 1 )); done

cat << EOF >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py

ReadMoleculeType('../Molecules/${name}_no_charge.xyz')
mol = GetMolecule('../Molecules/${name}_no_charge.xyz')
 
EOF

echo -en 'ModifyPolarizability(mol(),' >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py

atomtypecount=1
while [[ $atomtypecount -lt $natomtypes ]]; do echo -en "${atomtype[${atomtypecount}]}=m.values['${atomtype[${atomtypecount}]}']," >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py; atomtypecount=$(( $atomtypecount + 1 )); done
echo -en
echo -en "${atomtype[${natomtypes}]}=m.values['${atomtype[${natomtypes}]}'])" >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py



cat << EOF >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py

# Get energy, dips for 1,0,0 field

E0 = np.matrix([1.,0.,0.])
d = get_dipoles(E0=E0,cutoff=${cutset})
split_d = split_dipoles_onto_atoms(d)
tot = np.matrix([0.,0.,0.])
for dd in split_d:
	tot += dd
Uqq = np.multiply(get_U_qq(),27.211)
Uqd = np.multiply(get_U_qdip(),27.211)
Udd = np.multiply(get_U_dipdip(),27.211)
energyev = Uqq+Uqd+Udd
print 'energyev(E100), Uqq, Uqd, Udd'
print energyev,Uqq,Uqd,Udd
EOF

cat << EOF >> ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py

# Output data
sys.stdout = open('results/etafit_minuit_results.out', 'a')
print '\n${name}',$mvaluesout
sys.stdout = sys.__stdout__ 

EOF

#Make etafit w. scan

cp ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py ../Singlemol_Pol_Detemination/${name}_etafit_minuit_wscan.py
sed -i "s/m=minuit.Minuit(etacalc, ${minuitopts} cut=${cutset}, fix_cut=True)/m=minuit.Minuit(etacalc, ${minuitopts_wscan} cut=${cutset}, fix_cut=True)/" ../Singlemol_Pol_Detemination/${name}_etafit_minuit_wscan.py

#Make fits w. sum of abs_sq_dif in each dir rather than eta (favours getting larger components in mol. polarisability tensor right above small components)
cat ../Singlemol_Pol_Detemination/${name}_etafit_minuit.py | sed "s/eta('/abs_square_dif_fit('/g" | sed "s/etafit/abs_sq_dif_fit/g" > ../Singlemol_Pol_Detemination/${name}_abs_sq_dif_fit_minuit.py
cat ../Singlemol_Pol_Detemination/${name}_etafit_minuit_wscan.py | sed "s/eta('/abs_square_dif_fit('/g" | sed "s/etafit/abs_sq_dif_fit/g" > ../Singlemol_Pol_Detemination/${name}_abs_sq_dif_fit_minuit_wscan.py

rm coords.tmp

done

