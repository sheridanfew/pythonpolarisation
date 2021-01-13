arrowsizefactor="100"
xmin="-40"
xmax="40"
ymin="-40"
ymax="40"
zmin="-40"
zmax="40"

for file in $@
do
	NAtoms=`expr $(grep -n 'Molecule Ends.' "${file}" | awk -F":" 'NR==1{print $1}') - 2`
	linestartdip=`expr $(grep -n 'dipoles start.' "${file}" | awk -F":" 'NR==1{print $1}') + 1`
	base=$(basename $file | sed 's/\.dat//g')

	cat << EOF > dipoles_100_${base}.gnp
#set xrange [$xmin:$xmax]
#set yrange [$ymin:$ymax]
#set zrange [$zmin:$zmax]

#unset border
#unset ztics
#unset xtics
#unset ytics

EOF

	x=1
	while [[ $x -le $NAtoms ]]
		do
		echo "Recording dipole of atom ${x}"  
	 	dipline=$(( ${linestartdip} + ${x} - 1 ))

		xdip[$x]=$(cat ${file} | sed -n ${dipline}p |  sed 's/\[\[//g' |  sed 's/\]\]//g' | awk '{printf "%.7f",$1}')
		ydip[$x]=$(cat ${file} | sed -n ${dipline}p |  sed 's/\[\[//g' |  sed 's/\]\]//g' | awk '{printf "%.7f",$2}')
		zdip[$x]=$(cat ${file} | sed -n ${dipline}p |  sed 's/\[\[//g' |  sed 's/\]\]//g' | awk '{printf "%.7f",$3}')

		echo "dipole on atom ${x} = ${xdip[${x}]},${ydip[${x}]},${zdip[${x}]}"

		echo "Recording position of atom ${x}"        
		posnline=$(( ${x} + 1 ))

		atomtype[$x]=$(cat ${file} | sed -n ${posnline}p | awk '{printf "%s",$2}')
	
		if [[ ${atomtype[${x}]} == "C" ]]
		then
			col[${x}]="1"
		elif [[ ${atomtype[${x}]} == "S" ]]
		then
			col[${x}]="2"
		elif [[ ${atomtype[${x}]} == "H" ]]
		then
			col[${x}]="3"
		elif [[ ${atomtype[${x}]} == "O" ]]
		then
			col[${x}]="4"
		fi

		x[$x]=$(cat ${file} | sed -n ${posnline}p |  sed 's/\[\[//g' |  sed 's/\]\]//g' | awk '{printf "%.7f",$4}')
		y[$x]=$(cat ${file} | sed -n ${posnline}p |  sed 's/\[\[//g' |  sed 's/\]\]//g' | awk '{printf "%.7f",$5}')
		z[$x]=$(cat ${file} | sed -n ${posnline}p |  sed 's/\[\[//g' |  sed 's/\]\]//g' | awk '{printf "%.7f",$6}')

		echo "position of atom ${x} (${atomtype[${x}]}) = ${x[${x}]},${y[${x}]},${z[${x}]}"

		arrowstartx[$x]=$( echo "${x[${x}]} - (${xdip[${x}]}*${arrowsizefactor}/2)" | bc -l)
		arrowendx[$x]=$( echo "${x[${x}]} + (${xdip[${x}]}*${arrowsizefactor}/2)" | bc -l)
		arrowstarty[$x]=$( echo "${y[${x}]} - (${ydip[${x}]}*${arrowsizefactor}/2)" | bc -l)
		arrowendy[$x]=$( echo "${y[${x}]} + (${ydip[${x}]}*${arrowsizefactor}/2)" | bc -l)
		arrowstartz[$x]=$( echo "${z[${x}]} - (${zdip[${x}]}*${arrowsizefactor}/2)" | bc -l)
		arrowendz[$x]=$( echo "${z[${x}]} + (${zdip[${x}]}*${arrowsizefactor}/2)" | bc -l)

		echo "set arrow ${x} from ${arrowstartx[${x}]},${arrowstarty[${x}]},${arrowstartz[${x}]} to ${arrowendx[${x}]},${arrowendy[${x}]},${arrowendz[${x}]}  head filled size screen 0.01,16,20 lw 3 lc ${col[${x}]}" >> dipoles_100_${base}.gnp


head filled size screen 0.01,16,20 lw 3 lc 1


		x=$(( $x + 1 ))
	done

	echo -en "splot \"<echo '${x[1]} ${y[1]} ${z[1]}'\" with points ls 1, " >> dipoles_100_${base}.gnp

	x=2
	while [[ $x -le $(( $NAtoms - 1)) ]]
	do
		echo -en "\"<echo '${x[${x}]} ${y[${x}]} ${z[${x}]}'\" with points ls 1, " >> dipoles_100_${base}.gnp
		x=$(( $x + 1 ))
	done

	echo "\"<echo '${x[${NAtoms}]} ${y[${NAtoms}]} ${z[${NAtoms}]}'\" with points ls 1" >> dipoles_100_${base}.gnp

cat << EOF >> dipoles_100_${base}.gnp
set view equal xyz
unset key
replot
pause -1
EOF

done

