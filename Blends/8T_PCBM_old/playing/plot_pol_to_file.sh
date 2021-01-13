arrowsizefactor="100"
xmin="-40"
xmax="40"
ymin="-10"
ymax="10"
zmin="0"
zmax="40"

#Currently set up to ignore H atoms

for file in $@
do
	NAtoms=`expr $(grep -n 'Molecule Ends.' "${file}" | awk -F":" 'NR==1{print $1}') - 2`
	linestartdip=`expr $(grep -n 'dipoles start.' "${file}" | awk -F":" 'NR==1{print $1}') + 1`
	base=$(basename $file | sed 's/\.dat//g')

	cat << EOF > dipoles_100_${base}.gnp
set xrange [$xmin:$xmax]
set yrange [$ymin:$ymax]
set zrange [$zmin:$zmax]

unset border
unset ztics
unset xtics
unset ytics

unset key

set samples 600

set term png size 5000, 5000

set output "${base}.png"

set style line 1 lt 1 lw 5 pt 3 linecolor rgb "white"

set view equal xyz
unset key

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

		if [[ ${atomtype[${x}]} != "H" ]]; then
		  echo "set arrow ${x} from ${arrowstartx[${x}]},${arrowstarty[${x}]},${arrowstartz[${x}]} to ${arrowendx[${x}]},${arrowendy[${x}]},${arrowendz[${x}]}  head filled size screen 0.01,16,20 lw 3 lc ${col[${x}]}" >> dipoles_100_${base}.gnp
		fi

		x=$(( $x + 1 ))
	done


cat << EOF >> dipoles_100_${base}.gnp

splot "<echo '$xmin $ymin $zmin'" with points ls 1, "<echo '$xmax $ymin $zmin'" with points ls 1, "<echo '$xmin $ymax $zmin'" with points ls 1, "<echo '$xmax $ymax $zmin'" with points ls 1, "<echo '$xmin $ymin $zmax'" with points ls 1, "<echo '$xmax $ymin $zmax'" with points ls 1, "<echo '$xmin $ymax $zmax'" with points ls 1, "<echo '$xmax $ymax $zmax'" with points ls 1

EOF

done

