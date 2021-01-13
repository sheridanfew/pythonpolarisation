for mol in PTCDA Sexithiophene C8_BTBT Pc TIPS_Pc Rubrene
do
	mkdir Jobs/${mol}
	for state in cation neut anion
	do
		for insize in 0 1 2 3 4 5 6
		do
			for outsize in 0
			do
       		mkdir Jobs/${mol}/${mol}_${state}_neut_inner${insize}_outer${outsize}/
				outfile="Jobs/${mol}/${mol}_${state}_neut_inner${insize}_outer${outsize}/${mol}_${state}_neut_inner${insize}_outer${outsize}.py"
		 		cp header.py $outfile
		 		cat << EOF >> $outfile
name='${mol}_${state}_neut_inner${insize}_outer${outsize}'
#For crystals here, all cubic and centred at centre
insize=${insize}
#number of TVs in each dir central mol is from edge of inner region
outsize=${outsize}
state='${state}'
EOF
				cat mol_specific_parts/${mol}.py | sed "s/STATE/neut/g" >> $outfile
				if [[ $state == "neut" ]]
				then
					cat footer_diamond_coarseouter_noreorgs.py >> $outfile
				else
      				cat footer_diamond_coarseouter.py >> $outfile
				fi
			done
		done
	done
done
