for mol in Sexithiophene #C8_BTBT Pc TIPS_Pc Rubrene
do
	mkdir Jobs/${mol}
	for state in cation anion # neut
	do
		for insize in 2 #0 1 
		do
			for outsize in 0 # 1 2 3 4 5
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
EOF
				cat mol_specific_parts/${mol}.py | sed "s/STATE/${state}/g" >> $outfile
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
