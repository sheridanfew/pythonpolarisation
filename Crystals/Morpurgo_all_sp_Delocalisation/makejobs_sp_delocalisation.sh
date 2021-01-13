# Makes jobs to calculate polaron binding energy for a charge delocalised over a few molecules in the centre of the (sub)crystal

for mol in TIPS_Pc # point PTCDA Sexithiophene C8_BTBT Pc TIPS_Pc Rubrene
do
	mkdir Jobs/${mol}
	for state in neut cation anion
	do
		for insize in  4 # 6 7
		do
			for delocdist in 0 1 2
			do
       		mkdir Jobs/${mol}/${mol}_${state}_neut_delocdist${delocdist}_inner${insize}/
				outfile="Jobs/${mol}/${mol}_${state}_neut_delocdist${delocdist}_inner${insize}/${mol}_${state}_neut_delocdist${delocdist}_inner${insize}.py"
		 		cp header.py $outfile
		 		cat << EOF >> $outfile
name='${mol}_${state}_neut_delocdist${delocdist}_inner${insize}'
state='${state}'
# Max no of TVs from centre to delocalise the charge
delocdist=${delocdist}
#For crystals here, all diamond shape and centred at centre
insize=${insize}
outsize=0
#number of TVs in each dir central mol is from edge of inner region
EOF
				cat mol_specific_parts/${mol}.py | sed "s/STATE/neut/g" >> $outfile # N.B. Change neut here to ${state} if desire central molecule to have the polarisability tensor associated with a quantum chemical calculation of a molecule in state ${state}. Here we do not as we are looking at dif. charge delocalisations.
				cat footer_diamond_coarse_delocalisation.py >> $outfile
				echo "Jobs/${mol}/${mol}_${state}_neut_delocdist${delocdist}_inner${insize}/"
			done
		done
	done
done
