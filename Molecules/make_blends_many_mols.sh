for file in P3TI/cation/opt/33TI_cation_B3LYP_6-31gSTAR_opt/33TI_cation_B3LYP_6-31gSTAR_opt.log P4TI/cation/opt/44TI_cation_B3LYP_6-31gSTAR_opt/44TI_cation_B3LYP_6-31gSTAR_opt.log P5TI/cation/opt/35TI_cation_B3LYP_6-31gSTAR_opt/35TI_cation_B3LYP_6-31gSTAR_opt.log  PBTT-DPPb/cation/opt/3BTT-DPPb_cation_B3LYP_6-31gSTAR_opt_TD6_2013_04_12/3BTT-DPPb_cation_B3LYP_6-31gSTAR_opt_TD6_2013_04_12.log PDPP-TT-T-TT/cation/opt/3DPP-TT-T-TT_cation_B3LYP_6-31gSTAR_opt_TD4_2013_04_23/3DPP-TT-T-TT_cation_B3LYP_6-31gSTAR_opt_TD4_2013_04_23.log PDPP-TT-Tz-TT/cation/opt/3DPP-TT-Tz-TT_cation_B3LYP_6-31gSTAR_opt_TD4_2013_04_29/3DPP-TT-Tz-TT_cation_B3LYP_6-31gSTAR_opt_TD4_2013_04_29.log PFTBT/cation/opt/4FTBT_cation_B3LYP_6-31gSTAR_opt/4FTBT_cation_B3LYP_6-31gSTAR_opt.log PPV_MDMO/cation/opt/5PV_MDMO_cation_B3LYP_6-31gSTAR_opt_TD3/5PV_MDMO_cation_B3LYP_6-31gSTAR_opt_TD3.log PSi_CPDTBT/cation/opt/4Si_CPDTBT_cation_B3LYP_6-31gSTAR_opt/4Si_CPDTBT_cation_B3LYP_6-31gSTAR_opt.log PT/cation/opt/12T_cation_B3LYP_6-31gSTAR_opt/12T_cation_B3LYP_6-31gSTAR_opt.log PTI1/cation/opt/4TI1_cation_B3LYP_6-31gSTAR_opt_TD10/4TI1_cation_B3LYP_6-31gSTAR_opt_TD10.log PCDTBT/cation/opt/3CDTBT_cation_B3LYP_6-31gSTAR_opt/3CDTBT_cation_B3LYP_6-31gSTAR_opt.log PCPDT-BT/cation/opt/4CPDTBT_cation_B3LYP_6-31gSTAR_opt/4CPDTBT_cation_B3LYP_6-31gSTAR_opt.log
do
#	scp hpc:/work/spf310/Polymers/$file .
	log=$(basename $file)
	echo "log = $log"
#	rm $(echo "$log" | sed 's/.log//g' | sed 's/-/_/g' )*xyz
	echo "rm $(echo "$log" | sed 's/.log//g' | sed 's/-/_/g' )*xyz"
#	grep 'elpg' $log
	./gauss_log_to_py_xyz.sh $log
	xyzmul=$(echo "$log" | sed 's/.log/_mul.xyz/g' | sed 's/-/_/g' )
	echo "xyzmul=$xyzmul"
	mol=$(echo $log | sed 's/_/ /g' | awk '{printf "%s",$1}')
	echo "mol = $mol"
	echo "makeblends_pol_model.sh -c -m $mol -s ${mol}_PCBM_coulomb $xyzmul PCBM_EDITEDO_anion_mul.xyz"
	makeblends_pol_model_coulonly.sh -c -h 3.5 -i 3.8 -j 0.5 -m $mol -s ${mol}_PCBM_coulomb -u 15.0 $xyzmul PCBM_EDITEDO_anion_mul.xyz
	makeblends_pol_model_coulonly.sh -c -h 3.5 -i 3.8 -j 0.5 -m $mol -s ${mol}_PC71BM_coulomb -u 15.0 $xyzmul PC71BM_anion_B3LYP_6_31gSTAR_opt_mul.xyz
done


