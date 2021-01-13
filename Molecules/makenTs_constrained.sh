for thiophene in 2T 3T 4T 5T 6T 7T 8T
do
	./gauss_log_to_py_xyz.sh -h ${thiophene}_CONSTRAINED_cation_B3LYP_6-31gSTAR_opt_2013_04_15*HF*.log ${thiophene}_CONSTRAINED_cation_B3LYP_6-31gSTAR_opt_2013_04_15.log ${thiophene}_CONSTRAINED_cation
	./gauss_log_to_py_xyz.sh -h ${thiophene}_CONSTRAINED_neutral_B3LYP_6-31gSTAR_opt_2013_04_15*HF*.log ${thiophene}_CONSTRAINED_neutral_B3LYP_6-31gSTAR_opt_2013_04_15.log ${thiophene}_CONSTRAINED_neutral
done

./multifitpol.sh -n nTneut_CONSTRAINED 1T_neutral 2T_CONSTRAINED_neutral 3T_CONSTRAINED_neutral 4T_CONSTRAINED_neutral 5T_CONSTRAINED_neutral 6T_CONSTRAINED_neutral 7T_CONSTRAINED_neutral 8T_CONSTRAINED_neutral

./multifitpol.sh -n nTcat_CONSTRAINED 1T_cation 2T_CONSTRAINED_cation 3T_CONSTRAINED_cation 4T_CONSTRAINED_cation 6T_CONSTRAINED_cation 7T_CONSTRAINED_cation
