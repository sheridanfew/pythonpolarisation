for thiophene in 1T 2T 3T 4T 5T 6T 7T 8T 12T
do
	./gauss_log_to_py_xyz.sh -h ${thiophene}_cation_B3LYP_6-31gSTAR_opt*HF*.log ${thiophene}_cation_B3LYP_6-31gSTAR_opt.log ${thiophene}_cation
	./gauss_log_to_py_xyz.sh -h ${thiophene}_neutral_B3LYP_6-31gSTAR_opt*HF*.log ${thiophene}_neutral_B3LYP_6-31gSTAR_opt.log ${thiophene}_neutral
done

./multifitpol.sh -n nTneut 1T_neutral 2T_neutral 3T_neutral 4T_neutral 5T_neutral 6T_neutral 7T_neutral 8T_neutral 12T_neutral

./multifitpol.sh -n nTcat 1T_cation 2T_cation 4T_cation
