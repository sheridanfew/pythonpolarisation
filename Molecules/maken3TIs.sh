	./gauss_log_to_py_xyz.sh -h 13TI_cation_B3LYP_6-31gSTAR_opt_TD5_2013_04_08_HF_6-31gSTAR.log 13TI_cation_B3LYP_6-31gSTAR_opt_TD5_2013_04_08.log 13TI_cation
	./gauss_log_to_py_xyz.sh -h 33TI_cation_B3LYP_6-31gSTAR_opt_HF_6-31gSTAR.log 33TI_cation_B3LYP_6-31gSTAR_opt.log 33TI_cation
	./gauss_log_to_py_xyz.sh -h 13TI_neutral_B3LYP_6-31gSTAR_opt_TD5_2013_04_08_HF_6-31gSTAR.log 13TI_neutral_B3LYP_6-31gSTAR_opt_TD5_2013_04_08.log 13TI_neutral
	./gauss_log_to_py_xyz.sh -h 33TI_neutral_B3LYP_6-31gSTAR_opt_TD3_HF_6-31gSTAR.log 33TI_neutral_B3LYP_6-31gSTAR_opt_TD3.log 33TI_neutral


./multifitpol.sh -n n3TIneut 13TI_neutral 33TI_neutral

./multifitpol.sh -n n3TIcat 13TI_cation 33TI_cation
