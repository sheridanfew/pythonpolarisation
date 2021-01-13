CWD=$(pwd)

list=""
listpi=""

for struct in 1,3,5,7,9-undecapentaene 1,3,5,7,9-undecapentayne 1,3,5,7-nonatetraene 1,3,5,7-nonatetrayne 1,3,5-heptatriene 1,3,5-heptatriyne 1,3-pentadiene 1,3-pentadiyne 2,4,6,8,10-dodecapentaene 2,4,6,8,10-dodecapentayne 2,4,6,8-decatetraene 2,4,6,8-decatetrayne 2,4,6-octatriene 2,4,6-octatriyne 2,4-hexadiene 2,4-hexadiyne butane butene butyne decane dodecane ethane heptane hexane methane nonane octane pentane propane undecane
do
	name=chain_$(echo ${struct} | sed 's/,/_/g' | sed 's/-/_/g' )
	if [[ $struct == *ene* ]]; then
		echo "ALKENE: $struct"
		#./gauss_log_to_py_xyz_anisobonds_CORRECTEDBONDS_nosingfit.sh -n ${name}_pi -w -y ${struct}_connectivity.dat -z 8.0 ${struct}.log
		#./gauss_log_to_py_xyz_anisobonds_CORRECTEDBONDS_nosingfit.sh -n ${name} -y ${struct}_connectivity.dat -z 8.0 ${struct}.log  
		listpi=$(echo -en "$listpi ${name}_pi_fit_aniso" )
	elif [[ $struct == *yne* ]]; then
		echo "ALKYNE: $struct"
		#./gauss_log_to_py_xyz_anisobonds_CORRECTEDBONDS_nosingfit.sh -n ${name}_doublepi -x -y ${struct}_connectivity.dat -z 8.0 ${struct}.log
		#./gauss_log_to_py_xyz_anisobonds_CORRECTEDBONDS_nosingfit.sh -n ${name} -y ${struct}_connectivity.dat -z 8.0 ${struct}.log
		listpi=$(echo -en "$listpi ${name}_doublepi_fit_aniso" )
	else
		#./gauss_log_to_py_xyz_anisobonds_CORRECTEDBONDS_nosingfit.sh -n ${name} -y ${struct}_connectivity.dat -z 8.0 ${struct}.log 
		listpi=$(echo -en "$listpi ${name}_fit_aniso" )
	fi
	list=$(echo -en "$list ${name}_fit_aniso" )

done

./multifitpol.sh -n alkanes_enes_ynes $list
./multifitpol.sh -n alkanes_enes_ynes_wpi $listpi

./multifitpol.sh -n alkanes_enes_ynes_acenes_wpi $listpi BENZENE_picen_fit_aniso NAPHTA_picen_fit_aniso TETCEN_picen_fit_aniso ANTCEN_picen_fit_aniso Pc_picen_fit_aniso
./multifitpol.sh -n alkanes_enes_ynes_acenes $list BENZENE_singcen_fit_aniso NAPHTA_singcen_fit_aniso TETCEN_singcen_fit_aniso ANTCEN_singcen_fit_aniso Pc_singcen_fit_aniso

./multifitpol.sh -n alkanes_enes_ynes_acenes_nTs_wpi $listpi BENZENE_picen_fit_aniso NAPHTA_picen_fit_aniso TETCEN_picen_fit_aniso ANTCEN_picen_fit_aniso Pc_picen_fit_aniso thio_1T_neut_fit_aniso thio_2T_neut_fit_aniso thio_3T_neut_fit_aniso thio_4T_neut_fit_aniso thio_5T_neut_fit_aniso thio_6T_neut_fit_aniso thio_7T_neut_fit_aniso thio_8T_neut_fit_aniso
