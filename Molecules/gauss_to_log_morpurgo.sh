CWD=$(pwd)

for struct in Rubrene # sexithiophene_wH # TIPS-Pc C8_BTBT Rubrene
do
	for mol in a b
	do
	  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}_mol${mol}.xyz -n ${struct}_mol${mol}_neut -y ${struct}_connectivity.dat -z 8.0 ${struct}_neut.log 
	  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}_mol${mol}.xyz -n ${struct}_mol${mol}_cation -y ${struct}_connectivity.dat -z 8.0 ${struct}_cation.log 
	  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}_mol${mol}.xyz -n ${struct}_mol${mol}_anion -y ${struct}_connectivity.dat -z 8.0 ${struct}_anion.log 
	done
done

#for struct in BENZENE
#do
#  ./gauss_log_to_py_xyz_anisobonds.sh -a ${struct}_cif_no_gauss_orient.xyz -b ${struct}.xyz -y ${struct}_neut_HF631PLUSGSTAR_connectivity.dat -z 8.0 ${struct}_neut_HF631PLUSGSTAR.log 
#done

