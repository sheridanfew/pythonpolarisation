for struct in TIPS-Pc
do
# for mol in a b
# do
  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}.xyz -n TIPS_Pc_neut -y ${struct}_connectivity.dat -z 8.0 ${struct}_neut.log
  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}.xyz -n TIPS_Pc_cation -y ${struct}_connectivity.dat -z 8.0 ${struct}_cation.log   
  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}.xyz -n TIPS_Pc_anion -y ${struct}_connectivity.dat -z 8.0 ${struct}_anion.log
# done
done

#for struct in C8_BTBT
#do
# for mol in a b
# do
#  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}_mol${mol}.xyz -n ${struct}_mol${mol}_neut -y ${struct}_connectivity.dat -z 8.0 ${struct}_neut.log
#  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}_mol${mol}.xyz -n ${struct}_mol${mol}_cation -y ${struct}_connectivity.dat -z 8.0 ${struct}_cation.log   
#  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}_mol${mol}.xyz -n ${struct}_mol${mol}_anion -y ${struct}_connectivity.dat -z 8.0 ${struct}_anion.log
# done
#done

#for struct in Rubrene
#do
# for mol in a b
# do
#  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}_mol${mol}.xyz -n ${struct}_mol${mol}_neut -y ${struct}_connectivity.dat -z 8.0 ${struct}_neutral_HF_6-31PLUSgSTAR.log
#  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}_mol${mol}.xyz -n ${struct}_mol${mol}_cation -y ${struct}_connectivity.dat -z 8.0 ${struct}_cation_HF_6-31PLUSgSTAR.log   
#  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}_mol${mol}.xyz -n ${struct}_mol${mol}_anion -y ${struct}_connectivity.dat -z 8.0 ${struct}_anion_HF_6-31PLUSgSTAR.log
# done
#done
