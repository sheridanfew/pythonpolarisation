CWD=$(pwd)

for struct in 1 2 3 4 5 6 7 8
do
  ./gauss_log_to_py_xyz_anisobonds_CORRECTEDBONDS.sh -n thio_${struct}T_neut -y thio_${struct}T_connectivity.dat -z 8.0 thio_${struct}T_neut.log 
done


#	  ./gauss_log_to_py_xyz_anisobonds.sh -b ${struct}_mol${mol}.xyz -n ${struct}_mol${mol}_neut -y ${struct}_connectivity.dat -z 8.0 ${struct}_neut.log 

#for struct in BENZENE
#do
#  ./gauss_log_to_py_xyz_anisobonds.sh -a ${struct}_cif_no_gauss_orient.xyz -b ${struct}.xyz -y ${struct}_neut_HF631PLUSGSTAR_connectivity.dat -z 8.0 ${struct}_neut_HF631PLUSGSTAR.log 
#done

