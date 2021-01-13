CWD=$(pwd)

for struct in BENZENE # NAPHTA TETCEN ANTCEN Pc
do
  rm ${struct}*aniso*xyz
  ./gauss_log_to_py_xyz_anisobonds_CORRECTEDBONDS.sh -b ${struct}.xyz -n ${struct}_picen -y ${struct}_neut_HF631PLUSGSTAR_connectivity.dat -z 8.0 ${struct}_neut_HF631PLUSGSTAR.log 
  ./gauss_log_to_py_xyz_anisobonds_CORRECTEDBONDS.sh -b ${struct}.xyz -n ${struct}_singcen -y ${struct}_singcen_connectivity.dat -z 8.0 ${struct}_neut_HF631PLUSGSTAR.log 
done
