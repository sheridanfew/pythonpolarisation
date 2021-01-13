CWD=$(pwd)

for struct in C2H4 SiH4 ethene CH4 C3H8 C2H6  C2H2
do
  ./gauss_log_to_py_xyz_anisobonds.sh -y ${struct}_connectivity.dat -z 8.0 ${struct}.log 
done

cd /home/spf310/Dropbox/Phd/PythonPolarization/Singlemol_Pol_Detemination

for struct in C2H4 SiH4 ethene CH4 C3H8 C2H6 C2H2
do
  python ${struct}*ani*abs*it.py | tee ${struct}_absfit.out &
  python ${struct}*ani*eta*it.py | tee ${struct}_etafit.out &
done

cd $CWD
