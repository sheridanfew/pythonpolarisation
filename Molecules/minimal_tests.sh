for file in C[0-9]*log ethene*log SiH4.log ; do echo $file;
        confile=$(echo $file | sed 's/.log//g' )_connectivity.dat
	./gauss_log_to_py_xyz_newetafit.sh -y $confile -z 8 $file
done
