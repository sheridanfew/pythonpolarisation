for n in 1 2 3 4 5 6 7 8
do
	file="thio_${n}T_neut_aniso_chelpg.xyz"
	newfilename_lin_mean=$(echo $file | sed 's/.xyz/_thole_Lin_mean.xyz/')
	newfilename_exp_mean=$(echo $file | sed 's/.xyz/_thole_Exp_mean.xyz/')
	newfilename_lin_components=$(echo $file | sed 's/.xyz/_thole_Lin_components.xyz/')
	newfilename_exp_components=$(echo $file | sed 's/.xyz/_thole_Exp_components.xyz/')
	newfilename_lin_empirical=$(echo $file | sed 's/.xyz/_thole_Lin_empirical.xyz/')
	newfilename_exp_empirical=$(echo $file | sed 's/.xyz/_thole_Exp_empirical.xyz/')
	Natoms=$(echo "(7 * $n) + 2" | bc -l)
	head -n $Natoms $file > atoms.tmp

	#C	
	cat atoms.tmp | sed 's/iso=9.0799/iso=7.7764/g' > tmp; mv tmp $newfilename_lin_components
	cat atoms.tmp | sed 's/iso=9.0799/iso=6.578/g' > tmp; mv tmp $newfilename_exp_components
	cat atoms.tmp | sed 's/iso=9.0799/iso=8.7622/g' > tmp; mv tmp $newfilename_lin_mean
	cat atoms.tmp | sed 's/iso=9.0799/iso=7.9421/g' > tmp; mv tmp $newfilename_exp_mean
	cat atoms.tmp | sed 's/iso=9.0799/iso=10.1756/g' > tmp; mv tmp $newfilename_lin_empirical
	cat atoms.tmp | sed 's/iso=9.0799/iso=8.6959/g' > tmp; mv tmp $newfilename_exp_empirical
	#H
	cat atoms.tmp | sed 's/iso=0.001/iso=0.5974/g' > tmp; mv tmp $newfilename_lin_components
	cat atoms.tmp | sed 's/iso=0.001/iso=1.3849/g' > tmp; mv tmp $newfilename_exp_components
	cat atoms.tmp | sed 's/iso=0.001/iso=1.5350/g' > tmp; mv tmp $newfilename_lin_mean
	cat atoms.tmp | sed 's/iso=0.001/iso=1.3368/g' > tmp; mv tmp $newfilename_exp_mean
	cat atoms.tmp | sed 's/iso=0.001/iso=3.5020/g' > tmp; mv tmp $newfilename_lin_empirical
	cat atoms.tmp | sed 's/iso=0.001/iso=2.7927/g' > tmp; mv tmp $newfilename_exp_empirical
	#S
	cat atoms.tmp | sed 's/iso=15.4061/iso=11.0208/g' > tmp; mv tmp $newfilename_lin_components
	cat atoms.tmp | sed 's/iso=15.4061/iso=10.8651/g' > tmp; mv tmp $newfilename_exp_components
	cat atoms.tmp | sed 's/iso=15.4061/iso=15.0931/g' > tmp; mv tmp $newfilename_lin_mean
	cat atoms.tmp | sed 's/iso=15.4061/iso=13.9500/g' > tmp; mv tmp $newfilename_exp_mean
	cat atoms.tmp | sed 's/iso=15.4061/iso=19.7422/g' > tmp; mv tmp $newfilename_lin_empirical
	cat atoms.tmp | sed 's/iso=15.4061/iso=16.6984/g' > tmp; mv tmp $newfilename_exp_empirical

	rm atoms.tmp
done

: <<'END'
Thole ab initio components fit	lin	expon
H	0.5974	1.3849
C	7.7764	6.578
N	4.7575	2.9855
O	2.6442	2.3184
F	0.8149	1.4906
S	11.0208	10.8651
Cl	7.5682	6.7129
Br	12.7339	12.3384
I	20.2867	20.0081
		
a	1.6623	2.8990
END
