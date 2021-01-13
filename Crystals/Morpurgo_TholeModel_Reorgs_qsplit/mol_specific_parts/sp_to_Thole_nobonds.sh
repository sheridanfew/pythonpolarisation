#mols_cen=['sp_C8_BTBT_mola_STATE.xyz','sp_C8_BTBT_molb_neut.xyz']
#mols_sur=['sp_C8_BTBT_mola_neut.xyz','sp_C8_BTBT_molb_neut.xyz']
#mols_outer=['sp_C8_BTBT_mola_neut.xyz','sp_C8_BTBT_molb_neut.xyz']

#Rubrene_mola_anion_cifstruct_ESP.xyz

for file in *py; 
	do
	cat $file | sed 's/Pc/Pentacene/g' > tmp
	mv tmp $file
#	cat $file | sed 's/.xyz/_nobonds_cifstruct_ESP.xyz/g'  | sed 's/neut/neutral/g' | sed 's/sp_//g' | sed 's/sexithiophene_/sexithiophene_wH_/g' > tmp
#	mv tmp $file
#	cat $file | sed "s/\['/\['Morpurgo\//g"  | sed "s/,'/,'Morpurgo\//g" > tmp
#	mv tmp $file
done
