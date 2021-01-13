for file in Alkanes/*ane.gjf ConjugatedMolecules/*ene.gjf ConjugatedMolecules/*yne.gjf
do
	if $( echo "$file" | grep 'yclo') == ''
	then
		extract_connectivity_newmac.sh $file
		cp *onnectivity* (Molecules)
	fi
done


