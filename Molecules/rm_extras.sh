  for file in C8*xyz Rubre*xyz TIPS*xyz Pc*xyz sexi*xyz; do if [[ $file == *mola_* && $file != *cif* ]]; then mv $file $(echo $file | sed 's/mola_//'); fi; done
  for file in C8*xyz Rubre*xyz TIPS*xyz Pc*xyz sexi*xyz; do if [[ $file == *mol?_* && $file != *cif* ]]; then rm $file; fi; done

