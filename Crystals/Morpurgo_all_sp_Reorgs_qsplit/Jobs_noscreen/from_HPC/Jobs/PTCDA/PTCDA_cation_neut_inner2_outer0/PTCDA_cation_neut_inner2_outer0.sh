#!/bin/sh
#PBS -l walltime=71:58:02
#PBS -l select=1:ncpus=8:mem=11800mb

echo "Execution started:"
date

cd /work/spf310/PhD/Polarisation/PythonPolarization/Crystals/Morpurgo_all_sp_Reorgs_qsplit/Jobs/PTCDA/PTCDA_cation_neut_inner2_outer0//

pbsexec python PTCDA_cation_neut_inner2_outer0.py > PTCDA_cation_neut_inner2_outer0.out

echo -en "\n\ncsvs:\n"
cat *csv

echo -en "\n\ndat files:\n"
cat *dat

echo "Job finished:"
date

