#!/bin/sh
#PBS -l walltime=71:58:02
#PBS -l select=1:ncpus=8:mem=11800mb

echo "Execution started:"
date

cd /work/spf310/PhD/Polarisation/PythonPolarization/Crystals/Morpurgo_all_sp_Reorgs_qsplit/Jobs/TIPS_Pc/TIPS_Pc_cation_neut_inner5_outer0//

pbsexec python TIPS_Pc_cation_neut_inner5_outer0.py > TIPS_Pc_cation_neut_inner5_outer0.out

echo -en "\n\ncsvs:\n"
cat *csv

echo -en "\n\ndat files:\n"
cat *dat

echo "Job finished:"
date

