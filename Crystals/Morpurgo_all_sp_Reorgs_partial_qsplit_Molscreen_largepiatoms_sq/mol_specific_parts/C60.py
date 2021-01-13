mols_cen=['sp_C60_mola_STATE.xyz','sp_C60_molb_neut.xyz','sp_C60_molc_neut.xyz','sp_C60_mold_neut.xyz']
mols_sur=['sp_C60_mola_neut.xyz','sp_C60_molb_neut.xyz','sp_C60_molc_neut.xyz','sp_C60_mold_neut.xyz']
mols_outer=['sp_C60_mola_neut.xyz','sp_C60_molb_neut.xyz','sp_C60_molc_neut.xyz','sp_C60_mold_neut.xyz']

Natoms=60

#From cif:
'''
C60
_cell_length_a                   14.052(5)
_cell_length_b                   14.052(5)
_cell_length_c                   14.052(5)
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                90
_cell_volume                     2774.69
_cell_formula_units_Z            4
'''
#Get translation vectors:

a=14.0525/0.5291772109217
b=14.0525/0.5291772109217
c=14.0525/0.5291772109217

alpha=90*(pi/180)
beta=90*(pi/180)
gamma=90*(pi/180)

cif_unit_cell_volume=2774.69/(a*b*c*(0.5291772109217**3))
