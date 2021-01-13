mols_cen=['C60_STATE_aniso_cifstruct_chelpg.xyz']
mols_sur=['C60_neut_aniso_cifstruct_chelpg.xyz']
mols_outer=['sp_C60_neut.xyz']

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
