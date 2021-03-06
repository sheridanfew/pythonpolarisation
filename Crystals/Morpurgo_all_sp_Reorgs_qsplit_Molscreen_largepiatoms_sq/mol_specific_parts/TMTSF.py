mols_cen=['sp_TMTSF_mola_STATE.xyz']
mols_sur=['sp_TMTSF_mola_neut.xyz']
mols_outer=['sp_TMTSF_mola_neut.xyz']

Natoms=10

#From cif:
'''
TMTSF
_cell_length_a                   6.935(1)
_cell_length_b                   8.092(2)
_cell_length_c                   6.314(1)
_cell_angle_alpha                105.51(2)
_cell_angle_beta                 95.39(1)
_cell_angle_gamma                108.90(2)
_cell_volume                     316.619

'''
#Get translation vectors:

a=6.9351/0.5291772109217
b=8.0922/0.5291772109217
c=6.3141/0.5291772109217

alpha=105.512*(pi/180)
beta=95.391*(pi/180)
gamma=108.902*(pi/180)

cif_unit_cell_volume=316.619/(a*b*c*(0.5291772109217**3))
