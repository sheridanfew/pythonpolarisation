mols_cen=['C8_BTBT_mola_STATE_aniso_cifstruct_mul.xyz','C8_BTBT_molb_neut_aniso_cifstruct_mul.xyz']
mols_sur=['C8_BTBT_mola_neut_aniso_cifstruct_mul.xyz','C8_BTBT_molb_neut_aniso_cifstruct_mul.xyz']
mols_outer=['sp_C8_BTBT_mola_neut.xyz','sp_C8_BTBT_molb_neut.xyz']

screenradius=2.2302979232

#From cif:
'''
C8_BTBT
_cell_length_a                   5.927(7)
_cell_length_b                   7.88(1)
_cell_length_c                   29.18(4)
_cell_angle_alpha                90
_cell_angle_beta                 92.443(4)
_cell_angle_gamma                90
_cell_volume                     1361.61

'''
#Get translation vectors:

a=5.9277/0.5291772109217
b=7.881/0.5291772109217
c=29.184/0.5291772109217

alpha=90*(pi/180)
beta=92.4434*(pi/180)
gamma=90*(pi/180)

cif_unit_cell_volume=1361.61/(a*b*c*(0.5291772109217**3))
