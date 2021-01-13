mols_cen=['Pentacene_mola_STATE_aniso_cifstruct_mul.xyz','Pentacene_molb_neut_aniso_cifstruct_mul.xyz']
mols_sur=['Pentacene_mola_neut_aniso_cifstruct_mul.xyz','Pentacene_molb_neut_aniso_cifstruct_mul.xyz']
mols_outer=['sp_Pc_mola_neut.xyz','sp_Pc_molb_neut.xyz']

screenradius=2.6738915344

#From cif:
'''
Pc
_cell_length_a                   7.900
_cell_length_b                   6.060
_cell_length_c                   16.010
_cell_angle_alpha                101.90
_cell_angle_beta                 112.60
_cell_angle_gamma                85.80
_cell_volume                     692.384

'''
#Get translation vectors:

a=7.900/0.5291772109217
b=6.060/0.5291772109217
c=16.010/0.5291772109217

alpha=101.90*(pi/180)
beta=112.60*(pi/180)
gamma=90*(pi/180)

cif_unit_cell_volume=692.384/(a*b*c*(0.5291772109217**3))

