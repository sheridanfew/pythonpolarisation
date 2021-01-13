mols_cen=['sp_Pc_mola_STATE.xyz','sp_Pc_molb_neut.xyz']
mols_sur=['sp_Pc_mola_neut.xyz','sp_Pc_molb_neut.xyz']
mols_outer=['sp_Pc_mola_neut.xyz','sp_Pc_molb_neut.xyz']

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

