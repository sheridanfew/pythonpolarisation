mols_cen=['PDIF_CN2_mola_STATE_aniso_cifstruct_mul.xyz']
mols_sur=['PDIF_CN2_mola_neut_aniso_cifstruct_mul.xyz']
mols_outer=['sp_PDIFCN2_neut.xyz']

screenradius=2.5533199878

#From cif:
'''
PDIF-CN2
_cell_length_a                   5.2320(14)
_cell_length_b                   7.638(2)
_cell_length_c                   18.819(5)
_cell_angle_alpha                92.512(5)
_cell_angle_beta                 95.247(5)
_cell_angle_gamma                104.730(4)
_cell_volume                     722.5(3)

'''
#Get translation vectors:

a=5.232014/0.5291772109217
b=7.6382/0.5291772109217
c=18.8195/0.5291772109217

alpha=92.5125*(pi/180)
beta=95.2475*(pi/180)
gamma=104.7304*(pi/180)

cif_unit_cell_volume=722.53/(a*b*c*(0.5291772109217**3))
