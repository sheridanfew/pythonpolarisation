mols_cen=['sp_sexithiophene_mola_STATE.xyz','sp_sexithiophene_molb_neut.xyz','sp_sexithiophene_molc_neut.xyz','sp_sexithiophene_mold_neut.xyz']
mols_sur=['sp_sexithiophene_mola_neut.xyz','sp_sexithiophene_molb_neut.xyz','sp_sexithiophene_molc_neut.xyz','sp_sexithiophene_mold_neut.xyz']
mols_outer=['sp_sexithiophene_mola_neut.xyz','sp_sexithiophene_molb_neut.xyz','sp_sexithiophene_molc_neut.xyz','sp_sexithiophene_mold_neut.xyz']

Natoms=44

#From cif:
'''
Sexithiophene
_symmetry_equiv_pos_as_xyz
1 x,y,z
2 1/2-x,1/2+y,1/2-z
3 -x,-y,-z
4 1/2+x,1/2-y,1/2+z
_cell_length_a                   44.708(6)
_cell_length_b                   7.851(3)
_cell_length_c                   6.029(2)
_cell_angle_alpha                90
_cell_angle_beta                 90.76(2)
_cell_angle_gamma                90
_cell_volume                     2116.01
'''
#Get translation vectors:

a=44.7086/0.5291772109217
b=7.8513/0.5291772109217
c=6.0292/0.5291772109217

alpha=90*(pi/180)
beta=90.762*(pi/180)
gamma=90*(pi/180)

cif_unit_cell_volume=2116.01/(a*b*c*(0.5291772109217**3))
