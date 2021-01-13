mols_cen=['Morpurgo/sexithiophene_wH_mola_STATE_nobonds_cifstruct_ESP.xyz','Morpurgo/sexithiophene_wH_molb_neutral_nobonds_cifstruct_ESP.xyz','Morpurgo/sexithiophene_wH_molc_neutral_nobonds_cifstruct_ESP.xyz','Morpurgo/sexithiophene_wH_mold_neutral_nobonds_cifstruct_ESP.xyz']
mols_sur=['Morpurgo/sexithiophene_wH_mola_neutral_nobonds_cifstruct_ESP.xyz','Morpurgo/sexithiophene_wH_molb_neutral_nobonds_cifstruct_ESP.xyz','Morpurgo/sexithiophene_wH_molc_neutral_nobonds_cifstruct_ESP.xyz','Morpurgo/sexithiophene_wH_mold_neutral_nobonds_cifstruct_ESP.xyz']
mols_outer=['Morpurgo/sexithiophene_wH_mola_neutral_nobonds_cifstruct_ESP.xyz','Morpurgo/sexithiophene_wH_molb_neutral_nobonds_cifstruct_ESP.xyz','Morpurgo/sexithiophene_wH_molc_neutral_nobonds_cifstruct_ESP.xyz','Morpurgo/sexithiophene_wH_mold_neutral_nobonds_cifstruct_ESP.xyz']

#From cif:
'''
Sexithiophene
_symmetry_equiv_pos_as_nobonds_cifstruct_ESP.xyz
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
