mols_cen=['Morpurgo/C8_BTBT_mola_STATE_nobonds_cifstruct_ESP.xyz','Morpurgo/C8_BTBT_molb_neutral_nobonds_cifstruct_ESP.xyz']
mols_sur=['Morpurgo/C8_BTBT_mola_neutral_nobonds_cifstruct_ESP.xyz','Morpurgo/C8_BTBT_molb_neutral_nobonds_cifstruct_ESP.xyz']
mols_outer=['Morpurgo/C8_BTBT_mola_neutral_nobonds_cifstruct_ESP.xyz','Morpurgo/C8_BTBT_molb_neutral_nobonds_cifstruct_ESP.xyz']

#From cif:
'''
PTCDA
_cell_length_a                   3.740
_cell_length_b                   11.96
_cell_length_c                   17.34
_cell_angle_alpha                90
_cell_angle_beta                 98.8
_cell_angle_gamma                90
_cell_volume                     766.495

'''
#Get translation vectors:

a=3.740/0.5291772109217
b=11.96/0.5291772109217
c=17.34/0.5291772109217

alpha=90*(pi/180)
beta=98.8*(pi/180)
gamma=90*(pi/180)

cif_unit_cell_volume=766.495/(a*b*c*(0.5291772109217**3))
