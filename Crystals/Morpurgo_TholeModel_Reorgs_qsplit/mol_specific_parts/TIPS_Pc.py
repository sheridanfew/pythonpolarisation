mols_cen=['Morpurgo/TIPS_Pentacene_STATE_nobonds_cifstruct_ESP.xyz']
mols_sur=['Morpurgo/TIPS_Pentacene_neutral_nobonds_cifstruct_ESP.xyz']
mols_outer=['Morpurgo/TIPS_Pentacene_neutral_nobonds_cifstruct_ESP.xyz']

#From cif:
'''
TIPS

data_k01029 

_cell_length_a                    7.5650(15) 
_cell_length_b                    7.7500(15) 
_cell_length_c                    16.835(3) 
_cell_angle_alpha                 89.15(3) 
_cell_angle_beta                  78.42(3) 
_cell_angle_gamma                 83.63(3) 

_cell_volume                      960.9(3) 

'''
#Get translation vectors:

a=7.565015/0.5291772109217
b=7.750015/0.5291772109217
c=16.8353/0.5291772109217

alpha=89.153*(pi/180)
beta=78.423*(pi/180)
gamma=83.633*(pi/180)

cif_unit_cell_volume=960.9/(a*b*c*(0.5291772109217**3))
