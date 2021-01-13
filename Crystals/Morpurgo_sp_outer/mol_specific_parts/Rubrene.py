mols_cen=['Rubrene_mola_STATE_aniso_cifstruct_chelpg.xyz','Rubrene_molb_neut_aniso_cifstruct_chelpg.xyz']
mols_sur=['Rubrene_mola_neut_aniso_cifstruct_chelpg_edited.xyz','Rubrene_molb_neut_aniso_cifstruct_chelpg.xyz']
mols_outer=['sp_Rubrene_mola_neut.xyz','sp_Rubrene_molb_neut.xyz']

#centres=['Rubrene_mola_anion_aniso_cifstruct_chelpg_edited.xyz','Rubrene_molb_neut_aniso_cifstruct_chelpg_edited.xyz','Rubrene_mola_neut_aniso_cifstruct_chelpg_edited.xyz','Rubrene_molb_neut_aniso_cifstruct_chelpg_edited.xyz']

#surroundings=['Rubrene_mola_neut_aniso_cifstruct_chelpg_edited.xyz','Rubrene_molb_neut_aniso_cifstruct_chelpg_edited.xyz','Rubrene_mola_neut_aniso_cifstruct_chelpg_edited.xyz','Rubrene_molb_neut_aniso_cifstruct_chelpg_edited.xyz']

#From cif:
'''
Rubrene
_cell_length_a                   7.184(1)
_cell_length_b                   14.433(3)
_cell_length_c                   26.897(7)
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                90
_cell_volume                     2788.86
_cell_formula_units_Z            4
'''
#Get translation vectors:

a=7.1841/0.5291772109217
b= 14.4333/0.5291772109217
c= 26.8977/0.5291772109217

alpha=90*(pi/180)
beta=90*(pi/180)
gamma=90*(pi/180)

cif_unit_cell_volume=2788.86/(a*b*c*(0.5291772109217**3))
