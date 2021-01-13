mkdir C8_BTBT
mkdir dip
mkdir Pc
mkdir TIPS_Pc
mkdir Sexithiophene
mkdir Rubrene

for state in cation neut anion
do
  for innersize in 0 1 2
  do
	for outeradd in 1 2 3 4
     for shape in diamond cube
     do
       mkdir C8_BTBT/C8_BTBT_${state}_neut_${shape}_inner${size}_outer${outersize}/
		 cp header.py C8_BTBT/C8_BTBT_${state}_neut_${shape}_size${size}/C8_BTBT_${state}_neut_${shape}_size${size}.py
		 cat << EOF >> C8_BTBT/C8_BTBT_${state}_neut_${shape}_size${size}/C8_BTBT_${state}_neut_${shape}_size${size}.py
###################################
#START OF MOLECULE SPECIFIC SECTION
###################################

name='C8_BTBT_${state}_neut_${shape}'

inner_size=$inner_size
outeradd=$outer_add

centres=['C8_BTBT_mola_${state}_aniso_cifstruct_chelpg.xyz','C8_BTBT_molb_neut_aniso_cifstruct_chelpg.xyz']
surroundings=['C8_BTBT_mola_neut_aniso_cifstruct_chelpg.xyz','C8_BTBT_molb_neut_aniso_cifstruct_chelpg.xyz']
outer=['sp_C8_BTBT_mola_neut.xyz','sp_C8_BTBT_molb_neut.xyz']


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


###################################
#END OF MOLECULE SPECIFIC SECTION
###################################

EOF
      cat footer_${shape}.py >> C8_BTBT/C8_BTBT_${state}_neut_${shape}_size${size}/C8_BTBT_${state}_neut_${shape}_size${size}.py







###################################
### Pc
###################################
      mkdir Pc/Pc_${state}_neut_${shape}_size${size}/
      cp header.py Pc/Pc_${state}_neut_${shape}_size${size}/Pc_${state}_neut_${shape}_size${size}.py
      cat << EOF >> Pc/Pc_${state}_neut_${shape}_size${size}/Pc_${state}_neut_${shape}_size${size}.py
###################################
#START OF MOLECULE SPECIFIC SECTION
###################################

name='Pc_${state}_neut_diamond'

inner_size=$inner_size
outeradd=$outer_add

centres=['Pc_mola_${state}_aniso_cifstruct_chelpg.xyz','Pc_molb_neut_aniso_cifstruct_chelpg.xyz']

surroundings=['Pc_mola_neut_aniso_cifstruct_chelpg.xyz','Pc_molb_neut_aniso_cifstruct_chelpg.xyz']

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

###################################
#END OF MOLECULE SPECIFIC SECTION
###################################

EOF
      cat footer_${shape}.py >> Pc/Pc_${state}_neut_${shape}_size${size}/Pc_${state}_neut_${shape}_size${size}.py







###################################
### Rubrene
###################################

      mkdir Rubrene/Rubrene_${state}_neut_${shape}_size${size}/
      cp header.py Rubrene/Rubrene_${state}_neut_${shape}_size${size}/Rubrene_${state}_neut_${shape}_size${size}.py
      cat << EOF >> Rubrene/Rubrene_${state}_neut_${shape}_size${size}/Rubrene_${state}_neut_${shape}_size${size}.py
###################################
#START OF MOLECULE SPECIFIC SECTION
###################################

name='Rubrene_${state}_neut'

inner_size=$inner_size
outeradd=$outer_add

centres=['Rubrene_mola_${state}_aniso_cifstruct_chelpg_edited.xyz','Rubrene_molb_neut_aniso_cifstruct_chelpg_edited.xyz','Rubrene_mola_neut_aniso_cifstruct_chelpg_edited.xyz','Rubrene_molb_neut_aniso_cifstruct_chelpg_edited.xyz']

surroundings=['Rubrene_mola_neut_aniso_cifstruct_chelpg_edited.xyz','Rubrene_molb_neut_aniso_cifstruct_chelpg_edited.xyz','Rubrene_mola_neut_aniso_cifstruct_chelpg_edited.xyz','Rubrene_molb_neut_aniso_cifstruct_chelpg_edited.xyz']

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

###################################
#END OF MOLECULE SPECIFIC SECTION
###################################

EOF
      cat footer_${shape}.py >> Rubrene/Rubrene_${state}_neut_${shape}_size${size}/Rubrene_${state}_neut_${shape}_size${size}.py








###################################
### Sexithiophene
###################################

      mkdir Sexithiophene/Sexithiophene_${state}_neut_${shape}_size${size}/
      cp header.py Sexithiophene/Sexithiophene_${state}_neut_${shape}_size${size}/Sexithiophene_${state}_neut_${shape}_size${size}.py
      cat << EOF >> Sexithiophene/Sexithiophene_${state}_neut_${shape}_size${size}/Sexithiophene_${state}_neut_${shape}_size${size}.py
###################################
#START OF MOLECULE SPECIFIC SECTION
###################################
name='Sexithiophene_${state}_neut'

inner_size=$inner_size
outeradd=$outer_add

centres=['sexithiophene_mola_${state}_aniso_cifstruct_chelpg.xyz','sexithiophene_molb_neut_aniso_cifstruct_chelpg.xyz','sexithiophene_molc_neut_aniso_cifstruct_chelpg.xyz','sexithiophene_mold_neut_aniso_cifstruct_chelpg.xyz']

surroundings=['sexithiophene_mola_neut_aniso_cifstruct_chelpg.xyz','sexithiophene_molb_neut_aniso_cifstruct_chelpg.xyz','sexithiophene_molc_neut_aniso_cifstruct_chelpg.xyz','sexithiophene_mold_neut_aniso_cifstruct_chelpg.xyz']

#From cif:
'''
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

###################################
#END OF MOLECULE SPECIFIC SECTION
###################################

EOF
      cat footer_${shape}.py >> Sexithiophene/Sexithiophene_${state}_neut_${shape}_size${size}/Sexithiophene_${state}_neut_${shape}_size${size}.py









###################################
### TIPS
###################################

      mkdir TIPS_Pc/TIPS_Pc_${state}_neut_${shape}_size${size}/
      cp header.py TIPS_Pc/TIPS_Pc_${state}_neut_${shape}_size${size}/TIPS_Pc_${state}_neut_${shape}_size${size}.py
      cat << EOF >> TIPS_Pc/TIPS_Pc_${state}_neut_${shape}_size${size}/TIPS_Pc_${state}_neut_${shape}_size${size}.py
###################################
#START OF MOLECULE SPECIFIC SECTION
###################################

name='TIPS_Pc_${state}_neut'

inner_size=$inner_size
outeradd=$outer_add

centres=['TIPS_Pc_${state}_aniso_cifstruct_chelpg.xyz']
surroundings=['TIPS_Pc_neut_aniso_cifstruct_chelpg.xyz']

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

###################################
#END OF MOLECULE SPECIFIC SECTION
###################################

EOF
      cat footer_${shape}.py >> TIPS_Pc/TIPS_Pc_${state}_neut_${shape}_size${size}/TIPS_Pc_${state}_neut_${shape}_size${size}.py








###################################
### POINT DIPOLE TEST
###################################
      mkdir dip/dip_${state}_neut_${shape}_size${size}
      cp header.py dip/dip_${state}_neut_${shape}_size${size}/dip_${state}_neut_${shape}_size${size}.py
      cat << EOF >> dip/dip_${state}_neut_${shape}_size${size}/dip_${state}_neut_${shape}_size${size}.py
###################################
#START OF MOLECULE SPECIFIC SECTION
###################################

name='dip_${state}_neut_${shape}'

inner_size=$inner_size
outeradd=$outer_add

centres=['dip_mola_${state}.xyz','dip_molb_neut.xyz']
surroundings=['dip_mola_neut.xyz','dip_molb_neut.xyz']

#From cif:
'''
dip

'''
#Get translation vectors:

a=10/0.5291772109217
b=10/0.5291772109217
c=10/0.5291772109217

alpha=90*(pi/180)
beta=90*(pi/180)
gamma=90*(pi/180)

#cif_unit_cell_volume=1361.61/(a*b*c*(0.5291772109217**3))

###################################
#END OF MOLECULE SPECIFIC SECTION
###################################

EOF
      cat footer_${shape}.py >> dip/dip_${state}_neut_${shape}_size${size}/dip_${state}_neut_${shape}_size${size}.py
    done
  done
done


