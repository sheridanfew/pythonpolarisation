
cell_volume=sqrt(1 - (cos(alpha)**2) - (cos(beta)**2) - (cos(gamma)**2) + (2*cos(alpha)*cos(beta)*cos(gamma)))

print 'calc cell_volume, cif cell_volume:',cell_volume, ', ', cif_unit_cell_volume

if cell_volume > (cif_unit_cell_volume*1.1) or cell_volume < (cif_unit_cell_volume*0.9):
	print 'angles wrong?'
	exit (2)

#Converts frac coords to carts
matrix_to_cartesian=np.matrix( [[a, b*cos(gamma), c*cos(beta)],
[0, b*sin(gamma), c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)],
[0, 0, c*cell_volume/sin(gamma)]])

#carts to frac
matrix_to_fractional=matrix_to_cartesian.I

#TVs, TV[0,1,2] are the three translation vectors.
TV=matrix_to_cartesian.T

print 'TV:'
print TV

print 'TV[0,1,2] are the three translation vectors.'

f = open('%s_TVs.dat' % name, 'w')
for i in range(0,3,1):
	for j in range (0,3,1):
		f.write(str(TV[i,j]) + ' ')
	f.write('\n')
f.flush()
f.close()

for i in np.arange(0,len(centres),1):
	ReadMoleculeType('../../../Molecules/./%s' % centres[i])
	ReadMoleculeType('../../../Molecules/./%s' % surroundings[i])

#Reflection Properties (see http://mathworld.wolfram.com/Reflection.html):
''' No reflection in Pc
a_refl=-1
b_refl=1
c_refl=-1
d_refl=1.5

#in fractional coords:
#perp vector to refl plane in fractional coords
n_unnorm_frac=np.matrix( [[a_refl, b_refl, c_refl]] )
#perp vector to refl plane in fractional coords, normalised
n_refl_frac=np.matrix( [[a_refl, b_refl, c_refl]] )/ np.sqrt(np.dot(np.matrix( [[a_refl, b_refl, c_refl]] ),np.matrix( [[a_refl, b_refl, c_refl]] ).T))
#dist from origin to refl plane in fractional coords, will be in direction of normal to plane
dist_refl_norm_frac=d_refl/np.sqrt(np.dot(np.matrix( [[a_refl, b_refl, c_refl]] ),np.matrix( [[a_refl, b_refl, c_refl]] ).T))

print 'n_refl_frac', n_refl_frac

print 'd_refl', d_refl

#in cartesians
# Perpendicular cartesian vector, 
n_refl=n_refl_frac*TV
n_refl_norm=n_refl/sqrt(np.dot(n_refl,n_refl.T))

print 'n_refl', n_refl
print 'n_refl_norm', n_refl_norm

'''



cut=8.0

# Place Molecules
mols=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(len(centres))]
#mols contains all molecules.
#mols[0] contains a list of all molecules in position a, mols[1] all mols in pos'n b, etc.
#mols[0][x,y,z] contains molecule a in position x,y,z
#mols may as such be iterated over in a number of ways to consider different molecules.

print mols[0]
for x in np.arange(-size,size+1,1):
	xlist=x+size
	for y in np.arange(-size,size+1,1):
		ylist=y+size
		for z in np.arange(-size,size+1,1):
			zlist=z+size
			if x == 0 and y == 0 and z == 0:
				#central mol
				for i in np.arange(0,len(centres),1):
					mols[i][xlist][ylist].append(GetMolecule('../../../Molecules/./%s' % centres[i]))
				print 'xlist,ylist,zlist',xlist,ylist,zlist
			elif abs(x)+abs(y)+abs(z) <= size:
				# Places molecule of this type if not in centres[0], and (AS DIAMOND) max. 'size' translation vectors from centre unit cell.
				for i in np.arange(0,len(surroundings),1):
					mols[i][xlist][ylist].append(GetMolecule('../../../Molecules/./%s' % surroundings[i]))
					#Move to correct pos'n
					for atom in mols[i][xlist][ylist][zlist]():
						atom().move(x*TV[0]+y*TV[1]+z*TV[2])
				print 'xlist,ylist,zlist',xlist,ylist,zlist
			else:
				for molincell in np.arange(0,len(centres),1):
					mols[molincell][xlist][ylist].append('empty')


#Calculate Properties:

print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'size', size
E0 = np.matrix([0.,0.,0.])

print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'Calc jm'
jm = JMatrix(cutoff=cut)
print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'Calc dips:'
d = get_dipoles(E0=E0,jm=jm._m,cutoff=cut)
print strftime("%a, %d %b %Y %X +0000", gmtime())
Efield = get_electric_field(E0)
potential = get_potential()

print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'dips', d
split_d = split_dipoles_onto_atoms(d)
tot = np.matrix([0.,0.,0.])
for dd in split_d:
	tot += dd

print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'tot', tot
Uqq = np.multiply(get_U_qq(potential=potential),27.211)
print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'Uqq', Uqq
Uqd = np.multiply(get_U_qdip(jm=jm._m,dips=d,Efield=Efield),27.211)
print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'Uqd', Uqd
Udd = np.multiply(get_U_dipdip(jm=jm._m,dips=d.T),27.211)
print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'Udd', Udd
energyev = Udd+Uqd+Uqq
print 'energyev', energyev
energy=energyev/27.211

print strftime("%a, %d %b %Y %X +0000", gmtime())
print 'Making .dat cross sections for gnuplot'

# print dipoles
f = open('%s_size_%s_dipoles.dat' % (name,size), 'w')
for dd in split_d:
	dstr=str(dd)
	f.write(dstr)
	f.write('\n')
f.flush()
f.close()

# new.dats for gnuplot with individual molecules
for x in np.arange(-size,size+1,1):
	xlist=x+size
	for y in np.arange(-size,size+1,1):
		ylist=y+size
		for z in np.arange(-size,size+1,1):
			zlist=z+size
			for molincell in np.arange(0,len(centres),1):
				if mols[molincell][xlist][ylist][zlist] != 'empty':
					f = open('%s_size_%s_%s%s%s_mol%s.dat' % (name,size,xlist,ylist,zlist,molincell), 'w')
					f.write('Molecule Starts:\n')
					for molecule in mols[molincell][xlist][ylist][zlist]():
						  molstr=str(molecule())
						  for atom in Register.GetRegister('Atom'):
							    astr=str(atom)
							    if molstr in astr:
									f.write(astr)
									f.write('\n')
					f.write('Molecule Ends.')
					f.flush()
					f.close()


#REORG ENERGY with other mol in unit cells
reorg=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(len(centres))]
Uqdlist=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(len(centres))]
Uqqlist=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(len(centres))]
#reorg[i][x][y][z] contains this property for charge of central molecule moved to molecule in poition [i] in unit cell at [x][y][z].
#Other matrices same for each peroperty

print 'Moving charges...'
for x in np.arange(-size,size+1,1):
	xlist=x+size
	for y in np.arange(-size,size+1,1):
		ylist=y+size
		for z in np.arange(-size,size+1,1):
			zlist=z+size
			if mols[molincell][xlist][ylist][zlist] != 'empty':
	  			for molincell in np.arange(0,len(centres),1):
					print 'swapping charges of molecule at 0,0,0, (mols[0]([',size,'][',size,'][',size,']) and molecule',molincell,' at ', x, ',', y, ',', z, '.'
					for cen_atom, x_atom in zip(mols[0][size][size][size](), mols[molincell][xlist][ylist][zlist]()):
						  cen_initial=cen_atom()._crg
						  x_initial=x_atom()._crg
						  cen_atom()._crg=x_initial
						  x_atom()._crg=cen_initial

						#Calc energy following reorg

						Efield=get_electric_field(E0)
						newUqd = np.multiply(get_U_qdip(jm=jm._m,dips=d,Efield=Efield),27.211)

						print 'New Uqd (eV): ', newUqd
						NewEnergyev = Udd+newUqd+Uqq
						print 'New Utot (eV): ',NewEnergyev
						reorgeV=newUqd-Uqd
						print 'Reorginisation Energy (', x, ',', y, ',', z, ') (eV): ',reorgeV
						newpotential = get_potential()
						newUqq = np.multiply(get_U_qq(potential=potential),27.211)
						print 'New Uqq (eV) (not used for reorg energy calc as should be same in perfect cryst): ',newUqq
						print 'Change in Uqq (eV): ',(newUqq-Uqq)

						reorg[molincell][xlist][ylist].append(newUqd-Uqd)
						Uqdlist[molincell][xlist][ylist].append(newUqd)
						Uqqlist[molincell][xlist][ylist].append(newUqd-Uqd)

						time=strftime("%Y %b %d %X", gmtime())
						print time

						f = open('%s_size_%s_%s%s%s_mol%s.dat' % (name,size,xlist,ylist,zlist,molincell), 'a')
						f.write('\n\nReorganisation energy:\t%s eV\n' % str(reorgeV) )
						f.write('Polaron Binding\t%s eV\n' % str(Uqd+Udd) )
						f.write('Other Properties:\n')
						f.write ('Date\tname\tcentres[0]\tsurroundings[0]\tsize\tpos_a\tpos_b\tpos_c\tNewEnergyev\treorgeV\tUqd\tUqq\tUdd')
						f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,centres[0],surroundings[0],size,x,y,z,molincell,NewEnergyev,reorgeV,newUqd,newUqq,Udd))
						f.flush()
						f.close()

						f = open('../../../Crystals/reorg_energies.csv', 'a')
						f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,centres[0],surroundings[0],size,x,y,z,molincell,NewEnergyev,reorgeV,newUqd,newUqq,Udd))
						f.flush()
						f.close()

						f = open('./%s_size_%s_new.csv' % (name,size), 'w')
						f.write ('Date\tname\tcentres[0]\tsurroundings[0]\tsize\tpos_a\tpos_b\tpos_c\tNewEnergyev\treorgeV\tUqd\tUqq\tUdd')
						f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (time,name,centres[0],surroundings[0],size,x,y,z,molincell,energyev,reorgeV,Uqd,newUqq,Udd))
						f.flush()
						f.close()

						#Put charge back
						for cen_atom, x_atom in zip(mols[0][size][size][size](), mols[molincell][xlist][ylist][zlist]()):
							  cen_initial=cen_atom()._crg
							  x_initial=x_atom()._crg
							  cen_atom()._crg=x_initial
							  x_atom()._crg=cen_initial


#csv file



f = open('../../../Crystals/crystals.csv', 'a')
f.write ('\n%s\t%s\t%s\tN/A\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,centres[0],surroundings[0],size,energyev,Uqq,Uqd,Udd,tot[0,0],tot[0,1],tot[0,2]))
f.flush()
f.close()

#_new.dat for gnuplot
print 'making .dat'


f = open('%s_size_%s_new.dat' % (name, size), 'w')
f.write('Molecule Starts:\n')
for atom in Register.GetRegister('Atom'):
	astr=str(atom)
	f.write(astr)
	f.write('\n')
f.write('Molecule Ends.\n\n')
f.write ("Crystal Struct\tcentres[0]\tsurroundings[0]\tSize\tEnergy\tUqq\tUqd\tUdd\tDipole Moment")
f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (name,centres[0],surroundings[0],size,energyev,Uqq,Uqd,Udd,tot))
f.flush()
f.close()

'''#Gradual charge movement

f = open('../../../Crystals/reorg_energies_partial.csv', 'a')
f.write ('Date\tname\tcentres[0]\tsurroundings[0]\tsize\tpos_a\tpos_b\tpos_c\tReorgCT0.1\tReorgCT0.2\tReorgCT0.3\tReorgCT0.4\tReorgCT0.5\tReorgCT0.6\tReorgCT0.7\tReorgCT0.8\tReorgCT0.9\tReorgCT1.0\t')
f.flush()
f.close()

f = open('./%s_size_%s_partial.csv' % (name,size), 'w')
f.write ('Date\tname\tcentres[0]\tsurroundings[0]\tsize\tpos_a\tpos_b\tpos_c\tReorgCT0.1\tReorgCT0.2\tReorgCT0.3\tReorgCT0.4\tReorgCT0.5\tReorgCT0.6\tReorgCT0.7\tReorgCT0.8\tReorgCT0.9\tReorgCT1.0\t')
f.flush()
f.close()

print 'Moving charges...'
for x in np.arange(-size,size+1,1):
	xlist=x+size
	for y in np.arange(-size,size+1,1):
		ylist=y+size
		for z in np.arange(-size,size+1,1):
			zlist=z+size
				for molincell in np.arange(0,len(centres),1):
				if mols[molincell][zlist][ylist][xlist] != 'empty':
						f = open('../../../Crystals/reorg_energies_partial.csv', 'a')
						f.write ('Date\tname\tcentres[0]\tsurroundings[0]\tsize\tpos_a\tpos_b\tpos_c\tMolincell\tReorgCT0.1\tReorgCT0.2\tReorgCT0.3\tReorgCT0.4\tReorgCT0.5\tReorgCT0.6\tReorgCT0.7\tReorgCT0.8\tReorgCT0.9\tReorgCT1.0\t')
						f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,centres[0],surroundings[0],size,x,y,z,molincell))
						f.flush()
						f.close()

						f = open('./%s_size_%s_partial.csv' % (name,size), 'a')
						f.write ('Date\tname\tcentres[0]\tsurroundings[0]\tsize\tpos_a\tpos_b\tpos_c\tMolincell\tReorgCT0.1\tReorgCT0.2\tReorgCT0.3\tReorgCT0.4\tReorgCT0.5\tReorgCT0.6\tReorgCT0.7\tReorgCT0.8\tReorgCT0.9\tReorgCT1.0\t')
						f.write ('\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (time,name,centres[0],surroundings[0],size,x,y,z,molincell))
						f.flush()
						f.close()

						for qtrans in np.arange(0.1,1.1,0.1):
							  print 'swapping', qtrans, ' of charge of molecule at 0,0,0, (mols[0]([',size,'][',size,'][',size,']) and molecule',molincell,' at ', x, ',', y, ',', z, '.'
							  cen_ini=[]
							  x_ini=[]
							  for cen_atom, x_atom in zip(mols[0][size][size][size](), mols[molincell][xlist][ylist][zlist]()):
								    print 'cen atom q pre swap', cen_atom()._crg
								    print 'x atom q pre swap', x_atom()._crg
							   		cen_initial=np.float64(cen_atom()._crg)
							     		x_initial=np.float64(x_atom()._crg)
								    cen_ini.append(cen_initial)
								    x_ini.append(x_initial)
								    cen_atom()._crg=(qtrans*x_initial + (1-qtrans)*cen_initial)
								    x_atom()._crg=(qtrans*cen_initial + (1-qtrans)*x_initial)
								    print 'cen atom q after swap', cen_atom()._crg
								    print 'x atom q after swap', x_atom()._crg
							print 'len(cen_ini)',len(cen_ini), 'len(x_ini)', len(x_ini), 'len(mols[0][size][size][size]()),'len(mols[0][size][size][size]()),'len(mols[molincell][xlist][ylist][zlist]()),'len(mols[molincell][xlist][ylist][zlist]())
							  #Calc energy following reorg


							  Efield=get_electric_field(E0)
							  newUqd = np.multiply(get_U_qdip(jm=jm._m,dips=d,Efield=Efield),27.211)

							  print 'New Uqd (eV): ', newUqd
							  NewEnergyev = Udd+newUqd+Uqq
							  print 'New Utot (eV): ',NewEnergyev
							  reorgeV=newUqd-Uqd
							  print 'Reorginisation Energy (', x, ',', y, ',', z, ') (eV): ',reorgeV

							  time=strftime("%Y %b %d %X", gmtime())
							  print time

							  f = open('../../../Crystals/reorg_energies_partial.csv', 'a')
							  f.write ('\t%s' % (reorgeV))
							  f.flush()
							  f.close()

							  f = open('./%s_size_%s_partial.csv' % (name,size), 'a')
							  f.write ('\t%s' % (reorgeV))
							  f.flush()
							  f.close()

							  #Put charge back
							  for cen_atom, cen_initial in zip(mols[0][size][size][size](), cen_ini):
								    cen_atom()._crg=cen_initial
								    print 'cen atom q after swap back', cen_atom()._crg

							  for x_atom, x_initial in zip(mols[molincell][xlist][ylist][zlist](), x_ini):
								    x_atom()._crg=x_initial
								    print 'x atom q after swap back', x_atom()._crg'''
