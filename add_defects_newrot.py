import Basicelements.Rotation as rot

						#Using Rodrigues rotation formula for rotation about vector http://math.stackexchange.com/questions/142821/matrix-for-rotation-around-a-vector

						# ROT by alpha around a, gets ac face on ab, => want new crystal above in c => find max posn of atoms in bottom half crystal in c , min pos. in c of new crystal.
						# BUT only need to search through one unit cell at max disp. in c. mols in old crystal, most extreme in -b in new crystal
						# Move new crystal so centred in a and b (dot each posn with a, take mean, subtract mean * versor in a from each posn, same for b).
						# Move new crystal such that c.(mincposncry2) - c.(maxcposncry1) / root(c^2) = (7 Bohr?)
						# Move by c * (7 - [c.(mincposncry2) - c.(maxcposncry1)] )/ root(c^2)

NEED:

-Separate sections to build each crystal half
-Variables:
	-nperp Number of unit cells in (one) direction perp to face
	-npar Number of unit cells in (two) directions par to face

OR (actually probs better, as in some instances one TV much larger than another
	ntop[0,1,2] size in a, b, c of top cryst
	nbot[0,1,2] size in a, b, c of bot cryst

BUILD up at top, down at bottom

	perptop (=0,1,2) dimension (a,b, or c) perpendicular to face of top crystal which orients beside bottom crystal
	perpbot (=0,1,2)


perptop=1
perpbot=0

possindices=[0,1,2]

partop=[ v for v in possindices if not v == perptop ]
parbot[ v for v in possindices if not v == perpbot ]
# list of indices for vectors parallel to the face on top and bottom

#Crystal goes from -
for i in partop:
	mintop[i]=-ntop[i]

for i in parbot:
	minbot[i]=-nbot[i]

mintop[perptop]=0
minbot[perpbot]=0

if perptop != perpbot:
	anglerot = [ v for v in possindices if not v == perptop and not v == perpbot ][0]
	# angle of rotation to get from top to bottom orient is that joining the two vectors perp. fo face. ie. that with the other of the three indices 0,1,2
else:
	anglerot=perpbot
	# If both are same orientation, rotation about this angle will not affect this




Charge will always be on mol on bottom half posn 0,0,0 for dip calculation.

	



molposnslist = [None] * 3

if perptop != perpbot:
	Rotmat_forfaces=rot.Rotation(rodrigues=([TV[anglerot,0],TV[anglerot,1],TV[anglerot,2],angle[anglerot])
	# Rotation matrix which rotates top half of crystal to make two vectors we want perp. to face parallel to one another.
else:
	Rotmat_forfaces=rot.Rotation(rodrigues=([TV[anglerot,0],TV[anglerot,1],TV[anglerot,2],0.)
	# Don't do this rotation if crystals already aligned. In this instance only do manual rotation (if any) following

Rotmat_givenfaces=rot.Rotation(rodrigues=([TV[perpbot,0],TV[perpbot,1],TV[perpbot,2],ROTATION)
# Rotation matrix which rotates top half of crystal by angle ROTATION about vector perp to face joining crystals.

moltop=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(len(centres))]
molbot=[ [[[] for i in range((2*size)+1)] for i in range((2*size)+1)] for __ in xrange(len(centres))]

# MAKE BOTTOM HALF OF CRYSTAL (todo)
for molposnlist[0] in [0,nbot[0]+1,1]:
# looping over a posns, replace xlist
	for molposnlist[1] in [0,nbot[1]+1,1]:
		for molposnlist[2] in [0,nbot[2]+1,1]:
			if x == 0 and y == 0 and z == 0:
				#central mol

				for i in np.arange(0,len(centres),1):
					mols[i][xlist][ylist].append(GetMolecule('../../../Molecules/./%s' % centres[i]))
				print 'xlist,ylist,zlist',xlist,ylist,zlist

				#move to correct pos'n in unit cell
				#for atom in mols[2][xlist][ylist][zlist]():
				#	atom().move(0.5*TV[0]+0.5*TV[2])
				#for atom in mols[3][xlist][ylist][zlist]():
				#	atom().move(0.5*TV[0]+0.5*TV[2])

			elif abs(x)+abs(y)+abs(z) <= size:
				# Places molecule of this type if not in centres[0], and (AS DIAMOND) max. 'size' translation vectors from centre unit cell.
				for i in np.arange(0,len(surroundings),1):
					mols[i][xlist][ylist].append(GetMolecule('../../../Molecules/./%s' % surroundings[i]))
					#Move to correct pos'n
					for atom in mols[i][xlist][ylist][zlist]():
						atom().move(x*TV[0]+y*TV[1]+z*TV[2])
				print 'xlist,ylist,zlist',xlist,ylist,zlist

				#move to correct pos'n in unit cell
				#for atom in mols[2][xlist][ylist][zlist]():
				#	atom().move(0.5*TV[0]+0.5*TV[2])
				#for atom in mols[3][xlist][ylist][zlist]():
				#	atom().move(0.5*TV[0]+0.5*TV[2])



# MAKE TOP HALF OF CRYSTAL:
for molposnlist[0] in [0,ntop[0]+1,1]:
# looping over a posns, replace xlist
	for molposnlist[1] in [0,ntop[1]+1,1]:
		for molposnlist[2] in [0,ntop[2]+1,1]:
PLACE

# ROTATE TOP HALF OF CRYSTAL:
for molposnlist[0] in [0,ntop[0]+1,1]:
# looping over a posns, replace xlist
	for molposnlist[1] in [0,ntop[1]+1,1]:
		for molposnlist[2] in [0,ntop[2]+1,1]:
					moltop[molposnlist[0]][molposnlist[1]][molposnlist[2]].rotate(Rotmat_forfaces)
					Rotmat_givenfaces

# Get CENTRES of molecules wanting to align. Want mol in position nmols/2 in two directions perp to face, position 0 in direction par to face.
align_mol_top[perptop]=0
align_mol_top[parptopa]=(ntop[partopa]+1)/2
align_mol_top[parptopb]=(ntop[partopb]+1)/2

centop=molstop[align_mol_top[0],align_mol_top[1],align_mol_top[2]]._get_com

align_mol_bot[perptop]=0
align_mol_bot[parptop[0]]=(ntop[parbot[0]]+1)/2
align_mol_bot[parptop[1]]=(ntop[parbot[1]]+1)/2

cenbot=molsbot[align_mol_bot[0],align_mol_bot[1],align_mol_bot[2]]._get_com

# find min (top) and max (bot) atom posns along vector perpendicular to face in order to set a separation of SEP (3.5A?)

# top

posns_along_perp_top=[]

for atom in molstop[align_mol_top[0],align_mol_top[1],align_mol_top[2]]:
	posns_along_perp_top.append(np.dot(TVs[perpbot],atom()._pos))

minalongperp_aftercentre_top=(min(posns_along_perp_top)-np.dot(TVs[perpbot],centop)

# bot

posns_along_perp_bot=[]

for atom in molsbot[align_mol_bot[0],align_mol_bot[1],align_mol_bot[2]]:
	posns_along_perp_bot.append(np.dot(TVs[perpbot],atom()._pos))

maxalongperp_aftercentre_bot=(max(posns_along_perp_bot)-np.dot(TVs[perpbot],cenbot)

# MOVE all molecules in bottom half by -cenbot -maxalongperp_aftercentre_bot*TVs[perpbot]
# MOVE all molecules in top half by -centop -minalongperp_aftercentre_top*TVs[perpbot] + SEP*perpbot


# This is minimum posn value of posns in top crystal in direction perpendicular to face joining crystals (following centring this mol on 0,0,0. Translation will occur later so that translation making this happen for this mol will happen to whole crystal). Need to find this distance to get correct separation.


							meanpostop=mean(posns[])
							# Required to align centres of crystals
							minalongperp_aftercentre=(min(posns_along_perp)-np.dot(dimension_perp_face_top,meanpostop)
							# This is minimum posn value of posns in top crystal in direction perpendicular to face joining crystals (following centring this mol on 0,0,0. Translation will occur later so that translation making this happen for this mol will happen to whole crystal). Need to find this distance to get correct separation.

						if xlist==0 and	 ylist==0 and zlist==0:
							# EDIT such that 0 in dimension along which vector joining crys will point, 1 in other dims
							posns_along_perp=[]
							posns=[]
							for atom in mol:
								posns.append(atom()._pos)
								posns_along_perp.append(np.dot(dimension_perp_face_bot,atom()._pos)
							meanposbot=mean(posns[])
							maxalongperp=(max(posns_along_perp))-np.dot(dimension_perp_face_bot,meanposbot)
							# This is maximum posn value of posns in bottom crystal in direction perpendicular to face joining crystals (following centring this mol on 0,0,0. Translation will occur later so that translation making this happen for this mol will happen to whole crystal). Need to find this distance to get correct separation.


for mol in... 
	




a=5.9277/0.5291772109217
b=7.881/0.5291772109217
c=29.184/0.5291772109217

alpha=90*(pi/180)
beta=92.4434*(pi/180)
gamma=90*(pi/180)

#TVs, TV[0,1,2] are the three translation vectors.
TV=matrix_to_cartesian.T

print 'TV:'
print TV

print 'TV[0,1,2] are the three translation vectors.'


						# Rotate z to x

						atom().pos[0]=zpos
						atom().pos[2]=-xpos

						# Rotate z to y

						atom().pos[1]=zpos
						atom().pos[2]=-ypos


						zposns_rotated.append(atom().pos[2])


						print 'POLARISABILITY FLIP:'
						print 'atom()._pol', atom()._pol

						#Vector containing diagonal elements of polarisability tensor
						atom_pol_diags=np.matrix( [[ atom()._pol[0,0], atom()._pol[1,1], atom()._pol[2,2] ]])
						print 'Vector containing diagonal elements of polarisability tensor, atom_pol_diags', atom_pol_diags

						# Rotate z to x
						atom()._pol=Polarizability(noniso =[ [ np.sqrt(atom_pol_diags_flipped[2,0]**2),0.0,0.0], [0.0,np.sqrt(atom_pol_diags_flipped[1,0]**2),0.0], [0.0,0.0,np.sqrt(atom_pol_diags_flipped[0,0]**2)] ])		

						# Rotate z to y
						atom()._pol=Polarizability(noniso =[ [ np.sqrt(atom_pol_diags_flipped[0,0]**2),0.0,0.0], [0.0,np.sqrt(atom_pol_diags_flipped[2,0]**2),0.0], [0.0,0.0,np.sqrt(atom_pol_diags_flipped[1,0]**2)] ])	
	
						print 'NEW type(atom()._pol)', type(atom()._pol)

						print 'OLD type(atom()._pol)', type(atom()._pol)
						print 'type(atom()._pol[0,0])', type(atom()._pol[0,0])	
						print 'type(atom()._pol[0,1])', type(atom()._pol[0,1])	

				for molincell in np.arange(0,nmolsincell,1):
					for atom in mols[molincell][xlist][ylist][zlist]():
						atom().move(np.vector[0.,0.,-min(zposns_rotated) + max(z_posns_original) + (3.5//0.5291772109217))



