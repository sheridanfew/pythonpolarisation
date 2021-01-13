# BUILD crystal in orginal orientation down, flipped up. Ensure connecting surfaces of mols centred on same point.

molposnslist = [None] * 3

for molposnlist[0] in ...


			zposns_original=[]
						zposns_rotated.append(atom().pos[2])

			zposns_rotated=[]
				for molincell in np.arange(0,nmolsincell,1):
					for atom in mols[molincell][xlist][ylist][zlist]():
						print 'atom()._pos', atom()._pos

						xpos=atom()._pos[0]
						ypos=atom()._pos[1]
						zpos=atom()._pos[2]

						#Using Rodrigues rotation formula for rotation about vector http://math.stackexchange.com/questions/142821/matrix-for-rotation-around-a-vector

						# ROT by alpha around a, gets ac face on ab, => want new crystal above in c => find max posn of atoms in bottom half crystal in c , min pos. in c of new crystal.
						# BUT only need to search through one unit cell at max disp. in c. mols in old crystal, most extreme in -b in new crystal
						# Move new crystal so centred in a and b (dot each posn with a, take mean, subtract mean * versor in a from each posn, same for b).
						# Move new crystal such that c.(mincposncry2) - c.(maxcposncry1) / root(c^2) = (7 Bohr?)
						# Move by c * (7 - [c.(mincposncry2) - c.(maxcposncry1)] )/ root(c^2)

						# ac face on ab here
						orient=[aconab,...
						dimension_perp_face_bot=[2]
						dimension_perp_face_top=[1]

						#dimension_perp_face: 0,1,2 correspond to a,b,c. This refers to the vector perpendicular to the face joining the crystal in initial orietation (bottom) with the rotated (top) crystal


					   ux=TV[0,0]
						uy=TV[0,1]
						uz=TV[0,2]

						rotangle=alpha

					   W=np.matrix( [[0.0,-uz,uy],[uz,0.0,-ux][-uy,ux,0.0]] )

						I=np.matrix(np.eye(3))

					   ROTMAT=I + (np.sin(rotangle)*W) + (2*((npsin(rotangle/2)**2))*(W**2))

						atom()._pos=ROTMAT*atom()._pos

						print 'atom()._pos after flip', atom()._pos
						
						print 'POLARISABILITY FLIP:'
						print 'atom()._pol', atom()._pol

						#Vector containing diagonal elements of polarisability tensor
						atom_pol_diags=np.matrix( [[ atom()._pol[0,0], atom()._pol[1,1], atom()._pol[2,2] ]])
						print 'Vector containing diagonal elements of polarisability tensor, atom_pol_diags', atom_pol_diags
					
						atom_pol_diags_rotated=ROTMAT*atom_pol_diags

						atom()._pol=Polarizability(noniso =[ [ np.sqrt(atom_pol_diags_rotated[0,0]**2),0.0,0.0], [0.0,np.sqrt(atom_pol_diags_rotated[1,0]**2),0.0], [0.0,0.0,np.sqrt(atom_pol_diags_rotated[0,2]**2)] ])		
		
						if xlist==1 and	 ylist==1 and zlist==1:
							posns_along_perp=[]
							posns=[]
							for atom in mol:
								posns.append(atom()._pos)
								posns_along_perp.append(np.dot(dimension_perp_face_top,atom()._pos))
							meanpostop=mean(posns[])
							# Required to align centres of crystals
							minalongperp_aftercentre=(min(posns_along_perp)-np.dot(dimension_perp_face_top,meanpostop)
							# This is minimum posn value of posns in top crystal in direction perpendicular to face joining crystals (following centring this mol on 0,0,0. Translation will occur later so that translation making this happen for this mol will happen to whole crystal). Need to find this distance to get correct separation.

						if xlist==0 and	 ylist==0 and zlist==0:
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



