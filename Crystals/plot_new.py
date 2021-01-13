# -*- coding: utf-8 -*-
arrowsizefactor=10

from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from math import *
from itertools import product, combinations
#fig = plt.figure(figsize=plt.figaspect(4.0)*4.00) #Adjusts the aspect ratio and enlarges the figure (text does not enlarge)

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

#import importlib
#my_module = importlib.import_module("package.path.%s" % module_name)

class Arrow3D(FancyArrowPatch):
	zorder=1
	def __init__(self, xs, ys, zs, *args, **kwargs):
		FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
		self._verts3d = xs, ys, zs

	def draw(self, renderer):
		xs3d, ys3d, zs3d = self._verts3d
		xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M, )
		self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
		FancyArrowPatch.draw(self, renderer)

from protneuttest_properties import *

norm = mpl.colors.Normalize(vmin=min(reorgseV), vmax=max(reorgseV))
cmap = cm.jet
m = cm.ScalarMappable(norm=norm, cmap=cmap)

fig = plt.figure(str(name))
ax = fig.gca(projection='3d')
ax.set_aspect("equal")

arrow=[[] for i in range(len(posns)) ]

xlist=[]
ylist=[]
zlist=[]
molcentres=[]

for mol in range(0,len(posns),1):
	for i in range(0,len(posns[mol]),1):
		arrow[mol].append(Arrow3D([posns[mol][i][0]-(arrowsizefactor*dips[mol][i][0]), posns[mol][i][0]+(arrowsizefactor*dips[mol][i][0])],[posns[mol][i][1]-(arrowsizefactor*dips[mol][i][1]), posns[mol][i][1]+(arrowsizefactor*dips[mol][i][1])],
[posns[mol][i][2]-(arrowsizefactor*dips[mol][i][2]), posns[mol][i][2]+(arrowsizefactor*dips[mol][i][2])], mutation_scale=20, lw=1, arrowstyle="-|>", zorder=1, color=m.to_rgba(reorgseV[mol])))
		ax.add_artist(arrow[mol][i])
		#draw a point
		#ax.scatter(posns[molincell][xlist][ylist][zlist][i][0],posns[molincell][xlist][ylist][zlist][i][1],posns[molincell][xlist][ylist][zlist][i][2],color=m.to_rgba(reorgseV[molincell][xlist][ylist][zlist]),marker='.',s=0.02)
		xlist.append(posns[mol][i][0])
		ylist.append(posns[mol][i][1])
		zlist.append(posns[mol][i][2])
	molcentres.append([sum(x[0] for x in posns[mol])/len(posns[mol]),sum(x[1] for x in posns[mol])/len(posns[mol]),sum(x[2] for x in posns[mol])/len(posns[mol])])
	#molcentres[2][0,1,2] are centre of molecule 2 in x,y,z
	ax.text(molcentres[mol][0], molcentres[mol][1], molcentres[mol][2], str(reorgseV[mol]) + 'eV' , bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

meanpos=np.array([np.mean(xlist),np.mean(ylist),np.mean(zlist)])
maxpos=np.array([max(xlist),max(ylist),max(zlist)])
minpos=np.array([min(xlist),min(ylist),min(zlist)])

ax.scatter(maxpos[0],maxpos[1],maxpos[2],color=m.to_rgba(reorgseV[0]),marker='.',s=0.02)
ax.scatter(minpos[0],minpos[1],minpos[2],color=m.to_rgba(reorgseV[0]),marker='.',s=0.02)
#ax.text(1, 1, 1, "red", color='red')
#ax.text(0,0,0, 'Blah' , bbox={'facecolor':'white', 'alpha':0.5, 'pad':10,'zorder':10}, zorder=9)

#ax.plot(maxpos[0],maxpos[1],maxpos[2], 'w')
#ax.plot(minpos[0],minpos[1],minpos[2], 'w')

plt.axis('off')

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.normalize(vmin=-0.1, vmax=2.15))
# fake up the array of the scalar mappable. Urgh...
sm._A = []
plt.colorbar(sm)

#ax.view_init(elev=0,azim=0)

TVnames=['a','b','c']

#for TV in range(0,1,3):

TV=2

theta=(180/pi)*acos(TVs[TV,2]/sqrt(TVs[TV,0]**2+TVs[TV,1]**2+TVs[TV,2]**2))

phi=(180/pi)*atan(TVs[TV,1]/(TVs[TV,0]+0.001))

print 'elev=',(90.0 - theta), 'azim=', phi
ax.view_init(elev=(90.0 - theta), azim=phi)

	#img_name= str(name) + '_',TVnames[TV],'.png'
	#pl.savefig(img_name,dpi=200,bbox_inches='tight')
	#fig.clf()
	#ax.cla()
	#plotfigure()

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.normalize(vmin=-20, vmax=10))
# fake up the array of the scalar mappable. Urgh...
sm._A = []
plt.colorbar(sm)

plt.show()

