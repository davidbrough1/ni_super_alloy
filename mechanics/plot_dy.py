# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colorbar as clb
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.animation as animation
import time
from mpl_toolkits.mplot3d import Axes3D
import os

lab_s = ['$\sigma_{P}$, MPa']
p_alpha = 1

fig = plt.figure(figsize=(9, 6))

def animate(j):

	plt.clf()
	ax = fig.add_subplot(1,1,1)

	comp = []	
	for i in range(12,-1,-1):
		fname = 'SS'+`i`+'.dat'

		if os.path.isfile(fname):
			comp = np.loadtxt(fname)
#			print i,comp

			linecol = 'b'
			if i==0:
				linecol = 'r'

			plt.plot(comp[:,0], comp[:,1], c=linecol, linewidth=5.0)

	ax.set_xlabel(u'Strain, %', fontsize=30, labelpad=15)
	ax.set_ylabel(u'Stress, MPa', fontsize=30, labelpad=15)
	ax.tick_params(axis='x', labelsize=20, pad = 10)
	ax.tick_params(axis='y', labelsize=20, pad = 10)

	ax.grid(True)
	fig.tight_layout()

#	ax.set_xlim(0, 10)
#	ax.set_ylim(400, 700)


ani = animation.FuncAnimation(fig,animate, interval=1000)
plt.show()

