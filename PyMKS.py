import numpy as np
import matplotlib.pyplot as plt

import pymks
from pymks.tools import draw_microstructures
from pymks.datasets import make_delta_microstructures

from pymks.datasets import make_elastic_FE_strain_delta
from pymks.datasets import make_delta_microstructures
from pymks.tools import draw_microstructure_strain
from pymks.datasets.elastic_FE_simulation import ElasticFESimulation

from pymks import MKSRegressionModel
from pymks import DiscreteIndicatorBasis

import matplotlib.pyplot as plt
import math
import sys

from numpy import array
import random


elastic_modulus = (295., 289)
poissons_ratio = (0.3, 0.3)
volume_fraction = 0.2

macro_strain = 0.02

L = 151
size=(L, L)

fig, ax = plt.subplots(figsize=(12, 9))
ax.set_xlim(0, L-1)
ax.set_ylim(0, L-1)

#R_precipitate = math.sqrt(volume_fraction*(L**2)/math.pi)
R_precipitate = 10.

X_delta = make_delta_microstructures(n_phases=2, size=(L, L))

sim = ElasticFESimulation(elastic_modulus=elastic_modulus, poissons_ratio=poissons_ratio, macro_strain=macro_strain)
sim.run(X_delta)
y_stress_xx = sim.stress[..., 0]
y_strain_xx = sim.strain[..., 0]

basis = DiscreteIndicatorBasis(n_states=2)

model_strain = MKSRegressionModel(basis=basis)
model_stress = MKSRegressionModel(basis=basis)

model_strain.fit(X_delta, y_strain_xx)
model_stress.fit(X_delta, y_stress_xx)


# ---- square phase in the middle of the domain ----
X = []
Y1 = []

for j in range(L):
	Z = [] 
	for i in range(L):
		dis = 0.
		dis = math.sqrt((i-0.5*L)**2.+(j-0.5*L)**2.)
        
		if (dis < R_precipitate):
			Z.append(1) 
		else:
			Z.append(0) 
	Y1.append(Z) 

Y = []
Y.append(Y1)
X = array(Y)

stress_pred = model_stress.predict(X)
strain_pred = model_strain.predict(X)

total_result = []
R_deviation = 0.3


# ++++ Gaussian probability of precipitate radius ++++
Gau_R_LB = -int(R_precipitate*R_deviation)
Gau_R_UB = int(R_precipitate*R_deviation)+1
Gau_no_samples = Gau_R_UB - Gau_R_LB

Gau_R_dis = []
Gau_delta = 1.5
Gau_mu = 0

Gau_temp = []
prob = 0.
prob_tot = 0.
R_temp = 0.

for R_dis in range(Gau_R_LB,Gau_R_UB):
	R_temp = int(R_precipitate + R_dis)
	prob = np.exp(-np.power(R_dis+Gau_mu, 2.) / (2 * np.power(Gau_delta, 2.)) ) / math.sqrt(2*math.pi*np.power(Gau_delta, 2.))

	Gau_temp = np.array((R_temp,prob))
	Gau_R_dis=np.append(Gau_R_dis,Gau_temp)

Gau_R_dis = np.reshape(Gau_R_dis, (-1,2))
prob_tot = Gau_R_dis.sum(axis=0)


Gau_R_dis[0][1] = Gau_R_dis[0][1] / prob_tot[1]

for i in range(1,Gau_no_samples):
	Gau_R_dis[i][1] = Gau_R_dis[i-1][1] + Gau_R_dis[i][1] / prob_tot[1]

UB_Gau = Gau_R_dis[-1][1]


vf_in_domain = 0.
no_particles = 0

R_Gau = []
while (vf_in_domain <= volume_fraction):
	R_ind = random.uniform(0,UB_Gau)
    
	if ((R_ind > 0.) and (R_ind <= Gau_R_dis[0][1])):
		R_Gau.append(Gau_R_dis[0][0])
	else:
		for i in range(1,Gau_no_samples):
			if ((R_ind > Gau_R_dis[i-1][1]) and (R_ind <= Gau_R_dis[i][1])):
				R_Gau.append(Gau_R_dis[i][0])
				break

	vf_in_domain += math.pi*(R_Gau[-1]**2.)/(L**2.)
	no_particles += 1
# ++++++++++++

# **** Distribution of the particles ****
inter_particle_distance = 2.

ind_phase1_r = []

p_center = [[] for i in range(no_particles)]

p_x = random.randint(0,L-1)
p_y = random.randint(0,L-1)
p_center[0].append(p_x)
p_center[0].append(p_y)


iR=-int(R_Gau[0])
iL=int(R_Gau[0])+1


for i in range(iR,iL):
	for j in range(iR,iL):
		R_p_x = i + p_x
		R_p_y = j + p_y

		R_test = math.sqrt((R_p_x-p_x)**2.+(R_p_y-p_y)**2.)

		Rtemp = []
		if (R_test <= R_Gau[0]):

			if (R_p_x < 0):
				R_p_x += L
			elif (R_p_x > L-1):
				R_p_x -= L

			if (R_p_y < 0):
				R_p_y += L
			elif (R_p_y > L-1):
				R_p_y -= L


			Rtemp = (R_p_x,R_p_y)
			ind_phase1_r.append(Rtemp)

_iterations = 0
ind_var = 1
while (ind_var < no_particles):
	p_x = random.randint(0,L-1)
	p_y = random.randint(0,L-1)
	NA = 0
	_iterations += 1

	if(_iterations > 100):
		print "++++ Volume fraction of L12 is too high ++++"
		sys.exit()


    
	for i in range(0,ind_var):

		if( abs(p_x-p_center[i][0]) > (L/2+1) ): 
			disx1 = L-abs(p_x-p_center[i][0])
#			print "A"
		else:
			disx1 = p_x-p_center[i][0]
#			print "AA"

		if( abs(p_y-p_center[i][1]) > (L/2+1) ):
			disy1 = L-abs(p_y-p_center[i][1])
#			print "B"
		else:
			disy1 = p_y-p_center[i][1]
#			print "BB"


		R_test = math.sqrt(disx1**2.+disy1**2.)

#		print ind_var
#		print p_x,p_center[i][0]
#		print p_y,p_center[i][1]
#		print disx1,disy1
#		print R_test,(R_Gau[ind_var] + R_Gau[i] + inter_particle_distance)

		if (R_test < (R_Gau[ind_var] + R_Gau[i] + inter_particle_distance) ):
			NA = 1
			break

	if ( NA == 0 ):
		p_center[ind_var].append(p_x)
		p_center[ind_var].append(p_y)
        
		iR=-int(R_Gau[ind_var])
		iL=int(R_Gau[ind_var])+1
        
		for i in range(iR,iL):
			for j in range(iR,iL):
				R_p_x = i + p_x
				R_p_y = j + p_y

				R_test = math.sqrt((R_p_x-p_x)**2.+(R_p_y-p_y)**2.)

				Rtemp = []
				if (R_test <= R_Gau[ind_var]):

					if (R_p_x < 0):
						R_p_x += L
					elif (R_p_x > L-1):
						R_p_x -= L

					if (R_p_y < 0):
						R_p_y += L
					elif (R_p_y > L-1):
						R_p_y -= L


					Rtemp = (R_p_x,R_p_y)
					ind_phase1_r.append(Rtemp)

		ind_var += 1
		_iterations = 0

# ************

X = []
Y1 = [[0 for x in xrange(L)] for x in xrange(L)]

for i in range(len(ind_phase1_r)):
	Y1[ind_phase1_r[i][0]][ind_phase1_r[i][1]] = 1

Y = []
Y.append(Y1)
X = array(Y)

stress_pred = model_stress.predict(X)
strain_pred = model_strain.predict(X)

temp_total = (vf_in_domain,R_precipitate,stress_pred.mean(),strain_pred.mean())
total_result.append(temp_total)
        

total_result = array(total_result)

print total_result
print "Effective E: ",stress_pred.mean()/strain_pred.mean()
print elastic_modulus[0]*(1.-vf_in_domain)+elastic_modulus[1]*vf_in_domain


draw_microstructure_strain(X[0], strain_pred[0])
ax.imshow(X[0].swapaxes(0, 1), cmap=plt.cm.gray,interpolation='none')
plt.show()





