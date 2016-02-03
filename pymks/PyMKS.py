import os
import time
import numpy as np
import math
import numpy.ctypeslib as npct
import ctypes

# **** For data schema ****
import collections
import dicttoxml
import xmltodict
from xml.dom.minidom import parseString
#from ngc_nist import NGC

# ---- for PyMKS ----
import pymks
from pymks.datasets import make_delta_microstructures
from pymks.datasets import make_elastic_FE_strain_delta
from pymks.datasets import make_delta_microstructures
from pymks.datasets.elastic_FE_simulation import ElasticFESimulation
from pymks import MKSRegressionModel
from pymks import DiscreteIndicatorBasis


# ************************** PyMKS ****************************
def PyMKSElasticity(volumefraction, radius):

	inter_particle_distance = 1

# --- Unit: GPa ---
#       elastic_modulus = (295., 289.)
#       poissons_ratio = (0.3, 0.3)
#       macro_strain = 0.02

	R_deviation = 0.2
	L_PyMKS = 31

        # +++ Thomas et al., J Mat. Pro. Tech, 177, 2006, 469 +++
#       shear_modulus = 83100*(1.-0.5*(T-300.)/1673.)
        # +++ Huang et al., Mat. Sci. Tech., 23, 2007, 1105 +++
#       shear_modulus = 87416.4-32.73*T+0.00295*T*T
        # +++ Fisher, Scripta Meta., 20, 1986, 279 +++
#       Elastic_Moduli = 295 (gamma) (AVG),  GPa
#       Elastic_Moduli = 289 (gamma_prime) (AVG),  GPa

	elas_gamma = (298.*3./8.)*(1.-0.5*(T_service-300.)/1673.)
	elas_gammaprime = (289.*3./8.)*(1.-0.5*(T_service-300.)/1673.)


	elastic_modulus = (elas_gamma, elas_gammaprime)
	poissons_ratio = (0.3, 0.3)
	macro_strain = 0.02

	X_delta = make_delta_microstructures(n_phases=2, size=(L_PyMKS, L_PyMKS))

	sim = ElasticFESimulation(elastic_modulus=elastic_modulus, poissons_ratio=poissons_ratio, macro_strain=macro_strain)
	sim.run(X_delta)
	y_stress_xx = sim.stress[..., 0]
	y_strain_xx = sim.strain[..., 0]

	basis = DiscreteIndicatorBasis(n_states=2)

	model_strain = MKSRegressionModel(basis=basis)
	model_stress = MKSRegressionModel(basis=basis)

	model_strain.fit(X_delta, y_strain_xx)
	model_stress.fit(X_delta, y_stress_xx)

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	E_eff = np.array([0]*no_samples,dtype=float)
	print volumefraction

	if (no_samples > 1):
		for sam in range(no_samples):

			if(volumefraction[sam] > 0 ):
				volume_fraction = volumefraction[sam]

				ind_phase1_r = []

				p_x = int(L_PyMKS/2)
				p_y = int(L_PyMKS/2)

				iR=-int(L_PyMKS*math.sqrt(volume_fraction)/2)
				iL=int(L_PyMKS*math.sqrt(volume_fraction)/2)+1

	                        # +++ ignore the corner +++
				max_dis = 1.3*L_PyMKS/2.
				for i in range(iR,iL):
					for j in range(iR,iL):
						if ( math.sqrt(i**2.+j**2.) < max_dis) :
							R_p_x = i + p_x
							R_p_y = j + p_y

							Rtemp = []
							if (R_p_x < 0):
								R_p_x += L_PyMKS
							elif (R_p_x > L_PyMKS-1):
								R_p_x -= L_PyMKS

							if (R_p_y < 0):
								R_p_y += L_PyMKS
							elif (R_p_y > L_PyMKS-1):
								R_p_y -= L_PyMKS

							Rtemp = (R_p_x,R_p_y)
							ind_phase1_r.append(Rtemp)

				X = []
				Y1 = [[0 for x in xrange(L_PyMKS)] for x in xrange(L_PyMKS)]

				for i in range(len(ind_phase1_r)):
					Y1[ind_phase1_r[i][0]][ind_phase1_r[i][1]] = 1

				Y = []
				Y.append(Y1)
				X = np.array(Y)

				stress_pred = model_stress.predict(X)
				strain_pred = model_strain.predict(X)

				E_eff[sam] = stress_pred.mean()/strain_pred.mean()

			else:
				E_eff[sam] = 0.



	else:
		volume_fraction = volumefraction

		ind_phase1_r = []

		p_x = int(L_PyMKS/2)
		p_y = int(L_PyMKS/2)

		iR=-int(L_PyMKS*math.sqrt(volume_fraction)/2)
		iL=int(L_PyMKS*math.sqrt(volume_fraction)/2)+1

		# +++ ignore the corner +++
		max_dis = 1.3*L_PyMKS/2.
		for i in range(iR,iL):
			for j in range(iR,iL):
				if ( math.sqrt(i**2.+j**2.) < max_dis) :
					R_p_x = i + p_x
					R_p_y = j + p_y

					Rtemp = []
					if (R_p_x < 0):
						R_p_x += L_PyMKS
					elif (R_p_x > L_PyMKS-1):
						R_p_x -= L_PyMKS

					if (R_p_y < 0):
						R_p_y += L_PyMKS
					elif (R_p_y > L_PyMKS-1):
						R_p_y -= L_PyMKS

					Rtemp = (R_p_x,R_p_y)
					ind_phase1_r.append(Rtemp)

		X = []
		Y1 = [[0 for x in xrange(L_PyMKS)] for x in xrange(L_PyMKS)]

		for i in range(len(ind_phase1_r)):
			Y1[ind_phase1_r[i][0]][ind_phase1_r[i][1]] = 1

		Y = []
		Y.append(Y1)
		X = np.array(Y)

		stress_pred = model_stress.predict(X)
		strain_pred = model_strain.predict(X)

		E_eff[0] = stress_pred.mean()/strain_pred.mean()


	return E_eff



kinetics = []
with open('../xmls/Kinetics.xml') as inpread:
	kinetics = xmltodict.parse(inpread.read())

vf = float(kinetics['kinetics']['volumeFraction']['value'])
ra = float(kinetics['kinetics']['averageRadius']['value'])
T_service = float(kinetics['kinetics']['serviceTemperature']['value'])
no_samples = 1


#-------------------------------------------------------------- Knime starts
Elasticity = np.array([0]*no_samples,dtype=float)

print "PyMKS:",time.asctime(time.localtime(time.time()))
Elasticity = PyMKSElasticity(vf, ra)
Elasticity *= 1000

print Elasticity


ind_xml = int(0)
dictMKS = []
dictMKS.append(('optimumVolumeFraction',{'value':vf}))
dictMKS.append(('optimumRadius',({'value':ra},{'unit':'meter'})))
dictMKS.append(('youngsModulus',({'value':Elasticity[ind_xml]},{'unit':'MPa'})))
dictMKS = collections.OrderedDict(dictMKS)

MKSxml = dicttoxml.dicttoxml(dictMKS,custom_root='pyMks',attr_type=False)
PyMKSschema = open('PyMKS.xml', 'w')
PyMKSschema.write("%s\n" % (parseString(MKSxml).toprettyxml()))

