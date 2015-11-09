#@@@@@@@@@@@@@ Micro-GA @@@@@@@@@@@@@@@@@
# V1.0 07-31-2014
# Shengyen Li @ NIST
# shengyen.li@nist.gov
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import random
import numpy as np
from numpy import array
import ctypes

import operator
import math
import sys
import os
import random

import matplotlib.pyplot as plt

import time

# ---- for PyMKS ----
import pymks
from pymks.datasets import make_delta_microstructures
from pymks.datasets import make_elastic_FE_strain_delta
from pymks.datasets import make_delta_microstructures
from pymks.datasets.elastic_FE_simulation import ElasticFESimulation
from pymks import MKSRegressionModel
from pymks import DiscreteIndicatorBasis

# ++++ Fortran function for plastic deformation ++++
import irreverisble

# **** For data schema ****
import collections
import dicttoxml
import xmltodict
from xml.dom.minidom import parseString

# ************************ Micro-GA ************************
# -- Initialization --
def geno_pheno(i,temp_rand):
	return L_sample_range[i]+(H_sample_range[i]-L_sample_range[i])*temp_rand/d_size[i]

def GA_samp_ini():
	Pheno_ = []
	Geno_ = ''
	Geno_temp = []
	for i in range(no_variables):
		temp_rand = random.randint(0,d_size[i]-1)

		temp_bin = str(bin(temp_rand))
		geno_length=len(temp_bin)
		temp_bin2 = temp_bin[2:]
		Mem_diff = Memory[i] - geno_length + 2

		for j in range(Mem_diff):
			temp_bin2 = '0' + temp_bin2

		Geno_temp.append(temp_bin2)
		Geno_ = Geno_ + temp_bin2
		Pheno_.append(geno_pheno(i,temp_rand))

	return (Geno_,Pheno_)

# ---- Cross Over ----
def Cross_Over(RK_1,RK_2):

	while True:
		temp_rand_0 = random.randint(0,tot_Memory-1)
		temp_rand_1 = random.randint(0,tot_Memory-1)
		if temp_rand_0 != temp_rand_1:
			break

	rand_0 = min(temp_rand_0,temp_rand_1)
	rand_1 = max(temp_rand_0,temp_rand_1)

	New_Geno_0=RK_1[0:rand_0] + RK_2[rand_0:rand_1] + RK_1[rand_1:tot_Memory]
	New_Geno_1=RK_2[0:rand_0] + RK_1[rand_0:rand_1] + RK_2[rand_1:tot_Memory]

#	print RK_1,RK_2
#	print rand_0,rand_1
#	print New_Geno_0,New_Geno_1

	return (New_Geno_0,New_Geno_1)

# ---- Mutation ----
def Mutation(Geno_1,Geno_2,mu_rate):

	G_Mu = []
	G_Mu.append(Geno_1)
	G_Mu.append(Geno_2)

	temp_rand_0 = []
	temp_rand_1 = []

	for i in range(mu_rate):

		temp_rand_0 = random.randint(1,2) - 1         # list starts from 0
		temp_rand_1 = random.randint(0,tot_Memory-1)  # in range of 0,tot_Memory-1

#		print temp_rand_1
#		print "012345678901234567890123456789"

		temp_list_1 = G_Mu[temp_rand_0][0:temp_rand_1]

		mu_number = abs( int(G_Mu[temp_rand_0][temp_rand_1]) - 1 )
		temp_list = G_Mu[temp_rand_0][0:temp_rand_1] + str(mu_number)
		temp_list += G_Mu[temp_rand_0][temp_rand_1+1:]

#		print G_Mu[temp_rand_0]
		G_Mu[temp_rand_0] = temp_list
#		print G_Mu[temp_rand_0]

	Geno_1 = G_Mu[0]
	Geno_2 = G_Mu[1]
	
	return (Geno_1,Geno_2)

# ---- Convergence ----
# /// Only top two samples are considered in convergence examination ///
#
def Convergence(current_sample,generation):

	convergence = 0
	conv_rate = 0.
	no_reini = 0

	for k in range(tot_Memory):
		if(current_sample[1][k] == current_sample[3][k]):
			convergence += 1
#			Evolution.write("%s, %s, %s\n" % (current_sample[1][k],current_sample[3][k],convergence)) 

#	Evolution.write("\n") 
	conv_rate = (convergence*1.0 / tot_Memory)

	if(conv_rate > 0.98):
		current_sample = []
		current_sample = initialize()
		no_reini = generation
		print "rein",generation,no_reini

#		convergence = 0
#		for k in range(tot_Memory):
#			if(current_sample[1][k] == current_sample[3][k]):
#				convergence += 1
#				Evolution.write("%s, %s, %s\n" % (current_sample[1][k],current_sample[3][k],convergence)) 
#		sys.exit("re-initiate")

	return current_sample,no_reini

# **********************************************************
def initialize():
	current_sample = []
	if Calculating_model == 1:

		for i in range(no_variables):
			d_size.append(2**Memory[i]-1)

#		Pheno_string = []
#		Geno_string = []

		for i in range(no_samples):
			temp_G = ''
			temp_P = []
			temp_G, temp_P = GA_samp_ini()

#			Pheno_string.append(temp_P)
#			Geno_string.append(temp_G)
			current_sample.append(temp_P)
			current_sample.append(temp_G)
#		print Pheno_string,Geno_string
#		print current_sample
	else:
		print "NOT implemented yet"

	return current_sample


def gaussian(mode,x,mu,delta):
	dis = np.exp(-np.power(x - mu, 2.) / (2 * np.power(delta, 2.)) )
	probability = dis / math.sqrt(2*math.acos(-1.)*np.power(delta, 2.))

	# mode==1, return Gaussian distribution
	# else, return Gaussian probability
	if mode == 1:
		re_value = dis
	else:
		re_value = probability

	return re_value



def evaluation(gen,start_samples,variables):

	test_model = 0
# =0 full calculation (CALPHAD+PyMKS+IrTM)
# =1 Mechanistic models (PyMKS+IrTM)

# --------------------------------- initialize model parameters -------------------------------------
	ini_nudensity = 4E+26
	no_outputs = 10

# --------------------------------- initialize the arrays -------------------------------------
	fitness = np.array([0]*no_samples,dtype=float)
	variables = np.asarray(variables).reshape(no_samples,no_variables)
	Vf_gammaprime = np.array([0]*no_samples,dtype=float)
	csize_gammaprime = np.array([0]*no_samples,dtype=float)
	Elasticity = np.array([0]*no_samples,dtype=float)

#	print start_samples,len(variables[0])
	print variables
#----------------------------------------------------------------------------------------------

	if (test_model == 0):

		print "CALPHAD:",time.asctime(time.localtime(time.time()))
		Time.write("CALPHAD: %s\n" % (time.asctime(time.localtime(time.time()))) )
		microstructure = EQ_TCAPI(gen,no_samples,start_samples,len(variables[0]),variables,ini_nudensity,no_outputs)
		microstructure = np.fromiter(microstructure, dtype=float, count=no_samples*no_outputs).reshape(no_samples,no_outputs)

	# /*/*/*/ volume fraction and C.S. of gamma prime /*/*/*/
		vf = microstructure[:,0]
		ra = microstructure[:,1]

		SS_stress = microstructure[:,3]
		prec_stress = microstructure[:,4]

#		print microstructure
#		print "Solid Solution Strengthening",SS_stress
#		print "Precipitation",prec_stress

#		Ttest=variables[:,3]
		Ttest= [T_service for x in xrange(no_samples)]
#		print irreverisble.mechanics.__doc__

# --------- Output to XML ---------
		ind_xml = 0
		dictKinetics = []
		dictKinetics.append(('composition',({'element':({'chemicalElement':'Al'},{'value':variables[ind_xml,0]},{'unit':'molar fraction'})},{'element':({'chemicalElement':'Cr'},{'value':variables[ind_xml,1]},{'unit':'molar fraction'})})))
#		dictKinetics.append(('composition',('element':{'value':variables[ind_xml,0]})))
		dictKinetics.append(('processTemperature',({'value':variables[ind_xml,2]},{'unit':'kelvin'})))
		dictKinetics.append(('serviceTemperature',({'value':T_service},{'unit':'kelvin'})))
		dictKinetics.append(('initialNumberDensity',({'value':ini_nudensity},{'unit':'1/m^3'})))
		dictKinetics.append(('alloyDensity',({'value':microstructure[ind_xml,2]},{'unit':'kg/m^3'})))
		dictKinetics.append(('eqVolumeFraction',{'value':microstructure[ind_xml,6]}))
		dictKinetics.append(('interfaceEnergy',({'value':microstructure[ind_xml,7]},{'unit':'J/m^2'})))
		dictKinetics.append(('apbEnergy',({'value':microstructure[ind_xml,8]},{'unit':'J/m^2'})))
		dictKinetics.append(('processingTime',({'value':microstructure[ind_xml,9]},{'unit':'minutes'})))
		dictKinetics.append(('optimumVolumeFraction',{'value':microstructure[ind_xml,0]}))
		dictKinetics.append(('optimumRadius',({'value':microstructure[ind_xml,1]},{'unit':'meter'})))
		dictKinetics.append(('solidSolutionStress',({'value':microstructure[ind_xml,3]},{'unit':'MPa'})))
		dictKinetics.append(('precipitateStress',({'value':microstructure[ind_xml,4]},{'unit':'MPa'})))
		dictKinetics.append(('yieldStress',({'value':microstructure[ind_xml,5]},{'unit':'MPa'})))
		dictKinetics = collections.OrderedDict(dictKinetics)
#		print(dictKinetics)

		Kineticsxml = dicttoxml.dicttoxml(dictKinetics,custom_root='kinetics',attr_type=False)
		Kineticsschema = open('Kinetics.xml', 'w')
#		Kineticsschema.write("%s\n" % (Kineticsxml))
		Kineticsschema.write("%s\n" % (parseString(Kineticsxml).toprettyxml()))
#
# the present PyMKS model does not take radius into account
#
		print "PyMKS:",time.asctime(time.localtime(time.time()))
		Time.write("PyMKS: %s\n" % (time.asctime(time.localtime(time.time()))) )
		Elasticity = PyMKSElasticity(vf, ra)
		Elasticity *= 1000

# --------- Output to XML ---------
# *** Output PyMKS results as .XML ***
# --- ./envs/default/lib/python2.7/site-packages/dicttoxml.py has been modified (06/12/2015) ---
#
		ind_xml = 0
		dictMKS = []
		dictMKS.append(('optimumVolumeFraction',{'value':vf[ind_xml]}))
		dictMKS.append(('optimumRadius',({'value':ra[ind_xml]},{'unit':'meter'})))
		dictMKS.append(('youngsModulus',({'value':Elasticity[ind_xml]/1000.},{'unit':'MPa'})))
		dictMKS = collections.OrderedDict(dictMKS)
#		dictMKS = collections.OrderedDict([('optimumVolumeFraction',{'value':vf[ind_xml]}),('optimumRadius',({'value':ra[ind_xml]},{'unit':'meter'})),('youngsModulus',({'value':Elasticity[ind_xml]/1000.},{'unit':'MPa'}))])

		MKSxml = dicttoxml.dicttoxml(dictMKS,custom_root='pyMks',attr_type=False)
		PyMKSschema = open('PyMKS.xml', 'w')
		PyMKSschema.write("%s\n" % (parseString(MKSxml).toprettyxml()))

		print "Irreversible:",time.asctime(time.localtime(time.time()))
		Time.write("Irreversible: %s\n" % (time.asctime(time.localtime(time.time()))) )
		Stress_strain = irreverisble.mechanics(Elasticity,elastic_modulus,vf,prec_stress,SS_stress,Ttest,no_samples)
		Stress_strain = np.array(Stress_strain).reshape(no_samples,3)

# --------- Output to XML ---------
		ind_xml = 0
		dictPlastics = []
		dictPlastics.append(('solidSolutionStress',({'value':SS_stress[ind_xml]},{'unit':'MPa'})))
		dictPlastics.append(('precipitateStress',({'value':prec_stress[ind_xml]},{'unit':'MPa'})))
		dictPlastics.append(('yieldStress',({'value':microstructure[ind_xml,5]},{'unit':'MPa'})))
		dictPlastics.append(('youngsModulus',({'value':Elasticity[ind_xml]},{'unit':'MPa'})))
		dictPlastics.append(('optimumVolumeFraction',{'value':vf[ind_xml]}))
		dictPlastics.append(('serviceTemperature',({'value':T_service},{'unit':'Kelvin'})))
		dictPlastics.append(('ultimateTensileStress',({'value':Stress_strain[ind_xml,0]},{'unit':'MPa'})))
		dictPlastics.append(('ultimateTensileStrain',{'value':Stress_strain[ind_xml,1]}))
		dictPlastics.append(('workToNecking',({'value':Stress_strain[ind_xml,2]},{'unit':'MPa'})))
		dictPlastics = collections.OrderedDict(dictPlastics)

		Plasticsxml = dicttoxml.dicttoxml(dictPlastics,custom_root='plasticDeformation',attr_type=False)
		Plasticsschema = open('Plasticdeformation.xml', 'w')
		Plasticsschema.write("%s\n" % (parseString(Plasticsxml).toprettyxml()))

		fitness = Stress_strain[:,2]

		print "Finish:",time.asctime(time.localtime(time.time()))
		Time.write("Finish: %s\n" % (time.asctime(time.localtime(time.time()))) )

	else:
		vf = [0.158167]
		ra = [5]

		Elasticity = PyMKSElasticity(vf, ra)
		Elasticity *= 1000
		print "Results", Elasticity

		prec_stress = [220., 330., 140., 120.]
		SS_stress = [210., 100., 140., 200.]
		Ttest= [923., 923., 923., 923.]
#		Ttest= [300., 300., 350., 923.]

		Stress_strain = irreverisble.mechanics(Elasticity,elastic_modulus,vf,prec_stress,SS_stress,Ttest,no_samples)
		Stress_strain = np.array(Stress_strain).reshape(no_samples,3)
		print Stress_strain

		fitness = Stress_strain[:,0]


	return fitness

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class DoubleArrayType:
	def from_param(self, param):
		typename = type(param).__name__
		if hasattr(self, 'from_' + typename):
			return getattr(self, 'from_' + typename)(param)
		elif isinstance(param, ctypes.Array):
			return param
		else:
			raise TypeError("Can't convert %s" % typename)

    # Cast from array.array objects
	def from_array(self, param):
		if param.typecode != 'd':
			raise TypeError('must be an array of doubles')
			ptr, _ = param.buffer_info()
		return ctypes.cast(ptr, ctypes.POINTER(ctypes.c_double))

    # Cast from lists/tuples
	def from_list(self, param):
		val = ((ctypes.c_double)*len(param))(*param)
		return val

	from_tuple = from_list

    # Cast from a numpy array
	def from_ndarray(self, param):
		return param.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# ************************** PyMKS ****************************
def PyMKSElasticity(volumefraction, radius):

	inter_particle_distance = 1

# --- Unit: GPa ---
#	elastic_modulus = (295., 289.)
#	poissons_ratio = (0.3, 0.3)
#	macro_strain = 0.02

# the smaller R_deviation means the deviation is smaller
	R_deviation = 0.2

	L = L_PyMKS

#	L = 31
#	size=(L, L)
# --- plot the results ----
#	fig, ax = plt.subplots(figsize=(12, 9))
#	ax.set_xlim(0, L-1)
#	ax.set_ylim(0, L-1)
# R_precipitate = math.sqrt(volume_fraction*(L**2)/math.pi)
# -------------------------

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Initialize ^^^
#	X_delta = make_delta_microstructures(n_phases=2, size=(L, L))

#	sim = ElasticFESimulation(elastic_modulus=elastic_modulus, poissons_ratio=poissons_ratio, macro_strain=macro_strain)
#	sim.run(X_delta)
#	y_stress_xx = sim.stress[..., 0]
#	y_strain_xx = sim.strain[..., 0]

#	basis = DiscreteIndicatorBasis(n_states=2)

#	model_strain = MKSRegressionModel(basis=basis)
#	model_stress = MKSRegressionModel(basis=basis)

#	model_strain.fit(X_delta, y_strain_xx)
#	model_stress.fit(X_delta, y_stress_xx)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	E_eff = np.array([0]*no_samples,dtype=float)

	for sam in range(no_samples):
#		print volumefraction[sam]

		if(volumefraction[sam] > 0 ):
			volume_fraction = volumefraction[sam]

			ind_phase1_r = []

			p_x = int(L/2)
			p_y = int(L/2)

			iR=-int(L*math.sqrt(volume_fraction)/2)
			iL=int(L*math.sqrt(volume_fraction)/2)+1

			# +++ ignore the corner +++
			max_dis = 1.3*L/2.
			for i in range(iR,iL):
				for j in range(iR,iL):
					if ( math.sqrt(i**2.+j**2.) < max_dis) :
						R_p_x = i + p_x
						R_p_y = j + p_y

						Rtemp = []
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

			X = []
			Y1 = [[0 for x in xrange(L)] for x in xrange(L)]

			for i in range(len(ind_phase1_r)):
				Y1[ind_phase1_r[i][0]][ind_phase1_r[i][1]] = 1

			Y = []
			Y.append(Y1)
			X = array(Y)

			stress_pred = model_stress.predict(X)
			strain_pred = model_strain.predict(X)

			E_eff[sam] = stress_pred.mean()/strain_pred.mean()

		else:
			E_eff[sam] = 0.


#	ax.imshow(X[0].swapaxes(0, 1), cmap=plt.cm.gray,interpolation='none')
#	plt.show()

	return E_eff


# *************************************************************

if __name__ == '__main__':

# Two models: GA and enumerate (the method is not decided)
# can be chosen
# =1 - GA, for higher order system
# =2 - EM, for low dimensional issues

# -- Inputs --
	global no_samples
	global no_variables
	global L_sample_range
	global H_sample_range
	global Memory
	global tot_Memory
	global no_reini
	global Calculating_model

# -- Physical Inputs --
#	global samples_in_pop
#	global inp_phase_name
#	global inp_XT

	global L_PyMKS
	global T_service
	global elastic_modulus
	global poissons_ratio

	Calculating_model = 1
	no_reini = 0
	no_reini_p = 0

	L_PyMKS = 21


# *** Physical Inputs ***
#	no_samples = 1
#	no_variables = 3
#	L_sample_range = [ 0.10, 0.10, 1123.]
#	H_sample_range = [ 0.25, 0.20, 1473.]
#	T_service = 1123
#	generations = 1
#	Memory = [8 for x in xrange(no_variables)]
#	Memory = [6,6,8]

# --- read inputs from xml file ---
	GA_inputs = []
	with open('xsd/GA_inputs.xml') as inpread:
		GA_inputs = xmltodict.parse(inpread.read())

	no_samples = int(GA_inputs['gaInputs']['noSamples']['value'])         # number of samples in 1 generation
	no_variables = int(GA_inputs['gaInputs']['noVariables']['value'])     # number of variables in 1 sample
	generations = int(GA_inputs['gaInputs']['noGenerations']['value'])    # number of generations
	T_service = float(GA_inputs['gaInputs']['serviceTemperature']['value']) # service temperature

	L_sample_range = [0 for x in xrange(no_variables)]
	H_sample_range = [0 for x in xrange(no_variables)]
	Memory = [0 for x in xrange(no_variables)]

	for i in range(0,no_variables):
		L_sample_range[i] = float(GA_inputs['gaInputs']['lsampleRange']['value']['item'][i])
		H_sample_range[i] = float(GA_inputs['gaInputs']['hsampleRange']['value']['item'][i])
		Memory[i] = int(GA_inputs['gaInputs']['binMemory']['value']['item'][i])


# +++ Numerical Parameters +++
	History = []
	keep_history = 1
	History_Check = 0
                       # =1, while the sample has been evaluated, generate another one
                       # it benefits for the complicated objective calculations
                       # but, the simple one
	H_fitness = -1E10

	tot_Memory = 0
	for i in range(no_variables):
		tot_Memory += Memory[i]        # *** size of geno-string ***

	tot_sample = 2**tot_Memory
	mu_rate = int(tot_Memory/no_samples)  # *** mutation rate ***

	d_size = []                           # *** mash size ***

	current_sample = []
	evaluated_sample = []

	if keep_history == 1:
		results = open('results.dat', 'w')
		Evolution = open('fitness.dat', 'w')

	Best = open('best.dat', 'w')
	Time = open('time.dat', 'w')

# +++++++ Loading TC API function +++++++
	_file = 'TC_OPT.so'
	_path = os.path.join(*(os.path.split(os.path.realpath(__file__))[:-1]+(_file,)))
	C_subs = ctypes.cdll.LoadLibrary(_path)

	EQ_TCAPI = C_subs.API_TC
	DoubleArrayType = DoubleArrayType()
	EQ_TCAPI.argtypes = (ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,DoubleArrayType,ctypes.c_double,ctypes.c_int)

	EQ_TCAPI.restype = ctypes.POINTER(ctypes.c_double)


# $$$$ Initialize PyMKS $$$$
 # This is not the best way of using PyMKS:
 # the shear modulus and other properties are function of temperature and composition
 # but, in this function, in order to reduce the regression frequency
 # PyMKS is initialized in the main program

# --- Unit: GPa ---
	# +++ Thomas et al., J Mat. Pro. Tech, 177, 2006, 469 +++
#	shear_modulus = 83100*(1.-0.5*(T-300.)/1673.)
        # +++ Huang et al., Mat. Sci. Tech., 23, 2007, 1105 +++
#	shear_modulus = 87416.4-32.73*T+0.00295*T*T
	# +++ Fisher, Scripta Meta., 20, 1986, 279 +++
#	Elastic_Moduli = 295 (gamma) (AVG),  GPa
#	Elastic_Moduli = 289 (gamma_prime) (AVG),  GPa

	elas_gamma = (298.*3./8.)*(1.-0.5*(T_service-300.)/1673.)
	elas_gammaprime = (289.*3./8.)*(1.-0.5*(T_service-300.)/1673.)


	elastic_modulus = (elas_gamma, elas_gammaprime)
#	elastic_modulus = (298.,289.)
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


# $$$$$$$$$$$$$$$$$$$$$$$$$$

# ----- Initial Sample Preparation -----
	current_sample = initialize()

	tot_eav_samples = 0

	for gen in range(generations):

		print "==== Generation: ",gen
#		if gen%10000 == 1:
#			print "==== Generation: ",gen

# -- check the calculated samples
#		if no_reini-no_reini_p == 1:
#			Evolution.write("3**** %s\n" % (current_sample[1])) 
#			Evolution.write("4**** %s\n" % (current_sample[3])) 
#			sys.exit("re-initiate")
#		no_reini_p = no_reini

		# ----- Evaluation -----
		if (gen == 0) or (no_reini>0):
			start = 0
			evaluated_sample = []
		else:
			start = 2

			if Calculating_model == 1:
				l = 2
			else:
				l = 1

			k = len(evaluated_sample)-1

			for i in range(k,1,-1):
				del evaluated_sample[i]

# Evaluation ==========
		temp_sample = np.array(np.asarray(current_sample).reshape(no_samples,2)[:,0]).tolist()

		start_samples = start
		fitness_arr = [0 for x in xrange(no_samples)]
		fitness_arr = evaluation(gen,start_samples,temp_sample)

#		for i in range(no_samples):
#			print "results: ",fitness_arr[i]
# =====================
		for i in range(start,no_samples):

			if Calculating_model == 1:
				j=2*i
			else:
				j=i

			fitness = 0
			fitness = fitness_arr[i]


			temp_array = []
			if Calculating_model == 1:
				temp_array.append(current_sample[j])
				temp_array.append(current_sample[j+1])
				k=2
			else:
				temp_array.append(current_sample[j])
				k=1

			temp_array.append(fitness)
			evaluated_sample.append(temp_array)
			tot_eav_samples += 1

#		print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
#		print evaluated_sample

		if (History_Check == 1) :
			History.append(evaluated_sample)  # [Generation][no_samples][0-2]

		evaluated_sample = sorted(evaluated_sample, key=operator.itemgetter(k), reverse=True)

#		for i in range(no_samples):
		Evolution.write("%s, %s, " % (gen,evaluated_sample[0][2]))


		avg_fitness = 0.
		for i in range(no_samples):
			avg_fitness += evaluated_sample[i][2]
		avg_fitness /= no_samples
		Evolution.write("%s\n" % (avg_fitness))


		if keep_history == 1:
#			results.write("==== Generations: %s - # %s ==== %s ====\n" % (gen,tot_eav_samples,avg_fitness))
			for i in range(no_samples):
				if(evaluated_sample[-1]>0):
					results.write("%s\n" % evaluated_sample[i])

		if Calculating_model == 1:
			i=2
		else:
			i=1

		if (evaluated_sample[0][i] > H_fitness):
			Best.write("%s, %s\n" % (tot_eav_samples,evaluated_sample[0]) )
			H_fitness = evaluated_sample[0][i]
			print "Generation: ",gen,H_fitness
		# ----- Evolution -----
		if Calculating_model == 1:
			current_sample = []

			current_sample.append(evaluated_sample[0][0])  #RK_1, 0-Pheno
			current_sample.append(evaluated_sample[0][1])  #RK_1, 1-Geno
			current_sample.append(evaluated_sample[1][0])  #RK_2, 2-Pheno
			current_sample.append(evaluated_sample[1][1])  #RK_2, 3-Geno

			for sam in range(2,no_samples,2):
				# ----- Decoding Coding & History Check -----
#				brk_out = 1  revised on 11/19/2014
				brk_out = 0

				while True:
					New_Geno_1, New_Geno_2 = Cross_Over(current_sample[1],current_sample[3])
					New_Geno_1, New_Geno_2 = Mutation(New_Geno_1, New_Geno_2,mu_rate)
	
					if (History_Check == 1) :
						no_his_checked = 0
						for i in range(gen+1):
							for j in range(no_samples):
								if (New_Geno_1 == History[i][j][1]) or (New_Geno_2 == History[i][j][1]):
									brk_out = 0
									no_his_checked += 1
								else:
									brk_out = 1
					else:
						brk_out = 1

					if brk_out == 1:
						break

					if (no_his_checked > int(tot_sample*0.8)):
						# ----- Re_initialize -----
						no_his_checked = 0
						current_sample = []
						current_sample = initialize()
						break


				New_Pheno_1 = []
				New_Pheno_2 = []
				temp_New_Pheno_1 = 0.
				temp_New_Pheno_2 = 0.
				temp = 0

				for i in range(no_variables):
					if i == 0:
						temp_New_Pheno_1 = int(New_Geno_1[0:Memory[i]],2)
						temp_New_Pheno_2 = int(New_Geno_2[0:Memory[i]],2)
						temp = temp + Memory[i]
					else:
						temp_New_Pheno_1 = int(New_Geno_1[temp:(temp+Memory[i])],2)
						temp_New_Pheno_2 = int(New_Geno_2[temp:(temp+Memory[i])],2)
						temp = temp + Memory[i]

					temp_New_Pheno_1 = geno_pheno(i,temp_New_Pheno_1)
					temp_New_Pheno_2 = geno_pheno(i,temp_New_Pheno_2)

					New_Pheno_1.append(temp_New_Pheno_1)
					New_Pheno_2.append(temp_New_Pheno_2)

				current_sample.append(New_Pheno_1)
				current_sample.append(New_Geno_1)
				current_sample.append(New_Pheno_2)
				current_sample.append(New_Geno_2)

			# ----- Convergence -----

#			Evolution.write("1---- %s\n" % (current_sample[1]))
#			Evolution.write("2---- %s\n" % (current_sample[3]))

			current_sample, no_reini = Convergence(current_sample,gen)
#			print no_reini,gen

		else:
			print "Not Implemented!!"

# Test Gaussian Function
#	for i in range(30):
#		x = gaussian(2,i,0,2)
#		print i,x
