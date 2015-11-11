import os
import time
import numpy as np
import numpy.ctypeslib as npct
import ctypes

# **** For data schema ****
import collections
import dicttoxml
from xml.dom.minidom import parseString
import xmltodict

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


#def evaluation(gen,no_samples,start_samples,variables,no_variables,T_service):

if __name__ == '__main__':

	GA_select = []
	with open('xsd/GA_select.xml') as inpread:
		GA_select = xmltodict.parse(inpread.read())

	gen = int(GA_select['ngcInputs']['noGeneration']) 
	no_samples = int(GA_select['ngcInputs']['noSamples']) 
	start_samples = int(GA_select['ngcInputs']['startSamples']) 

	no_variables = len(GA_select['ngcInputs']['variables']['value'])
	variables = []

	for i in range (no_variables):
		variables.append(float(GA_select['ngcInputs']['variables']['value'][i]))

	no_variables = no_variables/no_samples

	T_service = float(GA_select['ngcInputs']['serviceTemperature']['value'])

# +++++++ Loading TC API function +++++++
	_file = 'TC_OPT.so'
	_path = os.path.join(*(os.path.split(os.path.realpath(__file__))[:-1]+(_file,)))
	C_subs = ctypes.cdll.LoadLibrary(_path)

	EQ_TCAPI = C_subs.API_TC

#	array_inp = npct.ndpointer(dtype=np.double, flags='C')
#	print array_inp
	DoubleArrayType = DoubleArrayType()
	EQ_TCAPI.argtypes = (ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,DoubleArrayType,ctypes.c_double,ctypes.c_int)
	EQ_TCAPI.restype = ctypes.POINTER(ctypes.c_double)


# --------------------------------- initialize model parameters -------------------------------------
	ini_nudensity = float(4E+26)
	no_outputs = int(10)
# --------------------------------- initialize the arrays -------------------------------------
	fitness = np.array([0]*no_samples,dtype=float)
	variables = np.asarray(variables).reshape(no_samples,no_variables)
	Vf_gammaprime = np.array([0]*no_samples,dtype=float)
	csize_gammaprime = np.array([0]*no_samples,dtype=float)
	Elasticity = np.array([0]*no_samples,dtype=float)

#----------------------------------------------------------------------------------------------
#	print "CALPHAD:",time.asctime(time.localtime(time.time()))
#	print gen,no_samples,start_samples,no_variables,variables,ini_nudensity,no_outputs


	microstructure = EQ_TCAPI(gen,no_samples,start_samples,no_variables,variables,ini_nudensity,no_outputs)
	microstructure = np.fromiter(microstructure, dtype=float, count=no_samples*no_outputs).reshape(no_samples,no_outputs)



	# /*/*/*/ volume fraction and C.S. of gamma prime /*/*/*/
	vf = microstructure[:,0]
	ra = microstructure[:,1]

	SS_stress = microstructure[:,3]
	prec_stress = microstructure[:,4]


	Ttest= [T_service for x in xrange(no_samples)]


# --------- Output to XML ---------
	ind_xml = int(0)
	dictKinetics = []
	dictKinetics.append(('composition',({'element':({'chemicalElement':'Al'},{'value':variables[ind_xml,0]},{'unit':'molar fraction'})},{'element':({'chemicalElement':'Cr'},{'value':variables[ind_xml,1]},{'unit':'molar fraction'})})))

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

	Kineticsxml = dicttoxml.dicttoxml(dictKinetics,custom_root='kinetics',attr_type=False)
	Kineticsschema = open('Kinetics.xml', 'w')

	Kineticsschema.write("%s\n" % (parseString(Kineticsxml).toprettyxml()))

	
