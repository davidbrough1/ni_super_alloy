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
GA_select = []
with open('../xmls/GA_select.xml') as inpread:
	GA_select = xmltodict.parse(inpread.read())

gen = int(GA_select['ngcInputs']['noGeneration']) 
no_samples = int(GA_select['ngcInputs']['noSamples']) 
start_samples = int(GA_select['ngcInputs']['startSamples']) 

no_variables = len(GA_select['ngcInputs']['variables']['value'])
variables = []

GA_inputs = []
with open('GA_inputs.xml') as inpread:
        GA_inputs = xmltodict.parse(inpread.read())

input_xml_ID = str(GA_inputs['gaInputs']['iD'])
Search_domain = [['' for x in xrange(2)] for j in range(no_variables)]
for i in range(0,no_variables):
	Search_domain[i][0] = GA_inputs['gaInputs']['searchDomain']['variable'][i]  # number after ['item']
	Search_domain[i][1] = GA_inputs['gaInputs']['searchDomain']['unit'][i]


for i in range (no_variables):
	variables.append(float(GA_select['ngcInputs']['variables']['value'][i]))

no_variables = no_variables/no_samples

T_service = float(GA_select['ngcInputs']['serviceTemperature']['value'])


#---------------------------------------------- Knime start from here
# +++++++ Loading TC API function +++++++
_file = 'TC_OPT.so'
_path = os.path.join(*(os.path.split(os.path.realpath(__file__))[:-1]+(_file,)))
C_subs = ctypes.cdll.LoadLibrary(_path)

EQ_TCAPI = C_subs.API_TC
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

microstructure = EQ_TCAPI(gen,no_samples,start_samples,no_variables,variables,ini_nudensity,no_outputs)

print "volume fraction:", microstructure[0]
print "mean radius:", microstructure[1]
print "yield stress:", microstructure[5]
print "precipitation stress:", microstructure[4]
print "processing time:", microstructure[9]

'''
ind_xml = 0
chem_comp = []
dictKinetics = []
Major_element = 1.
sampleName = str('Gen'+str(gen)+'Sam'+str(ind_xml))


for i in range(no_variables-1):
	chem_comp.append(({'chemicalElement':Search_domain[i][0]},{'value':variables[ind_xml,i]},{'unit':Search_domain[i][1]}))
	Major_element -= variables[ind_xml,i]

print microstructure[0,0] 

# This schema is designed for Ni-based alloy
dictKinetics.append(('gaInputID',input_xml_ID))
chem_comp.append(({'chemicalElement':'Ni'},{'value':Major_element},{'unit':Search_domain[i][1]}))
dictKinetics.append(('composition',chem_comp))
dictKinetics.append(('processTemperature',({'value':variables[ind_xml,2]},{'unit':'kelvin'})))
dictKinetics.append(('serviceTemperature',({'value':T_service},{'unit':'kelvin'})))
dictKinetics.append(('initialNumberDensity',({'value':ini_nudensity},{'unit':'1/m^3'})))
dictKinetics.append(('alloyDensity',({'value':microstructure[ind_xml,2]},{'unit':'kg/m^3'})))
dictKinetics.append(('eqVolumeFraction',{'value':microstructure[ind_xml,6]}))
dictKinetics.append(('interfaceEnergy',({'value':microstructure[ind_xml,7]},{'unit':'J/m^2'})))
dictKinetics.append(('apbEnergy',({'value':microstructure[ind_xml,8]},{'unit':'J/m^2'})))
dictKinetics.append(('processingTime',({'value':microstructure[ind_xml,9]},{'unit':'minutes'})))
dictKinetics.append(('volumeFraction',{'value':microstructure[ind_xml,0]}))
dictKinetics.append(('averageRadius',({'value':microstructure[ind_xml,1]},{'unit':'meter'})))
dictKinetics.append(('solidSolutionStress',({'value':microstructure[ind_xml,3]},{'unit':'MPa'})))
dictKinetics.append(('precipitateStress',({'value':microstructure[ind_xml,4]},{'unit':'MPa'})))
dictKinetics.append(('yieldStress',({'value':microstructure[ind_xml,5]},{'unit':'MPa'})))
dictKinetics = collections.OrderedDict(dictKinetics)

Kineticsxml = dicttoxml.dicttoxml(dictKinetics,custom_root='kinetics',attr_type=False)
filename = str('Kinetics.xml')
Kineticsschema = open(filename, 'w')
Kineticsschema.write("%s\n" % (parseString(Kineticsxml).toprettyxml()))
'''
