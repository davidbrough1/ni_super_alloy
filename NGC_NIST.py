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
with open('xmls/GA_select.xml') as inpread:
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


print variables,T_service


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
print "mean radius (m):", microstructure[1]
print "yield stress (MPa):", microstructure[5]
print "precipitation stress (MPa):", microstructure[4]
print "processing time (min):", microstructure[9]

