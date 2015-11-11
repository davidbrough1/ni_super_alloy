import os
import time
import numpy as np
import numpy.ctypeslib as npct
import ctypes

# **** For data schema ****
import collections
import dicttoxml
from xml.dom.minidom import parseString
#from ngc_nist import NGC

# ++++ Fortran function for plastic deformation ++++
import irreverisble



if __name__ == '__main__':

        kinetics = []
        with open('Kinetics.xml') as inpread:
                kinetics = xmltodict.parse(inpread.read())

        vf = int(GA_select['kinetics']['optimumVolumeFraction'])
        ra = int(GA_select['kinetics']['optimumRadius'])
        T_service = float(GA_select['kinetics']['serviceTemperature']['value'])
	prec_stress = 
	SS_stress =
	

        PyMKS = []
        with open('PyMKS.xml') as inpread:
                PyMKS = xmltodict.parse(inpread.read())

        Elasticity = float(PyMKS['pyMks']['youngsModulus'])


	print "Irreversible:",time.asctime(time.localtime(time.time()))
	Stress_strain = irreverisble.mechanics(Elasticity,elastic_modulus,vf,prec_stress,SS_stress,T_service,no_samples)
	Stress_strain = np.array(Stress_strain).reshape(no_samples,3)

# --------- Output to XML ---------
	ind_xml = int(0)
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

