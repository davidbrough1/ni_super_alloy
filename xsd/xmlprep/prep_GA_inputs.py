import collections
import dicttoxml
from xml.dom.minidom import parseString

#no_samples = 1
#no_variables = 3
#L_sample_range = [ 0.10, 0.10, 1123.]
#H_sample_range = [ 0.25, 0.20, 1473.]
#T_service = 1123
#generations = 1
#Memory = [8 for x in xrange(no_variables)]
#Memory = [6,6,8]

# --- read inputs from dictionary in .dat file ---
GA_inputs = []
with open('GA_inputs.dat','r') as inf:
	GA_inputs = eval(inf.read())

no_samples = GA_inputs['no_samples']      # number of samples in 1 generation
no_variables = GA_inputs['no_variables']  # number of variables in 1 sample
generations = GA_inputs['generations']    # number of generations
T_service = GA_inputs['T_service']        # service temperature

L_sample_range = [0 for x in xrange(no_variables)]
H_sample_range = [0 for x in xrange(no_variables)]
Memory = [0 for x in xrange(no_variables)]

for i in range(0,no_variables):
	L_sample_range[i] = GA_inputs['L_sample_range'][i]
	H_sample_range[i] = GA_inputs['H_sample_range'][i]
	Memory[i] = GA_inputs['Memory'][i]

# --- save GA_inputs to XML file ---
dictGAinputs = []
dictGAinputs.append(('noSamples',{'value':no_samples}))
dictGAinputs.append(('noVariables',{'value':no_variables}))
dictGAinputs.append(('noGenerations',{'value':generations}))
dictGAinputs.append(('serviceTemperature',({'value':T_service},{'unit':'kelvin'})))
dictGAinputs.append(('lsampleRange',{'value':L_sample_range}))
dictGAinputs.append(('hsampleRange',{'value':H_sample_range}))
dictGAinputs.append(('binMemory',{'value':Memory}))
dictGAinputs = collections.OrderedDict(dictGAinputs)

GAinputsxml = dicttoxml.dicttoxml(dictGAinputs,custom_root='gaInputs',attr_type=False)
GAinputsxmlschema = open('GA_inputs.xml', 'w')
GAinputsxmlschema.write("%s\n" % (parseString(GAinputsxml).toprettyxml()))

