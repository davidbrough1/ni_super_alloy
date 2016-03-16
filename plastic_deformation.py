import os
import time
import numpy as np
import numpy.ctypeslib as npct
import ctypes

# **** For data schema ****
import collections
import dicttoxml
import xmltodict
from xml.dom.minidom import parseString
#from ngc_nist import NGC

# ++++ Fortran function for plastic deformation ++++
import irreverisble


def stress_strain(elasticity, volume_fraction, prec_stress,
                  SS_stress, T_service):

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

    no_samples = 1
    print "Irreversible:", time.asctime(time.localtime(time.time()))
    Stress_strain, WTN = irreverisble.mechanics(elasticity, elastic_modulus,
                                                volume_fraction, prec_stress,
                                                SS_stress, T_service,
                                                no_samples)
    Stress_strain = np.array(np.trim_zeros(Stress_strain)).reshape(-1, 2)
    return Stress_strain, WTN

if __name__ == '__main__':
    kinetics = []
    with open('xmls/Kinetics.xml') as inpread:
        kinetics = xmltodict.parse(inpread.read())

    volume_fraction = float(kinetics['kinetics']['volumeFraction']['value'])
    ra = float(kinetics['kinetics']['averageRadius']['value'])
    T_service = float(kinetics['kinetics']['serviceTemperature']['value'])
    prec_stress = float(kinetics['kinetics']['precipitateStress']['value'])
    SS_stress = float(kinetics['kinetics']['precipitateStress']['value'])

    PyMKS = []
    with open('xmls/PyMKS.xml') as inpread:
    	PyMKS = xmltodict.parse(inpread.read())

    Elasticity = float(PyMKS['pyMks']['youngsModulus']['value'])