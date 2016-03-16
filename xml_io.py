import xmltodict
import dicttoxml
from xml.dom.minidom import parseString
import collections


def load_variables(xml_file_path, file_key, value_keys):
    _f = []
    with open(xml_file_path) as inpread:
        _f = xmltodict.parse(inpread.read())
    print tuple([_f[file_key][v]['value'] for v in value_keys])
    return tuple([float(_f[file_key][v]['value']) for v in value_keys])


def save_variables(xml_file_path, file_key, list_of_variables):
    list_of_dicts = [(v[0], ({'value': v[1]}, {'unit': v[2]}))
                     for v in list_of_variables]
    _dict = collections.OrderedDict(list_of_dicts)
    _xml = dicttoxml.dicttoxml(_dict, custom_root=file_key, attr_type=False)
    _schema = open(xml_file_path, 'w')
    _schema.write("%s\n" % (parseString(_xml).toprettyxml()))

if __name__ == '__main__':
    vf, ra, Elasticity = 1, 2, 3
    dictMKS = []
    dictMKS.append(('optimumVolumeFraction', {'value': vf}))
    dictMKS.append(('optimumRadius', ({'value': ra},
                   {'unit': 'meter'})))
    dictMKS.append(('youngsModulus', ({'value': Elasticity},
                   {'unit': 'MPa'})))
    print dictMKS
    dictMKS = collections.OrderedDict(dictMKS)
    print '\n', dictMKS
    MKSxml = dicttoxml.dicttoxml(dictMKS, custom_root='pyMks',
                                 attr_type=False)
    save_variables('PyMKS.xml', 'pyMks', dictMKS)
    # print '\n', MKSxml
    # PyMKSschema = open('PyMKS.xml', 'w')
    # PyMKSschema.write("%s\n" % (parseString(MKSxml).toprettyxml()))

    # xml_file_path = 'xmls/Kinetics.xml'
    # file_key = 'kinetics'
    # value_keys = ['volumeFraction', 'averageRadius', 'serviceTemperature',
    #               'precipitateStress', 'solidSolutionStress']
    # print load_variables(xml_file_path, file_key, value_keys)

    # kinetics = []
    # with open('xmls/Kinetics.xml') as inpread:
    #     kinetics = xmltodict.parse(inpread.read())
    # tmp = []
    # tmp.append(float(kinetics['kinetics']['volumeFraction']['value']))
    # tmp.append(float(kinetics['kinetics']['averageRadius']['value']))
    # tmp.append(float(kinetics['kinetics']['serviceTemperature']['value']))
    # tmp.append(float(kinetics['kinetics']['precipitateStress']['value']))
    # tmp.append(float(kinetics['kinetics']['precipitateStress']['value']))
    # print tuple(tmp)
