from wl1989stoich import *
import numpy as np


def kdCalc_original(components, T):
    """
    This is the original Weaver and Langmuir
    calculation of Kd's. Should double check the
    cpx CaSiO3 calculation.
    """
    cpx = {key:0 for key in components.keys()}
    plg = {key:0 for key in components.keys()}
    ol = {key:0 for key in components.keys()}
    anorthite=components['CaAl2O4']/(components['CaAl2O4']+1.5*components['NaAlO2'])
    plg['CaAl2O4'] = np.power(10.,((2446./T) - (1.122+0.2562*anorthite)))
    plg['NaAlO2'] = np.power(10.,(((3195.+(3283.*anorthite))/T) - (2.318 + (1.885*anorthite))))
    cpx['MgO'] = np.power(10.,((3798./T) - 2.28))
    cpx['FeO'] = 0.24*cpx['MgO']
    cpx['CaSiO3'] = np.power(10.,((1783./T)-0.753))
    cpx['CaAl2O4'] = np.power(10.,((2418./T)-2.3))
    cpx['NaAlO2'] = np.power(10.,((5087./T)-4.48))
    cpx['TiO2'] = np.power(10.,((1034./T)-1.27))
    ol['MgO'] = np.power(10.,((2715./T)-1.158))
    ol['FeO'] = np.power(10.,((4230./T)-2.741))
    kd = {'cpx':cpx, 'ol':ol, 'plg':plg}
    return kd

def kdCalc(components, T):
    """
    This uses Reynold's Thesis Calculations.
    """
    cpx = {key:0 for key in components.keys()}
    plg = {key:0 for key in components.keys()}
    ol = {key:0 for key in components.keys()}
    anorthite=components['CaAl2O4']/(components['CaAl2O4']+1.5*components['NaAlO2'])
    plg['CaAl2O4'] = 0.99*np.power(10.,((2446./T) - (1.122+0.2562*anorthite)))
    plg['NaAlO2'] = 0.92*np.power(10.,(((3195.+(3283.*anorthite))/T) - (2.318 + (1.885*anorthite))))
    cpx['MgO'] = np.power(10.,((3798./T) - 2.28))
    cpx['FeO'] = 0.24*cpx['MgO']
    cpx['CaSiO3'] = np.power(10.,((1783./T)-0.753))
    cpx['CaAl2O4'] = np.power(10.,((2418./T)-2.3))
    cpx['NaAlO2'] = np.power(10.,((5087./T)-4.48))
    cpx['TiO2'] = np.power(10.,((1034./T)-1.27))
    ol['MgO'] = np.power(10.,((3386./T)-1.599))
    ol['FeO'] = np.power(10.,((4230./T)-2.707))
    kd = {'cpx':cpx, 'ol':ol, 'plg':plg}
    return kd