from wl1989stoich import *
import numpy as np
import pandas as pd

def kdCalc_original(components, T, P = None):
    """
    This is the original Weaver and Langmuir
    calculation of Kd's. T in Kelvin.
    """
    keys = ['CaAl2O4', 'NaAlO2', 'MgO', 'FeO', 'CaSiO3', 'TiO2']
    cpx = {key:0 for key in keys}
    plg = {key:0 for key in keys}
    ol = {key:0 for key in keys}
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

def kdCalc_reynolds(components, T, P = None):
    """
    This uses Reynold's Thesis Calculations.
    """
    keys = ['CaAl2O4', 'NaAlO2', 'MgO', 'FeO', 'CaSiO3', 'TiO2']
    cpx = {key:0 for key in keys}
    plg = {key:0 for key in keys}
    ol = {key:0 for key in keys}
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
    
def kdCalc_langmuir1992(components, T, P):
    """
    This uses Langmuir Et Al. 1992 Calculations.
    P in bars, T in Kelvin.
    """
    keys = ['CaAl2O4', 'NaAlO2', 'MgO', 'FeO', 'CaSiO3', 'TiO2']
    cpx = {key:0 for key in keys}
    plg = {key:0 for key in keys}
    ol = {key:0 for key in keys}
    anorthite=components['CaAl2O4']/(components['CaAl2O4']+1.5*components['NaAlO2'])
    plg['CaAl2O4'] = np.power(10.,(2446./T) - 1.122 + (0.2562*anorthite))
    plg['NaAlO2'] = np.power(10.,((3195. + (3283*anorthite) + (0.0506*P))/T) - (1.885*anorthite) -2.3715)
    cpx['MgO'] = np.power(10.,(((3798. + (0.021*P))/T) - 2.28))
    cpx['FeO'] = cpx['MgO']*np.power(10.,-0.6198)
    cpx['CaSiO3'] = np.power(10.,((1783. + (0.0038*P))/T) -0.753)
    cpx['CaAl2O4'] = np.power(10.,((2418. + (0.068*P))/T) -2.3)
    cpx['NaAlO2'] = np.power(10.,((5087. + (0.073*P))/T) - 4.48)
    cpx['TiO2'] = np.power(10.,((1034. + (0.053*P))/T) - 1.27)
    ol['MgO'] = np.power(10.,((6921./T) + (3.4*components['NaAlO2']/2) + 
                (6.3*components['KO5']) + (0.00001154*P) - 3.27)/2.3)
    ol['FeO'] = ol['MgO']*(np.power(10., -0.523))
    kd = {'cpx':cpx, 'ol':ol, 'plg':plg}
    return kd