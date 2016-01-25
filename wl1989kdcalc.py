from wl1989stoich import *
import numpy as np
import pandas as pd

def kdCalc_original(components, T, P = None):
    """
    This is the original Weaver and Langmuir
    calculation of Kd's. T in Kelvin.
    """
    keys = ['MgO', 'FeO', 'TiO2', 'CaAl2O4', 'NaAlO2', 'CaSiO3', 'PO52', 'MnO', 'KAlO2']
    cpx = {key:0. for key in keys}
    plg = {key:0. for key in keys}
    ol = {key:0. for key in keys}
    anorthite=components['CaAl2O4']/(components['CaAl2O4']+1.5*components['NaAlO2'])
    plg['CaAl2O4'] = np.power(10.,((2446./T) - (1.122+0.2562*anorthite)))
    plg['NaAlO2'] = np.power(10.,(((3195.+(3283.*anorthite))/T) - (2.318 + (1.885*anorthite))))
    plg['KAlO2'] = 0.15
    plg['MnO'] = 0.031
    plg['PO52'] = 0.1
    cpx['MgO'] = np.power(10.,((3798./T) - 2.28))
    cpx['FeO'] = 0.24*cpx['MgO']
    cpx['CaSiO3'] = np.power(10.,((1783./T)-0.753))
    cpx['CaAl2O4'] = np.power(10.,((2418./T)-2.3))
    cpx['NaAlO2'] = np.power(10.,((5087./T)-4.48))
    cpx['TiO2'] = np.power(10.,((1034./T)-1.27))
    cpx['KAlO2'] = 0.007
    cpx['MnO'] = 1.02
    cpx['PO52'] = 0.05
    ol['MgO'] = np.power(10.,((2715./T)-1.158))
    ol['FeO'] = np.power(10.,((4230./T)-2.741))
    ol['KAlO2'] = 0.001
    ol['MnO'] = 1.63
    ol['PO52'] = 0.2
    kd = {'cpx':cpx, 'ol':ol, 'plg':plg}
    return kd

def kdCalc_reynolds(components, T, P = None):
    """
    This uses Reynold's Thesis Calculations.
    """
    keys = ['MgO', 'FeO', 'TiO2', 'PO52', 'MnO', 'CaAl2O4', 'NaAlO2', 'KAlO2', 'CaSiO3']
    cpx = {key:0. for key in keys}
    plg = {key:0. for key in keys}
    ol = {key:0. for key in keys}
    anorthite=components['CaAl2O4']/(components['CaAl2O4']+1.5*components['NaAlO2'])
    plg['CaAl2O4'] = 0.99*np.power(10.,((2446./T) - (1.122+0.2562*anorthite)))
    plg['NaAlO2'] = 0.92*np.power(10.,(((3195.+(3283.*anorthite))/T) - (2.318 + (1.885*anorthite))))
    plg['KAlO2'] = 0.15
    plg['MnO'] = 0.031
    plg['PO52'] = 0.1
    cpx['MgO'] = np.power(10.,((3798./T) - 2.28))
    cpx['FeO'] = 0.24*cpx['MgO']
    cpx['CaSiO3'] = np.power(10.,((1783./T)-0.753))
    cpx['CaAl2O4'] = np.power(10.,((2418./T)-2.3))
    cpx['NaAlO2'] = np.power(10.,((5087./T)-4.48))
    cpx['TiO2'] = np.power(10.,((1034./T)-1.27))
    cpx['KAlO2'] = 0.007
    cpx['MnO'] = 1.02
    cpx['PO52'] = 0.05
    ol['MgO'] = np.power(10.,((3386./T)-1.599))
    ol['FeO'] = np.power(10.,((4230./T)-2.707))
    ol['KAlO2'] = 0.001
    ol['MnO'] = 1.63
    ol['PO52'] = 0.2
    kd = {'cpx':cpx, 'ol':ol, 'plg':plg}
    return kd
    
def kdCalc_basalt4hiP(components, T, P):
    keys = ['MgO', 'FeO', 'TiO2', 'PO52', 'MnO', 'CaAl2O4', 'NaAlO2', 'KAlO2', 'CaSiO3']
    cpx = {key:0. for key in keys}
    plg = {key:0. for key in keys}
    ol = {key:0. for key in keys}
    anorthite=components['CaAl2O4']/(components['CaAl2O4']+1.5*components['NaAlO2'])
    plg['CaAl2O4'] = np.power(10.,((2446./T) - (1.122+0.2562*anorthite)))
    plg['NaAlO2'] = np.power(10.,((((0.0345*P) + (3283.*anorthite) + 3195.)/T) - (1.885*anorthite) - 2.349))
    plg['KAlO2'] = 0.15
    plg['MnO'] = 0.031
    plg['PO52'] = 0.1
    cpx['MgO'] = np.power(10.,((3798./(T-(0.008*P))) - 2.25))
    cpx['FeO'] = (0.24 + (0.000007*P))*cpx['MgO']
    cpx['CaSiO3'] = np.power(10.,((2357./(T-(0.01*P)))-1.22))
    cpx['CaAl2O4'] = np.power(10.,(14380./(T - (0.013*P))) - 10.75)
    cpx['NaAlO2'] = np.power(10.,((5087./(T-(0.016*P)))-4.48))
    cpx['TiO2'] = 0.3
    cpx['KAlO2'] = 0.007
    cpx['MnO'] = 1.02
    cpx['PO52'] = 0.05
    ol['MgO'] = np.power(10., (2715./(T - (0.004*P))) - 1.158)
    ol['FeO'] = ol['MgO']*0.3
    ol['KAlO2'] = 0.001
    ol['MnO'] = 1.63
    ol['PO52'] = 0.2
    kd = {'cpx':cpx, 'ol':ol, 'plg':plg}
    return kd
    
def kdCalc_langmuir1992(components, T, P):
    """
    This uses Langmuir Et Al. 1992 Calculations.
    P in bars, T in Kelvin.
    """
    keys = ['MgO', 'FeO', 'TiO2', 'PO52', 'MnO', 'CaAl2O4', 'NaAlO2', 'KAlO2', 'CaSiO3']
    cpx = {key:0. for key in keys}
    plg = {key:0. for key in keys}
    ol = {key:0. for key in keys}
    anorthite=components['CaAl2O4']/(components['CaAl2O4']+1.5*components['NaAlO2'])
    plg['CaAl2O4'] = np.power(10.,(2446./T) - (1.122  + 0.2562*anorthite))
    plg['NaAlO2'] = np.power(10.,((3195. + (3283.*anorthite) + (0.0506*P))/T) - (1.885*anorthite) -2.3715)
    plg['KAlO2'] = 0.15
    plg['MnO'] = 0.031
    plg['PO52'] = 0.1
    cpx['MgO'] = np.power(10.,(((3798. + (0.021*P))/T) - 2.28))
    cpx['FeO'] = cpx['MgO']*np.power(10.,-0.6198)
    cpx['CaSiO3'] = np.power(10.,((1783. + (0.0038*P))/T) -0.753)
    cpx['CaAl2O4'] = np.power(10.,((2418. + (0.068*P))/T) -2.3)
    cpx['NaAlO2'] = np.power(10.,((5087. + (0.073*P))/T) - 4.48)
    cpx['TiO2'] = np.power(10.,((1034. + (0.053*P))/T) - 1.27)
    cpx['KAlO2'] = 0.007
    cpx['MnO'] = 1.02
    cpx['PO52'] = 0.05
    ol['MgO'] = np.exp((6921./T) + (3.4*components['NaAlO2']/2.) + 
                (6.3*components['KAlO2']/2.) + (0.00001154*P) - 3.27)
    ol['FeO'] = ol['MgO']*(np.power(10., -0.523))
    ol['KAlO2'] = 0.001
    ol['MnO'] = 1.63
    ol['PO52'] = 0.2
    kd = {'cpx':cpx, 'ol':ol, 'plg':plg}
    return kd
    
    
def kdCalc_hbasalt(components, T, P):
    """
    This uses Langmuir Et Al. 1992 Calculations.
    P in bars, T in Kelvin.
    """
    H = 0.
    keys = ['MgO', 'FeO', 'TiO2', 'PO52', 'MnO', 'CaAl2O4', 'NaAlO2', 'KAlO2', 'CaSiO3']
    cpx = {key:0. for key in keys}
    plg = {key:0. for key in keys}
    ol = {key:0. for key in keys}
    anorthite=components['CaAl2O4']/(components['CaAl2O4']+1.5*components['NaAlO2'])
    plg['CaAl2O4'] = np.power(10.,(2405./T) - (1.12  + 0.2562*anorthite)+0.000005771*P - 0.00712275*H)
    plg['NaAlO2'] = np.power(10.,((3195. + (3283.*anorthite) + (0.0254*P) -(30.*H))/T) - (1.885*anorthite) -2.35)
    plg['KAlO2'] = 0.15
    plg['MnO'] = 0.031
    plg['PO52'] = 0.1
    cpx['MgO'] = np.power(10.,(((3798. + (0.017*P))/T) - 2.28 - (0.0036*H)))
    cpx['FeO'] = cpx['MgO']*0.24
    cpx['CaSiO3'] = np.power(10.,((1783. + (0.0038*P))/T) -0.753)
    cpx['CaAl2O4'] = np.power(10.,((2432. + (0.062*P))/T) -2.3 - (0.000095*H*components['KAlO2']/2.))
    cpx['NaAlO2'] = np.power(10.,((4829. + (0.079*P))/T) - 4.48-(0.0045*H))
    cpx['TiO2'] = np.power(10.,((1034. + (0.0531*P))/T) - 1.27)
    cpx['KAlO2'] = 0.007
    cpx['MnO'] = 1.02
    cpx['PO52'] = 0.05
    ol['MgO'] = np.exp((6277./T) + (3.3*components['NaAlO2']/2.) + 
                (6.3*components['KAlO2']/2.) + (0.0000134*P) - 2.8474 + 
                (0.000095*H*components['KAlO2']/2.) + (0.000538*components['KAlO2']/2.))
    ol['FeO'] = ol['MgO']*0.29
    ol['KAlO2'] = 0.001
    ol['MnO'] = 1.63
    ol['PO52'] = 0.2
    kd = {'cpx':cpx, 'ol':ol, 'plg':plg}
    return kd
    
