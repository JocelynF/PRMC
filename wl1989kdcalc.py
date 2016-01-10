from wl1989stoich import *
import numpy as np
import pandas as pd

def kdCalc_original(components, T):
    """
    This is the original Weaver and Langmuir
    calculation of Kd's. Should double check the                              
    cpx CaSiO3 calculation.
    """
    columns = ['ol', 'plg', 'cpx']
    index = ['CaAl2O4', 'NaAlO2', 'MgO', 'FeO', 'CaSiO3', 'CaAl2O4', 'TiO2']
    kd = pd.DataFrame(0,index = index, columns = columns)
    anorthite=components['CaAl2O4']/(components['CaAl2O4']+1.5*components['NaAlO2'])
    kd.loc['CaAl2O4', 'plg'] = np.power(10.,((2446./T) - (1.122+0.2562*anorthite)))
    kd.loc['NaAlO2', 'plg'] = np.power(10.,(((3195.+(3283.*anorthite))/T) - (2.318 + (1.885*anorthite))))
    kd.loc['MgO', 'cpx'] = np.power(10.,((3798./T) - 2.28))
    kd.loc['FeO', 'cpx'] = 0.24*kd.loc['MgO', 'cpx']
    kd.loc['CaSiO3', 'cpx'] = np.power(10.,((1783./T)-0.753))
    kd.loc['CaAl2O4', 'cpx'] = np.power(10.,((2418./T)-2.3))
    kd.loc['NaAlO2', 'cpx'] = np.power(10.,((5087./T)-4.48))
    kd.loc['TiO2', 'cpx'] = np.power(10.,((1034./T)-1.27))
    kd.loc['MgO', 'ol'] = np.power(10.,((2715./T)-1.158))
    kd.loc['FeO', 'ol'] = np.power(10.,((4230./T)-2.741))
    return kd

def kdCalc(components, T):
    """
    This uses Reynold's Thesis Calculations.
    """
    columns = ['ol', 'plg', 'cpx']
    index = ['CaAl2O4', 'NaAlO2', 'MgO', 'FeO', 'CaSiO3', 'CaAl2O4', 'TiO2']
    kd = pd.DataFrame(0,index = index, columns = columns)
    anorthite=components['CaAl2O4']/(components['CaAl2O4']+1.5*components['NaAlO2'])
    kd.loc['CaAl2O4', 'plg'] = 0.99*np.power(10.,((2446./T) - (1.122+0.2562*anorthite)))
    kd.loc['NaAlO2', 'plg'] = 0.92*np.power(10.,(((3195.+(3283.*anorthite))/T) - (2.318 + (1.885*anorthite))))
    kd.loc['MgO', 'cpx'] = np.power(10.,((3798./T) - 2.28))
    kd.loc['FeO', 'cpx'] = 0.24*kd.loc['MgO', 'cpx']
    kd.loc['CaSiO3', 'cpx'] = np.power(10.,((1783./T)-0.753))
    kd.loc['CaAl2O4', 'cpx'] = np.power(10.,((2418./T)-2.3))
    kd.loc['NaAlO2', 'cpx'] = np.power(10.,((5087./T)-4.48))
    kd.loc['TiO2', 'cpx'] = np.power(10.,((1034./T)-1.27))
    kd.loc['MgO', 'ol'] = np.power(10.,((3386./T)-1.599))
    kd.loc['FeO', 'ol'] = np.power(10.,((4230./T)-2.707))
    return kd