from wl1989stoich import *
from wl1989kdcalc import *
from wlState import *
from mixing import *
import numpy as np
import pandas as pd

kd_dict = {'cpx':{'La':0.07,   'Sr':0.1,    'Ba':0.0005,   'Sc':2,    'Th':0.0013,  'Ni':4.4,  'K': 0.007, 'V':0.6, 'P':0.05},
          'plg':{'La':0.08,    'Sr':1.5,    'Ba':0.3,      'Sc':0.08, 'Th':0.13,    'Ni':0.06, 'K':0.15, 'V':0.035, 'P':0.1},
          'ol':{'La':0.000001, 'Sr':0.0001, 'Ba':0.000002, 'Sc':0.3,  'Th':0.00001, 'Ni':10., 'K':0.001, 'V':0.09, 'P':0.2}}


ta = pd.Series([1.,1., 2./3.],index = ['ol', 'plg', 'cpx'])
uaj_plg = {'CaAl2O4':5./3., 'NaAlO2':2.5, 'MgO':0, 'FeO':0, 'CaSiO3':0, 'TiO2':0}#, 'HO5':0,'KAlO2':0}
uaj_ol = {'CaAl2O4':0, 'NaAlO2':0, 'MgO':1, 'FeO':1, 'CaSiO3':0, 'TiO2':0}#, 'HO5':0,'KAlO2':0}
uaj_cpx = {'CaAl2O4':4./3., 'NaAlO2':2., 'MgO':2., 'FeO':2., 'CaSiO3':1., 'TiO2':1.}#, 'HO5':0,'KAlO2':0}
uaj = {'ol':uaj_ol, 'plg':uaj_plg, 'cpx':uaj_cpx}
uaj = pd.DataFrame(uaj)


def get_first_T(system_components):
    firstT = 2000.
    deltaT = 100.
    qa, fa, major_liquid_components, num_iter = state(system_components,firstT,uaj, ta)
    fl = 1-np.um(fa.values)
    while (fl == 1.) or (deltaT > 1.):
        if fl == 1.:
            firstT = firstT-deltaT
        elif (fl < 1.) and (deltaT > 1.):
            firstT = firstT+deltaT
            deltaT = deltaT/10.
            firstT=firstT-deltaT
        qa, fa, major_liquid_components, num_iter = state(system_components,firstT,uaj, ta)
        fl = 1-sum(fa.values)
    return firstT


#Equilibrium Crystallization - System Components never Change 
def eq_model_trange(t_start, t_stop, major_start_comp, trace_start_comp):
    columns = ['T', 'fl', 'cpx', 'ol', 'plg']
    columns.extend(list(major_start_comp.index))
    columns.extend(list(trace_start_comp.index))
    tstep = 1
    bulk_d = pd.Series(0.,index = trace_start_comp.index)
    trange = np.arange(t_stop,t_start, tstep)
    trange[::-1].sort()
    results = pd.DataFrame(columns = columns, index = trange)
    system_components = cationMolFracComponent(major_start_comp)
    for i in xrange(len(trange)):
        qa, fa, major_liquid_components, num_iter = state(system_components,trange[i],uaj, ta)
        major_liquid_components['SiO2'] = 1. - np.sum(major_liquid_components.values)
        major_oxides = cationFracToWeight(major_liquid_components)
        liq = 1.-np.sum(fa.values)
        results.loc[trange[i], 'fl'] = liq
        results.loc[trange[i], 'cpx'] = fa['cpx']
        results.loc[trange[i], 'ol'] = fa['ol']
        results.loc[trange[i], 'plg'] = fa['plg']
        for element in major_oxides.index:
            results.loc[trange[i], element]
        #Trace Elements
        for elem in trace_start_comp.index:
            #Calculate Bulk D
            fa_tot = np.sum(fa.values)
            if fa_tot != 0:
                for phase in fa.index:
                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
            #Add erupted composition to eruption dictionary
            trace_elem = trace_start_comp[elem]/(liq +(1.-liq)*bulk_d[elem])  
            results.loc[trange[i], elem] = trace_elem            
    return results


##Fractional Crystallization - System Components change after each iteration
#def frac_model_trange(t_start, t_stop, major_start_comp, trace_start_comp):
#    columns = ['T', 'fl', 'cpx', 'ol', 'plg']
#    columns.append(major_start_comp.columns)
#    tstep = 1.
#    bulk_d = pd.Series(0.,index = trace_start_comp.index)
#    trange = np.arange(t_stop,t_start, tstep)
#    trange[::-1].sort()
#    results = pd.DataFrame(columns = columns, index = trange)
#    system_components = cationMolFracComponent(major_start_comp)
#    major_liquid_components = system_components.copy()
#    trace_liquid_comp = trace_start_comp.copy()
#    for i in xrange(len(trange)):
#        qa, fa, major_liquid_components, num_iter = state(major_liquid_components,trange[i],uaj, ta)
#        major_liquid_components['SiO2'] = 1. - np.sum(major_liquid_components.values)
#        major_oxides = cationFracToWeight(major_liquid_components)
#        liq = (1. - sum(fa.values()))
#        if i == 0:
#            results.loc[trange[i], 'fl'] = liq
#        else:
#            results.loc[trange[i], 'fl'] = liq*results.loc[trange[i-1], 'fl']        
#        for element in major_oxides.index:
#            results.loc[trange[i], major_oxides[element]]
#        #Trace Elements
#        for elem in trace_start_comp.index:
#            #Calculate Bulk D
#            fa_tot = np.sum(fa.values)
#            if fa_tot != 0:
#                for phase in fa.index:
#                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
#            #Add erupted composition to eruption dictionary
#            trace_liquid_comp[elem] = trace_liquid_comp[elem]/(liq +(1.-liq)*bulk_d[elem])
#            results.loc[trange[i], elem] = trace_liquid_comp[elem] 
#    return results
#
##Equilibrium Crystallization - System Components never Change 
#def eq_model_fstop(f_stop, major_start_comp, trace_start_comp):
#    columns = ['T', 'fl', 'cpx', 'ol', 'plg']
#    columns.append(major_start_comp.columns)
#    columns.append(trace_start_comp.columns)
#    tstep = 1.
#    bulk_d = pd.Series(0.,index = trace_start_comp.index)
#    results = pd.DataFrame(columns = columns)
#    i = 0
#    system_components = cationMolFracComponent(major_start_comp)
#    t = get_first_T(system_components)
#    qa, fa,major_liquid_components, num_iter = state(system_components,t,uaj, ta)
#    major_liquid_components['SiO2'] = 1.-np.sum(major_liquid_components.values)
#    major_oxides = cationFracToWeight(major_liquid_components)
#    results.loc[i, 'T'] = t
#    results.loc[i, 'cpx'] = fa['cpx']
#    results.loc[i, 'ol'] = fa['ol']
#    results.loc[i, 'plg'] = fa['plg']
#    for element in major_oxides.index:
#        results.loc[i,element] = major_oxides[element]
#    liq = 1.-np.sum(fa.values)
#    results.loc[i, 'fl'] = liq
#    #Trace Elements
#    for elem in trace_start_comp.index:
#        #Calculate Bulk D
#        fa_tot = np.sum(fa.values)
#        if fa_tot != 0:
#            for phase in fa.index:
#                bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
#        #Add erupted composition to eruption dictionary
#        trace_elem = trace_start_comp[elem]/(liq +(1.-liq)*bulk_d[elem])
#        results.loc[i,elem] = trace_elem
#    while liq > f_stop:
#        i += 1
#        t = t - tstep
#        qa, fa,major_liquid_components, num_iter = state(system_components,t,uaj, ta)
#        major_liquid_components['SiO2'] = 1. - np.sum(major_liquid_components.values)
#        results.loc[i, 'T'] = t
#        results.loc[i, 'cpx'] = fa['cpx']
#        results.loc[i, 'ol'] = fa['ol']
#        results.loc[i, 'plg'] = fa['plg']
#        major_oxides = cationFracToWeight(major_liquid_components)
#        for element in major_oxides.index:
#            results.loc[i,element] = major_oxides[element]
#        liq = 1.-np.sum(fa.values)
#        results.loc[i, 'fl'] = liq
#        #Trace Elements
#        for elem in trace_start_comp.keys():
#            #Calculate Bulk D
#            fa_tot = np.sum(fa.values)
#            if fa_tot != 0:
#                for phase in fa.index:
#                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
#            #Add erupted composition to eruption dictionary
#        trace_elem = trace_start_comp[elem]/(liq +(1.-liq)*bulk_d[elem])
#        results.loc[i,elem] = trace_elem
#        return fl, major_oxide_dict, trace_dict
#
#
##Fractional Crystallization - System Components change after each iteration
#def frac_model_fstop(f_stop, major_start_comp, trace_start_comp): 
#    i = 0    
#    columns = ['T', 'fl', 'cpx', 'ol', 'plg']
#    columns.append(major_start_comp.columns)
#    columns.append(trace_start_comp.columns)
#    tstep = 1.
#    bulk_d = pd.Series(0.,index = trace_start_comp.index)
#    results = pd.DataFrame(columns = columns)
#    system_components = cationMolFracComponent(major_start_comp)
#    t = get_first_T(system_components)
#    major_liquid_components = system_components.copy()
#    trace_liquid_comp = trace_start_comp.copy()
#    qa, fa, major_liquid_components, num_iter = state(major_liquid_components,t,uaj, ta)
#    major_liquid_components['SiO2']= 1.-np.sum(major_liquid_components.values)
#    major_oxides = cationFracToWeight(major_liquid_components)
#    results.loc[i, 'T'] = t
#    results.loc[i, 'cpx'] = fa['cpx']
#    results.loc[i, 'ol'] = fa['ol']
#    results.loc[i, 'plg'] = fa['plg']
#    for element in major_oxides.index:
#        results.loc[i,element] = major_oxides[element]
#    liq = 1.-np.sum(fa.values)
#    results.loc[i, 'fl'] = liq
#    #Trace Elements    
#    for elem in trace_start_comp.index:
#        #Calculate Bulk D
#        fa_tot = np.sum(fa.values)
#        if fa_tot != 0:
#            for phase in fa.index:
#                bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
#        #Add erupted composition to eruption dictionary
#	  trace_elem = trace_start_comp[elem]/(liq +(1.-liq)*bulk_d[elem])
#	  results.loc[i,elem] = trace_elem
#    while results.loc[i,'fl'] > f_stop:
#        i += 1
#        t = t-tstep
#        qa, fa, major_liquid_components, num_iter = state(major_liquid_components,t,uaj, ta)
#        major_liquid_components['SiO2'] = 1. - np.sum(major_liquid_components.value)
#        major_oxides = cationFracToWeight(major_liquid_components)
#        for element in major_oxides.index:
#            results.loc[i,element] = major_oxides[element]
#        liq = (1. - np.sum(fa.values))
#        liq = liq*results.loc[i-1, 'fl']
#        results.loc[i, 'fl'] = liq
#        #Trace Elements
#        for elem in trace_start_comp.index:
#            #Calculate Bulk D
#            fa_tot = np.sum(fa.values)
#            if fa_tot != 0:
#                for phase in fa.index:
#                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
#            #Add erupted composition to eruption dictionary
#            trace_liquid_comp[elem] = trace_liquid_comp[elem]/(liq +(1.-liq)*bulk_d[elem])
#            results.loc[i,elem] = trace_liquid_comp[elem]
#    return results
#
#    
#def insitu_model(f_stop, major_start_comp, trace_start_comp, small_f=0.5):
#    sz = 0.01 #percent going to solidification zone
#    i = 0    
#    columns = ['T', 'ml/mo', 'cpx', 'ol', 'plg']
#    columns.append(major_start_comp.columns)
#    columns.append(trace_start_comp.columns)
#    tstep = 1.
#    bulk_d = pd.Series(0.,index = trace_start_comp.columns)
#    insitu_results = pd.DataFrame(columns = columns)
#    major_sz_comp = major_start_comp.copy()
#    trace_sz_comp = trace_start_comp.copy()
#    #Crystallize solidification zone to small_f
#    eq_results = eq_model_fstop(small_f, major_sz_comp, trace_sz_comp)
#    #Final Oxide composition is the solidification zone liquid composition
#    for element in major_start_comp.columns:
#        major_sz_comp[element] = eq_results.loc[-1,element]
#    for elem in trace_start_comp.columns:
#        trace_sz_comp[key] = eq_results[-1,elem]
#    #Calculate percent liquid in solidification zone
#    liq = eq_results.loc[-1, 'fl']*sz
#    #Add Ml/Mo to list
#    big_f.append(1.-sz+liq)
#    #Mix solidification zone liquid with starting composition
#    major_liquid_comp = magma_mixing(major_sz_comp, major_start_comp, liq/big_f[-1])
#    trace_liquid_comp = magma_mixing(trace_sz_comp, trace_start_comp, liq/big_f[-1])
#    #Add this composition to the final liquid dictionary
#    for key in major_liquid_comp:
#        major_liquid_dict[key].append(major_liquid_comp[key])
#    for key in trace_liquid_comp:
#        trace_liquid_dict[key].append(trace_liquid_comp[key])    
#    #repeat this process until big_f is equal to f_stop   
#    while big_f[-1] > f_stop:  
#        #Crystallize solidification zone to small_f
#        fl, major_oxide_dict, trace_dict = eq_model_fstop(small_f, major_liquid_comp, trace_liquid_comp)
#        #Final Oxide composition is the solidification zone liquid composition
#        for key in major_oxide_dict:
#            major_sz_comp[key] = major_oxide_dict[key][-1]
#        for key in trace_dict:
#            trace_sz_comp[key] = trace_dict[key][-1]
#        #Calculate percent liquid in solidification zone
#        liq = fl[-1]*sz
#        #Add Ml/Mo to list
#        big_f.append(big_f[-1]-sz+liq)
#        #Mix solidification zone liquid with starting composition
#        major_liquid_comp = magma_mixing(major_sz_comp, major_liquid_comp, liq/big_f[-1])
#        trace_liquid_comp = magma_mixing(trace_sz_comp, trace_liquid_comp, liq/big_f[-1])
#        #Add this composition to the final liquid dictionary
#        for key in major_liquid_comp:
#            major_liquid_dict[key].append(major_liquid_comp[key])
#        for key in trace_liquid_comp:
#            trace_liquid_dict[key].append(trace_liquid_comp[key])
#    return big_f, major_liquid_dict, trace_liquid_dict
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    