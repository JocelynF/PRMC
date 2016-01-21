from wl1989stoich import *
from wl1989kdcalc import *
from wlState import *
from mixing import *
import numpy as np
import pandas as pd

kd_dict = {'cpx':{'La':0.07,   'Sr':0.1,    'Ba':0.0005,   'Sc':2,    'Th':0.0013,  'Ni':4.4,  'K2O': 0.007, 'V':0.6, 'P2O5':0.05, 'MnO':0.85},
          'plg':{'La':0.08,    'Sr':1.5,    'Ba':0.3,      'Sc':0.08, 'Th':0.13,    'Ni':0.06, 'K2O':0.15, 'V':0.035, 'P2O5':0.1, 'MnO':0.06},
          'ol':{'La':0.000001, 'Sr':0.0001, 'Ba':0.000002, 'Sc':0.3,  'Th':0.00001, 'Ni':10., 'K2O':0.001, 'V':0.09, 'P2O5':0.2, 'MnO':1.0}}


ta = {'cpx':1., 'plg':1., 'ol':2./3.}
uaj_plg = {'CaAl2O4':5./3., 'NaAlO2':2.5, 'MgO':0., 'FeO':0., 'CaSiO3':0., 'TiO2':0., 'KAlO2':0., 'PO52':0., 'MnO':0.}
uaj_ol = {'CaAl2O4':0., 'NaAlO2':0., 'MgO':1, 'FeO':1, 'CaSiO3':0., 'TiO2':0., 'KAlO2':0., 'PO52':0., 'MnO':0.}
uaj_cpx = {'CaAl2O4':4./3., 'NaAlO2':2., 'MgO':2., 'FeO':2., 'CaSiO3':1., 'TiO2':1., 'KAlO2':0., 'PO52':0., 'MnO':0.}
uaj = {'ol':uaj_ol, 'plg':uaj_plg, 'cpx':uaj_cpx}


def get_first_T(system_components, P = 1., kdCalc = kdCalc_original):
    firstT = 2000.
    deltaT = 100.
    qa, fa, major_liquid_components, num_iter = state(system_components,firstT,uaj, ta, P=P, kdCalc= kdCalc)
    fl = 1-sum(fa.values())
    if num_iter == 30000:
        print 'MAX ITERATION!'
    while (fl == 1.) or (deltaT > 1.):
        if fl == 1.:
            firstT = firstT-deltaT
        elif (fl < 1.) and (deltaT > 1.):
            firstT = firstT+deltaT
            deltaT = deltaT/10.
            firstT=firstT-deltaT
        qa, fa, major_liquid_components, num_iter = state(system_components,firstT,uaj, ta, P=P, kdCalc= kdCalc)
        fl = 1-sum(fa.values())
        if num_iter == 30000:
            print 'MAX ITERATION!'
            break
        print fl
        print deltaT
        print firstT
    return firstT


#Equilibrium Crystallization - System Components never Change 
def eq_model_trange(t_start, t_stop, major_start_comp, trace_start_comp, P = 1., kdCalc = kdCalc_original):
    tstep = 1.
    bulk_d = {key:0. for key in trace_start_comp}
    trange = np.arange(t_stop,t_start, tstep)
    system_components = oxideToComponent(major_start_comp)
    major_oxide_dict = {key:[] for key in major_start_comp}
    trace_dict = {key:[] for key in trace_start_comp}
    fl = []
    for i in xrange(len(trange)):
        # Major Elements
        #if i == 0:
        qa, fa, major_liquid_components, num_iter = state(system_components,trange[-i-1],uaj, ta, P = P, kdCalc = kdCalc)
        #else:
            #major_liquid_components = oxideToComponent(major_oxides)
            #qa, fa, major_liquid_components, num_iter = state(system_components,trange[-i-1],uaj, ta, P = P, fa = fa, liquid_components = major_liquid_components, kdCalc = kdCalc)
        major_oxides = cationFracToWeight(major_liquid_components)
        liq = 1.-sum(fa.values())
        fl.append(liq)
        # MnO, P2O5, K2O
        for elem in ['K2O', 'P2O5', 'MnO']:
            #Calculate Bulk D
            bulk_d[elem] = 0.
            fa_tot = sum(fa.values())
            if fa_tot != 0.:
                for phase in fa:
                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
            #Add erupted composition to eruption dictionary
            major_oxides[elem] = major_oxides[elem]/(liq +(1.-liq)*bulk_d[elem])
        for key in major_oxides:
            major_oxide_dict[key].append(major_oxides[key])
        #Trace Elements
        for elem in trace_start_comp:
            #Calculate Bulk D
            bulk_d[elem] = 0.
            fa_tot = sum(fa.values())
            if fa_tot != 0.:
                for phase in fa:
                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
            #Add erupted composition to eruption dictionary
            trace_dict[elem].append(trace_start_comp[elem]/(liq +(1.-liq)*bulk_d[elem]))
    return fl, major_oxide_dict, trace_dict


#Fractional Crystallization - System Components change after each iteration
def frac_model_trange(t_start, t_stop, major_start_comp, trace_start_comp, P=1., kdCalc = kdCalc_original):
    tstep = 1.
    bulk_d = {key:0. for key in trace_start_comp}
    trange = np.arange(t_stop,t_start, tstep)
    system_components = oxideToComponent(major_start_comp)
    major_liquid_components = system_components.copy()
    trace_liquid_comp = trace_start_comp.copy()
    major_oxide_dict = {key:[] for key in major_start_comp}
    trace_dict = {key:[] for key in trace_start_comp}
    fl = []
    for i in xrange(len(trange)):
        # Major Elements
        if i == 0:
            qa, fa, major_liquid_components, num_iter = state(major_liquid_components,trange[-i-1],uaj, ta, P=P, kdCalc = kdCalc)
        else:
            major_liquid_components = oxideToComponent(major_oxides)
            qa, fa, major_liquid_components, num_iter = state(major_liquid_components,trange[-i-1],uaj, ta, P=P,fa = fa, liquid_components = major_liquid_components, kdCalc = kdCalc)
        major_oxides = cationFracToWeight(major_liquid_components)
        liq = (1. - sum(fa.values()))
        if i == 0:
            fl.append(liq)
        else:
            fl.append(liq*fl[-1])
        # MnO, P2O5, K2O
        for elem in ['K2O', 'P2O5', 'MnO']:
            #Calculate Bulk D
            bulk_d[elem] = 0.
            fa_tot = sum(fa.values())
            if fa_tot != 0.:
                for phase in fa:
                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
            #Add erupted composition to eruption dictionary
            major_oxides[elem] = major_oxides[elem]/(liq +(1.-liq)*bulk_d[elem])
        for key in major_oxides:
            major_oxide_dict[key].append(major_oxides[key])
        #Trace Elements
        for elem in trace_start_comp:
            #Calculate Bulk D
            bulk_d[elem] = 0.
            fa_tot = sum(fa.values())
            if fa_tot != 0.:
                for phase in fa:
                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
            #Add erupted composition to eruption dictionary
            trace_liquid_comp[elem] = trace_liquid_comp[elem]/(liq +(1.-liq)*bulk_d[elem])
            trace_dict[elem].append(trace_liquid_comp[elem])
    return fl, major_oxide_dict, trace_dict


##Equilibrium Crystallization - System Components never Change 
def eq_model_fstop(f_stop, major_start_comp, trace_start_comp, P = 1., kdCalc = kdCalc_original):
    tstep = 1.
    bulk_d = {key:0. for key in trace_start_comp}
    system_components = oxideToComponent(major_start_comp)
    t = get_first_T(system_components)
    major_oxide_dict = {key:[] for key in major_start_comp}
    trace_dict = {key:[] for key in trace_start_comp}
    fl = []
    qa, fa,major_liquid_components, num_iter = state(system_components,t,uaj, ta, P = P, kdCalc = kdCalc)
    major_oxides = cationFracToWeight(major_liquid_components)
    liq = 1.-sum(fa.values())
    fl.append(liq)
    # MnO, P2O5, K2O
    for elem in ['K2O', 'P2O5', 'MnO']:
        #Calculate Bulk D
        bulk_d[elem] = 0.
        fa_tot = sum(fa.values())
        if fa_tot != 0.:
            for phase in fa:
                bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
        #Add erupted composition to eruption dictionary
        major_oxides[elem] = major_oxides[elem]/(liq +(1.-liq)*bulk_d[elem])
    for key in major_oxides:
        major_oxide_dict[key].append(major_oxides[key])
    #Trace Elements
    for elem in trace_start_comp:
        #Calculate Bulk D
        bulk_d[elem] = 0.
        fa_tot = sum(fa.values())
        if fa_tot != 0.:
            for key in fa:
                bulk_d[elem] += (fa[key]/fa_tot)*kd_dict[key][elem]
        #Add erupted composition to eruption dictionary
        trace_dict[elem].append(trace_start_comp[elem]/(liq +(1.-liq)*bulk_d[elem]))
    while fl[-1]>f_stop:
        t = t - tstep
        qa, fa,major_liquid_components, num_iter = state(system_components,t,uaj, ta, P = P, kdCalc = kdCalc)
        major_oxides = cationFracToWeight(major_liquid_components)
        liq = 1.-sum(fa.values())
        fl.append(liq)
        # MnO, P2O5, K2O
        for elem in ['K2O', 'P2O5', 'MnO']:
            #Calculate Bulk D
            bulk_d[elem] = 0.
            fa_tot = sum(fa.values())
            if fa_tot != 0.:
                for phase in fa:
                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
            #Add erupted composition to eruption dictionary
            major_oxides[elem] = major_oxides[elem]/(liq +(1.-liq)*bulk_d[elem])
        for key in major_oxides:
            major_oxide_dict[key].append(major_oxides[key])
        #Trace Elements
        for elem in trace_start_comp:
            #Calculate Bulk D
            bulk_d[elem] = 0.
            fa_tot = sum(fa.values())
            if fa_tot != 0.:
                for key in fa:
                    bulk_d[elem] += (fa[key]/fa_tot)*kd_dict[key][elem]
            #Add erupted composition to eruption dictionary
            trace_dict[elem].append(trace_start_comp[elem]/(liq +(1.-liq)*bulk_d[elem]))
    return fl, major_oxide_dict, trace_dict


##Fractional Crystallization - System Components change after each iteration
def frac_model_fstop(f_stop, major_start_comp, trace_start_comp, P = 1., kdCalc = kdCalc_original): 
    tstep = 1.
    bulk_d = {key:0. for key in trace_start_comp}
    system_components = oxideToComponent(major_start_comp)
    t = get_first_T(system_components)
    major_liquid_components = system_components.copy()
    trace_liquid_comp = trace_start_comp.copy()
    major_oxide_dict = {key:[] for key in major_start_comp}
    trace_dict = {key:[] for key in trace_start_comp}
    fl = []
    qa, fa, major_liquid_components, num_iter = state(major_liquid_components,t,uaj, ta, P = P, kdCalc = kdCalc)
    major_oxides = cationFracToWeight(major_liquid_components)
    liq = (1. - sum(fa.values()))
    fl.append(liq)
    # MnO, P2O5, K2O
    for elem in ['K2O', 'P2O5', 'MnO']:
        #Calculate Bulk D
        bulk_d[elem] = 0.
        fa_tot = sum(fa.values())
        if fa_tot != 0.:
            for phase in fa:
                bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
        #Add erupted composition to eruption dictionary
        major_oxides[elem] = major_oxides[elem]/(liq +(1.-liq)*bulk_d[elem])
    for key in major_oxides:
        major_oxide_dict[key].append(major_oxides[key])
    #Trace Elements
    for elem in trace_start_comp:
        #Calculate Bulk D
        bulk_d[elem] = 0.
        fa_tot = sum(fa.values())
        if fa_tot != 0.:
            for key in fa:
                bulk_d[elem] += (fa[key]/fa_tot)*kd_dict[key][elem]
        #Add erupted composition to eruption dictionary
        trace_liquid_comp[elem] = trace_liquid_comp[elem]/(liq +(1.-liq)*bulk_d[elem])
        trace_dict[elem].append(trace_liquid_comp[elem])
    while fl[-1]>f_stop:
        t = t-tstep
        major_liquid_components = oxideToComponent(major_oxides)
        qa, fa, major_liquid_components, num_iter = state(major_liquid_components,t,uaj, ta, P = P, kdCalc = kdCalc)
        major_oxides = cationFracToWeight(major_liquid_components)
        liq = (1. - sum(fa.values()))
        fl.append(liq*fl[-1])
        for elem in ['K2O', 'P2O5', 'MnO']:
            #Calculate Bulk D
            bulk_d[elem] = 0.
            fa_tot = sum(fa.values())
            if fa_tot != 0.:
                for phase in fa:
                    bulk_d[elem] += (fa[phase]/fa_tot)*kd_dict[phase][elem]
            #Add erupted composition to eruption dictionary
            major_oxides[elem] = major_oxides[elem]/(liq +(1.-liq)*bulk_d[elem])
        for key in major_oxides:
            major_oxide_dict[key].append(major_oxides[key])
        #Trace Elements
        for elem in trace_start_comp:
            #Calculate Bulk D
            bulk_d[elem] = 0
            fa_tot = sum(fa.values())
            if fa_tot != 0:
                for key in fa:
                    bulk_d[elem] += (fa[key]/fa_tot)*kd_dict[key][elem]
            #Add erupted composition to eruption dictionary
            trace_liquid_comp[elem] = trace_liquid_comp[elem]/(liq +(1.-liq)*bulk_d[elem])
            trace_dict[elem].append(trace_liquid_comp[elem])
    return fl, major_oxide_dict, trace_dict


   
#def insitu_model(f_stop, major_start_comp, trace_start_comp, small_f=0.5):
#    sz = 0.01 #percent going to solidification zone
#    major_sz_comp = major_start_comp.copy()
#    trace_sz_comp = trace_start_comp.copy()
#    #Initialize liquid dictionaries
#    major_liquid_dict = {key:[] for key in major_start_comp.keys()}
#    trace_liquid_dict = {key:[] for key in trace_start_comp.keys()}
#    #Big F is the total Ml/Mo from Langmuir 1989
#    big_f = []
#    #Crystallize solidification zone to small_f
#    fl, major_oxide_dict, trace_dict = eq_model_fstop(small_f, major_sz_comp, trace_sz_comp)
#    #Final Oxide composition is the solidification zone liquid composition
#    for key in major_oxide_dict:
#        major_sz_comp[key] = major_oxide_dict[key][-1]
#    for key in trace_dict:
#        trace_sz_comp[key] = trace_dict[key][-1]
#    #Calculate percent liquid in solidification zone
#    liq = fl[-1]*sz
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
#    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    