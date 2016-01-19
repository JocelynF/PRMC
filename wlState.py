import numpy as np
from wl1989stoich import *
from wl1989kdcalc import *

fa_guess = {'plg':0., 'ol':0., 'cpx':0.}
kdCalc = kdCalc_langmuir1992
    
def state(system_components, T, P, uaj, ta, fa = fa_guess):
    tolerance = np.power(10,-5.)
    max_iter = 3000
    qa = {'plg':0., 'ol':0., 'cpx':0.}
    dfa = {}
    # Variable a tells it whether or not to break the loop
    a = True
    # Liquid composition is initally the same as system since all fa's equal 0
    liquid_components = {key:system_components[key] for key in system_components.keys()}
    # Calculate new solid fractions in equilibrium with the current state
    phase_list = []
    kdaj = kdCalc(liquid_components, T)
    for p in fa.keys():
        #Calculate the initial saturation for each phase
        qa[p] = calculate_Qa(fa, kdaj, p, system_components, ta, uaj)
        if (qa[p]>0) or (fa[p]>0):
            # Add phase to list if it is oversaturated or it is undersaturated
            # and fa is greater than 0
            # DOUBLE CHECK IF I SHOULD USE STOICH
            phase_list.append(p)
    # If there is no phase that is saturated, no need to enter the loop
    if len(phase_list) == 0:
        a = False
    i = 0
    while (a == True) and (i<max_iter):
        i += 1
        # Use Newton Method to find new Fa if there are phases present  
        if len(phase_list) !=0:
#            fa_new = newton(qa, fa, kdaj, liquid_components, uaj)
            # THE FOLLOWING CODE USES THE MATRIX METHOD
            pab_dict = create_Pab_dict(fa, kdaj, system_components, uaj, phase_list)
            dfa = solve_matrix(pab_dict, qa, phase_list)
            if dfa == 'Singular':
                print 'Singular'
            fa_new = {}
            tst = 0.
            for phase in phase_list:
                # Check to make sure the new Fa is greater than 0 and less than 1
                fa_new[phase] = fa[phase] + dfa[phase]
                if fa_new[phase]<0:
                    fa_new[phase] = 0.1*fa[phase]
                elif fa_new[phase]>1:
                    fa_new[phase]= 0.9 + .1*fa[phase]
                tst += abs(fa[phase] - fa_new[phase])
                fa[phase] = fa_new[phase]
            # Recalculate Liquid Percent
            x = tst/len(phase_list)
            if x <= tolerance:
                a = False
            for component in liquid_components.keys():
                if component == 'SiO2':
                    pass
                else:                
                    liquid_components[component] = calculate_liquidComp(fa, kdaj, component, system_components)
            # Recalculate Saturation
            phase_list = []
            kdaj = kdCalc(liquid_components, T)
            for p in fa.keys():
                qa[p] = calculate_Qa(fa, kdaj, p, system_components, ta, uaj)
                if (qa[p]>0) or (fa[p]>0):
                    phase_list.append(p)

            if sum(fa.values())>=1:
                a = True
            elif sum(fa.values())<0:
                a = True
        else:
            a = False
    return qa, fa,liquid_components, i
            
            
# The following functions are used in State
            
def calculate_Rj(fa, kd, component):
    # Called by calculate_Pab and calculate_Qa
    temp = 0.
    for p in fa.keys():
        temp += fa[p]*(kd[p][component]-1)
    rj = 1./(1.+temp)
    return rj

def calculate_liquidComp(fa, kd, component, system_components):
    rj = calculate_Rj(fa, kd, component)
    clj = system_components[component]*rj
    return clj
    
def calculate_Qa(fa, kd, phase, system_components, ta, uaj):
    # Called by State
    # Initialize rj, clj, and caj
    caj = {'plg':{key:0 for key in system_components.keys()}, 
    'cpx':{key:0 for key in system_components.keys()}, 
    'ol':{key:0 for key in system_components.keys()}}
    # Given composition in component form and the Temperature,
    # calculate the saturation of a given phase.
    qa = -ta[phase]
    for component in system_components.keys():
        if component == 'SiO2':
            pass
        else:
        # Calculate the values for Rj
        # Calculate the liquid composition - doesn't need to be saved
        # Calculate the composition of the phases
            clj = calculate_liquidComp(fa,kd, component, system_components)
            caj[phase][component] = clj*kd[phase][component]
            qa += uaj[phase][component]*caj[phase][component]
    return qa

def calculate_Pab(fa, kd, phase1, phase2, system_components,uaj):
    # Called by create_pab_dict
    rj = {}
    pab = 0.
    for component in system_components.keys():
        if component == 'SiO2':
            pass
        else:
            rj = calculate_Rj(fa,kd,component)
            pab += uaj[phase1][component]*system_components[component]*kd[phase1][component]*(kd[phase2][component]-1)*rj*rj
    return pab

def create_Pab_dict(fa, kd, system_components, uaj, phase_list):
    # Called by State
    pab = {'plg':{}, 'cpx':{}, 'ol':{}}
    for phase1 in phase_list:
        for phase2 in phase_list:
            pab[phase1][phase2] = calculate_Pab(fa, kd, phase1, phase2, system_components,uaj)
    return pab
    
def solve_matrix(pab, qa, phase_list):
    # NEED TO ADD PART IN CASE OF SINGULAR MATRIX
    # Called by State
    # Convert dictionaries to arrays then use numpy to solve
    pab_array = np.zeros([len(phase_list),len(phase_list)])
    qa_array = np.zeros([len(phase_list),1])
    dfa = {}
    for i in xrange(0,len(phase_list)):
        qa_array[i] = qa[phase_list[i]]
        for j in xrange(0,len(phase_list)):
            pab_array[i,j] = pab[phase_list[i]][phase_list[j]]
    det = np.linalg.det(pab_array)
    if det == 0:
         dfa = 'Singular'
    # Solve Matrix
    else:
        pab_array = pab_array
        dfa_array = np.dot(np.linalg.inv(pab_array),(qa_array))
        # Convert back to Dicitonaries
        for k in xrange(0,len(phase_list)):
            dfa[phase_list[k]] = dfa_array[k][0]
    return dfa       
    
def newton(qa, fa, kd, system_components, uaj):
    # Don't want to divide by a number smaller than epsilon
    epsilon = np.power(10, -14)
    fa_new = {}
    for phase1 in fa.keys():
        qa_prime = 0
        for component in system_components.keys():
            if component == 'SiO2':
                pass
            else:
                kb = 0
                for phase2 in fa.keys():
                    kb+=kd[phase2][component]-1
                    rj = calculate_Rj(fa, kd, component)
                    qa_prime -= kd[phase1][component]*system_components[component]*rj*rj*kb
        #if np.abs(qa_prime) < epsilon:
            #return "Qa_prime Too Small"
        #else:
        fa_new[phase1] = fa[phase1]-(qa[phase1]/qa_prime)
    return fa_new       