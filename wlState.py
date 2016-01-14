import numpy as np
from wl1989stoich import *
from wl1989kdcalc import *
  
def state(system_components, T, uaj, ta):
    tolerance = np.power(10,-5.)
    max_iter = 3000
    # Initiate a first guess of Fa
    fa = pd.Series(0., index = ['ol', 'plg', 'cpx'])
    qa = pd.Series(0., index = ['ol', 'plg', 'cpx'])
    dfa = pd.Series(0., index = ['ol', 'plg', 'cpx'])
    # Variable a tells it whether or not to break the loop
    a = True
    # Liquid composition is initally the same as system since all fa's equal 0
    liquid_components = system_components.copy()
    # Calculate new solid fractions in equilibrium with the current state
    phase_list = []
    kdaj = kdCalc(liquid_components, T)
    for p in fa.index:
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
            pab = create_pab_matrix(fa, kdaj, system_components, uaj, phase_list)
            dfa = solve_matrix(pab, qa, phase_list)
            if dfa == 'Singular':
                print('Singular')
            fa_new = fa + dfa
            tst = 0.
		if ((fa_new > 0)&(fa_new < 1)).all() == False:
			for phase in phase_list:
			    # Check to make sure the new Fa is greater than 0 and less than 1
			    if fa_new[phase]<0:
				  fa_new[phase] = 0.1*fa[phase]
			    elif fa_new[phase]>1:
				  fa_new[phase]= 0.9 + .1*fa[phase]
		tst = np.sum(abs(fa-fa_new))
		fa = fa_new.copy()
            # Recalculate Liquid Percent
            x = tst/len(phase_list)
            if x <= tolerance:
                a = False
            liquid_components = calculate_liquidComp(fa, kdaj, component, system_components)
            # Recalculate Saturation
            phase_list = []
            kdaj = kdCalc(liquid_components, T)
            for p in fa.keys():
                qa[p] = calculate_Qa(fa, kdaj, p, system_components, ta, uaj)
                if (qa[p]>0) or (fa[p]>0):
                    phase_list.append(p)

            qa_new = [qa[x] for x in phase_list]
            #if all(np.abs(value) <= tolerance for value in qa_new):
                #a = False
            if sum(fa.values())>=1:
                a = True
            elif sum(fa.values())<0:
                a = True
        else:
            a = False
    return qa, fa,liquid_components, i
            
            
# The following functions are used in State
            
def calculate_Rj(fa, kd):
    # Called by calculate_Pab and calculate_Qa
    rj = pd.Series(0., index = kd.index) 
    #kd index does not include SiO2
    fadot = fa.dot
    for component in kd.index:
        temp = fadot(kd.loc[component,:]-1.)
        rj[component] = 1./(1.+temp)
    return rj

def calculate_liquidComp(fa, kd, system_components):
    rj = calculate_Rj(fa, kd)
    clj = system_components.multiply(rj)
    return clj
   
    
def calculate_Qa(fa, kd, system_components, ta, uaj):
    # Called by State
    # Initializecaj
    caj = pd.DataFrame(0,index = kd.index, columns = kd.columns)
    qa = pd.Series(0, index = kd.columns)
    #kd index does not include SiO2
    # Given composition in component form and the Temperature,
    # calculate the saturation of a given phase.
    # Calculate the values for Rj
    # Calculate the liquid composition - doesn't need to be saved
    # Calculate the composition of the phases
    clj = calculate_liquidComp(fa,kd, system_components)
    #The following line are just for speedin python. 
    #It avoids the use of dots in a loop.
    cljmultiply = clj.multiply
    for phase in kd.columns: 
        caj.loc[:,phase] = cljmultiply(kd.loc[:, phase]
        qa[phase] = -ta[phase] + caj.loc[:,phase].dot(uaj.loc[:,phase])
    return qa

def calculate_Pab(fa, kd, phase1, phase2, system_components,uaj):
    # Called by create_pab_dict
    for component in system_components.keys():
        rj = calculate_Rj(fa,kd)
	  clj = calculate_liquidComp(fa, kd, system_components)
	  cljrj = clj.multiply(rj)
	  caj = clj.multiply(kd.loc[:, phase1])
	  cljrjcaj = cljrj.multiply(caj)
	  cljrjcajuaj = cljrjcaj.multiply(uaj.loc[:, phase1]
        pab = np.sum(cljrjcajuaj.multiply(kd.loc[:,phase2]-1))
    return pab

def create_Pab_matrix(fa, kd, system_components, uaj, phase_list):
    # Called by State
    pab = pd.DataFrame(columns = kd.columns, index = kd.columns)
    for phase1 in phase_list:
        for phase2 in phase_list:
            pab.loc[phase1, phase2] = calculate_Pab(fa, kd, phase1, phase2, system_components,uaj)
    return pab
    
def solve_matrix(pab, qa, phase_list):
    # NEED TO ADD PART IN CASE OF SINGULAR MATRIX
    # Called by State
    # Convert dictionaries to arrays then use numpy to solve
    pab_array = pab.as_matrix()
    qa_array = qa.as_matrix()
    det = np.linalg.det(pab_array)
    if det == 0:
         dfa = 'Singular'
    # Solve Matrix
    else:
        dfa_array = np.dot(np.linalg.inv(pab_array),(qa_array))
        # Convert back to Dicitonaries
        for k in xrange(len(pab.index)):
            dfa[pab_index[k]] = dfa_array[k][0]
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