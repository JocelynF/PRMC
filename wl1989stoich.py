import numpy as np
import pandas as pd

<<<<<<< Updated upstream
mass = {'Si':28.0855, 'Ti':47.867, 'Al':26.9815, 'Fe':55.845, 'Mg':24.305, 
        'Ca':40.078, 'Na':22.98977,'O':15.999,'K':39.0987, 'P':30.973, 'Mn':54.938}

def cationMolFracElement(oxides):
    #Oxides is one row of a dataframe --- Should probably convert to series
    elem = ['SiO2', 'TiO2', 'AlO15', 'CaO', 'NaO5', 'MgO', 'FeO']#, 'KO5', 'PO52', 'MnO']
    mol = pd.Series(index = elem)
    oxide_norm = pd.Series(index = oxides.index)
    cationMolFrac = pd.Series(index = elem)
    #Normalize oxides to 100wt%
    for element in oxides.index:
        oxide_norm[element] = (oxides[element])/np.sum(oxides.values)
=======
mass = {'Si':28.0855, 'Ti':47.867, 'Al':26.9815, 'Fe':55.845, 'Mg':24.305, 'Ca':40.078, 'Na':22.98977,'O':15.999,'K':39.0987}#, 'H':1.0081 , }

def oxideToMolFracElement(oxides):
    mol = {}
    oxide_norm = {}
    cationMolFrac = {}
    #Normalize oxides to 100wt%
    tot = sum(oxides.values())
    for element in oxides:
        oxide_norm[element] = (oxides[element])/tot
>>>>>>> Stashed changes
        #Get moles of each oxides
    mol['SiO2'] = oxide_norm['SiO2']/(mass['Si']+2*mass['O']) 
    mol['TiO2'] = oxide_norm['TiO2']/(mass['Ti']+2*mass['O']) 
    mol['AlO15'] = (oxide_norm['Al2O3']*2.)/(2.*mass['Al']+3.*mass['O']) 
    mol['CaO'] = oxide_norm['CaO']/(mass['Ca']+mass['O']) 
    mol['NaO5'] = (oxide_norm['Na2O']*2.)/(2.*mass['Na']+mass['O']) 
    mol['MgO'] = oxide_norm['MgO']/(mass['Mg']+mass['O']) 
<<<<<<< Updated upstream
    mol['FeO'] = oxide_norm['FeO']/(mass['Fe']+mass['O'])
#    mol['KO5'] = (oxide_norm['K2O']*2.)/(2.*mass['K'] + mass['O'])
#    mol['PO52']= (oxide_norm['P2O5']*2.)/(2.*mass['P'] + 5.*mass['O'])
#    mol['MnO']= oxide_norm['MnO']/(mass['Mn'] + mass['O'])
    #Calculate mole fraction of each element
    tot = np.sum(mol.values)
    for element in mol.index:
        cationMolFrac[element] = mol[element]/tot
    return cationMolFrac
      
def cationMolFracComponent(oxides):
    components = 'CaAl2O4 NaAlO2 MgO FeO CaSiO3 TiO2'.split()
    compCationMolFrac = pd.Series(index = components)
    cationMolFrac = cationMolFracElement(oxides)
    #Components: CaAl2O4, NaAlO2, KAlO2, MgO, FeO, CaSiO3, TiO2, H2O, SiO2
=======
    mol['FeO'] = oxide_norm['FeO']/(mass['Fe']+mass['O']) 
    mol['KO5'] = (oxide_norm['K2O']*2.)/(2.*mass['K']+mass['O'])
    mol['PO52']=(oxide_norm['P2O5']*2.)/(2.*mass['P']+5.*mass['O'])
    mol['MnO'] = oxide_norm['MnO']/(mass['Mn']+mass['O']) 
    #mol['HO5'] = oxide_norm['H2O']/(mass['H']+0.5*mass['O'])
    #Calculate mole fraction of each element
    tot1 = sum(mol.values())
    for element in mol.keys():
        cationMolFrac[element] = mol[element]/tot1
    return cationMolFrac
      
def molFractoComponent(cationMolFrac):
    compCationMolFrac = {}
    #Components: CaAl2O4, NaAlO2, MgO, FeO, CaSiO3, TiO2, SiO2, KO5, PO52, MnO
>>>>>>> Stashed changes
    compCationMolFrac['MgO'] = cationMolFrac['MgO']
    compCationMolFrac['FeO'] = cationMolFrac['FeO']
    compCationMolFrac['KO5'] = cationMolFrac['KO5']
    compCationMolFrac['PO52'] = cationMolFrac['PO52']
    compCationMolFrac['MnO'] = cationMolFrac['MnO']
    compCationMolFrac['NaAlO2'] = 2.*cationMolFrac['NaO5']
    al_CaAl2O4 = (cationMolFrac['AlO15'] - cationMolFrac['NaO5'] - cationMolFrac['KO5'])
    compCationMolFrac['CaAl2O4'] = 1.5*al_CaAl2O4
    compCationMolFrac['CaSiO3'] = 2.*(cationMolFrac['CaO'] - (compCationMolFrac['CaAl2O4']/3.))
    compCationMolFrac['TiO2'] = cationMolFrac['TiO2']
    #compCationMolFrac['SiO2'] = cationMolFrac['SiO2'] - (compCationMolFrac['CaSiO3']/2.)
    #compCationMolFrac['HO5'] = cationMolFrac['HO5']
    return compCationMolFrac
 
<<<<<<< Updated upstream

=======
def oxideToComponent(oxides):
    cationMolFrac = oxideToMolFracElement(oxides)
    components = molFractoComponent(cationMolFrac)
    return components
 
>>>>>>> Stashed changes
def cationFracToWeight(components):
    index = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MgO', 'CaO', 'Na2O']
    oxide = pd.Series(index = index)
    oxide_final = pd.Series(index = index)
    oxide['Na2O'] = (components['NaAlO2']/2.)*(((2.*mass['Na'])+mass['O'])/2.)
    oxide['TiO2'] = components['TiO2']*(mass['Ti'] + (2*mass['O']))
    oxide['Al2O3'] = ((components['CaAl2O4']*(2./3.)) +(components['NaAlO2']/2.))*(((2.*mass['Al'])+(3.*mass['O']))/2.)
    oxide['FeO'] = components['FeO']*(mass['Fe']+mass['O'])
    oxide['MgO'] = components['MgO']*(mass['Mg']+mass['O'])
    oxide['CaO'] = ((components['CaAl2O4']/3.) + (components['CaSiO3']/2.))*(mass['Ca']+mass['O'])
    oxide['K2O'] = components['KO5']*(2.*mass['K']+mass['O'])/2.
    oxide['P2O5'] = components['PO5']*(2.*mass['P']+5.*mass['O'])/2.
    oxide['MnO'] = components['MnO']*(mass['Mn']+mass['O'])
    oxide['SiO2'] = (components['SiO2']+(components['CaSiO3']/2.))*(mass['Si']+(2*mass['O']))
    tot = np.sum(oxide.values)/100.
    for element in oxide.index:
            oxide_final[element] = oxide[element]/tot
    return oxide_final

