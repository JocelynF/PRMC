import numpy as np

mass = {'Si':28.0855, 'Ti':47.867, 'Al':26.9815, 'Fe':55.845, 'Mg':24.305, 
        'Ca':40.078, 'Na':22.98977,'O':15.999,'K':39.0987, 'P':30.973, 'Mn':54.938}

def oxideToMolFracElement(oxides):
    mol = {}
    oxide_norm = {}
    cationMolFrac = {}
    #Normalize oxides to 100wt%
    tot = sum(oxides.values())
    for element in oxides:
        oxide_norm[element] = (oxides[element])/tot
    #Get moles of each oxides
    mol['SiO2'] = oxide_norm['SiO2']/(mass['Si']+2*mass['O']) 
    mol['TiO2'] = oxide_norm['TiO2']/(mass['Ti']+2*mass['O']) 
    mol['AlO15'] = (oxide_norm['Al2O3']*2.)/(2.*mass['Al']+3.*mass['O']) 
    mol['CaO'] = oxide_norm['CaO']/(mass['Ca']+mass['O']) 
    mol['NaO5'] = (oxide_norm['Na2O']*2.)/(2.*mass['Na']+mass['O']) 
    mol['MgO'] = oxide_norm['MgO']/(mass['Mg']+mass['O']) 
    mol['FeO'] = oxide_norm['FeO']/(mass['Fe']+mass['O']) 
    mol['KO5'] = (oxide_norm['K2O']*2.)/(2.*mass['K']+mass['O'])
    mol['PO52']=(oxide_norm['P2O5']*2.)/(2.*mass['P']+5.*mass['O'])
    mol['MnO'] = oxide_norm['MnO']/(mass['Mn']+mass['O']) 
    #mol['HO5'] = oxide_norm['H2O']/(mass['H']+0.5*mass['O'])
    #Calculate mole fraction of each element
    tot1 = sum(mol.values())
    cationMolFrac = {element: mol[element]/tot1 for element in mol}
    return cationMolFrac
      
def molFractoComponent(cationMolFrac):
    compCationMolFrac = {}
    #Components: CaAl2O4, NaAlO2, MgO, FeO, CaSiO3, TiO2, SiO2, KO5, PO52, MnO
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
    compCationMolFrac['SiO2'] = cationMolFrac['SiO2'] - (compCationMolFrac['CaSiO3']/2.)
    #compCationMolFrac['HO5'] = cationMolFrac['HO5']
    return compCationMolFrac

def oxideToComponent(oxides):
    cationMolFrac = oxideToMolFracElement(oxides)
    components = molFractoComponent(cationMolFrac)
    return components
 
def cationFracToWeight(components):
    oxide = {}
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
    tot = sum(oxide.values())/100.
    oxide_dict = {element: oxide[element]/tot for element in oxide}
    return oxide_dict

