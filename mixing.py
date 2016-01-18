def magma_mixing(magma1, magma2, prop1):
    """
    Mixing of two compositions where
    the proportion of magma1 is prop1
    and the proportion of magma 2 is 
    (1-prop1).
    """
    new_comp = {}
    for component in magma1.keys():
        new_comp[component] = mixing(magma1[component], magma2[component], prop1)
    return new_comp
    
def mixing(comp1, comp2, prop1):
    new_comp = prop1*comp1 + (1-prop1)*comp2
    return new_comp
    