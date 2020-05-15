import numpy as np
from scipy.optimize import curve_fit 


def curve_fitting(x_array, y_array, restrict_MTO = False):
    ##Curve fit     
    fail = False
    if restrict_MTO == False:
        p0 = [8.798, 8.901, 0.640]
        para_bounds = ((8.0, 8.0, 0.0), (9.0, 9.5, 1.0))
        fit = mass_metal_fit
    else:
        p0 = [8.798, 0.640]
        para_bounds = ((8.0, 0.0), (9.0, 1.0))
        fit = mass_metal_fit_constMTO
        
    try:
        o11, o21 = curve_fit(fit, x_array, y_array, p0 = p0, bounds = para_bounds)
        print(o11)
    except ValueError:
        print('Failed curve fitting!')
        fail = True
        
    return o11, o21, fail
    
    
    
def mass_metal_fit(mass, a, b, g):
    '''
    MTO is NOT held constant.
    
    Andrews & Martini Mass-Metallicity Relation:
    8.798 - np.log10(1 + ((10**8.901)/(10**mass))**0.640)
    '''
    
    return a - np.log10(1 + ((10**b)/(10**mass))**g)   



def mass_metal_fit_constMTO(mass, a, g):
    '''
    MTO is held constant.
    
    Andrews & Martini Mass-Metallicity Relation:
    8.798 - np.log10(1 + ((10**8.901)/(10**mass))**0.640)
    '''
    
    return a - np.log10(1 + ((10**8.901)/(10**mass))**g) 
















