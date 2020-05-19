import numpy as np
from scipy.optimize import curve_fit 
from astropy.io import ascii as asc
from Metallicity_Stack_Commons.column_names import npz_filename_dict, filename_dict, temp_metal_names0, bin_names0


def curve_fitting(x_array, y_array, restrict_MTO = False):
    '''
    Purpose:
        This function uses SciPy's curve_fit function to produce a curve that fits the composite Metallicity
        vs Mass data points. It used values from Andrews & Martini (2013) as initial parameters.
           
    Params:
        x_array --> x values on plot (i.e. mass array).
        y_array --> y values on plot (i.e. metallicity array).
        restrict_MTO (OPTIONAL) --> if the mass turnover value should be held constant in the curve fit of
                                    Metallicity vs Mass, then restrict_MTO = True.
        
    Returns:
        o11 --> an array of the optimal values of the curve fit.
        o21 --> an array of the estimated covariance of o11.
        fail --> if there is a ValueError when running curve_fit, then fail = True. Otherwise, fail = False.
    '''
        
    fail = False
    if restrict_MTO == False:
        p0 = [8.798, 8.901, 0.640]
        para_bounds = ((8.0, 8.0, 0.0), (9.0, 9.5, 1.0))
    else:
        p0 = [8.798, 0.640]
        para_bounds = ((8.0, 0.0), (9.0, 1.0))
        
    try:
        o11, o21 = curve_fit(mass_metal_fit, x_array, y_array, p0 = p0, bounds = para_bounds)
        print(o11)
    except ValueError:
        print('Failed curve fitting!')
        fail = True
        
    return o11, o21, fail
    
    
    
def mass_metal_fit(mass, a, g, b = 8.901):
    '''
    Purpose:
        This function returns a curve calculated from best fit parameters from curve_fitting. The curve
        equation comes from Andrews & Martini (2013):
            8.798 - np.log10(1 + ((10**8.901)/(10**mass))**0.640)
           
    Params:
        mass --> array of log values of stellar masses.
        a --> log(O/H) asymptotic value.
        g --> gamma value.
        b (OPTIONAL) --> mass turnover value. For constant MTO, don't provide a value for b.
        
    Returns:
        o11 --> an array of the optimal values of the curve fit.
        o21 --> an array of the estimated covariance of o11.
        fail --> if there is a ValueError when running curve_fit, then fail = True. Otherwise, fail = False.
    '''
    
    return a - np.log10(1 + ((10**b)/(10**mass))**g)



def extract_error_bars(fitspath):
    '''
    Purpose:
        This function reads in the error bar arrays and reorders the low and high error arrays to plot.
           
    Params:
        fitspath --> string of the file path where data files are located.
        
    Returns:
        err_dict --> a dictionary containing temperature and metallicity low and high errors.
    '''
    
    der_prop_err = np.load(fitspath + npz_filename_dict['der_prop_errors'])
    der_prop_file = asc.read(fitspath + filename_dict['bin_derived_prop_rev'])
    valid_tbl = asc.read(fitspath + filename_dict['bin_valid_rev'], format = 'fixed_width_two_line')
        
    Te_err = der_prop_err['T_e_lowhigh_error']        
    metal_err = der_prop_err['12+log(O/H)_lowhigh_error']
    T_e = der_prop_file[temp_metal_names0[0]].data
    detect = np.where(valid_tbl[bin_names0[2]].data == 1.0)[0]
        
    Te_low_err = -1*np.log10(1 - Te_err[:,0]/T_e[detect])
    Te_high_err = np.log10(1 + Te_err[:,1]/T_e[detect])
    
    err_dict = {'T_e_lowhigh_error': [Te_low_err, Te_high_err], '12+log(O/H)_lowhigh_error': [metal_err[:,0], metal_err[:,1]]}
    
    der_prop_err.close()
    
    return err_dict















