'''
Purpose:
    This code runs all other codes and the entire binning and individual analyses.
'''

import sys
import os
from getpass import getuser
from os.path import exists
from datetime import date
from astropy.io import ascii as asc
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from chun_codes.cardelli import *
from Evolution_of_Galaxies import library, emission_line_fit, R_temp_calcul, valid_table, plots, indiv_gals
from Zcalbase_gal.Analysis.DEEP2_R23_O32 import error_prop
from Zcalbase_gal import histogram_plots
from Metallicity_Stack_Commons import get_user, dir_date, fitting_lines_dict, k_dict, attenuation

##New
path = get_user()
path_init = path + 'MZEvolve/'
path_init2 = path + 'Zcalbase_gal/'
##

##Old
'''
if getuser() == 'carol':
    path_init = 'C:/Users/carol/Google Drive/MZEvolve/'  
    path_init2 = 'C:/Users/carol/Google Drive/Zcalbase_gal/'                 
else:
    path_init = ''
'''
##   
    


def get_time(org_name, bin_pts = '', run_bin = False):
    '''
    Purpose: 
        This function finds and returns the path to a directory named after the current date (MMDDYYYY).
        If the directory doesn't exist yet, it creates a new directory named after the current date in the
        provided org_name directory.
        
    Usage:
        fitspath = general.get_time(org_name)
    
    Params:
        org_name --> a string of the directory that the date subdirectory will be in.
        bin_pts (OPTIONAL) --> a string of the number of galaxies in each bin. (e.g. '75_112_113_300_600_1444_1444')
        run_bin (OPTIONAL) --> a boolean that tells if the user is running get_time() for bin analysis or
            not. If True, it creates a subdirectory called bin_pts. Default is False.
        
        
    Returns:
        fitspath --> the path to the date/bin_pts directory.
        
    Outputs:    
        "Path already exists" --> prints this message if the current date directory already exists. 
        fitspath --> prints the path to the directory.
    
    '''
    
    today = date.today()
    fitspath = path_init + org_name + '/' + "%02i%02i%02i" % (today.month, today.day, today.year) + '/' 
    try:
        os.mkdir(fitspath)
        if run_bin:
            fitspath += bin_pts + '/'
            os.mkdir(fitspath)
    except FileExistsError:
        print("Path already exists")
    print(fitspath)
    
    return fitspath


def get_HB_luminosity():
    hdul = fits.open(path_init2 + 'DEEP2_Field_combined.fits')
    fits_table = hdul[1].data
    
    cosmo = FlatLambdaCDM(H0 = 70 * u.km / u.s / u.Mpc, Om0 = 0.3)
    lum_dist = np.log10(cosmo.luminosity_distance(fits_table['ZSPEC']).to_value(u.cm))
    lum = np.log10(4 * np.pi) + np.log10(fits_table['HB_FLUX_DATA']) + (2*lum_dist)
    
    hdul.close()
    
    return lum



def run_bin_analysis(err_prop = False):
    '''
    Purpose:
        This function runs the entire binning process: binning, emission line fitting, validation table,
        electron temperature and metallicity calculations, plotting results, and error propagation.
        
    Usage:
        general.run_bin_analysis()
        Requires user input: The user must specify which binning type (mass or masslhbeta). Spelling
                             must be exact, but case does not matter. Program exits if incorrect input.
        
    Params:
        err_prop (OPTIONAL) --> True if it is desired for the error propagation code to be run, False otherwise (by default).
        
    Returns:
        None
        
    Outputs:
        Calls other codes which produce output tables, pdfs, etc. (see function calls within code).    
    '''
    
    bin_type = input('Which binning type? mass or massLHbeta: ')
    if bin_type.lower() == 'mass':
        bin_type_str = 'massbin'
        HB_lum = []
        bool_hbeta_bin = False
        #fitspath = get_time('massbin', bin_pts_name, run_bin = True)
    elif bin_type.lower() == 'masslhbeta':
        bin_type_str = 'mass_LHbeta_bin'
        HB_lum = get_HB_luminosity()
        bool_hbeta_bin = True
        #fitspath = get_time('mass_LHbeta_bin', bin_pts_name, run_bin = True)
    else:
        print('Invalid binning type')
        sys.exit(0)
      
    bin_pts_input = [75, 112, 113, 300, 600, 1444, 1444]
    str_bin_pts_input = [str(val) for val in bin_pts_input]
    bin_pts_name = "_".join(str_bin_pts_input)
    
    ##New
    fitspath = dir_date(bin_type_str, path_init, year = True)
    fitspath += bin_pts_name + '/'
    try:
        os.mkdir(fitspath)
    except FileExistsError:
        print("Path already exists")
    ##
    
    
    #Run binning (in the case of adaptive binning)
    master_grid = path_init + 'Master_Grid.fits'
    master_mask_array = path_init + 'MastermaskArray.fits'
    result_revised = asc.read(path_init + 'results_deeps_revised.tbl')
    interp_file = path_init + 'cfht_I_band_revised_interp_data.npz'
    
    mass_revised = result_revised['best.stellar.m_star'].data

    
    plt.figure(figsize=(14,8))
    edge, flux = library.binning(mass_revised, result_revised['id'], bin_pts_input, interp_file,
                                 filename = master_grid, mname = master_mask_array, fitspath0 = fitspath,
                                 spectra_plot = True, adaptive = True, hbeta_bin = bool_hbeta_bin, lum = HB_lum)
    plt.tight_layout()
    plt.savefig(fitspath + 'composite_spectra_OHmasked_interp.pdf', bbox_inches = 'tight', pad_inches = 0)
    
    hdr = fits.getheader(master_grid)
    flux_fits_file = fitspath + 'flux.fits'
    fits.writeto(flux_fits_file, flux, hdr)
    
    
    
    #Run emission line fits     
    Spect_1D, header = fits.getdata(flux_fits_file, header = True)
    dispersion = header['CDELT1']
    wave = header['CRVAL1'] + dispersion*np.arange(header['NAXIS1'])
    
    ##New
    lambda0 = fitting_lines_dict['lambda0']
    line_type = fitting_lines_dict['line_type']
    line_name = fitting_lines_dict['line_name']
    ##
    
    ##Old
    '''
    lambda0 = [3726.18, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]
    line_type = ['Oxy2', 'Balmer', 'Balmer', 'Single', 'Balmer', 'Single', 'Single']
    line_name = ['OII_3727', 'HDELTA', 'HGAMMA', 'OIII_4363', 'HBETA', 'OIII_4958', 'OIII_5007']
    '''
    ##
    
    s = 1.0
    a = 1.0
    c = 2.0
    s1 = 1.3
    a1 = 1.5
    s2 = 5
    a2 = 1.8
    
    emission_line_fit.zm_general(fitspath, Spect_1D, dispersion, wave, lambda0, line_type,
                                 line_name, s, a, c, s1, a1, s2, a2, hbeta_bin = bool_hbeta_bin)
    
    
    
    #Run validation table
    valid_table.make_validation_table(fitspath, bin_type_str)
    
    
    #Run dust attenuation
    em_file = fitspath + 'emission_lines.tbl'
    combine_asc = asc.read(em_file)
    attenuation.compute_EBV(fitspath[:-1], combine_asc)
    
    
    #Run R, Te, and Metallicity calculations 
    k_4363 = k_dict['OIII_4363']
    k_5007 = k_dict['OIII_5007']
    EBV_file = asc.read(fitspath + 'dust_attenuation_values.tbl') 
    EBV = EBV_file['E(B-V)'].data
    metal_file = fitspath + 'derived_properties_metallicity'
    R_temp_calcul.run_function(em_file, metal_file, EBV, k_4363, k_5007)
    
    
    #Run plots
    out_fname = fitspath + 'derived_properties_metallicity.pdf'
    plots.bin_derived_props_plots(metal_file + '.tbl', em_file, out_fname, bool_hbeta_bin)
    
    
    #Run error propagation, histograms, and revised data plots
    if err_prop == True:
        valid_tbl = fitspath + 'validation.tbl'
        error_prop.error_prop_chuncodes(fitspath, em_file, metal_file + '.tbl', valid_tbl)
        TM_dict_list = [fitspath + 'Te_propdist_dict.npz', fitspath + 'Te_xpeaks.npz', fitspath + 'metal_errors.npz',
                        fitspath + 'metal_xpeaks.npz', fitspath + 'metallicity_pdf.npz', fitspath + 'Te_errors.npz']
        FR_dict_list = [fitspath + 'flux_propdist.npz', fitspath + 'flux_errors.npz', fitspath + 'flux_xpeak.npz']
        histogram_plots.run_histogram_TM(fitspath, metal_file + '.tbl', TM_dict_list, valid_tbl, sharex=False)
        histogram_plots.run_histogram_FR(fitspath, metal_file + '.tbl', FR_dict_list, valid_tbl, sharex=False)
        metal_file = fitspath + 'derived_properties_metallicityrevised'
        em_file = fitspath + 'emission_linesrevised.tbl'
        out_fname = fitspath + 'derived_properties_metallicityrevised.pdf'
        plots.bin_derived_props_plots(metal_file + '.tbl', em_file, out_fname, bool_hbeta_bin)


    
 
    
#date is a string of the date in MMDDYYYY format   
def run_indiv_analysis(date_mass, date_HB):
    '''
    Purpose:
        This function runs the entire individual galaxy analysis process, which is based off both binning
        results. It calls codes to consolidate binning and individual data, to calculate individual 
        metallicities based on bin temperatures, and to plot useful relationships.
        
    Usage:
        general.run_indiv_analysis()
        
    Params:
        date_mass --> a string of the date directory that the mass bin results are in.
        date_HB --> a string of the date directory that the mass-LHbeta bin results are in.
        
    Returns:
        None
        
    Outputs:
        Calls other codes which produce output tables, pdfs, etc. (see function calls within code).     
    '''
    
    ##New
    fitspath = dir_date('individual', path_init, year = True)
    ##
    
    ##Old
    #fitspath = get_time('individual')
    ##

    
    line_file = path_init2 + 'dataset/DEEP2_all_line_fit.fits'
    mass_bin_npz = path_init + 'massbin/' + date_mass + '/75_112_113_300_600_1444_1444/binning.npz'
    mass_bin_file = path_init + 'massbin/' + date_mass + '/75_112_113_300_600_1444_1444/binning.tbl'
    mass_Te_file = path_init + 'massbin/' + date_mass + '/75_112_113_300_600_1444_1444/derived_properties_metallicityrevised.tbl'
    HB_bin_npz = path_init + 'mass_LHbeta_bin/' + date_HB + '/75_112_113_300_600_1444_1444/binning.npz'  
    HB_bin_file = path_init + 'mass_LHbeta_bin/' + date_HB + '/75_112_113_300_600_1444_1444/binning.tbl'  
    HB_Te_file = path_init + 'mass_LHbeta_bin/' + date_HB + '/75_112_113_300_600_1444_1444/derived_properties_metallicityrevised.tbl'
    
    #Create individual Te and lines table
    indiv_gals.create_Te_line_table(fitspath, line_file, mass_bin_npz, HB_bin_npz, mass_Te_file, HB_Te_file)
    
    
    #Calculate individual metallicities 
    ##Need to add EBV values and k values to large table
    EBV = np.zeros(4088)
    k_4363 = np.zeros(4088)
    k_5007 = np.zeros(4088)
    em_file = fitspath + 'individual_Te_emLines.tbl'
    metal_file = fitspath + 'individual_derived_properties_metallicity'
    R_temp_calcul.run_function(em_file, metal_file, EBV, k_4363, k_5007)
    
    
    #Make individual galaxy plots
    MTO = ''
    plots.indiv_derived_props_plots(fitspath, metal_file + '.tbl', mass_bin_file, HB_bin_file, mass_Te_file,
                                    HB_Te_file, MTO)
    MTO = '_constantMTO'
    plots.indiv_derived_props_plots(fitspath, metal_file + '.tbl', mass_bin_file, HB_bin_file, mass_Te_file,
                                    HB_Te_file, MTO, restrict_MTO = True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    