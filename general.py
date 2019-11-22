'''
This is the general code that runs all codes.
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


if getuser() == 'carol':
    path_init = 'C:\\Users\\carol\\Google Drive\\MZEvolve\\'  
    path_init2 = 'C:\\Users\\carol\\Google Drive\\Zcalbase_gal\\'                 
else:
    path_init = ''



def get_time(org_name):
    today = date.today()
    fitspath = path_init + org_name + '\\' + "%02i%02i%02i" % (today.month, today.day, today.year) + '\\'
    try:
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



def run_bin_analysis():
    bin_type = input('Which binning type? mass or massLHbeta: ')
    if bin_type.lower() == 'mass':
        bin_type_str = 'massbin'
        HB_lum = []
        bool_hbeta_bin = False
        fitspath = get_time('massbin')
    elif bin_type.lower() == 'masslhbeta':
        bin_type_str = 'massLHbetabin'
        HB_lum = get_HB_luminosity()
        bool_hbeta_bin = True
        fitspath = get_time('mass_LHbeta_bin')
    else:
        print('Invalid binning type')
        sys.exit(0)
    
    
    
    #Run binning (in the case of adaptive binning)
    master_grid = path_init + 'Master_Grid.fits'
    master_mask_array = path_init + 'MastermaskArray.fits'
    result_revised = asc.read(path_init + 'results_deeps_revised.tbl')
    interp_file = path_init + 'cfht_I_band_revised_interp_data.npz'
    
    mass_revised = result_revised['best.stellar.m_star'].data
    
    bin_pts_input = [75, 112, 113, 300, 600, 1444, 1444]
    str_bin_pts_input = [str(val) for val in bin_pts_input]
    bin_pts_fname = "_".join(str_bin_pts_input)
    bin_pts_fname = bin_type_str + '_revised_' + bin_pts_fname

    
    plt.figure(figsize=(14,8))
    edge, flux = library.binning(mass_revised, result_revised['id'], bin_pts_input, interp_file, bin_pts_fname,
                                 filename = master_grid, mname = master_mask_array, fitspath0 = fitspath,
                                 spectra_plot = True, adaptive = True, hbeta_bin = bool_hbeta_bin, lum = HB_lum)
    plt.tight_layout()
    plt.savefig(fitspath + bin_pts_fname + '_composite_spectra_OHmasked_interp.pdf', bbox_inches = 'tight', pad_inches = 0)
    
    hdr = fits.getheader(master_grid)
    flux_fits_file = fitspath + bin_pts_fname + '_flux.fits'
    fits.writeto(flux_fits_file, flux, hdr)
    
    
    
    #Run emission line fits     
    Spect_1D, header = fits.getdata(flux_fits_file, header = True)
    dispersion = header['CDELT1']
    wave = header['CRVAL1'] + dispersion*np.arange(header['NAXIS1'])
    
    lambda0 = [3726.18, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]
    line_type = ['Oxy2', 'Balmer', 'Balmer', 'Single', 'Balmer', 'Single', 'Single']
    line_name = ['OII_3727', 'HDELTA', 'HGAMMA', 'OIII_4363', 'HBETA', 'OIII_4958', 'OIII_5007']
    
    s = 1.0
    a = 1.0
    c = 2.0
    s1 = 1.3
    a1 = 1.5
    s2 = 5
    a2 = 1.8
    
    emission_line_fit.zm_general(fitspath, bin_pts_fname, Spect_1D, header, dispersion, wave, lambda0, line_type,
                                 line_name, s, a, c, s1, a1, s2, a2, hbeta_bin = bool_hbeta_bin)
    
    
    
    #Run validation table
    valid_table.make_validation_table(fitspath, bin_pts_fname)
    
    
    #Run dust attenuation
    #change later once function is implemented and get actual values
    EBV = np.zeros(len(edge))
    k_4363 = np.zeros(len(edge))
    k_5007 = np.zeros(len(edge))
    '''
    em_file = fitspath + bin_pts_fname + '_emission_lines.tbl'
    dust_attenuation(fitspath, bin_pts_fname, em_file)
    '''
    
    
    #Run R, Te, and Metallicity calculations        
    em_file = fitspath + bin_pts_fname + '_emission_lines.tbl'
    metal_file = fitspath + bin_pts_fname + '_derived_properties_metallicity'
    R_temp_calcul.run_function(em_file, metal_file, EBV, k_4363, k_5007)
    
    
    #Run plots
    out_fname = fitspath + bin_pts_fname + '_derived_properties_metallicity.pdf'
    plots.bin_derived_props_plots(metal_file + '.tbl', em_file, out_fname)
    
    
    #Run error propagation
    error_prop.error_prop_chuncodes(fitspath, em_file, metal_file + '.tbl')
    dict_list = [fitspath + 'Te_propdist_dict.npz', fitspath + 'Te_xpeaks.npz', fitspath + 'metal_errors.npz',
                 fitspath + 'metal_xpeaks.npz', fitspath + 'metallicity_pdf.npz', fitspath + 'flux_propdist.npz',
                 fitspath + 'flux_errors.npz', fitspath + 'Te_errors.npz']
    histogram_plots.run_histogram(fitspath, metal_file + '.tbl', dict_list)
    
    
 
    
   
def run_indiv_analysis():
    fitspath = get_time('individual')
    
    bin_pts_input = [75, 112, 113, 300, 600, 1444, 1444]
    str_bin_pts_input = [str(val) for val in bin_pts_input]
    bin_pts_fname = "_".join(str_bin_pts_input)
    bin_pts_fname = 'revised_' + bin_pts_fname
    
    #fitspath += bin_pts_fname + '\\' + bin_type_str
    
    line_file = path_init2 + 'dataset\\DEEP2_all_line_fit.fits'
    mass_bin_npz = path_init + 'massbin\\11212019\\massbin_revised_75_112_113_300_600_1444_1444.npz'
    mass_bin_file = path_init + 'massbin\\11212019\\massbin_revised_75_112_113_300_600_1444_1444_binning.tbl'
    mass_Te_file = path_init + 'massbin\\11212019\\massbin_revised_75_112_113_300_600_1444_1444_derived_properties_metallicity.tbl'
    HB_bin_npz = path_init + 'mass_LHbeta_bin\\11212019\\massLHbetabin_revised_75_112_113_300_600_1444_1444.npz'
    HB_bin_file = path_init + 'mass_LHbeta_bin\\11212019\\massLHbetabin_revised_75_112_113_300_600_1444_1444_binning.tbl'
    HB_Te_file = path_init + 'mass_LHbeta_bin\\11212019\\massLHbetabin_revised_75_112_113_300_600_1444_1444_derived_properties_metallicity.tbl'
    
    #Create individual Te and lines table
    indiv_gals.create_Te_lineratio_table(fitspath, line_file, mass_bin_npz, HB_bin_npz, mass_Te_file, HB_Te_file)
    
    
    #Calculate individual metallicities 
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    