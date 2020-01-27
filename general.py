'''
Purpose:
    This code runs all other codes and the entire binning and individual analyses.
'''

import sys
import os
#from getpass import getuser
#from os.path import exists
from datetime import date
from astropy.io import ascii as asc
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.table import Table
#from chun_codes.cardelli import *
from Evolution_of_Galaxies import library, emission_line_fit, valid_table, plots, indiv_gals
from Zcalbase_gal.Analysis.DEEP2_R23_O32 import error_prop
from Zcalbase_gal import histogram_plots
from Metallicity_Stack_Commons import get_user, dir_date, fitting_lines_dict, k_dict, attenuation
from Metallicity_Stack_Commons.temp_metallicity_calc import R_calculation, temp_calculation, metallicity_calculation

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




def R_temp_metal_calc(line_file, outfile, EBV, k_4363, k_5007):
    '''
    Purpose:
        This function runs the R calculation, temperature calculation, and metallicity calculation
        functions for the case of individual and stacked spectra.
        
    Usage:
        R_temp_calcul.run_function(line_file, outfile, EBV, k_4363, k_5007)
        
    Params:
        line_file --> ascii table that has all the data for each emission line.
        outfile --> a string naming the tables that are produced by this function.
        EBV --> an array of E(B - V) values. The array is the size of the number of sources (or bins).
        k_4363 --> an array of dust extinction values at the OIII4363 emission line for each bin/individual source.
        k_5007 --> an array of dust extinction values at the OIII5007 emission line for each bin/individual source.
        
    Returns:
        None
        
    Outputs:
        out_ascii --> an ascii table containing the bin temperatures, bin or individual metallicities,
            and bin or individual line fluxes and S/N.
        out_fits --> a fits table containing the bin temperatures, bin or individual metallicities,
            and bin or individual line fluxes and S/N.
    '''
    
    line_table = asc.read(line_file)
    
    if 'Log10(Mass)' in line_table.keys():
        #Case for individual spectra
        out_ascii = outfile + '.tbl'
        out_fits = outfile + '.fits'
        
        OII = line_table['OII_Flux'].data
        SN_OII = line_table['OII_SN'].data
        OIII4959 = line_table['OIII4959_Flux'].data
        SN_4959 = line_table['OIII4959_SN'].data
        OIII5007 = line_table['OIII5007_Flux'].data
        SN_5007 = line_table['OIII5007_SN'].data
        HBETA = line_table['HBETA_Flux'].data
        SN_HBETA = line_table['HBETA_SN'].data
        
        source_ID = line_table['OBJNO'].data
        mass_bin_ID = line_table['Mass_Bin_ID'].data
        HB_bin_ID = line_table['Mass_LHBeta_Bin_ID'].data
        log_mass = line_table['Log10(Mass)'].data
        LHbeta = line_table['HBeta_Luminosity'].data
        mass_T_e = line_table['Mass_Bin_Te'].data
        HB_T_e = line_table['Mass_LHBeta_Bin_Te'].data
        HB_bin_detect = line_table['Mass_LHBeta_Bin_Detections'].data
        mass_bin_detect = line_table['Mass_Bin_Detections'].data
        mass_indiv_detect = line_table['Mass_Individual_Detections'].data 
        HB_indiv_detect = line_table['MassLHB_Individual_Detections'].data 
        ebv = line_table['E(B-V)'].data
        
        #Detection variables here include the non-detections with reliable limits so that the metallicity
        #is calculated for those sources.
        mass_detect = np.where(((mass_bin_detect == 1.0) | (mass_bin_detect == 0.5)) & (np.isfinite(OIII5007) == True) &
                               (OIII5007 >= 1e-18) & (OIII5007 <= 1e-15) & (np.isfinite(OII) == True) &
                               (OII >= 1e-18) & (OII <= 1e-15) & (np.isfinite(HBETA) == True) & 
                               (HBETA >= 1e-18) & (HBETA <= 1e-15))[0]
        HB_detect = np.where(((HB_bin_detect == 1.0) | (HB_bin_detect == 0.5)) & (np.isfinite(OIII5007) == True) & 
                             (OIII5007 >= 1e-18) & (OIII5007 <= 1e-15) & (np.isfinite(OII) == True) & 
                             (OII >= 1e-18) & (OII <= 1e-15) & (np.isfinite(HBETA) == True) & 
                             (HBETA >= 1e-18) & (HBETA <= 1e-15) & (LHbeta > 0))[0]
        mass_nondetect = np.where((mass_bin_detect == 0.5) & (np.isfinite(OIII5007) == True) & 
                                  (OIII5007 >= 1e-18) & (OIII5007 <= 1e-15) & (np.isfinite(OII) == True) & 
                                  (OII >= 1e-18) & (OII <= 1e-15) & (np.isfinite(HBETA) == True) & 
                                  (HBETA >= 1e-18) & (HBETA <= 1e-15))[0]
        HB_nondetect = np.where((HB_bin_detect == 0.5) & (np.isfinite(OIII5007) == True) & 
                                (OIII5007 >= 1e-18) & (OIII5007 <= 1e-15) & (np.isfinite(OII) == True) & 
                                (OII >= 1e-18) & (OII <= 1e-15) & (np.isfinite(HBETA) == True) & 
                                (HBETA >= 1e-18) & (HBETA <= 1e-15) & (LHbeta > 0))[0]
        
        #Correct detection markings (detection vs non-detection w/ reliable limits) are applied here
        mass_indiv_detect[mass_detect] = 1.0
        HB_indiv_detect[HB_detect] = 1.0
        mass_indiv_detect[mass_nondetect] = 0.5
        HB_indiv_detect[HB_nondetect] = 0.5
        
        #create zero arrays all same length
        two_beta = np.zeros(len(source_ID))
        three_beta = np.zeros(len(source_ID))
        R23 = np.zeros(len(source_ID))
        O32 = np.zeros(len(source_ID))  
        
        two_beta[mass_detect] = OII[mass_detect] / HBETA[mass_detect]
        two_beta[HB_detect] = OII[HB_detect] / HBETA[HB_detect]
        three_beta[mass_detect] = (OIII5007[mass_detect] * (1 + 1/3.1)) / HBETA[mass_detect]
        three_beta[HB_detect] = (OIII5007[HB_detect] * (1 + 1/3.1)) / HBETA[HB_detect]
        
        #Calculate R23 and O32
        R23[mass_detect] = np.log10((OII[mass_detect] + ((1 + 1/3.1) * OIII5007[mass_detect])) / HBETA[mass_detect])
        R23[HB_detect] = np.log10((OII[HB_detect] + ((1 + 1/3.1) * OIII5007[HB_detect])) / HBETA[HB_detect])
        O32[mass_detect] = np.log10(((1 + 1/3.1) * OIII5007[mass_detect]) / OII[mass_detect])
        O32[HB_detect] = np.log10(((1 + 1/3.1) * OIII5007[HB_detect]) / OII[HB_detect])
        
        
        mass_O_s_ion, mass_O_d_ion, mass_com_O_log, mass_log_O_s, mass_log_O_d = metallicity_calculation(mass_T_e, two_beta, three_beta)
        HB_O_s_ion, HB_O_d_ion, HB_com_O_log, HB_log_O_s, HB_log_O_d = metallicity_calculation(HB_T_e, two_beta, three_beta)
        
        n = ('Source_ID', 'Mass_Bin_ID', 'Mass_LHBeta_Bin_ID', 'Mass_Bin_Detections', 'Mass_LHBeta_Bin_Detections',
             'Mass_Individual_Detections', 'MassLHB_Individual_Detections', 'Log10(Mass)', 'HBeta_Luminosity', 'OIII_5007_Flux_Observed', 
             'OIII_4958_Flux_Observed', 'OII_3727_Flux_Observed', 'HBETA_Flux_Observed', 'Mass_Bin_Te',
             'Mass_LHBeta_Bin_Te', 'R23', 'O32', 'OII/HBeta', 'OIII/HBeta', 'Mass_Bin_O_s_ion',
             'Mass_Bin_O_d_ion', 'Mass_Bin_com_O_log', 'Mass_LHBeta_Bin_O_s_ion', 'Mass_LHBeta_Bin_O_d_ion',
             'Mass_LHBeta_Bin_com_O_log', 'E(B-V)')
        tab0 = Table([source_ID, mass_bin_ID, HB_bin_ID, mass_bin_detect, HB_bin_detect,
                      mass_indiv_detect, HB_indiv_detect, log_mass, LHbeta, OIII5007, OIII4959, OII, HBETA,
                      mass_T_e, HB_T_e, R23, O32, two_beta, three_beta, mass_O_s_ion, mass_O_d_ion, 
                      mass_com_O_log, HB_O_s_ion, HB_O_d_ion, HB_com_O_log, ebv], names = n)
        
    else:
        #Case for stacked spectra
        out_ascii = outfile + '.tbl'
        out_fits = outfile + '.fits'
        
        OII = line_table['OII_3727_Flux_Observed'].data
        SN_OII = line_table['OII_3727_S/N'].data
        OIII4363 = line_table['OIII_4363_Flux_Observed'].data
        SN_4363 = line_table['OIII_4363_S/N'].data
        OIII4959 = line_table['OIII_4958_Flux_Observed'].data
        SN_4959 = line_table['OIII_4958_S/N'].data
        OIII5007 = line_table['OIII_5007_Flux_Observed'].data
        SN_5007 = line_table['OIII_5007_S/N'].data
        HBETA = line_table['HBETA_Flux_Observed'].data
        SN_HBETA = line_table['HBETA_S/N'].data
        
        N_Galaxy = line_table['Number of Galaxies'].data
        avg_mass = line_table['mass_avg'].data 
        detection = line_table['Detection'].data
        ID = line_table['ID'].data
        log_mass = avg_mass
        
        two_beta = OII / HBETA
        three_beta = (OIII5007 * (1 + 1/3.1)) / HBETA
        
        #Calculate R23 composite and O32 composite
        R23_composite = np.log10((OII + ((1 + 1/3.1) * OIII5007)) / HBETA)
        O32_composite = np.log10(((1 + 1/3.1) * OIII5007) / OII)
        
        #R, Te, and metallicity calculations
        R_value = R_calculation(OIII4363, OIII5007, EBV, k_4363, k_5007)
        T_e = temp_calculation(R_value)
        O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metallicity_calculation(T_e, two_beta, three_beta)        
        
        n = ('ID', 'Detection', 'R23_Composite', 'O32_Composite', 'N_Galaxies', 'OIII_5007_Flux_Observed', 
             'OIII_5007_S/N', 'OIII_4958_Flux_Observed', 'OIII_4958_S/N', 'OIII_4363_Flux_Observed', 'OIII_4363_S/N',
             'HBETA_Flux_Observed', 'HBETA_S/N', 'OII_3727_Flux_Observed', 'OII_3727_S/N', 'Temperature',
             'log_O_s', 'log_O_d', 'O_s_ion', 'O_d_ion', 'com_O_log')
        tab0 = Table([ID, detection, R23_composite, O32_composite, N_Galaxy, OIII5007, SN_5007, OIII4959,
                      SN_4959, OIII4363, SN_4363, HBETA, SN_HBETA, OII, SN_OII, T_e, log_O_s, log_O_d,
                      O_s_ion, O_d_ion, com_O_log], names = n)

    
    asc.write(tab0, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    tab0.write(out_fits, format = 'fits', overwrite = True)    




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
    ##Old
    '''
    bin_pts_input = [75, 112, 113, 300, 600, 1444, 1444]
    str_bin_pts_input = [str(val) for val in bin_pts_input]
    bin_pts_name = "_".join(str_bin_pts_input)
    '''
    ##
    
    bin_type = input('Which binning type? mass or massLHbeta: ')
    if bin_type.lower() == 'mass':
        bin_type_str = 'massbin'
        HB_lum = []
        bool_hbeta_bin = False
        #fitspath = get_time('massbin', bin_pts_name, run_bin = True)  ###
    elif bin_type.lower() == 'masslhbeta':
        bin_type_str = 'mass_LHbeta_bin'
        HB_lum = get_HB_luminosity()
        bool_hbeta_bin = True
        #fitspath = get_time('mass_LHbeta_bin', bin_pts_name, run_bin = True)  ###
    else:
        print('Invalid binning type')
        sys.exit(0)
        
    ##New
    bin_pts_input = [75, 112, 113, 300, 600, 1444, 1444]
    str_bin_pts_input = [str(val) for val in bin_pts_input]
    bin_pts_name = "_".join(str_bin_pts_input)
    
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
    ##Old
    EBV = np.zeros(len(edge))
    k_4363 = np.zeros(len(edge))
    k_5007 = np.zeros(len(edge))
    ##
    em_file = fitspath + 'emission_lines.tbl'
    ##New
    '''
    combine_asc = asc.read(em_file)
    attenuation.compute_EBV(fitspath[:-1], combine_asc)
    '''
    ##
    
    
    #Run R, Te, and Metallicity calculations 
    ##New
    '''
    k_4363 = k_dict['OIII_4363']
    k_5007 = k_dict['OIII_5007']
    EBV_file = asc.read(fitspath + 'dust_attenuation_values.tbl') 
    EBV = EBV_file['E(B-V)'].data
    '''
    ##
    metal_file = fitspath + 'derived_properties_metallicity'
    ##Old  R_temp_calcul.run_function(em_file, metal_file, EBV, k_4363, k_5007)
    R_temp_metal_calc(em_file, metal_file, EBV, k_4363, k_5007)
    
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
    ##Old  R_temp_calcul.run_function(em_file, metal_file, EBV, k_4363, k_5007)
    R_temp_metal_calc(em_file, metal_file, EBV, k_4363, k_5007)
    
    
    #Make individual galaxy plots
    MTO = ''
    plots.indiv_derived_props_plots(fitspath, metal_file + '.tbl', mass_bin_file, HB_bin_file, mass_Te_file,
                                    HB_Te_file, MTO)
    MTO = '_constantMTO'
    plots.indiv_derived_props_plots(fitspath, metal_file + '.tbl', mass_bin_file, HB_bin_file, mass_Te_file,
                                    HB_Te_file, MTO, restrict_MTO = True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    