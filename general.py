import sys
import os
from astropy.io import ascii as asc
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

from .Analysis import library, emission_line_fit, indiv_gals, R_temp_calcul
from .Plotting.composite_plots import bin_derived_props_plots
from .Plotting.individual_plots import indiv_derived_props_plots, indiv_metal_mass_plot, get_indiv_detect
from .Plotting.relation_fitting import extract_error_bars

from Metallicity_Stack_Commons import get_user, dir_date, fitting_lines_dict
from Metallicity_Stack_Commons.column_names import filename_dict, indv_names0, temp_metal_names0
from Metallicity_Stack_Commons.column_names import bin_mzevolve_names0, bin_names0
from Metallicity_Stack_Commons.analysis.error_prop import fluxes_derived_prop
from Metallicity_Stack_Commons.analysis.composite_indv_detect import main
from Metallicity_Stack_Commons.plotting.balmer import HbHgHd_fits
from Metallicity_Stack_Commons import valid_table


path = get_user()
path_init = path + 'MZEvolve/'
path_init2 = path + 'Zcalbase_gal/'


def get_HB_luminosity():
    hdul = fits.open(path_init2 + 'DEEP2_Field_combined.fits')
    fits_table = hdul[1].data
    
    cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
    lum_dist = np.log10(cosmo.luminosity_distance(fits_table['ZSPEC']).to_value(u.cm))
    lum = np.log10(4 * np.pi) + np.log10(fits_table['HB_FLUX_DATA']) + (2*lum_dist)
    
    hdul.close()
    
    return lum   


def run_bin_analysis(valid_rev=False, dust_atten=False, err_prop=False, indiv=False):
    '''
    Purpose:
        This function runs the entire binning process: binning, emission line fitting, validation table,
        electron temperature and metallicity calculations, plotting results, and error propagation.
        
    Usage:
        general.run_bin_analysis()
        Requires user input: The user must specify which binning type (mass or masslhbeta). Spelling
                             must be exact, but case does not matter. Program exits if incorrect input.
        
    Params:
        err_prop (OPTIONAL) --> True if it is desired for the error propagation code to be run, False 
                                otherwise (by default).
        
    Returns:
        None
        
    Outputs:
        Calls other codes which produce output tables, pdfs, etc. (see function calls within code).  
    '''
    
    bin_type = input('Which binning type? mass or massLHbeta: ')
    if bin_type.lower() == 'mass':
        bin_type_str = 'massbin'
        bool_hbeta_bin = False
    elif bin_type.lower() == 'masslhbeta':
        bin_type_str = 'mass_LHbeta_bin'
        bool_hbeta_bin = True
    else:
        print('Invalid binning type')
        sys.exit(0)
        
    HB_lum = get_HB_luminosity()
        
    # Make working directory/get fitspath    
    bin_pts_input = [75, 112, 113, 300, 600, 1444, 1443]
    str_bin_pts_input = [str(val) for val in bin_pts_input]
    bin_pts_name = "_".join(str_bin_pts_input)
    
    fitspath = dir_date(bin_type_str, path_init, year=True)
    fitspath += bin_pts_name + '/'
    try:
        os.mkdir(fitspath)
    except FileExistsError:
        print("Path already exists")
    
    
    # Run binning (in the case of adaptive binning)
    master_grid = path_init + 'Master_Grid.fits'
    master_mask_array = path_init + 'MastermaskArray.fits'
    result_revised = asc.read(path_init + 'results_deeps_revised.tbl')
    interp_file = path_init + 'cfht_I_band_revised_interp_data.npz'
    
    mass_revised = result_revised['best.stellar.m_star'].data
 
    plt.figure(figsize=(14,8))
    flux = library.binning(mass_revised, result_revised['id'], bin_pts_input, interp_file,
                           filename=master_grid, mname=master_mask_array, fitspath0=fitspath,
                           spectra_plot=True, adaptive=True, hbeta_bin=bool_hbeta_bin, lum=HB_lum)
    plt.tight_layout()
    plt.savefig(fitspath + 'composite_spectra_OHmasked_interp.pdf', bbox_inches='tight', pad_inches=0)
    
    hdr = fits.getheader(master_grid)
    flux_fits_file = fitspath + filename_dict['comp_spec']
    fits.writeto(flux_fits_file, flux, hdr)
    
    
    
    # Run emission line fits     
    Spect_1D, header = fits.getdata(flux_fits_file, header=True)
    dispersion = header['CDELT1']
    wave = header['CRVAL1'] + dispersion*np.arange(header['NAXIS1'])
    lambda0 = fitting_lines_dict['lambda0']
    line_type = fitting_lines_dict['line_type']
    line_name = fitting_lines_dict['line_name']
    s = 1.0
    a = 1.0
    c = 2.0
    s1 = 1.3
    a1 = 1.5
    s2 = 5
    a2 = 1.8
    
    emission_line_fit.zm_general(fitspath, Spect_1D, dispersion, wave, lambda0, line_type,
                                 line_name, s, a, c, s1, a1, s2, a2, hbeta_bin=bool_hbeta_bin)
    
    
    
    # Run validation table
    valid_table.make_validation_table(fitspath)
    valid_file = fitspath + filename_dict['bin_valid']
    valid_rev_file = fitspath + filename_dict['bin_valid_rev']
    if bool_hbeta_bin:
        vtbl = asc.read(valid_file, format='fixed_width_two_line')
        detect = vtbl['Detection'].data
        detect[11] = 0.5
        vtbl.replace_column('Detection', detect)
        asc.write(vtbl, valid_rev_file, format='fixed_width_two_line', overwrite=True)
    else:
        vtbl = asc.read(valid_file, format='fixed_width_two_line')
        asc.write(vtbl, valid_rev_file, format='fixed_width_two_line', overwrite=True)
        
    
    # Run raw data derived properties calculations (option to apply dust correction)
    fluxes_derived_prop(fitspath, raw=True, binned_data=True, apply_dust=False, revised=False)
    fluxes_derived_prop(fitspath, raw=True, binned_data=True, apply_dust=False, revised=True)
    if dust_atten:
        fluxes_derived_prop(fitspath, raw=True, binned_data=True, apply_dust=True, revised=False)
        fluxes_derived_prop(fitspath, raw=True, binned_data=True, apply_dust=True, revised=True)
        
        
    # Run Monte Carlo randomization calculations (option to apply dust correction)
    fluxes_derived_prop(fitspath, raw=False, binned_data=True, apply_dust=False, revised=False)
    fluxes_derived_prop(fitspath, raw=False, binned_data=True, apply_dust=False, revised=True)
    if dust_atten:
        fluxes_derived_prop(fitspath, raw=False, binned_data=True, apply_dust=True, revised=False)
        fluxes_derived_prop(fitspath, raw=False, binned_data=True, apply_dust=True, revised=True)
        
    '''   
    # Run plots
    bin_derived_props_plots(fitspath, hbeta_bin=bool_hbeta_bin, err_bars=False, revised=False)  #raw data
    if dust_atten:
        HbHgHd_fits(fitspath)'''
        
    
    '''
    # Run dust attenuation
    bin_file = fitspath + filename_dict['bin_info']
    em_file = fitspath + filename_dict['bin_fit']
    HbHgHd_fits(fitspath)
    '''
    
    '''
    metal_file = fitspath + filename_dict['bin_derived_prop']
    R_temp_calcul.run_function(em_file, bin_file, metal_file)
    '''

    
    '''
    #Run error propagation and revised data plots
    if err_prop:
        fluxes_derived_prop(fitspath, binned_data=True)
        bin_derived_props_plots(fitspath, hbeta_bin=bool_hbeta_bin, err_bars=True, revised=True)
        HbHgHd_fits(fitspath, use_revised=True)
    '''
      
    '''    
    if indiv:
        # Create individual_bin_info table
        line_file = path_init2 + 'All Datasets/DEEP2_all_line_fit.fits'
        indiv_gals.indiv_bin_info_table(fitspath, line_file, use_revised=True)
 
        # Create individual_properties table
        indiv_gals.indiv_em_table(fitspath, line_file)
        
        # Create individual_derived_properties table
        main(fitspath, '', revised=False, det3=True)
        
        # Run individual plots
        indiv_derived_props_plots(fitspath, restrictMTO=True, revised=False, err_bars=False, 
                                  hbeta_bin=bool_hbeta_bin)
        indiv_derived_props_plots(fitspath, restrictMTO=False, revised=False, err_bars=False, 
                                  hbeta_bin=bool_hbeta_bin)
    '''
        
        

def run_indiv_metal_mass_plot(fitspathM, fitspathMLHb, restrictMTO=False, revised_files=False, 
                              error_bars=False):
    '''
    Purpose:
        This function extracts the necessary data from files and runs the individual Metallicity vs Mass
        plotting function, which produces a two-paneled Metallicity vs Mass plot containing data from mass
        bins and mass-LHbeta bins.
           
    Params:
        fitspathM --> a string containing the partial path of the location of the mass bin data.
                      (e.g., 'massbin/05182020/75_112_113_300_600_1444_1443/')
        fitspathMLHb --> a string containing the partial path of the location of the mass-LHbeta bin data.
                         (e.g., 'mass_LHbeta_bin/05182020/75_112_113_300_600_1444_1443/')
        restrictMTO (OPTIONAL) --> if the mass turnover value should be held constant in the curve fit of
                                   Metallicity vs Mass, then restrictMTO=True.
        revised_files (OPTIONAL) --> if the revised data tables should be used, then revised_files=True.
        error_bars (OPTIONAL) --> if error bars for metallicity and temperature should be plotted, then 
                                  error_bars=True.
        
    Returns:
        None
    '''
    
    # Read in individual data tables
    Mbin_indiv_derivedprops_tbl = asc.read(path_init + fitspathM + filename_dict['indv_derived_prop'], 
                                           format='fixed_width_two_line')
    Mbin_indiv_props_tbl = asc.read(path_init + fitspathM + filename_dict['indv_prop'], 
                                    format='fixed_width_two_line')
    Mbin_indiv_bininfo_tbl = asc.read(path_init + fitspathM + filename_dict['indv_bin_info'], 
                                      format='fixed_width_two_line')
    
    MLHbbin_indiv_derivedprops_tbl = asc.read(path_init + fitspathMLHb + filename_dict['indv_derived_prop'], 
                                              format='fixed_width_two_line')
    MLHbbin_indiv_props_tbl = asc.read(path_init + fitspathMLHb + filename_dict['indv_prop'], 
                                       format='fixed_width_two_line')
    MLHbbin_indiv_bininfo_tbl = asc.read(path_init + fitspathMLHb + filename_dict['indv_bin_info'], 
                                         format='fixed_width_two_line')
    
    
    # Read in composite data tables
    Mbin_valid_tbl = asc.read(path_init + fitspathM + filename_dict['bin_valid'], 
                              format='fixed_width_two_line')
    Mbininfo_tbl = asc.read(path_init + fitspathM + filename_dict['bin_info'], format='fixed_width_two_line')
    
    MLHbbin_valid_tbl = asc.read(path_init + fitspathMLHb + filename_dict['bin_valid'], 
                                 format='fixed_width_two_line')
    MLHbbininfo_tbl = asc.read(path_init + fitspathMLHb + filename_dict['bin_info'], 
                               format='fixed_width_two_line')
    if revised_files:
        Mbin_derivedprops_tbl = asc.read(path_init + fitspathM + filename_dict['bin_derived_prop_rev'], 
                                         format='fixed_width_two_line')
        MLHbbin_derivedprops_tbl = asc.read(path_init + fitspathMLHb + filename_dict['bin_derived_prop_rev'], 
                                            format='fixed_width_two_line')
    else:    
        Mbin_derivedprops_tbl = asc.read(path_init + fitspathM + filename_dict['bin_derived_prop'], 
                                         format='fixed_width_two_line')
        MLHbbin_derivedprops_tbl = asc.read(path_init + fitspathMLHb + filename_dict['bin_derived_prop'], 
                                            format='fixed_width_two_line')


    # Read in individual data
    Mbin_indiv_logM = Mbin_indiv_props_tbl[indv_names0[3]].data
    MLHbbin_indiv_logM = MLHbbin_indiv_props_tbl[indv_names0[3]].data
    
    Mbin_indiv_metal = Mbin_indiv_derivedprops_tbl[temp_metal_names0[1]].data
    Mbin_indiv_bin_detect = Mbin_indiv_bininfo_tbl[bin_names0[2]].data
    
    MLHbbin_indiv_metal = MLHbbin_indiv_derivedprops_tbl[temp_metal_names0[1]].data
    MLHbbin_indiv_bin_detect = MLHbbin_indiv_bininfo_tbl[bin_names0[2]].data
    
    
    # Read in composite data
    Mbin_logM = Mbininfo_tbl[bin_mzevolve_names0[2]].data      
    Mbin_metal = Mbin_derivedprops_tbl[temp_metal_names0[1]].data
    Mbin_detect_col = Mbin_valid_tbl[bin_names0[2]].data
    
    MLHbbin_logM = MLHbbininfo_tbl[bin_mzevolve_names0[2]].data      
    MLHbbin_metal = MLHbbin_derivedprops_tbl[temp_metal_names0[1]].data
    MLHbbin_detect_col = MLHbbin_valid_tbl[bin_names0[2]].data
    
    
    # Define detection and non-detection (with reliable 5007) arrays for bins
    Mbin_detect = np.where(Mbin_detect_col == 1.0)[0]
    Mbin_nondetect = np.where(Mbin_detect_col == 0.5)[0]
    
    MLHbbin_detect = np.where(MLHbbin_detect_col == 1.0)[0]
    MLHbbin_nondetect = np.where(MLHbbin_detect_col == 0.5)[0]
    
    # Define detection and non-detection (with reliable 5007) arrays for individual galaxies
    Mbin_indiv_detect, Mbin_indiv_nondetect = get_indiv_detect(Mbin_indiv_props_tbl, Mbin_indiv_bin_detect)
    MLHbbin_indiv_detect, MLHbbin_indiv_nondetect = get_indiv_detect(MLHbbin_indiv_props_tbl, 
                                                                     MLHbbin_indiv_bin_detect, 
                                                                     LHbeta_bins=True)
    
    
    # Define mass bin and mass-LHbeta bin dictionaries
    Mbin_dict = {'path': path_init + fitspathM, 'composite_logM': Mbin_logM, 'composite_detect': Mbin_detect, 
                 'composite_nondetect': Mbin_nondetect, 'composite_metallicity': Mbin_metal, 
                 'indiv_logM': Mbin_indiv_logM, 'indiv_detect': Mbin_indiv_detect, 
                 'indiv_nondetect': Mbin_indiv_nondetect, 'indiv_metallicity': Mbin_indiv_metal}
    MLHbbin_dict = {'path': path_init + fitspathMLHb, 'composite_logM': MLHbbin_logM, 
                    'composite_detect': MLHbbin_detect, 'composite_nondetect': MLHbbin_nondetect, 
                    'composite_metallicity': MLHbbin_metal, 'indiv_logM': MLHbbin_indiv_logM, 
                    'indiv_detect': MLHbbin_indiv_detect, 'indiv_nondetect': MLHbbin_indiv_nondetect, 
                    'indiv_metallicity': MLHbbin_indiv_metal}
    if error_bars:
        err_dictM = extract_error_bars(fitspathM)
        err_dictMLHb = extract_error_bars(fitspathM)
        Mbin_dict['composite_metal_errors'] = err_dictM['12+log(O/H)_lowhigh_error']
        MLHbbin_dict['composite_metal_errors'] = err_dictMLHb['12+log(O/H)_lowhigh_error']
    
    
    indiv_metal_mass_plot(Mbin_dict, MLHbbin_dict, restrictMTO=restrictMTO, revised=revised_files, 
                          err_bars=error_bars)
