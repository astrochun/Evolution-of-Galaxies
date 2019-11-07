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
from Evolution_of_Galaxies import library, emission_line_fit, R_temp_calcul, valid_table, plots #, indiv_gals
#from Zcalbase_gal import error_prop


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
    EBV = np.zeros(len(line_name))
    k_4363 = np.zeros(len(line_name))
    k_5007 = np.zeros(len(line_name))
    
    
    #Run R, Te, and Metallicity calculations        
    em_file = fitspath + bin_pts_fname + '_emission_lines.tbl'
    outfile = fitspath + bin_pts_fname + '_derived_properties_metallicity'
    R_temp_calcul.run_function(em_file, outfile, EBV, k_4363, k_5007)
    
    
    #Run plots
    plots.bin_derived_props_plots(fitspath, bin_pts_fname)
    
    
    '''
    #Run error propagation
    error_prop.error_prop_chuncodes(fitspath + bin_pts_fname, em_file, outfile + '.tbl')
    '''
    

    
    
'''    
#combine_ascii is emission_line.tbl for me
def dust_attenuation(fitspath, bin_pts_fname, combine_ascii):
    line_name = ['OII_3727', 'HDELTA', 'HGAMMA', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']
    
    combine_asc = asc.read(combine_ascii)
    ini_con = 0.468
    ID = combine_asc['ID']
    HBeta = combine_asc['HBETA_Flux_Observed'].data
    HGamma = combine_asc['HGAMMA_Flux_Observed'].data
    
    lam0_OII = combine_asc['OII_3727_X_bar'].data
    lam0_HDELTA = combine_asc['HDELTA_X_bar'].data
    lam0_Hgamma = combine_asc['HGAMMA_X_bar'].data
    lam0_HBETA = combine_asc['HBETA_X_bar'].data
    lam0_4363 = combine_asc['OIII_4363_X_bar'].data
    lam0_4958 = combine_asc['OIII_4958_X_bar'].data
    lam0_5007 = combine_asc['OIII_5007_X_bar'].data

    k_3727 = call_cardelli(lam0_OII)
    k_HDELTA = call_cardelli(lam0_HDELTA)
    k_Hgamma = call_cardelli(lam0_Hgamma)
    k_HBETA = call_cardelli(lam0_HBETA)
    k_4363 = call_cardelli(lam0_4363)
    k_4958 = call_cardelli(lam0_4958)
    k_5007 = call_cardelli(lam0_5007)

    
    
    EBV= np.log10((HBeta/HGamma)*(ini_con))*2.5*(1/(k_Hgamma-k_HBETA))
    for nn in range(len(HGamma)):
        if EBV[nn] <= 0: EBV[nn] = 0
    
    #print EBV
    A_3727 = EBV*k_3727
    A_HDELTA = EBV*k_HDELTA
    A_Hgamma = EBV*k_Hgamma
    A_HBETA = EBV*k_HBETA
    A_4363 = EBV*k_4363
    A_4958 = EBV*k_4958
    A_5007 = EBV*k_5007
    #print "A_3727:", A_3727

    out_ascii = fitspath+'/dust_attentuation_values.tbl'
    #if not exists(out_ascii_single):
    n2= ('ID','k_3727', 'k_HDELTA', 'k_Hgamma', 'k_HBETA', 'k_4363', 'k_4958', 'k_5007', 'E(B-V)')
    tab1 = Table([ID,k_3727, k_HDELTA, k_Hgamma, k_HBETA , k_4363, k_4958, k_5007, EBV], names=n2)
    asc.write(tab1, out_ascii, format='fixed_width_two_line')
    
    
def call_cardelli(lam0): #, extrapolate=False):
    #lambda0 =[3726.16, 3868.74, 3888.65, 3967.51, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]* u.angstrom
    line_name = ['OII_3727','NeIII','HeI','3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']
    lambda0 = lam0*u.angstrom
    k_values= cardelli(lambda0,R=3.1)
    return k_values
'''    
    