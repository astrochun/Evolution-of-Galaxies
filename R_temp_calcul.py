#Calculates the R value, electron temperature, and metallicity from the flux table
#produced by the zoom_and_gauss_general functions
#Currently running: Grid
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from scipy.optimize import curve_fit 
from getpass import getuser
    

#For generalizing for several users
if getuser() == 'carol':
    fitspath = "C:\\Users\\carol\\Google Drive\\MZEvolve\\"
    fitspath2 = fitspath + "massbin\\"
else:
    fitspath = "../DEEP2/" 
    fitspath2 = "../"
    
bin_pts_input = [75, 112, 113, 300, 600, 1444, 1444]
str_bin_pts_input = [str(val) for val in bin_pts_input]
bin_pts_fname = "_".join(str_bin_pts_input)
bin_pts_fname = 'hbeta_revised_' + bin_pts_fname
bin_pts_fname2 = 'individual'
 
#don't need --> just replace with bin_pts_fname   
#N_in_bin = bin_pts_fname

mark_nondet = False
if mark_nondet:
    updated = '_updated'
else:
    updated = ''


#Constants
a = 13205
b = 0.92506
c = 0.98062
 
    

def R_calculation(OIII4363, OIII5007):
    R_value = OIII4363 / (OIII5007 * (4/3))
    return R_value


def temp_calculation(R):
    T_e = a*(-np.log10(R) - b)**(-1*c)      
    return T_e


def metallicity_calculation(T_e, two_beta, three_beta): 
    O_s_ion = np.zeros(len(T_e))
    O_d_ion = np.zeros(len(T_e))
    com_O = np.zeros(len(T_e))
    com_O_log = np.zeros(len(T_e))
    O_s_ion_log = np.zeros(len(T_e))
    O_d_ion_log = np.zeros(len(T_e))
    t_3 = np.zeros(len(T_e))
    t_2 = np.zeros(len(T_e))
    x2 = np.zeros(len(T_e))
    
    detect = np.where((two_beta != 0) & (three_beta != 0))[0]
    
    t_3[detect] = T_e[detect] * 1e-4
    t_2[detect] = 0.7 * t_3[detect] + 0.17
    x2[detect] = 1e-4 * 1e3 * t_2[detect]**(-0.5)

    O_s_ion_log[detect] = np.log10(two_beta[detect]) + 5.961 + 1.676 / t_2[detect] - 0.4 * np.log10(t_2[detect]) - 0.034 * t_2[detect] + np.log10(1 + 1.35 * x2[detect]) - 12
    O_d_ion_log[detect] = np.log10(three_beta[detect]) + 6.200 + 1.251 / t_3[detect] - 0.55 * np.log10(t_3[detect]) - 0.014 * (t_3[detect]) - 12

    O_s_ion[detect] = 10**(O_s_ion_log[detect])
    O_d_ion[detect] = 10**(O_d_ion_log[detect])
    com_O[detect] = O_s_ion[detect] + O_d_ion[detect]
    com_O_log[detect] = np.log10(com_O[detect]) + 12

    return O_s_ion, O_d_ion, com_O_log, O_s_ion_log, O_d_ion_log



def run_function():
       
    line_file = fitspath2 + 'indivgals_Te_lineRatios.tbl'
    #line_file = fitspath2 + 'hbeta_revised_75_112_113_300_600_1444_1444_updated_massbin_emission_lines.tbl'

    line_table = asc.read(line_file)
    
    if 'two_beta' in line_table.keys():
        #Case for individual spectra 
        
        out_ascii = fitspath2 + bin_pts_fname2 + '_derived_properties_metallicity.tbl'
        out_fits = fitspath2 + bin_pts_fname2 + '_derived_properties_metallicity.fits'
        
        OIII4959 = line_table['OIII4959'].data
        OIII5007 = line_table['OIII5007'].data
        HBETA = line_table['HBeta'].data
        HGAMMA = line_table['HGamma'].data
        SNR_HG = line_table['SNR_HG'].data
        raw_OIII4363 = line_table['OIII4363'].data
        SNR_4363 = line_table['SNR_4363'].data
        R23_individual = line_table['Individual_R23'].data
        O32_individual = line_table['Individual_O32'].data
        detections = line_table['Detection'].data
        
        two_beta = line_table['two_beta'].data
        three_beta = line_table['three_beta'].data
        T_e = line_table['Temperature'].data
        source_ID = line_table['Source_ID'].data
        
        
        O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metallicity_calculation(T_e, two_beta, three_beta)
        
        n = ('Source_ID', 'R23', 'O32', 'Observed_Flux_5007', 'Observed_Flux_4959',
             'Observed_Flux_HBeta', 'Temperature', 'Detections', 'O_s_ion', 'O_d_ion', 'com_O_log')
        tab0 = Table([source_ID, R23_individual, O32_individual, OIII5007, OIII4959, HBETA, T_e, detections,
                      O_s_ion, O_d_ion, com_O_log], names = n)
    
    
    elif 'Log10(Mass)' in line_table.keys():
        #Case for individual spectra
    
        out_ascii = fitspath2 + bin_pts_fname2 + '_derived_properties_metallicity.tbl'
        out_fits = fitspath2 + bin_pts_fname2 + '_derived_properties_metallicity.fits'
        
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
        HB_bin_ID = line_table['LHBeta_Bin_ID'].data
        log_mass = line_table['Log10(Mass)'].data
        LHbeta = line_table['HBeta_Luminosity'].data
        mass_T_e = line_table['Mass Bin Te'].data
        HB_T_e = line_table['LHBeta Bin Te'].data
        HB_bin_detect = line_table['LHBeta Bin Detections'].data
        mass_bin_detect = line_table['Mass Bin Detections'].data
        indiv_detect = line_table['Individual Detections'].data 
        ebv = line_table['EBV'].data
        
        
        mass_detect = np.where((mass_bin_detect == 1.0) & (np.isfinite(OIII5007) == True) & (OIII5007 >= 1e-18) & 
                               (OIII5007 <= 1e-15) & (np.isfinite(OIII4959) == True) & (OIII4959 >= 1e-18) & 
                               (OIII4959 <= 1e-15) & (np.isfinite(OII) == True) & (OII >= 1e-18) & (OII <= 1e-15) &
                               (np.isfinite(HBETA) == True) & (HBETA >= 1e-18) & (HBETA <= 1e-15))[0]
        HB_detect = np.where((HB_bin_detect == 1.0) & (np.isfinite(OIII5007) == True) & (OIII5007 >= 1e-18) & 
                             (OIII5007 <= 1e-15) & (np.isfinite(OIII4959) == True) & (OIII4959 >= 1e-18) & 
                             (OIII4959 <= 1e-15) & (np.isfinite(OII) == True) & (OII >= 1e-18) & (OII <= 1e-15) &
                             (np.isfinite(HBETA) == True) & (HBETA >= 1e-18) & (HBETA <= 1e-15))[0]
        indiv_detect[mass_detect] = 1.0
        indiv_detect[HB_detect] = 1.0
        
        #create zero arrays all same length
        two_beta = np.zeros(len(source_ID))
        three_beta = np.zeros(len(source_ID))
        R23 = np.zeros(len(source_ID))
        O32 = np.zeros(len(source_ID))  
        
        two_beta[mass_detect] = OII[mass_detect] / HBETA[mass_detect]
        two_beta[HB_detect] = OII[HB_detect] / HBETA[HB_detect]
        three_beta[mass_detect] = (OIII5007[mass_detect] + OIII4959[mass_detect]) / HBETA[mass_detect]
        three_beta[HB_detect] = (OIII5007[HB_detect] + OIII4959[HB_detect]) / HBETA[HB_detect]
        
        #Calculate R23 and O32
        R23[mass_detect] = np.log10((OII[mass_detect] + ((4/3) * OIII5007[mass_detect])) / HBETA[mass_detect])
        R23[HB_detect] = np.log10((OII[HB_detect] + ((4/3) * OIII5007[HB_detect])) / HBETA[HB_detect])
        O32[mass_detect] = np.log10(((4/3) * OIII5007[mass_detect]) / OII[mass_detect])
        O32[HB_detect] = np.log10(((4/3) * OIII5007[HB_detect]) / OII[HB_detect])
        
        
        mass_O_s_ion, mass_O_d_ion, mass_com_O_log, mass_log_O_s, mass_log_O_d = metallicity_calculation(mass_T_e, two_beta, three_beta)
        HB_O_s_ion, HB_O_d_ion, HB_com_O_log, HB_log_O_s, HB_log_O_d = metallicity_calculation(HB_T_e, two_beta, three_beta)
        
        
        n = ('Source_ID', 'Mass Bin ID', 'LHBeta Bin ID', 'Mass Bin Detections', 'LHBeta Bin Detections',
             'Individual Detections', 'Log10(Mass)', 'HBeta_Luminosity', 'Observed_Flux_5007',
             'Observed_Flux_4959', 'Observed_Flux_3727', 'Observed_Flux_HBeta', 'Mass Bin Te', 'LHBeta Bin Te',
             'R23', 'O32', 'Mass Bin O_s_ion', 'Mass Bin O_d_ion', 'Mass Bin com_O_log', 'LHBeta Bin O_s_ion',
             'LHBeta Bin O_d_ion', 'LHBeta Bin com_O_log', 'EBV')
        tab0 = Table([source_ID, mass_bin_ID, HB_bin_ID, mass_bin_detect, HB_bin_detect, indiv_detect, log_mass,
                      LHbeta, OIII5007, OIII4959, OII, HBETA, mass_T_e, HB_T_e, R23, O32, mass_O_s_ion,
                      mass_O_d_ion, mass_com_O_log, HB_O_s_ion, HB_O_d_ion, HB_com_O_log, ebv], names = n)
        
    else:
        #Case for stacked spectra
        
        out_ascii = fitspath2 + bin_pts_fname + updated + '_derived_properties_metallicity.tbl'
        out_fits = fitspath2 + bin_pts_fname + updated + '_derived_properties_metallicity.fits'
        
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
        log_mass = avg_mass
        
        two_beta = OII / HBETA
        three_beta = (OIII5007 + OIII4959) / HBETA
        
        #Calculate R23 composite and O32 composite
        R23_composite = np.log10((OII + ((4/3) * OIII5007)) / HBETA)
        O32_composite = np.log10(((4/3) * OIII5007) / OII)
        
        #R, Te, and metallicity calculations
        R_value = R_calculation(OIII4363, OIII5007)
        T_e = temp_calculation(R_value)
        O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metallicity_calculation(T_e, two_beta, three_beta)
        
        
        n = ('R23_Composite', 'O32_Composite', 'N_Galaxies', 'Observed_Flux_5007', 'S/N_5007', 'Observed_Flux_4959',
             'S/N_4959', 'Observed_Flux_4363', 'S/N_4363', 'Observed_Flux_HBETA', 'S/N_HBETA', 'Observed_Flux_3727',
             'S/N_3727', 'Temperature', 'O_s_ion', 'O_d_ion', 'com_O_log')
        tab0 = Table([R23_composite, O32_composite, N_Galaxy, OIII5007, SN_5007, OIII4959, SN_4959,
                      OIII4363, SN_4363, HBETA, SN_HBETA, OII, SN_OII, T_e, O_s_ion, O_d_ion, com_O_log], names = n)

        

    
    asc.write(tab0, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    tab0.write(out_fits, format = 'fits', overwrite = True)
    