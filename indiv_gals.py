#Creates table including arrays for T_e and line ratios for individual galaxies.
import numpy as np
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import Table, Column
from getpass import getuser


if getuser() == 'carol':
    path = "C:\\Users\\carol\\Google Drive\\"
    path2 = path + "MZEvolve\\massbin\\"
    path3 = path + "Zcalbase_gal\\dataset\\"
else:
    path = "../DEEP2/" 
    path2 = "../"
    
    
HB_bin_type = 'hbeta_revised_75_112_113_300_600_1444_1444'
mass_bin_type = 'revised_75_112_113_300_600_1444_1444'
    
###INPUT FILES 
#contains bin indices
HB_bin_file = 'mass_bin_hbeta_revised_75_112_113_300_600_1444_1444.npz'

#contains electron temperatures
HB_derived_prop_file = 'hbeta_revised_75_112_113_300_600_1444_1444_updated_massbin_derived_properties_metallicity.tbl'
mass_derived_prop_file = 'revised_75_112_113_300_600_1444_1444_updated_massbin_derived_properties_metallicity.tbl'

#contains mass and luminosity values
mass_LHbeta_file = 'revised_75_112_113_300_600_1444_1444_mass_SFR_data.npz'

#contains combined visual detections and S/N > 3 detections
HB_valid_file = 'hbeta_revised_75_112_113_300_600_1444_1444_massbin_validation.tbl'
mass_valid_file = 'revised_75_112_113_300_600_1444_1444_massbin_validation.tbl'




def create_Te_lineratio_table():
    HB_bin_info = np.load(path2 + HB_bin_file)
    HB_bin_Te = asc.read(path2 + HB_derived_prop_file)
    mass_bin_Te = asc.read(path2 + mass_derived_prop_file)
    mass_LHbeta = np.load(path2 + mass_LHbeta_file)
    HB_valid_table = asc.read(path2 + HB_valid_file)
    mass_valid_table = asc.read(path2 + mass_valid_file)
    
    
    HB_bin_ind = HB_bin_info['bin_ind']                   #valid HBeta bin indices relative to unsorted data table
    HB_Te = HB_bin_Te['Temperature'].data                 #HBeta bin electron temperatures
    mass_Te = mass_bin_Te['Temperature'].data             #Mass bin electron temperatures
    mass = mass_LHbeta['mass']
    LHbeta = mass_LHbeta['lum']
    HB_detections = HB_valid_table['Detection'].data      #HBeta bin detections
    mass_detections = mass_valid_table['Detection'].data  #Mass bin detections
    
    hdu = fits.open(path3 + 'DEEP2_all_line_fit.fits')
    line_table = hdu[1].data                  
    
    
    #HBeta bin valid indices
    HB_valid_idx = []
    for ii in range(len(HB_bin_ind)):                      #creates 1D array of all indices
        HB_valid_idx += list(HB_bin_ind[ii])               
    idx_array = np.array(HB_valid_idx)                     #len = 4088
    
    EBV_array = np.zeros(len(idx_array))             
    indiv_detect_array = np.zeros(len(idx_array))        
    mass_array = np.log10(mass[idx_array])
    LHbeta_array = LHbeta[idx_array]        
    
    line_table = line_table[idx_array]              
    objno = line_table['OBJNO']
    OII_Flux = line_table['OII_FLUX_DATA']
    OII_SN = line_table['OII_SNR']
    OIII4959_Flux = line_table['OIIIB_FLUX_DATA']
    OIII4959_SN = line_table['OIIIB_SNR']
    OIII5007_Flux = line_table['OIIIR_FLUX_DATA']
    OIII5007_SN = line_table['OIIIR_SNR']
    HBETA_Flux = line_table['HB_FLUX_DATA']
    HBETA_SN = line_table['HB_SNR'] 
    
    #only works if mass bins are split into two SFR bins --> fix library.py and adapt for more general case
    mass_bin_ID = []
    LHbeta_bin_ID = []
    HB_Te_array = []
    mass_Te_array = []
    HB_detect_array = []
    mass_detect_array = []
    kk = 1    
    for ii in range(len(HB_bin_ind)):
        for jj in range(len(HB_bin_ind[ii])):
            mass_bin_ID.append(kk)
            LHbeta_bin_ID.append(ii + 1)
            HB_Te_array.append(HB_Te[ii])
            mass_Te_array.append(mass_Te[kk - 1])
        
            #Excluding the last bin --> S/N > 3, but possible AGN contribution
            if HB_bin_type == 'hbeta_revised_75_112_113_300_600_1444_1444':
                if (ii + 1 == 14):
                    HB_detect_array.append(0.0)
                else:
                    HB_detect_array.append(HB_detections[ii])    
                
            if mass_bin_type == 'revised_75_112_113_300_600_1444_1444':
                mass_detect_array.append(mass_detections[kk - 1])
            
        if (ii + 1) % 2 == 0:
            kk += 1
    
    out_ascii = path2 + 'indivgals_Te_lineRatios.tbl'
    n = ('OBJNO', 'Mass_Bin_ID', 'LHBeta_Bin_ID', 'Log10(Mass)', 'HBeta_Luminosity', 'Mass Bin Te', 'LHBeta Bin Te',
         'OII_Flux', 'OII_SN', 'OIII4959_Flux', 'OIII4959_SN', 'OIII5007_Flux', 'OIII5007_SN', 'HBETA_Flux',
         'HBETA_SN', 'Mass Bin Detections', 'LHBeta Bin Detections', 'Individual Detections', 'EBV')
    Te_ratio_table = Table([objno, mass_bin_ID, LHbeta_bin_ID, mass_array, LHbeta_array, mass_Te_array, HB_Te_array,
                            OII_Flux, OII_SN, OIII4959_Flux, OIII4959_SN, OIII5007_Flux, OIII5007_SN, HBETA_Flux,
                            HBETA_SN, mass_detect_array, HB_detect_array, indiv_detect_array, EBV_array], names = n)
    asc.write(Te_ratio_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    
    
    hdu.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    