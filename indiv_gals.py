#Creates table including arrays for T_e and line ratios for individual galaxies.
import numpy as np
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import Table, Column
from getpass import getuser

'''
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
'''



def create_Te_lineratio_table(fitspath, line_file, mass_bin_file, HB_bin_file, mass_Te_file, HB_Te_file):
    mass_bin_npz = np.load(mass_bin_file)
    mass_Te_tbl = asc.read(mass_Te_file)
    HB_bin_npz = np.load(HB_bin_file)
    HB_Te_tbl = asc.read(HB_Te_file)
    hdu = fits.open(line_file)
    
    
    line_table = hdu[1].data
    mass_bin_ind = mass_bin_npz['bin_ind']            #valid mass bin indices relative to unsorted data table
    HB_bin_ind = HB_bin_npz['bin_ind']                #valid mass-LHBeta bin indices relative to unsorted data table
    mass_Te = mass_Te_tbl['Temperature'].data         #Mass bin electron temperatures
    HB_Te = HB_Te_tbl['Temperature'].data             #HBeta bin electron temperatures
    mass = mass_bin_npz['mass']
    LHbeta = HB_bin_npz['lum']
    mass_detect = mass_Te_tbl['Detection'].data       #Mass bin detections
    HB_detect = HB_Te_tbl['Detection'].data           #HBeta bin detections            
    
    
    #Get all mass bin valid indices in 1D array
    mass_valid_idx = []
    for ii in range(len(mass_bin_ind)):
        mass_valid_idx += list(mass_bin_ind[ii])               
    idx_array = np.array(mass_valid_idx)                     
    
    
    EBV_array = np.zeros(len(idx_array))             
    indiv_detect_array = np.zeros(len(idx_array))        
    mass_array = np.log10(mass[idx_array])
    LHbeta_array = LHbeta[idx_array]        
    line_table = line_table[idx_array]   


    #Get individual lines fluxes and S/N           
    objno = line_table['OBJNO']
    OII_Flux = line_table['OII_FLUX_DATA']
    OII_SN = line_table['OII_SNR']
    OIII4959_Flux = line_table['OIIIB_FLUX_DATA']
    OIII4959_SN = line_table['OIIIB_SNR']
    OIII5007_Flux = line_table['OIIIR_FLUX_DATA']
    OIII5007_SN = line_table['OIIIR_SNR']
    HBETA_Flux = line_table['HB_FLUX_DATA']
    HBETA_SN = line_table['HB_SNR'] 
    
    
    #only works if mass bins are split into two MLHbeta bins --> fix library.py and adapt for more general case
    mass_bin_ID = []
    HB_bin_ID = []
    mass_Te_array = []
    HB_Te_array = []
    mass_detect_array = []
    HB_detect_array = []
    kk = 1    
    for ii in range(len(HB_bin_ind)):
        for jj in range(len(HB_bin_ind[ii])):
            mass_bin_ID.append(kk)
            HB_bin_ID.append(ii + 1)
            HB_Te_array.append(HB_Te[ii])
            mass_Te_array.append(mass_Te[kk - 1])
            
            HB_detect_array.append(HB_detect[ii]) 
            mass_detect_array.append(mass_detect[kk - 1])
            
        if (ii + 1) % 2 == 0:
            kk += 1
    
    out_ascii = fitspath + 'individual_Te_emLines.tbl'  
    n = ('OBJNO', 'Mass_Bin_ID', 'Mass_LHBeta_Bin_ID', 'Log10(Mass)', 'HBeta_Luminosity', 'Mass_Bin_Detections',
         'Mass_LHBeta_Bin_Detections', 'Individual_Detections', 'E(B-V)', 'Mass_Bin_Te', 'Mass_LHBeta_Bin_Te',
         'OII_Flux', 'OII_SN', 'OIII4959_Flux', 'OIII4959_SN', 'OIII5007_Flux', 'OIII5007_SN', 'HBETA_Flux',
         'HBETA_SN')
    Te_line_table = Table([objno, mass_bin_ID, HB_bin_ID, mass_array, LHbeta_array, mass_detect_array,
                            HB_detect_array, indiv_detect_array, EBV_array, mass_Te_array,
                            HB_Te_array, OII_Flux, OII_SN, OIII4959_Flux, OIII4959_SN, OIII5007_Flux,
                            OIII5007_SN, HBETA_Flux, HBETA_SN], names = n)
    asc.write(Te_line_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    
    
    hdu.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    