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
    
#make as inputs eventually  
bin_info_file = 'mass_bin_hbeta_revised_75_112_113_300_600_1444_1444.npz'
derived_prop_file = 'hbeta_revised_75_112_113_300_600_1444_1444_updated_massbin_derived_properties_metallicity.tbl'
mass_LHbeta_file = 'revised_75_112_113_300_600_1444_1444_mass_SFR_data.npz'
valid_file = 'hbeta_revised_75_112_113_300_600_1444_1444_massbin_validation.tbl'


def create_Te_lineratio_table():
    bin_info = np.load(path2 + bin_info_file)
    bin_Te = asc.read(path2 + derived_prop_file)      #from the stacks
    mass_LHbeta = np.load(path2 + mass_LHbeta_file)
    valid_table = asc.read(path2 + valid_file)
    
    bin_ind = bin_info['bin_ind']                   #valid mass indices relative to unsorted data table
    Te = bin_Te['Temperature'].data
    mass = mass_LHbeta['mass']
    LHbeta = mass_LHbeta['lum']
    detections = valid_table['Detection'].data
    
    hdu = fits.open(path3 + 'DEEP2_all_line_fit.fits')
    line_table = hdu[1].data                  
    
    valid_idx = []
    for i in range(len(bin_ind)):                   #bin_ind is 2D array of indices (an array for each bin)
        valid_idx += list(bin_ind[i])               #creates 1D array of all indices
    idx_array = np.array(valid_idx)
    
    EBV_array = np.zeros(len(idx_array))
    indiv_detect_array = np.zeros(len(idx_array))
    mass_array = np.log10(mass[valid_idx])
    LHbeta_array = LHbeta[valid_idx]        
    
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
    Te_array = []
    detect_array = []
    kk = 1
    for ii in range(len(bin_ind)):
        for jj in range(len(bin_ind[ii])):
            mass_bin_ID.append(kk)
            LHbeta_bin_ID.append(ii + 1)
            Te_array.append(Te[ii])
            detect_array.append(detections[ii])
        if (ii + 1) % 2 == 0:
            kk += 1
    
    out_ascii = path2 + 'indivgals_Te_lineRatios.tbl'
    n = ('OBJNO', 'Mass_Bin_ID', 'LHBeta_Bin_ID', 'Log10(Mass)', 'HBeta_Luminosity', 'Te', 'OII_Flux', 'OII_SN',
         'OIII4959_Flux', 'OIII4959_SN', 'OIII5007_Flux', 'OIII5007_SN', 'HBETA_Flux', 'HBETA_SN', 'Bin Detections',
         'Individual Detections', 'EBV')
    Te_ratio_table = Table([objno, mass_bin_ID, LHbeta_bin_ID, mass_array, LHbeta_array, Te_array, OII_Flux, OII_SN,
                            OIII4959_Flux, OIII4959_SN, OIII5007_Flux, OIII5007_SN, HBETA_Flux, HBETA_SN, detect_array,
                            indiv_detect_array, EBV_array], names = n)
    asc.write(Te_ratio_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    
    
    hdu.close()
    
    
    
    
    
    
    
    
    
    
    
    
    #'_updated_massbin_emission_lines.tbl' contains whether detection or not, bin sizes,
    #and average measurements for each bin
    
    #Make Te a 2D array?
    
    #Call R calculation function here (on bins with detections) to calculate R using average measurements
    #then call Te calculation to get Te and put value in inner array for each bin.
    
    #If non-detection, then set to some value
    
    #Call metallicity calculation, then create array of metallicity and add to table?
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    