#Creates table including arrays for T_e, line fluxes, and ratios for individual galaxies.
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
    
    

def create_Te_lineratio_table():
    bin_info = np.load(path2 + 'mass_bin_hbeta_revised_75_112_113_300_600_1444_1444.npz')
    bin_ind = bin_info['bin_ind']                   #valid indices relative to unsorted data table
    
    hdu = fits.open(path3 + 'DEEP2_all_line_fit.fits')
    line_table = hdu[1].data               
    
    valid_idx = []
    for i in range(len(bin_ind)):                   #bin_ind is 2D array of indices (an array for each bin)
        valid_idx += list(bin_ind[i])               #creates 1D array of all indices
    idx_array = np.array(valid_idx)        
    
    line_table = line_table[idx_array]              #get only valid indices
    Te_array = np.zeros(len(idx_array))
    
    out_ascii = path2 + 'indivgals_Te_lineRatios.tbl'
    Te_ratio_table = Table(line_table)
    Te_col = Column(Te_array, name = 'T_e')
    Te_ratio_table.add_column(Te_col, index = 0)
    asc.write(Te_ratio_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    
    hdu.close()