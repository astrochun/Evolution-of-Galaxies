import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from pylab import subplots_adjust
from getpass import getuser

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

N_in_bin = bin_pts_fname

hbeta_bin = True


    
def main():
    table0 = asc.read(fitspath2 + N_in_bin + '_massbin_emission_lines.tbl', format = 'fixed_width_two_line') 
    HG_flux = table0['HGAMMA_Flux_Observed'].data 
    HB_flux = table0['HBETA_Flux_Observed'].data
    avg_mass = table0['mass_avg'].data
    HG_to_HB = HG_flux / HB_flux        
            
    idx_low = np.arange(0, len(avg_mass), 2)
    idx_high = np.arange(1, len(avg_mass), 2)
    
    plt.scatter(avg_mass[idx_low], HG_to_HB[idx_low], color = 'blue')
    plt.scatter(avg_mass[idx_high], HG_to_HB[idx_high], color = 'orange')
    plt.axhline(y = 0.468)
    