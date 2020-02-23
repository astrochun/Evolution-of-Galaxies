from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt

from getpass import getuser
if getuser() == 'carol':
    path = "C:\\Users\\carol\\Google Drive\\MZEvolve\\"
    path2 = path + 'massbin\\'
else:
    path = "../DEEP2/" 
    path2 = "../"

result = ascii.read(path2 + 'results_deeps_revised.tbl')
mass = result['best.stellar.m_star'].data
deep_fields = ascii.read(path + 'magfiles/deep_fields.mag')

phot_cols = deep_fields.colnames[2:5]
phot_cols = [str0 for str0 in phot_cols if '_err' not in str0]

no_mass_idx = np.where((mass <= 0) | (np.isfinite(mass) == False))[0]
mag_no_mass = -2.5*np.log10(deep_fields['cfht_I'][no_mass_idx].data)+16.4

def binning(mass, deep_fields, phot_cols):
    fig, ax = plt.subplots(3, 1, figsize = (15, 15))
    plt.subplots_adjust(hspace = 0.7)
    valid_mass_idx = np.where((mass < 1e14) & (mass > 0))[0]
    valid_mass = np.log10(mass[valid_mass_idx])
    
    for ii in range(len(phot_cols)):
        mag = -2.5*np.log10(deep_fields[phot_cols[ii]].data)+16.4
        valid_mag = mag[valid_mass_idx]
        ax[ii].scatter(valid_mag, valid_mass, s=10, alpha=0.5)
        ax[ii].annotate(phot_cols[ii], [0.95, 0.95], xycoords = "axes fraction", ha = "right", va = "top")
        #galaxy_number = np.where(np.isfinite(mag) == True)[0] 
        
        #binning 
        grid = np.arange(np.nanmin(valid_mag), np.nanmax(valid_mag), 0.25)
        avg_mass = np.zeros(len(grid))
        standard_dev = np.zeros(len(grid))
        N_gals = np.zeros(len(grid))
        for jj in range(len(grid)):
            in_bin = np.where((valid_mag >= grid[jj]) & (valid_mag < grid[jj] + 0.25))[0]
            if len(in_bin) > 0:
                avg_mass[jj] = np.nanmean(valid_mass[in_bin])
                standard_dev[jj] = np.std(valid_mass[in_bin])
                N_gals[jj] = len(in_bin)
        nonzero = np.where(avg_mass != 0)[0]
        ax[ii].scatter(grid[nonzero] + 0.125, avg_mass[nonzero], s = 50, color = "black")
        ax[ii].errorbar(grid[nonzero] + 0.125, avg_mass[nonzero], yerr = standard_dev[nonzero], color = "black", fmt = "none")
        
        np.savez(path2 + 'revised_mag_mass_' + str(phot_cols[ii]) + '.npz', valid_mag = valid_mag,
                 valid_mass = valid_mass, grid = grid, average_mass = avg_mass,
                 standard_deviation = standard_dev, N_gals = N_gals)


def interpolation(filename, band, mag_no_mass, no_mass_idx):
    npz_file = np.load(filename)
    grid, avg_mass, valid_mag, valid_mass, standard_dev, N_gals = npz_file['grid'], npz_file['average_mass'], npz_file['valid_mag'], npz_file['valid_mass'], npz_file['standard_deviation'], npz_file['N_gals']
    gals_in_bin = np.where(N_gals > 0)[0]
    interp_data = interp1d(grid[gals_in_bin], avg_mass[gals_in_bin], fill_value="extrapolate")
    nonzero = np.where(avg_mass != 0)[0]
    
    plt.figure(figsize = (15, 5))
    plt.scatter(valid_mag, valid_mass, s = 10, alpha = 0.5)
    plt.annotate(band, [0.95, 0.95], xycoords = "axes fraction", ha = "right", va = "top")
    plt.scatter(grid[nonzero] + 0.125, avg_mass[nonzero], s = 50, color = "black")
    plt.errorbar(grid[nonzero] + 0.125, avg_mass[nonzero], yerr = standard_dev[nonzero], color = "black", fmt = "none")
    plt.plot(grid + 0.125, avg_mass, grid + 0.125, interp_data(grid), 'r-')
    
    np.savez(path2 + band + '_band_revised_interp_data.npz', interp_data = np.array([interp_data]),
             mag_no_mass = mag_no_mass, no_mass_idx = no_mass_idx)
    
    return interp_data
    
    
    
#binning(mass, deep_fields, phot_cols)    
#interpolation(path2 + 'mag_mass_cfht_B.npz', phot_cols[0])
#interpolation(path2 + 'mag_mass_cfht_R.npz', phot_cols[1])
interpolation(path2 + 'revised_mag_mass_cfht_I.npz', phot_cols[2], mag_no_mass, no_mass_idx)
