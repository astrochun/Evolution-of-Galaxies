from astropy.io import ascii as asc
from astropy.io import fits
import numpy as np
from astropy.table import Table, Column
from getpass import getuser
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


if getuser() == 'carol':
    path = "C:\\Users\\carol\\Google Drive\\MZEvolve\\"
    path2 = path + "massbin\\"
else:
    path = "../DEEP2/" 
    path2 = "../"

bin_pts_input = [75, 112, 113, 300, 600, 1444, 1444]
str_bin_pts_input = [str(val) for val in bin_pts_input]
bin_pts_fname = "_".join(str_bin_pts_input)
bin_pts_fname = 'hbeta_revised_' + bin_pts_fname

N_in_bin = bin_pts_fname

mark_nondet = True
if mark_nondet:
    updated = '_updated'
else:
    updated = ''


#revisit
'''def set_luminosity():
    #binning_file = np.load(path2 + 'mass_bin_' + bin_pts_fname + '.npz')
    #lum = binning_file['lum']
    em_table = asc.read(path2 + N_in_bin + '_massbin_emission_lines.tbl')
    O_5007_rms = em_table['OIII_5007_RMS'].data
    
    invalid_lum_idx = np.where(lum == -1)[0]
    updated_lum = lum
    updated_lum[invalid_lum_idx] = 3 * O_5007_rms[invalid_lum_idx]'''
    


def invalid_limit():
    em_table = asc.read(path2 + N_in_bin + '_massbin_emission_lines.tbl')
    valid_table = asc.read(path2 + N_in_bin + '_massbin_validation.tbl')
    
    O_4363_flux = em_table['OIII_4363_Flux_Observed'].data
    O_4363_SN = em_table['OIII_4363_S/N'].data
    HGamma_rms = em_table['HGAMMA_RMS'].data
    valid_O_4363 = valid_table['Valid_OIII_4363'].data
    
    invalid_stacks_idx = np.where((O_4363_SN < 3) | (valid_O_4363 == 0))[0]  
    detections = np.ones(len(valid_O_4363))
    detections[invalid_stacks_idx] = 0
    
    updated_O_4363_flux = O_4363_flux
    updated_O_4363_flux[invalid_stacks_idx] = 3 * HGamma_rms[invalid_stacks_idx]
    
    out_ascii = path2 + N_in_bin + '_updated_massbin_emission_lines.tbl'
    updated_em_table = Table(em_table)
    updated_em_table['OIII_4363_Flux_Observed'] = updated_O_4363_flux
    detection_col = Column(detections, name = 'Detection')
    updated_em_table.add_column(detection_col)
    asc.write(updated_em_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    


def make_validation_table():
    massbin_table = asc.read(path2 + N_in_bin + '_massbin.tbl', format = 'fixed_width_two_line')
    ID = massbin_table['ID'].data
    mass_min = massbin_table['mass_min'].data
    mass_max = massbin_table['mass_max'].data
    mass_avg = massbin_table['mass_avg'].data
    N = massbin_table['Number of Galaxies'].data
    valid_OIII_4363 = np.ones(len(ID))
    
    if N_in_bin == '800':
        valid_OIII_4363[1] = 0
        valid_OIII_4363[4] = 0
        
    if N_in_bin == '75_112_113_300_600_1444_1444':
        valid_OIII_4363[1] = 0
        valid_OIII_4363[5] = 0
        valid_OIII_4363[6] = 0   
        
    if N_in_bin == 'revised_75_112_113_300_600_1444_1444':
        valid_OIII_4363[3] = 0
        valid_OIII_4363[5] = 0
        valid_OIII_4363[6] = 0
        
    if N_in_bin == 'hbeta_revised_75_112_113_300_600_1444_1444':
        valid_OIII_4363[2] = 0
        valid_OIII_4363[8] = 0
        valid_OIII_4363[10] = 0
        valid_OIII_4363[12] = 0       
        
    
    out_ascii = path2 + N_in_bin + '_massbin_validation.tbl'
    n = ('ID', 'mass_min', 'mass_max', 'mass_avg', 'Number of Galaxies','Valid_OIII_4363')
    valid_table = Table([ID, mass_min, mass_max, mass_avg, N, valid_OIII_4363], names = n)
    asc.write(valid_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    
