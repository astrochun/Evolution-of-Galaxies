from astropy.io import ascii as asc
import numpy as np
from astropy.table import Table, Column
from getpass import getuser


if getuser() == 'carol':
    path = "C:\\Users\\carol\\Google Drive\\MZEvolve\\"
    path2 = path
else:
    path = "../DEEP2/" 
    path2 = "../"


N_in_bin = '800'

def invalid_limit():
    em_table = asc.read(path + N_in_bin + '_massbin_emission_lines.tbl')
    valid_table = asc.read(path + N_in_bin + '_massbin_validation.tbl')
    
    O_4363_flux = em_table['OIII_4363_Flux_Observed'].data
    O_4363_SN = em_table['OIII_4363_S/N'].data
    HGamma_rms = em_table['HGAMMA_RMS'].data
    valid_O_4363 = valid_table['Valid_OIII_4363'].data
    
    invalid_stacks_idx = np.where((O_4363_SN < 3) | (valid_O_4363 == 0))[0]  
    detections = np.ones(len(valid_O_4363))
    detections[invalid_stacks_idx] = 0
    
    updated_O_4363_flux = O_4363_flux
    updated_O_4363_flux[invalid_stacks_idx] = 2 * HGamma_rms[invalid_stacks_idx]
    
    out_ascii = path + N_in_bin + '_updated_massbin_emission_lines.tbl'
    updated_em_table = Table(em_table)
    updated_em_table['OIII_4363_Flux_Observed'] = updated_O_4363_flux
    detection_col = Column(detections, name = 'Detection')
    updated_em_table.add_column(detection_col)
    asc.write(updated_em_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    


def make_validation_table():
    massbin_table = asc.read(path + N_in_bin + '_massbin.tbl', format = 'fixed_width_two_line')
    ID = massbin_table['ID'].data
    mass_min = massbin_table['mass_min'].data
    mass_max = massbin_table['mass_max'].data
    mass_avg = massbin_table['mass_avg'].data
    N = massbin_table['Number of Galaxies'].data
    valid_OIII_4363 = np.ones(len(ID))
    
    if N_in_bin == '800':
        valid_OIII_4363[1] = 0
        valid_OIII_4363[4] = 0
    
    out_ascii = path + N_in_bin + '_massbin_validation.tbl'
    n = ('ID', 'mass_min', 'mass_max', 'mass_avg', 'Number of Galaxies','Valid_OIII_4363')
    valid_table = Table([ID, mass_min, mass_max, mass_avg, N, valid_OIII_4363], names = n)
    asc.write(valid_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    
