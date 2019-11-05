from astropy.io import ascii as asc
import numpy as np
from astropy.table import Table, Column
from getpass import getuser

'''
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
#bin_pts_fname = 'revised_75_112_113_300_600_1444_1444'

N_in_bin = bin_pts_fname

mark_nondet = True
if mark_nondet:
    updated = '_updated'
else:
    updated = ''
 '''
   

'''
def invalid_limit(fitspath, bin_pts_fname):
    em_table = asc.read(fitspath + bin_pts_fname + '_massbin_emission_lines.tbl')
    valid_table = asc.read(fitspath + bin_pts_fname + '_massbin_validation.tbl')
    
    O_4363_flux = em_table['OIII_4363_Flux_Observed'].data
    O_4363_SN = em_table['OIII_4363_S/N'].data
    O_5007_SN = em_table['OIII_5007_S/N'].data
    HGamma_rms = em_table['HGAMMA_RMS'].data
    valid_O_4363 = valid_table['Valid_OIII_4363'].data
    
    #correct?
    invalid_stacks_idx = np.where((O_4363_SN < 3) | (valid_O_4363 == 0))[0] 
    #implement 5007 cut?
    detections = np.ones(len(valid_O_4363))
    detections[invalid_stacks_idx] = 0
    
    updated_O_4363_flux = O_4363_flux
    updated_O_4363_flux[invalid_stacks_idx] = 3 * HGamma_rms[invalid_stacks_idx]
    
    out_ascii2 = path2 + N_in_bin + '_massbin_validation.tbl'
    updated_valid_table = Table(valid_table)
    valid_detect_col = Column(detections, name = 'Detection')
    updated_valid_table.add_column(valid_detect_col)
    asc.write(updated_valid_table, out_ascii2, format = 'fixed_width_two_line', overwrite = True)
    
    out_ascii = path2 + N_in_bin + '_updated_massbin_emission_lines.tbl'
    updated_em_table = Table(em_table)
    updated_em_table['OIII_4363_Flux_Observed'] = updated_O_4363_flux
    detection_col = Column(detections, name = 'Detection')
    updated_em_table.add_column(detection_col)
    asc.write(updated_em_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
 '''

'''
def temp_valid_table():
    path = 'C:\\Users\\carol\\Google Drive\\MZEvolve\\massbin\\10252019\\TESTrevised_75_112_113_300_600_1444_1444\\'
    bin_pts_fname = 'massbin_revised_75_112_113_300_600_1444_1444'
    massbin_table = asc.read(path + bin_pts_fname + '_binning.tbl', format = 'fixed_width_two_line')
    em_table = asc.read(path + bin_pts_fname + '_emission_lines.tbl')
    ID = massbin_table['ID'].data
    mass_min = massbin_table['mass_min'].data
    mass_max = massbin_table['mass_max'].data
    mass_avg = massbin_table['mass_avg'].data
    N = massbin_table['Number of Galaxies'].data
    
    O_4363_flux = em_table['OIII_4363_Flux_Observed'].data
    O_4363_SN = em_table['OIII_4363_S/N'].data
    O_5007_SN = em_table['OIII_5007_S/N'].data
    HGamma_rms = em_table['HGAMMA_RMS'].data
    
    detections = np.zeros(len(ID))
    detections[0] = 1
    detections[2] = 1
    detections[4] = 1
    
    invalid_stacks_idx = np.where(detections != 1)[0]
    
    updated_O_4363_flux = O_4363_flux
    updated_O_4363_flux[invalid_stacks_idx] = 3 * HGamma_rms[invalid_stacks_idx]
    
    out_ascii = path + bin_pts_fname + '_updated_emission_lines.tbl'
    updated_em_table = Table(em_table)
    updated_em_table['OIII_4363_Flux_Observed'] = updated_O_4363_flux
    detection_col = Column(detections, name = 'Detection')
    updated_em_table.add_column(detection_col)
    asc.write(updated_em_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
        
    
    out_ascii = path + bin_pts_fname + '_massbin_validation.tbl'
    n = ('ID', 'mass_min', 'mass_max', 'mass_avg', 'Number of Galaxies', 'Detection')
    valid_table = Table([ID, mass_min, mass_max, mass_avg, N, detections], names = n)
    asc.write(valid_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
  '''  


def make_validation_table(fitspath, bin_pts_fname):
    massbin_table = asc.read(fitspath + bin_pts_fname + '_binning.tbl', format = 'fixed_width_two_line')
    em_table = asc.read(fitspath + bin_pts_fname + '_emission_lines.tbl')
    
    ID = massbin_table['ID'].data
    mass_min = massbin_table['mass_min'].data
    mass_max = massbin_table['mass_max'].data
    mass_avg = massbin_table['mass_avg'].data
    N = massbin_table['Number of Galaxies'].data
    
    O_4363_flux = em_table['OIII_4363_Flux_Observed'].data
    O_4363_SN = em_table['OIII_4363_S/N'].data
    O_5007_SN = em_table['OIII_5007_S/N'].data
    HGamma_rms = em_table['HGAMMA_RMS'].data
    
    
    #Update OIII4363 values for non-detections
    detections = np.zeros(len(ID))            #one --> detection, zero --> non-detection
    valid_stacks_idx = np.where((O_4363_SN >= 3) & (O_5007_SN > 100))[0] 
    detections[valid_stacks_idx] = 1
    invalid_stacks_idx = np.where(detections != 1)[0]
    #updated_O_4363_flux = O_4363_flux
    updated_O_4363_flux = np.zeros(len(O_4363_flux))
    for i in range(len(O_4363_flux)):
        if detections[i] == 0:
            updated_O_4363_flux[i] = 3 * HGamma_rms[i]
        else:
            updated_O_4363_flux[i] = O_4363_flux[i] 
    #updated_O_4363_flux[invalid_stacks_idx] = 3 * HGamma_rms[invalid_stacks_idx]
    
    
    #Add updated OIII4363 flux column and detection column to emission_line.tbl
    out_ascii_em_table = fitspath + bin_pts_fname + '_emission_lines.tbl'
    updated_O_4363_col = Column(updated_O_4363_flux, name = 'Updated_OIII_4363_Flux_Observed')
    detection_col = Column(detections, name = 'Detection')
    em_table.add_column(updated_O_4363_col, index=30) 
    em_table.add_column(detection_col, index=7)
    asc.write(em_table, out_ascii_em_table, format = 'fixed_width_two_line', overwrite = True)
        
    #Create validation table
    out_ascii_valid_table = fitspath + bin_pts_fname + '_validation.tbl'
    n = ('ID', 'mass_min', 'mass_max', 'mass_avg', 'Number of Galaxies', 'Detection')
    valid_table = Table([ID, mass_min, mass_max, mass_avg, N, detections], names = n)
    asc.write(valid_table, out_ascii_valid_table, format = 'fixed_width_two_line', overwrite = True)
    
