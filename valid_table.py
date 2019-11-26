from astropy.io import ascii as asc
import numpy as np
from astropy.table import Table, Column
 


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
    
    
    ##Update OIII4363 values for non-detections
    #1.0 --> detection, 0.0 --> non-detection, 0.5 --> non-detection with reliable 5007
    detections = np.zeros(len(ID))
    valid_stacks_idx = np.where((O_4363_SN >= 3) & (O_5007_SN > 100))[0] 
    reliable_5007_stacks = np.where((O_4363_SN < 3) & (O_5007_SN > 100))[0]
    detections[valid_stacks_idx] = 1
    detections[reliable_5007_stacks] = 0.5
    QA_flag = np.zeros(len(O_4363_flux))
    QA_flag = quality_assurance(bin_pts_fname, QA_flag)
    updated_O_4363_flux = np.zeros(len(O_4363_flux))
    for i in range(len(O_4363_flux)):
        if detections[i] == 0.0 or detections[i] == 0.5 or QA_flag[i] == 1.0:
            updated_O_4363_flux[i] = 3 * HGamma_rms[i]
        elif detections[i] == 1.0:
            updated_O_4363_flux[i] = O_4363_flux[i] 
    
    
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
    
    
    
    
def quality_assurance(bin_pts_fname, QA_flag):
    #Limits 4363 line width
    if bin_pts_fname == 'massLHbetabin_revised_75_112_113_300_600_1444_1444':
        QA_flag[10] = 1.0
        QA_flag[11] = 1.0
    elif bin_pts_fname == 'massbin_revised_75_112_113_300_600_1444_1444':   
        QA_flag[5] = 1.0
        
    return QA_flag
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
