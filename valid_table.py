from astropy.io import ascii as asc
import numpy as np
from astropy.table import Table
from Metallicity_Stack_Commons.column_names import filename_dict, bin_names0
 


def make_validation_table(fitspath, bin_type_str):
    '''
    Purpose: 
        This function creates a validation table for a given binning set. The validation table
        contains a OIII4363 detection column where 1.0 means detection, 0.5 means non-detection with
        reliable OIII5007, and 0.0 means unreliable non-detection.
        
    Usage:
        valid_table.make_validation_table(fitspath)
        
    Params:
        fitspath --> a string of the file path where the input file is and where the output file
            will be placed.
        bin_type_str --> a string describing the binning type. (e.g. 'massLHbetabin' or 'massbin') 
        
    Returns: 
        None
        
    Outputs:
        fitspath + 'validation.tbl' --> a validation table containing bin IDs;
            minimum, maximum, and average masses in each bin; number of galaxies in each bin; and 
            a column indicating if the bin has an OIII4363 detection or non-detection.
    '''
    
    bin_table = asc.read(fitspath + filename_dict['bin_info'], format = 'fixed_width_two_line')
    em_table = asc.read(fitspath + filename_dict['bin_fit'])
    
    ID = bin_table[bin_names0[0]].data
    Nstack = bin_table[bin_names0[1]].data
    
    O_4363_flux = em_table['OIII_4363_Flux_Observed'].data
    O_4363_SN = em_table['OIII_4363_S/N'].data
    O_4363_sigma = em_table['OIII_4363_Sigma'].data
    O_5007_SN = em_table['OIII_5007_S/N'].data
    #HGamma_rms = em_table['HGAMMA_RMS'].data  ## temp fix
    
    
    ##Update OIII4363 values for non-detections
    #1.0 --> detection, 0.0 --> non-detection, 0.5 --> non-detection with reliable 5007 or detection with wide line
    detections = np.zeros(len(ID))
    valid_stacks_idx = np.where((O_4363_SN >= 3) & (O_5007_SN > 100) & (O_4363_sigma < 2))[0] 
    reliable_5007_stacks = np.where((O_4363_SN < 3) & (O_5007_SN > 100))[0]
    wide_lines_valid = np.where((O_4363_SN >= 3) & (O_5007_SN > 100) & (O_4363_sigma >= 2))[0]
    detections[valid_stacks_idx] = 1
    detections[reliable_5007_stacks] = 0.5
    detections[wide_lines_valid] = 0.5
    #QA_flag = np.zeros(len(O_4363_flux))
    #QA_flag = quality_assurance(bin_type_str, QA_flag)
    updated_O_4363_flux = np.zeros(len(O_4363_flux))
    for i in range(len(O_4363_flux)):
        #Temp fix to run
        '''
        if detections[i] == 0.0 or detections[i] == 0.5:  # or QA_flag[i] == 1.0:
            updated_O_4363_flux[i] = 3 * HGamma_rms[i]
        elif detections[i] == 1.0:
        '''
        updated_O_4363_flux[i] = O_4363_flux[i] 
    
    '''
    #Add updated OIII4363 flux column and detection column to emission_line.tbl
    out_ascii_em_table = fitspath + 'emission_lines.tbl'
    em_table['OIII_4363_Flux_Observed'].name = 'Original_OIII_4363_Flux_Observed' 
    updated_O_4363_col = Column(updated_O_4363_flux, name = 'OIII_4363_Flux_Observed')
    detection_col = Column(detections, name = 'Detection')
    em_table.add_column(updated_O_4363_col, index=31) 
    em_table.add_column(detection_col, index=7)
    asc.write(em_table, out_ascii_em_table, format = 'fixed_width_two_line', overwrite = True)
    '''
        
    #Create validation table
    out_ascii_valid_table = fitspath + filename_dict['bin_valid'] 
    n = tuple(bin_names0 + ['OIII_4363_Flux_Observed'])
    valid_table = Table([ID, Nstack, detections, O_4363_flux], names = n)
    asc.write(valid_table, out_ascii_valid_table, format = 'fixed_width_two_line', overwrite = True)
    
    
    
    
def quality_assurance(bin_type_str, QA_flag):
    '''
    Purpose:
        This function allows for manual flagging of sources for quality assurance of OIII4363 detections.
        Based on the bin_type_str keyword, the user can override the detection flag by setting a specific
        bin index equal to the desired flag (1.0, 0.5, or 0.0).
        
    Usage:
        valid_table.quality_assurance(bin_type_str, QA_flag)
        
    Params:
        bin_type_str --> a string describing the binning type. (e.g. 'massLHbetabin' or 'massbin')
        QA_flag --> a numpy zeros array the size of the number of bins. This is used to flag sources by
            changing the value at a specific index to the desired flag.
        
    Returns: 
        QA_flag --> the updated flag array.
        
    Outputs:
        None
    '''
    if bin_type_str == 'mass_LHbeta_bin':
        QA_flag[10] = 1.0    #has large line width on OIII4363
        QA_flag[11] = 1.0    #has large line width on OIII4363
    elif bin_type_str == 'massbin':   
        QA_flag[5] = 1.0     #has large line width on OIII4363
        
    return QA_flag
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
