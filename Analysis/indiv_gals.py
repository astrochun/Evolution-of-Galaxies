import numpy as np
from astropy.io import fits
from astropy.io import ascii as asc
from Metallicity_Stack_Commons.column_names import filename_dict, bin_names0, bin_mzevolve_names0, indv_names0, line_fit_suffix_add
from Metallicity_Stack_Commons import line_name, line_type
from Evolution_of_Galaxies.general import table_to_dict
    



def indiv_bin_info_table(fitspath, line_file, use_revised = False):
    '''
    Purpose:
        This function creates a table containing bin information such as bin ID numbers, minimum, 
        maximum, median, and average mass and/or HBeta luminosity measurements, bin detection markings,
        and indiviudal galaxy ID numbers.
        
    Params:
        fitspath --> file path where bin files are located and where the output file will be located
        line_file --> a fits table containing individual galaxy ID numbers
        use_revised (OPTIONAL) --> boolean to use the revised validation table or not (default is False)
        
    Returns:
        None
        
    Outputs:
        'individual_bin_info.tbl'
    '''
    
    bin_npz = np.load(fitspath + filename_dict['bin_info'].replace('.tbl', '.npz'))
    hdu = fits.open(line_file)
    if use_revised == True:
        valid_tbl = asc.read(fitspath + filename_dict['bin_valid_rev'])
    else:
        valid_tbl = asc.read(fitspath + filename_dict['bin_valid'])
        
    #Convert bin info table into a dictionary
    bin_dict = table_to_dict(fitspath + filename_dict['bin_info'])
    
    line_table = hdu[1].data
    objno = line_table['OBJNO']
    bin_idx = bin_npz['bin_ind']     #valid mass bin indices relative to unsorted data table 
    detect = valid_tbl[bin_names0[2]].data                   
    
    
    #Get all bin valid values in dictionary
    bin_cols = [bin_names0[0], bin_names0[2], indv_names0[0]] + bin_mzevolve_names0
    indiv_dict = {}
    for name in bin_cols:
        indiv_dict[name] = []
    
    valid_idx = []
    for ii in range(len(bin_idx)):
        valid_idx += list(bin_idx[ii])
        for jj in range(len(bin_idx[ii])):
            indiv_dict[bin_names0[0]] += [bin_dict[bin_names0[0]][ii]]   #bin_ID
            indiv_dict[bin_names0[2]] += [detect[ii]]                    #Detection
            indiv_dict[bin_mzevolve_names0[0]] += [bin_dict[bin_mzevolve_names0[0]][ii]]   #logM_min
            indiv_dict[bin_mzevolve_names0[1]] += [bin_dict[bin_mzevolve_names0[1]][ii]]   #logM_max
            indiv_dict[bin_mzevolve_names0[2]] += [bin_dict[bin_mzevolve_names0[2]][ii]]   #logM_avg
            indiv_dict[bin_mzevolve_names0[3]] += [bin_dict[bin_mzevolve_names0[3]][ii]]   #logM_median
            indiv_dict[bin_mzevolve_names0[4]] += [bin_dict[bin_mzevolve_names0[4]][ii]]   #logLHb_min
            indiv_dict[bin_mzevolve_names0[5]] += [bin_dict[bin_mzevolve_names0[5]][ii]]   #logLHb_max
            indiv_dict[bin_mzevolve_names0[6]] += [bin_dict[bin_mzevolve_names0[6]][ii]]   #logLHb_avg
            indiv_dict[bin_mzevolve_names0[7]] += [bin_dict[bin_mzevolve_names0[7]][ii]]   #logLHb_median
    idx_array = np.array(valid_idx)
    indiv_dict[indv_names0[0]] = objno[idx_array]      #ID     

    
    out_ascii = fitspath + filename_dict['indv_bin_info']
    asc.write(indiv_dict, names = tuple(bin_cols), output = out_ascii, format = 'fixed_width_two_line', overwrite = True)
    
    
    hdu.close()
  
    
    
    
    
def indiv_em_table(fitspath, line_file): 
    '''
    Purpose:
        This function creates a table containing individual galaxy ID numbers, masses, HBeta luminosities,
        and line measurements for OII_3727, HBeta, and OIII_5007.
    
    Params:
        fitspath --> file path where bin files are located and where the output file will be located
        line_file --> a fits table containing the individual line measurements
        
    Returns:
        None
        
    Outputs:
        'individual_properties.tbl'
    '''
    
    #Read in individual line measurement file and bin information file
    hdu = fits.open(line_file)
    bin_npz = np.load(fitspath + filename_dict['bin_info'].replace('.tbl', '.npz'))
    
    #Get indices valid for binning
    valid_idx = []
    for ii in range(len(bin_npz['bin_ind'])):
        valid_idx += list(bin_npz['bin_ind'][ii])
    idx_array = np.array(valid_idx)
    
    #Get column names from individual line measurement file
    line_tbl_names = ['_FLUX_MOD', '_FLUX_DATA', '_SNR', '_LAMBDA', '_PEAK', '_Y0', '_SIGMA', '_NOISE']
    line_tbl_prefixes = ['OIIB', 'HB', 'OIIIR']
    line_tbl_col = []
    for mm in range(len(line_tbl_prefixes)):
        for nn in range(len(line_tbl_names)):
            line_tbl_col += [line_tbl_prefixes[mm] + line_tbl_names[nn]]
    
    #Get column names for output table
    cols = []
    line_names = [line_name[0], line_name[4], line_name[-1]]
    line_types = [line_type[0], line_type[4], line_type[-1]]
    for ii in range(len(line_names)):
        cols += line_fit_suffix_add(line_names[ii], line_types[ii])
    cols.remove('HBETA_Abs_Norm')
    cols.remove('HBETA_Abs_Sigma')
    
    #Collect data in dictionary and write to output table
    line_table = hdu[1].data 
    line_table = line_table[idx_array]
    data = {indv_names0[0]:line_table['OBJNO'], indv_names0[3]:bin_npz['mass'][idx_array], 
            indv_names0[4]:bin_npz['lum'][idx_array]}
    for jj in range(len(cols)):
        data[cols[jj]] = line_table[line_tbl_col[jj]]

    out_ascii = fitspath + filename_dict['indv_prop']
    asc.write(data, names = tuple([indv_names0[0]] + indv_names0[3:5] + cols), output = out_ascii, format = 'fixed_width_two_line', overwrite = True) 