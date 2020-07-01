import numpy as np
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import Table
from Metallicity_Stack_Commons.column_names import filename_dict, bin_names0, bin_mzevolve_names0, indv_names0, line_fit_suffix_add
from Metallicity_Stack_Commons import line_name, line_type
    



def indiv_bin_info_table(fitspath, line_file, LHb_bin = False, use_revised = False):
    '''
    Purpose:
        This function creates a table containing bin information such as bin ID numbers, minimum, 
        maximum, median, and average mass and/or HBeta luminosity measurements, bin detection markings,
        and indiviudal galaxy ID numbers.
        
    Params:
        fitspath --> file path where bin files are located and where the output file will be located
        line_file --> a fits table containing individual galaxy ID numbers
        LHb_bin (OPTIONAL) --> boolean to determine which binning type was used (default is False)
        use_revised (OPTIONAL) --> boolean to use the revised validation table or not (default is False)
        
    Returns:
        None
        
    Outputs:
        'individual_bin_info.tbl'
    '''
    
    if LHb_bin == True:
        bin_type = 'logLHb'
    else:
        bin_type = 'logM'
    
    bin_npz = np.load(fitspath + filename_dict['bin_info'].replace('.tbl', '.npz'))
    bin_tbl = asc.read(fitspath + filename_dict['bin_info'])
    hdu = fits.open(line_file)
    if use_revised == True:
        valid_tbl = asc.read(fitspath + filename_dict['bin_valid_rev'])
    else:
        valid_tbl = asc.read(fitspath + filename_dict['bin_valid'])
    
    
    line_table = hdu[1].data
    objno = line_table['OBJNO']
    bin_ID = bin_tbl[bin_names0[0]].data
    bin_min = bin_tbl[bin_type + '_min'].data
    bin_max = bin_tbl[bin_type + '_max'].data
    bin_avg = bin_tbl[bin_type + '_avg'].data
    bin_median = bin_tbl[bin_type + '_median'].data
    bin_idx = bin_npz['bin_ind']     #valid mass bin indices relative to unsorted data table 
    detect = valid_tbl[bin_names0[2]].data                   
    
    
    #Get all bin valid values in 1D array
    valid_idx = []
    bin_ID_array = []
    bin_min_array = []
    bin_max_array = []
    bin_avg_array = []
    bin_median_array = []
    detect_array = []
    for ii in range(len(bin_idx)):
        valid_idx += list(bin_idx[ii])
        for jj in range(len(bin_idx[ii])):
            bin_ID_array.append(bin_ID[ii]) 
            bin_min_array.append(bin_min[ii])
            bin_max_array.append(bin_max[ii])
            bin_avg_array.append(bin_avg[ii])
            bin_median_array.append(bin_median[ii]) 
            detect_array.append(detect[ii])            
    idx_array = np.array(valid_idx)                     

    
    out_ascii = fitspath + filename_dict['indv_bin_info']  
    bin_cols = [bin_names0[0], bin_names0[2], indv_names0[0]]
    bin_cols += [col for col in bin_mzevolve_names0 if col.startswith(bin_type)]
    indiv_bin_info_tbl = Table([bin_ID_array, detect_array, objno[idx_array], bin_min_array, bin_max_array, bin_avg_array,
                                bin_median_array], names = tuple(bin_cols))
    asc.write(indiv_bin_info_tbl, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    
    
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
    indiv_props_tbl = Table(data, names = tuple([indv_names0[0]] + indv_names0[3:5] + cols))
    asc.write(indiv_props_tbl, out_ascii, format = 'fixed_width_two_line', overwrite = True) 