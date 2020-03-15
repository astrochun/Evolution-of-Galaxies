import numpy as np
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import Table
from Metallicity_Stack_Commons.column_names import filename_dict, bin_names0, bin_mzevolve_names0, indv_names0, gauss_names0
from Metallicity_Stack_Commons import line_name
#from Metallicity_Stack_Commons.analysis.composite_indv_detect import main
    



def indiv_bin_info_table(fitspath, line_file, bin_npz_file, valid_file, LHb_bin = False):
    #This function makes individual_bin_info.tbl   

    '''
    REDO DOCUMENTATION
    Purpose:
        This function creates a table the size of the number of individual galaxies with valid mass, 
        and it holds individual galaxy line measurements, bin temperatures (for both cases) applied to the
        galaxies within each bin, detection markings, and individual and bin ID numbers.
        
    Usage:
        indiv_gals.create_Te_line_table(fitspath, line_file, mass_bin_file, HB_bin_file, mass_Te_file, HB_Te_file)
        
    Params:
        fitspath --> a string of the file path where the input file is and where the output file will be placed.
        line_file --> a fits table containing the individual line measurements.
        mass_bin_file --> an npz file containing the indices of galaxies, relative to line_file, contained in 
            each mass bin. It also contains the mass values of all sources (indexed the same as line_file).
        HB_bin_file --> an npz file containing the indices of galaxies, relative to line_file, contained in 
            each mass-LHbeta bin. It also contains the log of the HBeta luminosity values of all sources.
        mass_Te_file --> an ascii table containing the electron temperatures and detection markings for each
            mass bin.
        HB_Te_file --> an ascii table containing the electron temperatures and detection markings for each
            mass-LHbeta bin.
        
    Returns:
        None
        
    Outputs:
        fitspath + 'individual_Te_emLines.tbl' --> an ascii table containing individual and bin ID numbers,
            bin temperatures correlated to individual galaxies in each bin, detection markings, and individual
            line fluxes and signal to noise.
    '''
    
    if LHb_bin == True:
        bin_type = 'logLHb'
    else:
        bin_type = 'logM'
    
    bin_npz = np.load(bin_npz_file)
    bin_tbl = asc.read(fitspath + filename_dict['bin_info'])
    valid_tbl = asc.read(fitspath + filename_dict['bin_valid'])
    hdu = fits.open(line_file)
    
    
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
  
    
    
    
    
def indiv_em_table(fitspath, line_file, bin_npz_file): 
    #This function makes individual_properties.tbl
    
    hdu = fits.open(line_file)
    bin_npz = np.load(bin_npz_file)
    
    #Get logMass and logLHbeta
    all_logM = bin_npz['mass']
    all_logLHb = bin_npz['lum']
    bin_idx = bin_npz['bin_ind']
    valid_idx = []
    for ii in range(len(bin_idx)):
        valid_idx += list(bin_idx[ii])
    idx_array = np.array(valid_idx)
    logM = all_logM[idx_array]
    logLHb = all_logLHb[idx_array]
    
    #Get individual lines fluxes and S/N 
    line_table = hdu[1].data 
    line_table = line_table[idx_array]
    ID = line_table['OBJNO']      
    OII_Flux = line_table['OII_FLUX_DATA']
    OII_SN = line_table['OII_SNR']
    HBETA_Flux = line_table['HB_FLUX_DATA']
    HBETA_SN = line_table['HB_SNR']
    OIII5007_Flux = line_table['OIIIR_FLUX_DATA']
    OIII5007_SN = line_table['OIIIR_SNR']
    
    line_cols = [gauss_names0[0], gauss_names0[2]]
    line_names = [line_name[0], line_name[4], line_name[-1]] 
    cols = []
    for ii in range(len(line_names)):
        cols += [line_names[ii] + '_' + line_cols[0]]
        cols += [line_names[ii] + '_' + line_cols[1]]

    out_ascii = fitspath + filename_dict['indv_prop']
    n = [indv_names0[0]] + indv_names0[3:5] + cols 
    indiv_props_tbl = Table([ID, logM, logLHb, OII_Flux, OII_SN, HBETA_Flux, HBETA_SN, OIII5007_Flux, 
                             OIII5007_SN], names = n)
    asc.write(indiv_props_tbl, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    
    

'''
def indiv_derived_props(fitspath, bin_derpropsrev_file, indiv_props_file, indiv_bin_file):
    #This function makes individual_derived_properties.tbl
    
    #outfile = fitspath + filename_dict['indv_derived_prop']
    outfile = fitspath + 'individual_derived_properties.tbl'
    
    #Create individual_derived_props
    main(fitspath, '', bin_derpropsrev_file, indiv_props_file, indiv_bin_file, outfile, det3 = False)  
'''   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    