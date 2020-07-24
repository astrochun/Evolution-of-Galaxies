from astropy.io import ascii as asc

from Metallicity_Stack_Commons.analysis.temp_metallicity_calc import temp_calculation, metallicity_calculation 
from Metallicity_Stack_Commons.column_names import temp_metal_names0, bin_ratios0, bin_names0
from Metallicity_Stack_Commons.analysis.ratios import flux_ratios
from Metallicity_Stack_Commons import line_name_short
    
    
def run_function(line_file, bin_file, outfile, EBV = None):
    '''
    Purpose:
        This function gets line ratios and runs the temperature and metallicity calculation functions
        for the case of stacked spectra.
        
    Usage:
        R_temp_calcul.run_function(line_file, outfile)
        
    Params:
        line_file --> ascii table that has emission line measurements.
        bin_file --> ascii table that contains bin information (e.g., bin IDs).
        outfile --> a string naming the temp/metal table produced by this function.
        EBV (OPTIONAL) --> an array of E(B-V) values.
        
    Returns:
        None
        
    Outputs:
        out_ascii --> ascii table containing bin line ratios, temperatures, and metallicities.
        out_fits --> fits table containing bin line ratioes, temperatures, and metallicities.
    '''
    
    print("Reading : ", line_file)
    line_table = asc.read(line_file)
    
    # Case for stacked spectra
    print("Reading : ", bin_file)
    bin_table = asc.read(bin_file, format='fixed_width_two_line')
    bin_IDs = bin_table[bin_names0[0]].data
    
    # Get flux ratios
    flux_dict = {line_name_short['HB']:line_table['HBETA_Flux_Observed'].data,
                 line_name_short['OII']:line_table['OII_3727_Flux_Observed'].data,
                 line_name_short['OIII']:line_table['OIII_5007_Flux_Observed'].data,
                 line_name_short['4363']:line_table['OIII_4363_Flux_Observed'].data}  
    flux_ratios_dict = flux_ratios(flux_dict, flux_type = 'composite')
  
    # Calculate composite electron temperature and metallicity
    Te = temp_calculation(flux_ratios_dict[bin_ratios0[-1]])
    metal_dict = metallicity_calculation(Te, flux_ratios_dict[bin_ratios0[2]], flux_ratios_dict[bin_ratios0[3]])
    
    # Get bin IDs and composite measurements in one dicitonary and write to table
    tbl_dict = {bin_names0[0]:bin_IDs}
    tbl_dict.update(flux_ratios_dict)
    tbl_dict.update(T_e = Te)
    tbl_dict.update(metal_dict)
    n = tuple([bin_names0[0]] + bin_ratios0 + temp_metal_names0)
    
    print("Writing : ", outfile)
    asc.write(tbl_dict, names=n, output=outfile, format='fixed_width_two_line', overwrite=True)
    