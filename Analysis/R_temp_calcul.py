from astropy.io import ascii as asc
from astropy.table import Table
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
    
    line_table = asc.read(line_file)
    
    #Case for stacked spectra
    bin_table = asc.read(bin_file, format = 'fixed_width_two_line')
    bin_IDs = bin_table[bin_names0[0]].data
    
    flux_dict = {line_name_short['HB']:line_table['HBETA_Flux_Observed'].data,
                 line_name_short['OII']:line_table['OII_3727_Flux_Observed'].data,
                 line_name_short['OIII']:line_table['OIII_5007_Flux_Observed'].data,
                 line_name_short['4363']:line_table['OIII_4363_Flux_Observed'].data}  
    flux_ratios_dict = flux_ratios(flux_dict)
    
    logR23 = flux_ratios_dict['logR23']
    logO32 = flux_ratios_dict['logO32']
    two_beta = flux_ratios_dict['two_beta']
    three_beta = flux_ratios_dict['three_beta']
    R = flux_ratios_dict['R']
    
    T_e = temp_calculation(R)
    metal_dict = metallicity_calculation(T_e, two_beta, three_beta)
    
    n = tuple([bin_names0[0]] + bin_ratios0 + temp_metal_names0)
    tab0 = Table([bin_IDs, logR23, logO32, two_beta, three_beta, T_e, metal_dict['12+log(O/H)'], 
                  metal_dict['log(O+/H)'], metal_dict['log(O++/H)'], metal_dict['O+/H'],
                  metal_dict['O++/H']], names = n) 

    
    asc.write(tab0, outfile, format = 'fixed_width_two_line', overwrite = True)
    tab0.write(outfile.replace('.tbl', '.fits'), format = 'fits', overwrite = True) 