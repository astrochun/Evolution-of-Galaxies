from astropy.io import ascii as asc
from astropy.table import Table
from Metallicity_Stack_Commons.analysis.temp_metallicity_calc import temp_calculation, metallicity_calculation 
from Metallicity_Stack_Commons.column_names import temp_metal_names0, bin_ratios0, bin_names0
from Metallicity_Stack_Commons.analysis.ratios import flux_ratios
from Metallicity_Stack_Commons import line_name_short
    
    
def run_function(line_file, bin_file, outfile):
    '''
    Purpose:
        This function runs the R calculation, temperature calculation, and metallicity calculation
        functions for the case of individual and stacked spectra.
        
    Usage:
        R_temp_calcul.run_function(line_file, outfile, EBV, k_4363, k_5007)
        
    Params:
        line_file --> ascii table that has all the data for each emission line.
        outfile --> a string naming the tables that are produced by this function.
        EBV --> an array of E(B - V) values. The array is the size of the number of sources (or bins).
        k_4363 --> an array of dust extinction values at the OIII4363 emission line for each bin/individual source.
        k_5007 --> an array of dust extinction values at the OIII5007 emission line for each bin/individual source.
        
    Returns:
        None
        
    Outputs:
        out_ascii --> an ascii table containing the bin temperatures, bin or individual metallicities,
            and bin or individual line fluxes and S/N.
        out_fits --> a fits table containing the bin temperatures, bin or individual metallicities,
            and bin or individual line fluxes and S/N.
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
    
    
    
    
    
    
    
    