import numpy as np
from astropy.io import ascii as asc
from astropy.table import Table
from Metallicity_Stack_Commons.analysis.temp_metallicity_calc import R_calculation, temp_calculation, metallicity_calculation 
from Metallicity_Stack_Commons.column_names import temp_metal_names0, bin_ratios0, bin_names0
from Metallicity_Stack_Commons import OIII_r
    
    
def run_function(line_file, bin_file, outfile, EBV, k_4363, k_5007):
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
    
    if 'Log10(Mass)' in line_table.keys():
        #Case for individual spectra
        '''
        OII = line_table['OII_3727_Flux_Observed'].data
        SN_OII = line_table['OII_3727_S/N'].data
        HBETA = line_table['HBETA_Flux_Observed'].data
        SN_HBETA = line_table['HBETA_S/N'].data
        OIII5007 = line_table['OIII_5007_Flux_Observed'].data
        SN_5007 = line_table['OIII_5007_S/N'].data
        
        source_ID = line_table['ID'].data
        mass_bin_ID = line_table['Mass_Bin_ID'].data
        HB_bin_ID = line_table['Mass_LHBeta_Bin_ID'].data
        log_mass = line_table['Log10(Mass)'].data
        LHbeta = line_table['HBeta_Luminosity'].data
        mass_T_e = line_table['Mass_Bin_Te'].data
        HB_T_e = line_table['Mass_LHBeta_Bin_Te'].data
        #HB_bin_detect = line_table['Mass_LHBeta_Bin_Detections'].data
        #mass_bin_detect = line_table['Mass_Bin_Detections'].data
        #mass_indiv_detect = line_table['Mass_Individual_Detections'].data 
        #HB_indiv_detect = line_table['MassLHB_Individual_Detections'].data 
        #ebv = line_table['E(B-V)'].data
        
        #Detection variables here include the non-detections with reliable limits so that the metallicity
        #is calculated for those sources.
        
        mass_detect = np.where(((mass_bin_detect == 1.0) | (mass_bin_detect == 0.5)) & (np.isfinite(OIII5007) == True) &
                               (OIII5007 >= 1e-18) & (OIII5007 <= 1e-15) & (np.isfinite(OII) == True) &
                               (OII >= 1e-18) & (OII <= 1e-15) & (np.isfinite(HBETA) == True) & 
                               (HBETA >= 1e-18) & (HBETA <= 1e-15))[0]
        HB_detect = np.where(((HB_bin_detect == 1.0) | (HB_bin_detect == 0.5)) & (np.isfinite(OIII5007) == True) & 
                             (OIII5007 >= 1e-18) & (OIII5007 <= 1e-15) & (np.isfinite(OII) == True) & 
                             (OII >= 1e-18) & (OII <= 1e-15) & (np.isfinite(HBETA) == True) & 
                             (HBETA >= 1e-18) & (HBETA <= 1e-15) & (LHbeta > 0))[0]
        mass_nondetect = np.where((mass_bin_detect == 0.5) & (np.isfinite(OIII5007) == True) & 
                                  (OIII5007 >= 1e-18) & (OIII5007 <= 1e-15) & (np.isfinite(OII) == True) & 
                                  (OII >= 1e-18) & (OII <= 1e-15) & (np.isfinite(HBETA) == True) & 
                                  (HBETA >= 1e-18) & (HBETA <= 1e-15))[0]
        HB_nondetect = np.where((HB_bin_detect == 0.5) & (np.isfinite(OIII5007) == True) & 
                                (OIII5007 >= 1e-18) & (OIII5007 <= 1e-15) & (np.isfinite(OII) == True) & 
                                (OII >= 1e-18) & (OII <= 1e-15) & (np.isfinite(HBETA) == True) & 
                                (HBETA >= 1e-18) & (HBETA <= 1e-15) & (LHbeta > 0))[0]
        
        
        #Correct detection markings (detection vs non-detection w/ reliable limits) are applied here
        #mass_indiv_detect[mass_detect] = 1.0
        #HB_indiv_detect[HB_detect] = 1.0
        #mass_indiv_detect[mass_nondetect] = 0.5
        #HB_indiv_detect[HB_nondetect] = 0.5
        
        #create zero arrays all same length
        two_beta = np.zeros(len(source_ID))
        three_beta = np.zeros(len(source_ID))
        R23 = np.zeros(len(source_ID))
        O32 = np.zeros(len(source_ID))  
        
        
        two_beta[mass_detect] = OII[mass_detect] / HBETA[mass_detect]
        two_beta[HB_detect] = OII[HB_detect] / HBETA[HB_detect]
        three_beta[mass_detect] = (OIII5007[mass_detect] * (1 + 1/OIII_r)) / HBETA[mass_detect]
        three_beta[HB_detect] = (OIII5007[HB_detect] * (1 + 1/OIII_r)) / HBETA[HB_detect]
        
        #Calculate R23 and O32
        R23[mass_detect] = np.log10((OII[mass_detect] + ((1 + 1/OIII_r) * OIII5007[mass_detect])) / HBETA[mass_detect])
        R23[HB_detect] = np.log10((OII[HB_detect] + ((1 + 1/OIII_r) * OIII5007[HB_detect])) / HBETA[HB_detect])
        O32[mass_detect] = np.log10(((1 + 1/OIII_r) * OIII5007[mass_detect]) / OII[mass_detect])
        O32[HB_detect] = np.log10(((1 + 1/OIII_r) * OIII5007[HB_detect]) / OII[HB_detect])
        
        
        
        mass_com_O_log, mass_metal_dict = metallicity_calculation(mass_T_e, two_beta, three_beta)
        HB_com_O_log, HB_metal_dict = metallicity_calculation(HB_T_e, two_beta, three_beta)
        mass_O_s_ion = mass_metal_dict["O_s_ion"]
        mass_O_d_ion = mass_metal_dict["O_d_ion"]
        HB_O_s_ion = HB_metal_dict["O_s_ion"]
        HB_O_d_ion = HB_metal_dict["O_d_ion"]
        
        n = ('Source_ID', 'Mass_Bin_ID', 'Mass_LHBeta_Bin_ID', 'Mass_Bin_Detections', 'Mass_LHBeta_Bin_Detections',
             'Mass_Individual_Detections', 'MassLHB_Individual_Detections', 'Log10(Mass)', 'HBeta_Luminosity', 'OIII_5007_Flux_Observed', 
             'OIII_4958_Flux_Observed', 'OII_3727_Flux_Observed', 'HBETA_Flux_Observed', 'Mass_Bin_Te',
             'Mass_LHBeta_Bin_Te', 'R23', 'O32', 'OII/HBeta', 'OIII/HBeta', 'Mass_Bin_O_s_ion',
             'Mass_Bin_O_d_ion', 'Mass_Bin_com_O_log', 'Mass_LHBeta_Bin_O_s_ion', 'Mass_LHBeta_Bin_O_d_ion',
             'Mass_LHBeta_Bin_com_O_log', 'E(B-V)')
        
        tab0 = Table([source_ID, mass_bin_ID, HB_bin_ID, mass_bin_detect, HB_bin_detect,
                      mass_indiv_detect, HB_indiv_detect, log_mass, LHbeta, OIII5007, OIII4959, OII, HBETA,
                      mass_T_e, HB_T_e, R23, O32, two_beta, three_beta, mass_O_s_ion, mass_O_d_ion, 
                      mass_com_O_log, HB_O_s_ion, HB_O_d_ion, HB_com_O_log, ebv], names = n)
        '''
        
    else:
        #Case for stacked spectra
        bin_table = asc.read(bin_file, format = 'fixed_width_two_line')
        bin_IDs = bin_table[bin_names0[0]].data
        
        OII = line_table['OII_3727_Flux_Observed'].data
        OIII4363 = line_table['OIII_4363_Flux_Observed'].data
        OIII5007 = line_table['OIII_5007_Flux_Observed'].data
        HBETA = line_table['HBETA_Flux_Observed'].data
        
        two_beta = OII / HBETA
        three_beta = (OIII5007 * (1 + 1/OIII_r)) / HBETA
        
        #Calculate R23 composite and O32 composite
        logR23_comp = np.log10((OII + ((1 + 1/OIII_r) * OIII5007)) / HBETA)
        logO32_comp = np.log10(((1 + 1/OIII_r) * OIII5007) / OII)
        
        #R, Te, and metallicity calculations
        R_value = R_calculation(OIII4363, OIII5007, EBV)
        T_e = temp_calculation(R_value)
        metal_dict = metallicity_calculation(T_e, two_beta, three_beta)
        
        n = tuple([bin_names0[0]] + bin_ratios0 + temp_metal_names0)
        tab0 = Table([bin_IDs, logR23_comp, logO32_comp, two_beta, three_beta, T_e, metal_dict['12+log(O/H)'], 
                      metal_dict['log(O+/H)'], metal_dict['log(O++/H)'], metal_dict['O+/H'],
                      metal_dict['O++/H']], names = n)        

    
    asc.write(tab0, outfile, format = 'fixed_width_two_line', overwrite = True)
    tab0.write(outfile.replace('.tbl', '.fits'), format = 'fits', overwrite = True) 
    
    
    
    
    
    
    
    