#Calculates the R value, electron temperature, and metallicity from the flux table
#produced by the zoom_and_gauss_general functions
#Currently running: Grid
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from getpass import getuser
    

#For generalizing for several users
if getuser() == 'carol':
    fitspath = "C:\\Users\\carol\\Google Drive\\MZEvolve\\"
    fitspath2 = fitspath + "massbin\\"
else:
    fitspath = "../DEEP2/" 
    fitspath2 = "../"
    
bin_pts_input = [75, 112, 113, 300, 600, 1444, 1444]
str_bin_pts_input = [str(val) for val in bin_pts_input]
bin_pts_fname = "_".join(str_bin_pts_input)
bin_pts_fname = 'hbeta_revised_' + bin_pts_fname
bin_pts_fname2 = 'individual'
#two/three components to bin_pts_fname: bin type (individual or stacks --> if stacks then size of stacks) and
#revised data or not

#make a marker to determine which case the file is & update file naming convention
 
#don't need --> just replace with bin_pts_fname   
#N_in_bin = bin_pts_fname

mark_nondet = False
if mark_nondet:
    updated = '_updated'
else:
    updated = ''


#Constants
a = 13205
b = 0.92506
c = 0.98062
 
    

def R_calculation(OIII4363, OIII5007):
    R_value = OIII4363 / (OIII5007 * (4/3))
    return R_value


def temp_calculation(R):
    T_e = a*(-np.log10(R) - b)**(-1*c)      
    return T_e


def metallicity_calculation(T_e, two_beta, three_beta): 
    O_s_ion = np.zeros(len(T_e))
    O_d_ion = np.zeros(len(T_e))
    com_O = np.zeros(len(T_e))
    com_O_log = np.zeros(len(T_e))
    O_s_ion_log = np.zeros(len(T_e))
    O_d_ion_log = np.zeros(len(T_e))
    t_3 = np.zeros(len(T_e))
    t_2 = np.zeros(len(T_e))
    x2 = np.zeros(len(T_e))
    
    detect = np.where((two_beta != 0) & (three_beta != 0))[0]
    
    t_3[detect] = T_e[detect] * 1e-4
    t_2[detect] = 0.7 * t_3[detect] + 0.17
    x2[detect] = 1e-4 * 1e3 * t_2[detect]**(-0.5)

    O_s_ion_log[detect] = np.log10(two_beta[detect]) + 5.961 + 1.676 / t_2[detect] - 0.4 * np.log10(t_2[detect]) - 0.034 * t_2[detect] + np.log10(1 + 1.35 * x2[detect]) - 12
    O_d_ion_log[detect] = np.log10(three_beta[detect]) + 6.200 + 1.251 / t_3[detect] - 0.55 * np.log10(t_3[detect]) - 0.014 * (t_3[detect]) - 12

    O_s_ion[detect] = 10**(O_s_ion_log[detect])
    O_d_ion[detect] = 10**(O_d_ion_log[detect])
    com_O[detect] = O_s_ion[detect] + O_d_ion[detect]
    com_O_log[detect] = np.log10(com_O[detect]) + 12

    return O_s_ion, O_d_ion, com_O_log, O_s_ion_log, O_d_ion_log



def derived_properties_plots():    
    metal_Te_file = fitspath2 + 'individual_derived_properties_metallicity.tbl'
    MT_ascii = asc.read(metal_Te_file)

    log_mass = MT_ascii['Log10(Mass)'].data
    LHbeta = MT_ascii['HBeta_Luminosity'].data
    ind_metal = MT_ascii['com_O_log'].data
    R23 = MT_ascii['R23'].data
    O32 = MT_ascii['O32'].data
    indiv_detect = MT_ascii['Individual Detections'].data
    
    detection = np.where((indiv_detect == 1.0) & (LHbeta > 0))[0]

    pdf_pages = PdfPages(fitspath2 + 'individual_metal_plots.pdf')

    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot1 = ax1.scatter(log_mass[detection], LHbeta[detection], 0.8, c=ind_metal[detection], marker='*')
    cb = fig1.colorbar(plot1)
    cb.set_label('Metallicity')
    ax1.set_xlabel('Mass')
    ax1.set_ylabel('LHBeta')
    ax1.set_title('Mass vs. Luminosity Colormap=Metallicity')
    fig1.set_size_inches(8, 8)
    fig1.savefig(pdf_pages, format='pdf')
    
    fig2, ax2 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot2 = ax2.scatter(R23[detection], O32[detection], 0.8, c=ind_metal[detection], marker='*')
    cb = fig2.colorbar(plot2)
    cb.set_label('Metallicity')
    ax2.set_xlabel('R23')
    ax2.set_ylabel('O32')
    ax2.set_title('O32 vs. R23 Colormap=Metallicity')
    fig2.set_size_inches(8, 8)
    fig2.savefig(pdf_pages, format='pdf')
    
    fig3, ax3 = plt.subplots()
    ax3.scatter(R23[detection], ind_metal[detection], marker='*', label = 'R23')
    ax3.scatter(O32[detection], ind_metal[detection], marker='*', label = 'O32')
    ax3.legend(loc = 'best')
    ax3.set_ylabel('Metallicity')
    ax3.set_title('Metallicity vs. R23 and O32')
    fig3.set_size_inches(8, 8)
    fig3.savefig(pdf_pages, format='pdf')
    
    pdf_pages.close()



def run_function():
       
    line_file = fitspath2 + 'indivgals_Te_lineRatios.tbl'
    #line_file = fitspath2 + 'hbeta_revised_75_112_113_300_600_1444_1444_updated_massbin_emission_lines.tbl'

    line_table = asc.read(line_file)
    
    if 'two_beta' in line_table.keys():
        #Case for individual spectra 
        
        out_ascii = fitspath2 + bin_pts_fname2 + '_derived_properties_metallicity.tbl'
        out_fits = fitspath2 + bin_pts_fname2 + '_derived_properties_metallicity.fits'
        
        OIII4959 = line_table['OIII4959'].data
        OIII5007 = line_table['OIII5007'].data
        HBETA = line_table['HBeta'].data
        HGAMMA = line_table['HGamma'].data
        SNR_HG = line_table['SNR_HG'].data
        raw_OIII4363 = line_table['OIII4363'].data
        SNR_4363 = line_table['SNR_4363'].data
        R23_individual = line_table['Individual_R23'].data
        O32_individual = line_table['Individual_O32'].data
        detections = line_table['Detection'].data
        
        two_beta = line_table['two_beta'].data
        three_beta = line_table['three_beta'].data
        T_e = line_table['Temperature'].data
        source_ID = line_table['Source_ID'].data
        
        
        O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metallicity_calculation(T_e, two_beta, three_beta)
        
        n = ('Source_ID', 'R23', 'O32', 'Observed_Flux_5007', 'Observed_Flux_4959',
             'Observed_Flux_HBeta', 'Temperature', 'Detections', 'O_s_ion', 'O_d_ion', 'com_O_log')
        tab0 = Table([source_ID, R23_individual, O32_individual, OIII5007, OIII4959, HBETA, T_e, detections,
                      O_s_ion, O_d_ion, com_O_log], names = n)
    
    
    elif 'Log10(Mass)' in line_table.keys():
        #Case for individual spectra
    
        out_ascii = fitspath2 + bin_pts_fname2 + '_derived_properties_metallicity.tbl'
        out_fits = fitspath2 + bin_pts_fname2 + '_derived_properties_metallicity.fits'
        
        OII = line_table['OII_Flux'].data
        SN_OII = line_table['OII_SN'].data
        OIII4959 = line_table['OIII4959_Flux'].data
        SN_4959 = line_table['OIII4959_SN'].data
        OIII5007 = line_table['OIII5007_Flux'].data
        SN_5007 = line_table['OIII5007_SN'].data
        HBETA = line_table['HBETA_Flux'].data
        SN_HBETA = line_table['HBETA_SN'].data
        
        source_ID = line_table['OBJNO'].data
        log_mass = line_table['Log10(Mass)'].data
        LHbeta = line_table['HBeta_Luminosity'].data
        T_e = line_table['Te'].data
        bin_detect = line_table['Bin Detections'].data
        indiv_detect = line_table['Individual Detections'].data 
        
        detect = np.where((bin_detect == 1.0) & (np.isfinite(OIII5007) == True) & (OIII5007 >= 1e-18) &
                          (OIII5007 <= 1e-15) & (np.isfinite(OIII4959) == True) & (OIII4959 >= 1e-18) & 
                          (OIII4959 <= 1e-15) & (np.isfinite(OII) == True) & (OII >= 1e-18) & (OII <= 1e-15) &
                          (np.isfinite(HBETA) == True) & (HBETA >= 1e-18) & (HBETA <= 1e-15))[0]
        indiv_detect[detect] = 1.0
        
        #create zero arrays all same length
        two_beta = np.zeros(len(T_e))
        three_beta = np.zeros(len(T_e))
        R23 = np.zeros(len(T_e))
        O32 = np.zeros(len(T_e))        
        
        two_beta[detect] = OII[detect] / HBETA[detect]
        three_beta[detect] = (OIII5007[detect] + OIII4959[detect]) / HBETA[detect]
        
        #Calculate R23 and O32
        R23[detect] = np.log10((OII[detect] + ((4/3) * OIII5007[detect])) / HBETA[detect])
        O32[detect] = np.log10(((4/3) * OIII5007[detect]) / OII[detect])
        
        
        O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metallicity_calculation(T_e, two_beta, three_beta)
        
        
        n = ('Source_ID', 'Log10(Mass)', 'HBeta_Luminosity', 'Observed_Flux_5007', 'Observed_Flux_4959',
             'Observed_Flux_HBeta', 'Temperature', 'Bin Detections','Individual Detections', 'R23',
             'O32', 'O_s_ion', 'O_d_ion', 'com_O_log')
        tab0 = Table([source_ID, log_mass, LHbeta, OIII5007, OIII4959, HBETA, T_e, bin_detect, indiv_detect,
                      R23, O32, O_s_ion, O_d_ion, com_O_log], names = n)
        
    else:
        #Case for stacked spectra
        
        out_ascii = fitspath2 + bin_pts_fname + updated + '_derived_properties_metallicity.tbl'
        out_fits = fitspath2 + bin_pts_fname + updated + '_derived_properties_metallicity.fits'
        
        OII = line_table['OII_3727_Flux_Observed'].data
        SN_OII = line_table['OII_3727_S/N'].data
        OIII4363 = line_table['OIII_4363_Flux_Observed'].data
        SN_4363 = line_table['OIII_4363_S/N'].data
        OIII4959 = line_table['OIII_4958_Flux_Observed'].data
        SN_4959 = line_table['OIII_4958_S/N'].data
        OIII5007 = line_table['OIII_5007_Flux_Observed'].data
        SN_5007 = line_table['OIII_5007_S/N'].data
        HBETA = line_table['HBETA_Flux_Observed'].data
        SN_HBETA = line_table['HBETA_S/N'].data
        
        N_Galaxy = line_table['Number of Galaxies'].data
        avg_mass = line_table['mass_avg'].data 
        log_mass = avg_mass
        
        two_beta = OII / HBETA
        three_beta = (OIII5007 + OIII4959) / HBETA
        
        #Calculate R23 composite and O32 composite
        R23_composite = np.log10((OII + ((4/3) * OIII5007)) / HBETA)
        O32_composite = np.log10(((4/3) * OIII5007) / OII)
        
        #R, Te, and metallicity calculations
        R_value = R_calculation(OIII4363, OIII5007)
        T_e = temp_calculation(R_value)
        O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metallicity_calculation(T_e, two_beta, three_beta)
        
        
        n = ('R23_Composite', 'O32_Composite', 'N_Galaxies', 'Observed_Flux_5007', 'S/N_5007', 'Observed_Flux_4959',
             'S/N_4959', 'Observed_Flux_4363', 'S/N_4363', 'Observed_Flux_HBETA', 'S/N_HBETA', 'Observed_Flux_3727',
             'S/N_3727', 'Temperature', 'O_s_ion', 'O_d_ion', 'com_O_log')
        tab0 = Table([R23_composite, O32_composite, N_Galaxy, OIII5007, SN_5007, OIII4959, SN_4959,
                      OIII4363, SN_4363, HBETA, SN_HBETA, OII, SN_OII, T_e, O_s_ion, O_d_ion, com_O_log], names = n)

        
    

    
    asc.write(tab0, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    tab0.write(out_fits, format = 'fits', overwrite = True)
    
    
    '''   
    #Plots
    #pdf_pages = PdfPages(fitspath2 + bin_pts_fname + updated + '_massbin_derived_properties_metallicity.pdf')
    pdf_pages = PdfPages(fitspath2 + bin_pts_fname2 + '_derived_properties_metallicity.pdf')

    if mark_nondet:
        #detections = valid_table['Valid_OIII_4363'].data
        non_detection_mark = np.where(detections == 0)[0]
        detection_mark = np.where(detections == 1)[0]
    else:
        detection_mark = np.arange(len(T_e))
    

    fig1, ax1 = plt.subplots()
    ax1.scatter(np.log10(T_e[detection_mark]), R23_composite[detection_mark], marker = '.')
    if mark_nondet:
        ax1.scatter(np.log10(T_e[non_detection_mark]), R23_composite[non_detection_mark], marker = '^')
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('R23')
    ax1.set_title('R23 vs. Temperatures')
    pdf_pages.savefig()
     
    fig2, ax2 = plt.subplots()
    ax2.scatter(np.log10(T_e[detection_mark]), O32_composite[detection_mark], marker = '.')
    if mark_nondet:
        ax2.scatter(np.log10(T_e[non_detection_mark]), O32_composite[non_detection_mark], marker = '^')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('O32')
    ax2.set_title('O32 vs. Temperatures')
    pdf_pages.savefig()

    fig3, ax3 = plt.subplots()
    ax3.scatter(R23_composite[detection_mark], com_O_log[detection_mark], marker = '.')
    if mark_nondet:
        ax3.scatter(R23_composite[non_detection_mark], com_O_log[non_detection_mark], marker = '^')
    ax3.set_xlabel('R23')
    ax3.set_ylabel('12+log(O/H) Te')
    ax3.set_title('Composite Metallicity vs. R23')
    pdf_pages.savefig()

    fig4, ax4 = plt.subplots()
    ax4.scatter(O32_composite[detection_mark], com_O_log[detection_mark], marker = '.')
    if mark_nondet:
        ax4.scatter(O32_composite[non_detection_mark], com_O_log[non_detection_mark], marker = '^')
    ax4.set_xlabel('O32')
    ax4.set_ylabel('12+log(O/H) Te')
    ax4.set_title('Composite Metallicity vs. O32')
    pdf_pages.savefig()
    
    fig5, ax5 = plt.subplots()
    ax5.scatter(avg_mass[detection_mark], np.log10(T_e[detection_mark]), marker = '.')
    if mark_nondet:
        ax5.scatter(avg_mass[non_detection_mark], np.log10(T_e[non_detection_mark]), marker = '^')
    ax5.set_xlabel('Avg Mass')
    ax5.set_ylabel('Temperature (K)')
    ax5.set_title('Temperatures vs. Avg Mass')
    pdf_pages.savefig()
    
    fig6, ax6 = plt.subplots()
    ax6.scatter(avg_mass[detection_mark], com_O_log[detection_mark], marker = '.')
    y = 8.798 - np.log10(1 + ((10**8.901)/(10**avg_mass))**0.640)
    ax6.plot(avg_mass, y, color='g', linestyle = '-', marker = '*')
    if mark_nondet:
        ax6.scatter(avg_mass[non_detection_mark], com_O_log[non_detection_mark], marker = '^')
    ax6.set_xlabel('Avg Mass')
    ax6.set_ylabel('12+log(O/H) Te')
    ax6.set_title('Composite Metallicity vs. Avg Mass')
    pdf_pages.savefig()
    
    pdf_pages.close()    

    #3D plots
    fig_3d = plt.figure(figsize = (10,8))
    ax_3d = plt.axes(projection = '3d')
    ax_3d.set_xlabel('R23')
    ax_3d.set_ylabel('O32')
    ax_3d.set_zlabel('Temperatures')
    ax_3d.scatter(R23_composite[detection_mark], O32_composite[detection_mark], T_e[detection_mark],
                  marker = '.', linewidths = None)
    if mark_nondet:
        ax_3d.scatter(R23_composite[non_detection_mark], O32_composite[non_detection_mark],
                      T_e[non_detection_mark], marker = '^')
    plt.show() 
    '''