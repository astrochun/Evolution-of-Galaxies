import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from Metallicity_Stack_Commons import line_name_short
from Metallicity_Stack_Commons.column_names import temp_metal_names0, bin_ratios0, bin_mzevolve_names0
from Metallicity_Stack_Commons.column_names import bin_names0, filename_dict
from Metallicity_Stack_Commons.analysis.ratios import flux_ratios
from .relation_fitting import extract_error_bars



def bin_derived_props_plots(fitspath, hbeta_bin = False, err_bars = False, revised = False):
    '''
    Purpose: 
        This function creates plots for the stacked spectra (bins): OII/HBeta vs Stellar Mass,
        OIII/HBeta vs Stellar Mass, OIII4363/OIII5007 vs Stellar Mass, R23 vs Temperature, O32 vs Temperature,
        Metallicity vs R23, Metallicity vs O32, Temperature vs Stellar Mass, and Metallicity vs Stellar Mass.
        
    Params:
        fitspath --> string of the file path where data files are located. 
        hbeta_bin (OPTIONAL) --> if the binning type is mass-LHbeta bins, then hbeta_bin = True.
        err_bars (OPTIONAL) --> if error bars for metallicity and temperature should be plotted, then 
                                err_bars = True.
        revised (OPTIONAL) --> if the revised data tables should be used, then revised = True.
        
    Returns:
        None
        
    Outputs:
        pdf_pages --> a pdf containing all of the plots for the stacked spectra on separate pages.
    '''
    
    #Read in data tables
    if revised:
        pdf_pages = PdfPages(fitspath + filename_dict['bin_derived_prop_rev'].replace('.tbl', '.pdf'))
        metal_table = asc.read(fitspath + filename_dict['bin_derived_prop_rev'])
        em_table = asc.read(fitspath + filename_dict['bin_fit_rev'])
        valid_table = asc.read(fitspath + filename_dict['bin_valid_rev'])
    else:
        pdf_pages = PdfPages(fitspath + filename_dict['bin_derived_prop'].replace('.tbl', '.pdf'))
        metal_table = asc.read(fitspath + filename_dict['bin_derived_prop'])
        em_table = asc.read(fitspath + filename_dict['bin_fit'])
        valid_table = asc.read(fitspath + filename_dict['bin_valid'])
    bin_table = asc.read(fitspath + filename_dict['bin_info'])
    
    #Read in composite values
    com_O_log = metal_table[temp_metal_names0[1]].data
    T_e = metal_table[temp_metal_names0[0]].data
    logR23 = metal_table[bin_ratios0[0]].data
    logO32 = metal_table[bin_ratios0[1]].data
    log_twoBeta = np.log10(metal_table[bin_ratios0[2]].data)
    log_threeBeta = np.log10(metal_table[bin_ratios0[3]].data)
    logM_avg = bin_table[bin_mzevolve_names0[2]].data
    logLHb_avg = bin_table[bin_mzevolve_names0[6]].data
    detection = valid_table[bin_names0[2]].data
    
    #Get flux ratios 
    flux_dict = {line_name_short['HB']:em_table['HBETA_Flux_Observed'].data,
                 line_name_short['OII']:em_table['OII_3727_Flux_Observed'].data,
                 line_name_short['OIII']:em_table['OIII_5007_Flux_Observed'].data,
                 line_name_short['4363']:em_table['OIII_4363_Flux_Observed'].data}  
    flux_ratios_dict = flux_ratios(flux_dict)
    logR = np.log10(flux_ratios_dict['R'])
    
    #Define detection and non-detection (with reliable 5007) arrays
    detect = np.where(detection == 1)[0]
    non_detect = np.where(detection == 0.5)[0]
    
    if err_bars == True:
        err_dict = extract_error_bars(fitspath)
    
    #Line Ratios vs Mass
    fig1, ax1 = plt.subplots(2, 3, sharex = True)
    
    for ii in range(len(detect)):
        if hbeta_bin == True:
            if detect[ii] % 2 == 0:
                pt_color = 'cyan'
            else:
                pt_color = 'blue'
            ax1[1, 2].scatter(logM_avg[detect[ii]], logLHb_avg[detect[ii]], marker = '.', color = pt_color)
        else:
            pt_color = 'blue'
        ax1[0, 0].scatter(logM_avg[detect[ii]], logR23[detect[ii]], marker = '.', color = pt_color)
        ax1[0, 1].scatter(logM_avg[detect[ii]], logO32[detect[ii]], marker = '.', color = pt_color)
        ax1[0, 2].scatter(logM_avg[detect[ii]], log_twoBeta[detect[ii]], marker = '.', color = pt_color)
        ax1[1, 0].scatter(logM_avg[detect[ii]], log_threeBeta[detect[ii]], marker = '.', color = pt_color)
        ax1[1, 1].scatter(logM_avg[detect[ii]], logR[detect[ii]], marker = '.', color = pt_color)
    for ii in range(len(non_detect)):
        if hbeta_bin == True:
            if non_detect[ii] % 2 == 0:
                pt_color = 'pink'
            else:
                pt_color = 'red'
            ax1[1, 2].scatter(logM_avg[non_detect[ii]], logLHb_avg[non_detect[ii]], marker = '^', color = pt_color)
        else:
            pt_color = 'orange'
        ax1[0, 0].scatter(logM_avg[non_detect[ii]], logR23[non_detect[ii]], marker = '^', color = pt_color)
        ax1[0, 1].scatter(logM_avg[non_detect[ii]], logO32[non_detect[ii]], marker = '^', color = pt_color)
        ax1[0, 2].scatter(logM_avg[non_detect[ii]], log_twoBeta[non_detect[ii]], marker = '^', color = pt_color)
        ax1[1, 0].scatter(logM_avg[non_detect[ii]], log_threeBeta[non_detect[ii]], marker = '^', color = pt_color)
        ax1[1, 1].scatter(logM_avg[non_detect[ii]], logR[non_detect[ii]], marker = '^', color = pt_color)
    
    ax1[0, 0].annotate('$R_{23}$', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')
    ax1[0, 1].annotate('$O_{32}$', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')
    ax1[0, 2].annotate('[OII]/H$\\beta$', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')
    ax1[1, 0].annotate('[OIII]/H$\\beta$', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')
    ax1[1, 1].annotate('R', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')
    ax1[1, 2].annotate('H$\\beta$ Luminosity', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')  
    
    for ii in range(0, 3):    
        ax1[1, ii].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')     
    for t_ax in ax1:
        for tt in range(len(t_ax)):
            t_ax[tt].tick_params(axis = 'x', labelbottom=True)
    plt.subplots_adjust(left = 0.075, right = 0.97, bottom = 0.14, top = 0.99, wspace = 0.225, hspace = 0.05)
    pdf_pages.savefig()  
    
    
    #log(R23) vs log(Temperature)
    fig2, ax2 = plt.subplots()
    ax2.scatter(np.log10(T_e[detect]), logR23[detect], marker = '.')
    ax2.scatter(np.log10(T_e[non_detect]), logR23[non_detect], marker = '<')
    if err_bars == True:
        ax2.errorbar(np.log10(T_e[detect]), logR23[detect], xerr = err_dict['T_e_lowhigh_error'], fmt = '.')
    ax2.set_xlabel('log($T_e$) (K)')
    ax2.set_ylabel('log($R_{23}$)')
    ax2.set_title('$R_{23}$ vs. Temperature')
    pdf_pages.savefig()
    
     
    #log(O32) vs log(Temperature)
    fig3, ax3 = plt.subplots()
    ax3.scatter(np.log10(T_e[detect]), logO32[detect], marker = '.')
    ax3.scatter(np.log10(T_e[non_detect]), logO32[non_detect], marker = '<')
    if err_bars == True:
        ax3.errorbar(np.log10(T_e[detect]), logO32[detect], xerr = err_dict['T_e_lowhigh_error'], fmt = '.')
    ax3.set_xlabel('log($T_e$) (K)')
    ax3.set_ylabel('log($O_{32}$)')
    ax3.set_title('$O_{32}$ vs. Temperature')
    pdf_pages.savefig()
    
    
    #log(Composite Metallicity) vs log(R23)
    fig4, ax4 = plt.subplots()
    ax4.scatter(logR23[detect], com_O_log[detect], marker = '.')
    ax4.scatter(logR23[non_detect], com_O_log[non_detect], marker = '^')
    if err_bars == True:
        ax4.errorbar(logR23[detect], com_O_log[detect], yerr = err_dict['12+log(O/H)_lowhigh_error'], fmt = '.')
    ax4.set_xlabel('log($R_{23}$)')
    ax4.set_ylabel('12+log(O/H) $T_e$')
    ax4.set_title('Composite Metallicity vs. $R_{23}$')
    pdf_pages.savefig()
    
    
    #log(Composite Metallicity) vs log(O32)
    fig5, ax5 = plt.subplots()
    ax5.scatter(logO32[detect], com_O_log[detect], marker = '.')
    ax5.scatter(logO32[non_detect], com_O_log[non_detect], marker = '^')
    if err_bars == True:
        ax5.errorbar(logO32[detect], com_O_log[detect], yerr = err_dict['12+log(O/H)_lowhigh_error'], fmt = '.')
    ax5.set_xlabel('log($O_{32}$)')
    ax5.set_ylabel('12+log(O/H) $T_e$')
    ax5.set_title('Composite Metallicity vs. $O_{32}$')
    pdf_pages.savefig()
    
    
    #log(Temperature) vs log(Avg Mass)
    fig6, ax6 = plt.subplots()
    ax6.scatter(logM_avg[detect], np.log10(T_e[detect]), marker = '.')
    ax6.scatter(logM_avg[non_detect], np.log10(T_e[non_detect]), marker = 'v')
    if err_bars == True:
        ax6.errorbar(logM_avg[detect], np.log10(T_e[detect]), yerr = err_dict['T_e_lowhigh_error'], fmt = '.')
    ax6.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax6.set_ylabel('log($T_e$) (K)')
    ax6.set_title('Temperature vs. Avg Mass')
    pdf_pages.savefig()
    
    
    #log(Composite Metallicity) vs log(Avg Mass)
    fig7, ax7 = plt.subplots()
    ax7.scatter(logM_avg[detect], com_O_log[detect], marker = '.')
    ax7.scatter(logM_avg[non_detect], com_O_log[non_detect], marker = '^')
    if err_bars == True:
        ax7.errorbar(logM_avg[detect], com_O_log[detect], yerr = err_dict['12+log(O/H)_lowhigh_error'], fmt = '.')
    mass_range = np.arange(8.2, 9.9, 0.05)
    AM_relation = 8.798 - np.log10(1 + ((10**8.901)/(10**mass_range))**0.640)
    ax7.plot(mass_range, AM_relation, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    ax7.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax7.set_ylabel('12+log(O/H) $T_e$')
    ax7.set_title('Composite Metallicity vs. Avg Mass')
    pdf_pages.savefig()
    
    
    pdf_pages.close()