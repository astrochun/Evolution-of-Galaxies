import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit 
from Metallicity_Stack_Commons import OIII_r
from Metallicity_Stack_Commons.column_names import indv_names0, temp_metal_names0, bin_ratios0, bin_mzevolve_names0, bin_names0, filename_dict
      
    

def bin_derived_props_plots(metal_file, em_file, bin_file, valid_file, out_file, hbeta_bin = False):
    '''
    Purpose: 
        This function creates plots for the stacked spectra (bins): OII/HBeta vs Stellar Mass,
        OIII/HBeta vs Stellar Mass, OIII4363/OIII5007 vs Stellar Mass, R23 vs Temperature, O32 vs Temperature,
        Metallicity vs R23, Metallicity vs O32, Temperature vs Stellar Mass, and Metallicity vs Stellar Mass.
        
    Usage:
        plots.bin_derived_props_plots(metal_file, em_file, out_file)
        
    Params:
        metal_file --> file containing bin metallicities, electron temperatures, R23, and O32 values.
        em_file --> file containing each bin's average stellar mass, detection marking, OII observed flux,
            OIII4363 observed flux, OIII5007 observed flux, and HBeta observed flux values.
        
    Returns:
        None
        
    Outputs:
        out_file --> a pdf containing plots for the stacked spectra: R23 vs Mass, O32 vs Mass,
            OII/HBeta vs Mass, OIII/HBeta vs Mass, OIII4363/OIII5007 vs Mass, R23 vs Temperature,
            O32 vs Temperature, Metallicity vs R23, Metallicity vs O32, Temperature vs Mass, 
            and Metallicity vs Mass.
    '''
    
    metal_table = asc.read(metal_file)
    em_table = asc.read(em_file)
    bin_table = asc.read(bin_file)
    valid_table = asc.read(valid_file)
    pdf_pages = PdfPages(out_file)
    
    com_O_log = metal_table[temp_metal_names0[1]].data
    T_e = metal_table[temp_metal_names0[0]].data
    R23_composite = metal_table[bin_ratios0[0]].data
    O32_composite = metal_table[bin_ratios0[1]].data
    
    avg_mass = bin_table[bin_mzevolve_names0[2]].data
    detection = valid_table[bin_names0[2]].data
    OII = em_table['OII_3727_Flux_Observed'].data
    OIII4363 = em_table['OIII_4363_Flux_Observed'].data
    OIII5007 = em_table['OIII_5007_Flux_Observed'].data
    HBeta = em_table['HBETA_Flux_Observed'].data
    
    detect = np.where(detection == 1)[0]
    non_detect = np.where(detection == 0.5)[0]   #non-detection with reliable 5007
        
    #Line ratios        
    OII_HBeta = OII / HBeta
    OIII_HBeta = (OIII5007 * (1 + 1/OIII_r)) / HBeta            
    OIII_ratio = OIII4363 / (OIII5007 * (1 + 1/OIII_r))            
    
    #Line Ratios vs Mass
    fig1, ax1 = plt.subplots(2, 3, sharex = True)
    
    if hbeta_bin == True:
        lum_avg = bin_table[bin_mzevolve_names0[6]].data
        for ii in range(len(detect)):
            if detect[ii] % 2 == 0:
                ax1[0, 0].scatter(avg_mass[detect[ii]], R23_composite[detect[ii]], marker = '.', color = 'cyan')
                ax1[0, 1].scatter(avg_mass[detect[ii]], O32_composite[detect[ii]], marker = '.', color = 'cyan')
                ax1[0, 2].scatter(avg_mass[detect[ii]], np.log10(OII_HBeta[detect[ii]]), marker = '.', color = 'cyan')
                ax1[1, 0].scatter(avg_mass[detect[ii]], np.log10(OIII_HBeta[detect[ii]]), marker = '.', color = 'cyan')
                ax1[1, 1].scatter(avg_mass[detect[ii]], np.log10(OIII_ratio[detect[ii]]), marker = '.', color = 'cyan')
                ax1[1, 2].scatter(avg_mass[detect[ii]], lum_avg[detect[ii]], marker = '.', color = 'cyan')
            else:
                ax1[0, 0].scatter(avg_mass[detect[ii]], R23_composite[detect[ii]], marker = '.', color = 'blue')
                ax1[0, 1].scatter(avg_mass[detect[ii]], O32_composite[detect[ii]], marker = '.', color = 'blue')
                ax1[0, 2].scatter(avg_mass[detect[ii]], np.log10(OII_HBeta[detect[ii]]), marker = '.', color = 'blue')
                ax1[1, 0].scatter(avg_mass[detect[ii]], np.log10(OIII_HBeta[detect[ii]]), marker = '.', color = 'blue')
                ax1[1, 1].scatter(avg_mass[detect[ii]], np.log10(OIII_ratio[detect[ii]]), marker = '.', color = 'blue')
                ax1[1, 2].scatter(avg_mass[detect[ii]], lum_avg[detect[ii]], marker = '.', color = 'blue')
        for ii in range(len(non_detect)):
            if non_detect[ii] % 2 == 0:
                ax1[0, 0].scatter(avg_mass[non_detect[ii]], R23_composite[non_detect[ii]], marker = '^', color = 'pink')
                ax1[0, 1].scatter(avg_mass[non_detect[ii]], O32_composite[non_detect[ii]], marker = '^', color = 'pink')
                ax1[0, 2].scatter(avg_mass[non_detect[ii]], np.log10(OII_HBeta[non_detect[ii]]), marker = '^', color = 'pink')
                ax1[1, 0].scatter(avg_mass[non_detect[ii]], np.log10(OIII_HBeta[non_detect[ii]]), marker = '^', color = 'pink')
                ax1[1, 1].scatter(avg_mass[non_detect[ii]], np.log10(OIII_ratio[non_detect[ii]]), marker = '^', color = 'pink')
                ax1[1, 2].scatter(avg_mass[non_detect[ii]], lum_avg[non_detect[ii]], marker = '^', color = 'pink')
            else:
                ax1[0, 0].scatter(avg_mass[non_detect[ii]], R23_composite[non_detect[ii]], marker = '^', color = 'red')
                ax1[0, 1].scatter(avg_mass[non_detect[ii]], O32_composite[non_detect[ii]], marker = '^', color = 'red')
                ax1[0, 2].scatter(avg_mass[non_detect[ii]], np.log10(OII_HBeta[non_detect[ii]]), marker = '^', color = 'red')
                ax1[1, 0].scatter(avg_mass[non_detect[ii]], np.log10(OIII_HBeta[non_detect[ii]]), marker = '^', color = 'red')
                ax1[1, 1].scatter(avg_mass[non_detect[ii]], np.log10(OIII_ratio[non_detect[ii]]), marker = '^', color = 'red')
                ax1[1, 2].scatter(avg_mass[non_detect[ii]], lum_avg[non_detect[ii]], marker = '^', color = 'red')
    else:
        ax1[0, 0].scatter(avg_mass[detect], R23_composite[detect], label = 'R_23')
        ax1[0, 0].scatter(avg_mass[non_detect], R23_composite[non_detect], marker = '^', label = 'R_23')
        ax1[0, 1].scatter(avg_mass[detect], O32_composite[detect], label = 'O_32')
        ax1[0, 1].scatter(avg_mass[non_detect], O32_composite[non_detect], marker = '^', label = 'O_32')
        ax1[0, 2].scatter(avg_mass[detect], np.log10(OII_HBeta[detect]), label = 'OII/HBeta')
        ax1[0, 2].scatter(avg_mass[non_detect], np.log10(OII_HBeta[non_detect]), marker = '^', label = 'OII/HBeta')
        ax1[1, 0].scatter(avg_mass[detect], np.log10(OIII_HBeta[detect]), label = 'OIII/HBeta')
        ax1[1, 0].scatter(avg_mass[non_detect], np.log10(OIII_HBeta[non_detect]), marker = '^', label = 'OIII/HBeta')
        ax1[1, 1].scatter(avg_mass[detect], np.log10(OIII_ratio[detect]), label = '4363/(5007 * (1 + 1/3.1))')
        ax1[1, 1].scatter(avg_mass[non_detect], np.log10(OIII_ratio[non_detect]), marker = '^',  label = '4363/(5007 * (1 + 1/3.1))')
        
        
    ax1[0, 0].annotate('R23', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')
    ax1[0, 1].annotate('O32', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')
    ax1[0, 2].annotate('OII/HBeta', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')
    ax1[1, 0].annotate('OIII/HBeta', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')
    ax1[1, 1].annotate('4363/(5007 * (1 + 1/3.1))', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')
    ax1[1, 2].annotate('HBeta Luminosity', [0.05, 0.97], xycoords = 'axes fraction', va = 'top', ha = 'left', fontsize = '7')  
        
    ax1[1, 0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1[1, 1].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1[1, 2].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')     
    for t_ax in ax1:
        for tt in range(len(t_ax)):
            t_ax[tt].tick_params(axis = 'x', labelbottom=True)
    plt.subplots_adjust(left = 0.075, right = 0.97, bottom = 0.075, top = 0.99, wspace = 0.225, hspace = 0.05)
    pdf_pages.savefig()    
    
    #R23 vs Temperature
    fig2, ax2 = plt.subplots()
    ax2.scatter(np.log10(T_e[detect]), R23_composite[detect], marker = '.')
    ax2.scatter(np.log10(T_e[non_detect]), R23_composite[non_detect], marker = '<')
    ax2.set_xlabel('$T_e$ (K)')
    ax2.set_ylabel('$R_{23}$')
    ax2.set_title('$R_{23}$ vs. Temperatures')
    pdf_pages.savefig()
     
    #O32 vs Temperature
    fig3, ax3 = plt.subplots()
    ax3.scatter(np.log10(T_e[detect]), O32_composite[detect], marker = '.')
    ax3.scatter(np.log10(T_e[non_detect]), O32_composite[non_detect], marker = '<')
    ax3.set_xlabel('$T_e$ (K)')
    ax3.set_ylabel('$O_{32}$')
    ax3.set_title('$O_{32}$ vs. Temperatures')
    pdf_pages.savefig()
    
    #Metallicity vs R23
    fig4, ax4 = plt.subplots()
    ax4.scatter(R23_composite[detect], com_O_log[detect], marker = '.')
    ax4.scatter(R23_composite[non_detect], com_O_log[non_detect], marker = '^')
    ax4.set_xlabel('$R_{23}$')
    ax4.set_ylabel('12+log(O/H) $T_e$')
    ax4.set_title('Composite Metallicity vs. $R_{23}$')
    pdf_pages.savefig()
    
    #Metallicity vs O32
    fig5, ax5 = plt.subplots()
    ax5.scatter(O32_composite[detect], com_O_log[detect], marker = '.')
    ax5.scatter(O32_composite[non_detect], com_O_log[non_detect], marker = '^')
    ax5.set_xlabel('$O_{32}$')
    ax5.set_ylabel('12+log(O/H) $T_e$')
    ax5.set_title('Composite Metallicity vs. O32')
    pdf_pages.savefig()
    
    #Temperature vs Mass
    fig6, ax6 = plt.subplots()
    ax6.scatter(avg_mass[detect], np.log10(T_e[detect]), marker = '.')
    ax6.scatter(avg_mass[non_detect], np.log10(T_e[non_detect]), marker = 'v')
    ax6.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax6.set_ylabel('$T_e$ (K)')
    ax6.set_title('Temperatures vs. Avg Mass')
    pdf_pages.savefig()
    
    #Metallicity vs Mass
    fig7, ax7 = plt.subplots()
    ax7.scatter(avg_mass[detect], com_O_log[detect], marker = '.')
    ax7.scatter(avg_mass[non_detect], com_O_log[non_detect], marker = '^')
    mass = np.arange(8.2, 9.9, 0.05)
    y = 8.798 - np.log10(1 + ((10**8.901)/(10**mass))**0.640)
    ax7.plot(mass, y, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    ax7.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax7.set_ylabel('12+log(O/H) $T_e$')
    ax7.set_title('Composite Metallicity vs. Avg Mass')
    pdf_pages.savefig()
    
    pdf_pages.close()
    





def indiv_derived_props_plots(fitspath, mass_indiv_bin_file, LHb_indiv_bin_file, mass_indiv_props_file, 
                              mass_indiv_metal_file, LHb_indiv_props_file, LHb_indiv_metal_file, 
                              mass_valid_file, LHb_valid_file, mass_bin_file, LHb_bin_file, 
                              mass_bin_metal_file, LHb_bin_metal_file, MTO, restrict_MTO = False):    
    '''
    REDO DOCUMENTATION
    Purpose:
        This function creates plots for the individual measurements within each bin: HBeta Luminosity vs
        Stellar Mass (color Metallicity), O32 vs R23 (color Metallicity), Metallicity vs R23 and O32,
        Metallicity vs Stellar Mass. (Each plot type has two plots - one for Mass Bins and one for Mass-LHBeta Bins)
         
        
    Usage:
        plots.indiv_derived_props_plots(fitspath, metal_Te_file, mass_bin_file, HB_bin_file, mass_metal_file,
                                        HB_metal_file, MTO, restrict_MTO)
        
    Params:
        fitspath --> a string of the file path where the input file is and where the output file will be placed.
        metal_Te_file --> table containing the calculations for individual metallicities, R23, O32, and 
            detection markings.
        mass_bin_file --> table containing the average mass for each mass bin.
        HB_bin_file --> table containing the average mass for each mass-LHbeta bin. 
        mass_metal_file --> table containing the detection markings and bin metallicities for each mass bin.
        HB_metal_file --> table containing the detection markings and bin metallicities for each mass-LHbeta bin.
        MTO --> a string that indicates in the pdf name whether or not the MTO parameter is held constant when
            a curve fit is applied to the Metallicity vs Mass plots.
        restrict_MTO (OPTIONAL) --> if you want to hold the MTO parameter constant for the curve fit of
            Metallicity vs Mass, then restrict_MTO = True. Its default value is False.
        
    Returns:
        None
        
    Outputs:
        fitspath + 'individual_metal_plots' + MTO + '.pdf' --> a pdf containing plots for the individual 
            measurements: HBeta Luminosity vs Stellar Mass (color Metallicity), O32 vs R23 (color Metallicity),
            Metallicity vs R23 and O32, Metallicity vs Stellar Mass.
    '''
    
    mass_indiv_MT_tbl = asc.read(mass_indiv_metal_file)
    mass_indiv_props_tbl = asc.read(mass_indiv_props_file)
    LHb_indiv_MT_tbl = asc.read(LHb_indiv_metal_file)
    LHb_indiv_props_tbl = asc.read(LHb_indiv_props_file)
    mass_indiv_bin_tbl = asc.read(mass_indiv_bin_file)
    LHb_indiv_bin_tbl = asc.read(LHb_indiv_bin_file)
    
    mass_valid_tbl = asc.read(mass_valid_file)
    LHb_valid_tbl = asc.read(LHb_valid_file)
    mass_bin_tbl = asc.read(mass_bin_file)
    LHb_bin_tbl = asc.read(LHb_bin_file)
    mass_metal_tbl = asc.read(mass_bin_metal_file)
    LHb_metal_tbl = asc.read(LHb_bin_metal_file)


    ##individual galaxy data, i.e. comes from MT_ascii
    Mbin_log_mass = mass_indiv_props_tbl[indv_names0[3]].data
    Mbin_log_LHb = mass_indiv_props_tbl[indv_names0[4]].data
    LHbbin_log_mass = LHb_indiv_props_tbl[indv_names0[3]].data
    LHbbin_log_LHb = LHb_indiv_props_tbl[indv_names0[4]].data
    
    Mbin_indiv_metal = mass_indiv_MT_tbl[temp_metal_names0[1]].data
    LHbbin_indiv_metal = LHb_indiv_MT_tbl[temp_metal_names0[1]].data
    Mbin_indiv_R23 = mass_indiv_MT_tbl[indv_names0[1]].data
    Mbin_indiv_O32 = mass_indiv_MT_tbl[indv_names0[2]].data
    LHbbin_indiv_R23 = LHb_indiv_MT_tbl[indv_names0[1]].data
    LHbbin_indiv_O32 = LHb_indiv_MT_tbl[indv_names0[2]].data
    Mbin_indiv_two_beta = mass_indiv_MT_tbl[indv_names0[5]].data
    Mbin_indiv_three_beta = mass_indiv_MT_tbl[indv_names0[6]].data
    LHbbin_indiv_two_beta = LHb_indiv_MT_tbl[indv_names0[5]].data
    LHbbin_indiv_three_beta = LHb_indiv_MT_tbl[indv_names0[6]].data
    
    mass_indiv_detect = mass_indiv_bin_tbl[bin_names0[2]].data
    LHb_indiv_detect = LHb_indiv_bin_tbl[bin_names0[2]].data 
    
    
    ##bin data, i.e. comes from either mass or massLHb specific files
    mass_avg_mass = mass_bin_tbl[bin_mzevolve_names0[2]].data
    LHb_avg_mass = LHb_bin_tbl[bin_mzevolve_names0[2]].data 
    mass_detect_col = mass_valid_tbl[bin_names0[2]].data
    LHb_detect_col = LHb_valid_tbl[bin_names0[2]].data        
    mass_metal = mass_metal_tbl[temp_metal_names0[1]].data
    LHb_metal = LHb_metal_tbl[temp_metal_names0[1]].data 
    
    
    ##detection determinations  
    #bins
    mass_detect = np.where(mass_detect_col == 1.0)[0]
    mass_nondetect = np.where(mass_detect_col == 0.5)[0]
    LHb_detect = np.where(LHb_detect_col == 1.0)[0]
    LHb_nondetect = np.where(LHb_detect_col == 0.5)[0]
    
    #individual    
    
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



    ###REDO by excluding those with bad lines
    mass_ind_detect = np.where(mass_indiv_detect == 1.0)[0]
    mass_ind_nondetect = np.where(mass_indiv_detect == 0.5)[0]
    LHb_ind_detect = np.where(LHb_indiv_detect == 1.0)[0]
    LHb_ind_nondetect = np.where(LHb_indiv_detect == 0.5)[0]
    ###

    if MTO == True:
        pdf_pages = PdfPages(fitspath + filename_dict['indv_derived_prop'].replace('.tbl', 'constMTO.pdf'))
    else:
        pdf_pages = PdfPages(fitspath + filename_dict['indv_derived_prop'].replace('.tbl', '.pdf'))
    
    
    ##HBeta Luminosity vs Mass ColorMap=Metallicity##
    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot1 = ax1.scatter(Mbin_log_mass[mass_ind_detect], Mbin_log_LHb[mass_ind_detect], 5.0,
                        c=Mbin_indiv_metal[mass_ind_detect], marker='*')
    plot1 = ax1.scatter(Mbin_log_mass[mass_ind_nondetect], Mbin_log_LHb[mass_ind_nondetect], 5.0, 
                        facecolors = 'None', c=Mbin_indiv_metal[mass_ind_nondetect], marker='^')
    cb = fig1.colorbar(plot1)
    cb.set_label('Metallicity')
    ax1.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1.set_ylabel('log(H$\\beta$ Luminosity)')
    ax1.set_title('$M_\star$ Bins: H$\\beta$ Luminosity vs. Mass Colormap=Metallicity')
    fig1.set_size_inches(8, 8)
    fig1.savefig(pdf_pages, format='pdf')
    
    fig2, ax2 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot2 = ax2.scatter(LHbbin_log_mass[LHb_ind_detect], LHbbin_log_LHb[LHb_ind_detect], 5.0, 
                        c=LHbbin_indiv_metal[LHb_ind_detect], marker='*')
    plot2 = ax2.scatter(LHbbin_log_mass[LHb_ind_nondetect], LHbbin_log_LHb[LHb_ind_nondetect], 5.0, 
                        facecolors = 'None', c=LHbbin_indiv_metal[LHb_ind_nondetect], marker='^')
    cb = fig2.colorbar(plot2)
    cb.set_label('Metallicity')
    ax2.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax2.set_ylabel('log(H$\\beta$ Luminosity)')
    ax2.set_title('$M_\star$-LH$\\beta$ Bins: H$\\beta$ Luminosity vs. Mass Colormap=Metallicity')
    fig2.set_size_inches(8, 8)
    fig2.savefig(pdf_pages, format='pdf')
    
    
    
    
    ##O32 vs R23 ColorMap=Metallicity##
    fig3, ax3 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot3 = ax3.scatter(Mbin_indiv_R23[mass_ind_detect], Mbin_indiv_O32[mass_ind_detect], 5.0, 
                        c=Mbin_indiv_metal[mass_ind_detect], marker='*')
    plot3 = ax3.scatter(Mbin_indiv_R23[mass_ind_nondetect], Mbin_indiv_O32[mass_ind_nondetect], 5.0, 
                        facecolors = 'None', c=Mbin_indiv_metal[mass_ind_nondetect], marker='^')
    cb = fig3.colorbar(plot3)
    cb.set_label('Metallicity')
    ax3.set_xlabel('$R_{23}$')
    ax3.set_ylabel('$O_{32}$')
    ax3.set_title('$M_\star$ Bins: $O_{32}$ vs. $R_{23}$ Colormap=Metallicity')
    fig3.set_size_inches(8, 8)
    fig3.savefig(pdf_pages, format='pdf')
    
    fig4, ax4 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot4 = ax4.scatter(LHbbin_indiv_R23[LHb_ind_detect], LHbbin_indiv_O32[LHb_ind_detect], 5.0, 
                        c=LHbbin_indiv_metal[LHb_ind_detect], marker='*')
    plot4 = ax4.scatter(LHbbin_indiv_R23[LHb_ind_nondetect], LHbbin_indiv_O32[LHb_ind_nondetect], 5.0, 
                        facecolors = 'None', c=LHbbin_indiv_metal[LHb_ind_nondetect], marker='^')
    cb = fig4.colorbar(plot4)
    cb.set_label('Metallicity')
    ax4.set_xlabel('$R_{23}$')
    ax4.set_ylabel('$O_{32}$')
    ax4.set_title('$M_\star$-LH$\\beta$ Bins: $O_{32}$ vs. $R_{23}$ Colormap=Metallicity')
    fig4.set_size_inches(8, 8)
    fig4.savefig(pdf_pages, format='pdf')
    
    
    
    
    ##Metallicity vs R23 and O32##
    fig5, ax5 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax5[0].scatter(Mbin_indiv_R23[mass_ind_detect], Mbin_indiv_metal[mass_ind_detect], facecolors = 'None',
                   edgecolors = 'blue', label = '$R_{23}$')
    ax5[1].scatter(Mbin_indiv_O32[mass_ind_detect], Mbin_indiv_metal[mass_ind_detect], facecolors = 'None', 
                   edgecolors = 'red', label = '$O_{32}$')
    ax5[0].scatter(Mbin_indiv_R23[mass_ind_nondetect], Mbin_indiv_metal[mass_ind_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'blue', label = '$R_{23}$', alpha = 0.5)
    ax5[1].scatter(Mbin_indiv_O32[mass_ind_nondetect], Mbin_indiv_metal[mass_ind_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'red', label = '$O_{32}$', alpha = 0.5)
    ax5[0].legend(loc = 'best')
    ax5[1].legend(loc = 'best')
    ax5[0].set_ylabel('Metallicity')
    ax5[0].set_title('$M_\star$ Bins: Metallicity vs. $R_{23}$')
    ax5[1].set_title('$M_\star$ Bins: Metallicity vs. $O_{32}$')
    fig5.set_size_inches(8, 8)
    fig5.savefig(pdf_pages, format='pdf')
    
    fig6, ax6 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax6[0].scatter(LHbbin_indiv_R23[LHb_ind_detect], LHbbin_indiv_metal[LHb_ind_detect], facecolors = 'None',
                   edgecolors = 'blue', label = '$R_{23}$')
    ax6[1].scatter(LHbbin_indiv_O32[LHb_ind_detect], LHbbin_indiv_metal[LHb_ind_detect], facecolors = 'None',
                   edgecolors = 'red', label = '$O_{32}$')
    ax6[0].scatter(LHbbin_indiv_R23[LHb_ind_nondetect], LHbbin_indiv_metal[LHb_ind_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'blue', label = '$R_{23}$', alpha = 0.5)
    ax6[1].scatter(LHbbin_indiv_O32[LHb_ind_nondetect], LHbbin_indiv_metal[LHb_ind_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'red', label = '$O_{32}$', alpha = 0.5)
    ax6[0].legend(loc = 'best')
    ax6[1].legend(loc = 'best')
    ax6[0].set_ylabel('Metallicity')
    ax6[0].set_title('$M_\star$-LH$\\beta$ Bins: Metallicity vs. $R_{23}$')
    ax6[1].set_title('$M_\star$-LH$\\beta$ Bins: Metallicity vs. $O_{32}$')
    fig6.set_size_inches(8, 8)
    fig6.savefig(pdf_pages, format='pdf')
    
    
    
    ##R23 and O32 vs Mass##
    fig7, ax7 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax7[0].scatter(Mbin_log_mass[mass_ind_detect], Mbin_indiv_R23[mass_ind_detect], facecolors = 'None',
                   edgecolors = 'blue', label = '$R_{23}$')    
    ax7[1].scatter(Mbin_log_mass[mass_ind_detect], Mbin_indiv_O32[mass_ind_detect], facecolors = 'None', 
                   edgecolors = 'red', label = '$O_{32}$')
    ax7[0].scatter(Mbin_log_mass[mass_ind_nondetect], Mbin_indiv_R23[mass_ind_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'blue', label = '$R_{23}$')
    ax7[1].scatter(Mbin_log_mass[mass_ind_nondetect], Mbin_indiv_O32[mass_ind_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'red', label = '$O_{32}$')
    ax7[0].legend(loc = 'best')
    ax7[1].legend(loc = 'best')
    ax7[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax7[0].set_title('$M_\star$ Bins: $R_{23}$ vs. Mass')
    ax7[1].set_title('$M_\star$ Bins: $O_{32}$ vs. Mass')
    fig7.set_size_inches(8, 8)
    fig7.savefig(pdf_pages, format='pdf')
    
    fig8, ax8 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax8[0].scatter(LHbbin_log_mass[LHb_ind_detect], LHbbin_indiv_R23[LHb_ind_detect], facecolors = 'None', 
                   edgecolors = 'blue', label = '$R_{23}$')    
    ax8[1].scatter(LHbbin_log_mass[LHb_ind_detect], LHbbin_indiv_O32[LHb_ind_detect], facecolors = 'None', 
                   edgecolors = 'red', label = '$O_{32}$')
    ax8[0].scatter(LHbbin_log_mass[LHb_ind_nondetect], LHbbin_indiv_R23[LHb_ind_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'blue', label = '$R_{23}$')
    ax8[1].scatter(LHbbin_log_mass[LHb_ind_nondetect], LHbbin_indiv_O32[LHb_ind_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'red', label = '$O_{32}$')
    ax8[0].legend(loc = 'best')
    ax8[1].legend(loc = 'best')
    ax8[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax8[0].set_title('$M_\star$-LH$\\beta$ Bins: $R_{23}$ vs. Mass')
    ax8[1].set_title('$M_\star$-LH$\\beta$ Bins: $O_{32}$ vs. Mass')
    fig8.set_size_inches(8, 8)
    fig8.savefig(pdf_pages, format='pdf')
    
    
    
    ##OII/HBeta and OIII/HBeta vs Mass##
    fig9, ax9 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax9[0].scatter(Mbin_log_mass[mass_ind_detect], np.log10(Mbin_indiv_two_beta[mass_ind_detect]), 
                   facecolors = 'None', edgecolors = 'cyan', label = '$\\frac{OII}{H\\beta}$')    
    ax9[1].scatter(Mbin_log_mass[mass_ind_detect], np.log10(Mbin_indiv_three_beta[mass_ind_detect]), 
                   facecolors = 'None', edgecolors = 'orange', label = '$\\frac{OIII}{H\\beta}$')
    ax9[0].scatter(Mbin_log_mass[mass_ind_nondetect], np.log10(Mbin_indiv_two_beta[mass_ind_nondetect]), 
                   marker='^', facecolors = 'None', edgecolors = 'cyan', label = '$\\frac{OII}{H\\beta}$')
    ax9[1].scatter(Mbin_log_mass[mass_ind_nondetect], np.log10(Mbin_indiv_three_beta[mass_ind_nondetect]), 
                   marker='^', facecolors = 'None', edgecolors = 'orange', label = '$\\frac{OIII}{H\\beta}$')
    ax9[0].legend(loc = 'best')
    ax9[1].legend(loc = 'best')
    ax9[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax9[0].set_title('$M_\star$ Bins: $\\frac{OII}{H\\beta}$ vs. Mass')
    ax9[1].set_title('$M_\star$ Bins: $\\frac{OIII}{H\\beta}$ vs. Mass')
    fig9.set_size_inches(8, 8)
    fig9.savefig(pdf_pages, format='pdf')
    
    fig10, ax10 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax10[0].scatter(LHbbin_log_mass[LHb_ind_detect], np.log10(LHbbin_indiv_two_beta[LHb_ind_detect]), 
                    facecolors = 'None', edgecolors = 'cyan', label = '$\\frac{OII}{H\\beta}$')    
    ax10[1].scatter(LHbbin_log_mass[LHb_ind_detect], np.log10(LHbbin_indiv_three_beta[LHb_ind_detect]), 
                    facecolors = 'None', edgecolors = 'orange', label = '$\\frac{OIII}{H\\beta}$')
    ax10[0].scatter(LHbbin_log_mass[LHb_ind_nondetect], np.log10(LHbbin_indiv_two_beta[LHb_ind_nondetect]),
                    marker='^', facecolors = 'None', edgecolors = 'cyan', label = '$\\frac{OII}{H\\beta}$')
    ax10[1].scatter(LHbbin_log_mass[LHb_ind_nondetect], np.log10(LHbbin_indiv_three_beta[LHb_ind_nondetect]),
                    marker='^', facecolors = 'None', edgecolors = 'orange', label = '$\\frac{OIII}{H\\beta}$')
    ax10[0].legend(loc = 'best')
    ax10[1].legend(loc = 'best')
    ax10[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax10[0].set_title('$M_\star$-LH$\\beta$ Bins: $\\frac{OII}{H\\beta}$ vs. Mass')
    ax10[1].set_title('$M_\star$-LH$\\beta$ Bins: $\\frac{OIII}{H\\beta}$ vs. Mass')
    fig10.set_size_inches(8, 8)
    fig10.savefig(pdf_pages, format='pdf')
    
    
    
    ##Metallcity vs Mass##
    fig11, ax11 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.98, wspace = 0.0)
    
    
    #Andrews&Martini fit
    mass = np.arange(7.5, 10, 0.05)
    y = 8.798 - np.log10(1 + ((10**8.901)/(10**mass))**0.640)
    ax11[0].plot(mass, y, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    ax11[1].plot(mass, y, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    
    
    #Individual detections and non-detections
    ax11[0].scatter(Mbin_log_mass[mass_ind_detect], Mbin_indiv_metal[mass_ind_detect], s = 15, 
                    facecolors = 'None', edgecolors = 'blue', label = 'Individual Detections')
    ax11[1].scatter(LHbbin_log_mass[LHb_ind_detect], LHbbin_indiv_metal[LHb_ind_detect], s = 15, 
                    facecolors = 'None', edgecolors = 'blue', label = 'Individual Detections')
    ax11[0].scatter(Mbin_log_mass[mass_ind_nondetect], Mbin_indiv_metal[mass_ind_nondetect], s = 15, 
                    marker='^', facecolors = 'None', edgecolors = 'blue', label = 'Individual Non-Detections', 
                    alpha = 0.5)
    ax11[1].scatter(LHbbin_log_mass[LHb_ind_nondetect], LHbbin_indiv_metal[LHb_ind_nondetect], s = 15, 
                    marker='^', facecolors = 'None', edgecolors = 'blue', label = 'Individual Non-Detections',
                    alpha = 0.5)
    
    print('Number of mass bin individual sources plotted:', len(Mbin_log_mass[mass_ind_detect]))
    print('Number of mass-LHbeta bin individual sources plotted:', len(LHbbin_log_mass[LHb_ind_detect]))
    
    #Mass bin detections and non-detections
    ax11[0].scatter(mass_avg_mass[mass_detect], mass_metal[mass_detect], s = 25, color = 'red',
                   label = 'Bin Detections')
    ax11[0].scatter(mass_avg_mass[mass_nondetect], mass_metal[mass_nondetect], s = 25, color = 'red',
                   marker = '^', label = 'Bin Non-Detections', alpha = 0.5)
    #HBeta bin detections and non-detections
    ax11[1].scatter(LHb_avg_mass[LHb_detect], LHb_metal[LHb_detect], s = 25, color = 'red',
                   label = 'Bin Detections')
    ax11[1].scatter(LHb_avg_mass[LHb_nondetect], LHb_metal[LHb_nondetect], s = 25, color = 'red',
                   marker = '^', label = 'Bin Non-Detections', alpha = 0.5)
    
    
    ##Curve fit     
    fail = False
    if restrict_MTO == False:
        p0 = [8.798, 8.901, 0.640]
        para_bounds = ((8.0, 8.0, 0.0), (9.0, 9.5, 1.0))
        fit = mass_metal_fit
    else:
        p0 = [8.798, 0.640]
        para_bounds = ((8.0, 0.0), (9.0, 1.0))
        fit = mass_metal_fit2
        
    try:
        o11, o21 = curve_fit(fit, mass_avg_mass[mass_detect], mass_metal[mass_detect], p0 = p0, bounds = para_bounds)
        print(o11)
    except ValueError:
        print('fail Mass Bins')
        fail = True
    try:
        o12, o22 = curve_fit(fit, LHb_avg_mass[LHb_detect], LHb_metal[LHb_detect], p0 = p0, bounds = para_bounds)
        print(o12)
    except ValueError:
        print('fail HBeta Bins')
        fail = True
        
    if not fail:
        if restrict_MTO == False:
            ax11[0].plot(mass, mass_metal_fit(mass, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
            ax11[1].plot(mass, mass_metal_fit(mass, *o12), alpha = 0.5, color = 'red', label = 'Our Fit')
        else:
            ax11[0].plot(mass, mass_metal_fit2(mass, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
            ax11[1].plot(mass, mass_metal_fit2(mass, *o12), alpha = 0.5, color = 'red', label = 'Our Fit')
          
    ax11[0].legend(title = '$M_\star$ Bins', fontsize = 5, loc = 'upper left')    
    ax11[1].legend(title = '$M_\star$-LH$\\beta$ Bins', fontsize = 5, loc = 'upper left')
    ax11[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax11[0].set_ylabel('12+log(O/H) $T_e$')
    ax11[1].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
   
        
    pdf_pages.savefig()
    
    pdf_pages.close()
    
    
    
def mass_metal_fit(mass, a, b, g):
    '''
    MTO is NOT held constant.
    
    Andrews & Martini Mass-Metallicity Relation:
    8.798 - np.log10(1 + ((10**8.901)/(10**mass))**0.640)
    '''
    
    return a - np.log10(1 + ((10**b)/(10**mass))**g)   



def mass_metal_fit2(mass, a, g):
    '''
    MTO is held constant.
    
    Andrews & Martini Mass-Metallicity Relation:
    8.798 - np.log10(1 + ((10**8.901)/(10**mass))**0.640)
    '''
    
    return a - np.log10(1 + ((10**8.901)/(10**mass))**g) 

 



def lum_vs_mass(fitspath, binning_file):
    '''
    Purpose:
        This function creates a plot of HBeta Luminosity vs Stellar Mass for individual galaxy values.
        
    Usage:
        plots.lum_vs_mass(fitspath, binning_file)
        
    Params:
        fitspath --> fitspath --> a string of the file path where the input file is and where the output file
            will be placed.
        binning_file --> npz file containing the mass, HBeta luminosity, and bin edge values for each bin.
            (Valid for either binning type.)
        
    Returns:
        None
        
    Outputs:
        fitspath + 'lum_vs_mass.pdf' --> a pdf containing a HBeta Luminosity vs Stellar Mass plot with each
            bin's values separated by dashed lines to show the galaxies that fall within each bin.
    '''
    
    binning = np.load(fitspath + binning_file)
    lum = binning['lum']
    mass = binning['mass']
    bin_idx = binning['bin_ind']
    bin_low = binning['lowest_hbeta']
    bin_high = binning['highest_hbeta']
    mass_max = binning['bin_redge']
    mass_min = binning['bin_edge']
    
    pdf_pages = PdfPages(fitspath + 'lum_vs_mass.pdf')
    
    
    idx = []
    #odd-numbered bins are 
    low_bins_high_lum = []           #odd numbered bins
    low_bins_low_lum = []
    high_bins_high_lum = []          #even numbered bins
    high_bins_low_lum = []
    low_bins_high_mass = []           #odd numbered bins
    low_bins_low_mass = []
    high_bins_high_mass = []          #even numbered bins
    high_bins_low_mass = []
    for ii in range(len(bin_idx)):
        for jj in range(len(bin_idx[ii])):
            idx.append(bin_idx[ii][jj])
        if (ii + 1) % 2 == 0:
            high_bins_high_lum.append(bin_high[ii])
            high_bins_low_lum.append(bin_low[ii])
            high_bins_high_mass.append(mass_max[ii])
            high_bins_low_mass.append(mass_min[ii])
        else:
            low_bins_high_lum.append(bin_high[ii])
            low_bins_low_lum.append(bin_low[ii])
            low_bins_high_mass.append(mass_max[ii])
            low_bins_low_mass.append(mass_min[ii])


    fig1, ax1 = plt.subplots()
    ax1.scatter(np.log10(mass[idx]), lum[idx], color = 'orange', s = 0.5, marker='.')
    ax1.scatter(low_bins_low_mass, low_bins_low_lum, color = 'black', s = 15, marker='.',
                label = 'Lowest HBeta Luminosity (Lower Bins)')
    ax1.scatter(low_bins_high_mass, low_bins_high_lum, color = 'green', s = 15, marker='.',
                label = 'Median HBeta Luminosity')
    #ax1.scatter(high_bins_low_mass, high_bins_low_lum, color = 'green', s = 15, marker='.')
    ax1.scatter(high_bins_high_mass, high_bins_high_lum, color = 'blue', s = 15, marker='.',
                label = 'Highest HBeta Luminosity (Upper Bins)')
    
    for ii in range(len(low_bins_low_mass)):
        ax1.axvline(x=low_bins_low_mass[ii], linestyle='--', color = 'black', linewidth = 0.5)
        ax1.axvline(x=high_bins_high_mass[ii], linestyle='--', color = 'black', linewidth = 0.5)
    
        if ii == (len(low_bins_low_mass) - 1): 
            ax1.hlines(y=low_bins_high_lum[ii], xmin = low_bins_low_mass[ii], xmax = high_bins_high_mass[ii],
                       linestyle='--', color = 'black', linewidth = 0.5)
        else:
            ax1.hlines(y=low_bins_high_lum[ii], xmin = low_bins_low_mass[ii], xmax = low_bins_low_mass[ii+1],
                       linestyle='--', color = 'black', linewidth = 0.5)
    ax1.set_ylim(36,43)
    ax1.set_xlim(6,12.5)
    ax1.set_xlabel("log($\\frac{M_\star}{M_{\odot}}$)")
    ax1.set_ylabel("log(H$\\beta$ Luminosity)")
    ax1.legend(loc = 'lower left', fontsize = 'xx-small')
    
    pdf_pages.savefig()
    
    pdf_pages.close()





def metallicity_vs_R23_O32():
    fitspath = 'C:/Users/carol/Google Drive/MZEvolve/'
    
    metal_file = fitspath + 'individual/11262019/individual_derived_properties_metallicity.tbl'
    mass_bin_file = fitspath + 'massbin/11222019/massbin_revised_75_112_113_300_600_1444_1444_binning.tbl'
    mass_Te_file = fitspath + 'massbin/11222019/massbin_revised_75_112_113_300_600_1444_1444_derived_properties_metallicity.tbl'    
    HB_bin_file = fitspath + 'mass_LHbeta_bin/11222019/massLHbetabin_revised_75_112_113_300_600_1444_1444_binning.tbl'
    HB_Te_file = fitspath + 'mass_LHbeta_bin/11222019/massLHbetabin_revised_75_112_113_300_600_1444_1444_derived_properties_metallicity.tbl'
    
    
    MT_ascii = asc.read(metal_file)
    mass_bin_tbl = asc.read(mass_bin_file)
    HB_bin_tbl = asc.read(HB_bin_file)
    mass_metal_tbl = asc.read(mass_Te_file)
    HB_metal_tbl = asc.read(HB_Te_file)


    ##individual galaxy data, i.e. comes from MT_ascii
    log_mass = MT_ascii['Log10(Mass)'].data
    LHbeta = MT_ascii['HBeta_Luminosity'].data
    mass_ind_metal = MT_ascii['Mass_Bin_com_O_log'].data
    HB_ind_metal = MT_ascii['Mass_LHBeta_Bin_com_O_log'].data
    R23 = MT_ascii['R23'].data
    O32 = MT_ascii['O32'].data
    indiv_detect = MT_ascii['Individual_Detections'].data
    mass_bin_detect = MT_ascii['Mass_Bin_Detections'].data    
    HB_bin_detect = MT_ascii['Mass_LHBeta_Bin_Detections'].data  
    mass_bin_ID = MT_ascii['Mass_Bin_ID'].data
    HB_bin_ID = MT_ascii['Mass_LHBeta_Bin_ID'].data
    
    
    ##bin data, i.e. comes from either mass or HB specific files
    mass_avg_mass = mass_bin_tbl['mass_avg'].data
    HB_avg_mass = HB_bin_tbl['mass_avg'].data 
    mass_detect_col = mass_metal_tbl['Detection'].data
    HB_detect_col = HB_metal_tbl['Detection'].data        
    mass_metal = mass_metal_tbl['com_O_log'].data
    HB_metal = HB_metal_tbl['com_O_log'].data 
    
    
    ##detection determinations  
    #bins
    mass_detect = np.where(mass_detect_col == 1.0)[0]
    mass_nondetect = np.where(mass_detect_col == 0.5)[0]
    HB_detect = np.where(HB_detect_col == 1.0)[0]
    HB_nondetect = np.where(HB_detect_col == 0.5)[0]
    
    #individual
    mass_ind_detect = np.where((indiv_detect == 1.0) & (mass_bin_detect == 1.0))[0]
    mass_ind_nondetect = np.where((indiv_detect == 1.0) & (mass_bin_detect == 0.5))[0]
    HB_ind_detect = np.where((indiv_detect == 1.0) & (LHbeta > 0) & (HB_bin_detect == 1.0))[0]
    HB_ind_nondetect = np.where((indiv_detect == 1.0) & (LHbeta > 0) & (HB_bin_detect == 0.5))[0]
    
    #detections
    mass_bin_one = np.where((indiv_detect == 1.0) & (mass_bin_detect == 1.0) & (mass_bin_ID == 1))[0]
    mass_bin_three = np.where((indiv_detect == 1.0) & (mass_bin_detect == 1.0) & (mass_bin_ID == 3))[0] 
    mass_bin_four = np.where((indiv_detect == 1.0) & (mass_bin_detect == 1.0) & (mass_bin_ID == 4))[0]
    mass_bin_five = np.where((indiv_detect == 1.0) & (mass_bin_detect == 1.0) & (mass_bin_ID == 5))[0]
    mass_bin_six = np.where((indiv_detect == 1.0) & (mass_bin_detect == 1.0) & (mass_bin_ID == 6))[0]
    #non-detections with reliable limits
    mass_bin_two = np.where((indiv_detect == 1.0) & (mass_bin_detect == 0.5) & (mass_bin_ID == 2))[0]

    #detections
    HB_bin_two = np.where((indiv_detect == 1.0) & (LHbeta > 0) & (HB_bin_detect == 1.0) & (HB_bin_ID == 2))[0]
    HB_bin_four = np.where((indiv_detect == 1.0) & (LHbeta > 0) & (HB_bin_detect == 1.0) & (HB_bin_ID == 4))[0] 
    HB_bin_six = np.where((indiv_detect == 1.0) & (LHbeta > 0) & (HB_bin_detect == 1.0) & (HB_bin_ID == 6))[0]
    HB_bin_ten = np.where((indiv_detect == 1.0) & (LHbeta > 0) & (HB_bin_detect == 1.0) & (HB_bin_ID == 10))[0]
    #non-detections with reliable limits
    HB_bin_three = np.where((indiv_detect == 1.0) & (LHbeta > 0) & (HB_bin_detect == 0.5) & (HB_bin_ID == 3))[0]
    HB_bin_seven = np.where((indiv_detect == 1.0) & (LHbeta > 0) & (HB_bin_detect == 0.5) & (HB_bin_ID == 7))[0]
    HB_bin_eight = np.where((indiv_detect == 1.0) & (LHbeta > 0) & (HB_bin_detect == 0.5) & (HB_bin_ID == 8))[0]
    HB_bin_twelve = np.where((indiv_detect == 1.0) & (LHbeta > 0) & (HB_bin_detect == 0.5) & (HB_bin_ID == 12))[0]


    pdf_pages = PdfPages(fitspath + 'individual/12072019/metallicity_vs_R23_O32.pdf')
    
    
    #####Metallicity vs R23 and O32#####
    #Mass bin Metallicity vs R23 (one plot)
    fig1, ax1 = plt.subplots()
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax1.scatter(R23[mass_bin_one], mass_ind_metal[mass_bin_one], facecolors = 'None', edgecolors = 'red', label = 'Bin 1')
    ax1.scatter(R23[mass_bin_three], mass_ind_metal[mass_bin_three], facecolors = 'None', edgecolors = 'orange', label = 'Bin 3')
    ax1.scatter(R23[mass_bin_four], mass_ind_metal[mass_bin_four], facecolors = 'None', edgecolors = 'magenta', label = 'Bin 4')
    ax1.scatter(R23[mass_bin_five], mass_ind_metal[mass_bin_five], facecolors = 'None', edgecolors = 'green', label = 'Bin 5')
    ax1.scatter(R23[mass_bin_six], mass_ind_metal[mass_bin_six], facecolors = 'None', edgecolors = 'blue', label = 'Bin 6')
    
    ax1.scatter(R23[mass_bin_two], mass_ind_metal[mass_bin_two], marker='^', facecolors = 'None', edgecolors = 'purple', label = 'Bin 2')
    
    ax1.legend(loc = 'best')
    ax1.set_ylabel('Metallicity')
    ax1.set_xlabel('$R_{23}$')
    ax1.set_title('$M_\star$ Bins: Metallicity vs. $R_{23}$')
    fig1.set_size_inches(8, 8)
    fig1.savefig(pdf_pages, format='pdf')
    
    
    #Mass bin Metallicity vs R23 (per bin)
    fig2, ax2 = plt.subplots(nrows = 3, ncols = 2, sharex = True, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)    
    ax2[0, 0].scatter(R23[mass_bin_one], mass_ind_metal[mass_bin_one], facecolors = 'None', edgecolors = 'red', label = 'Bin 1')
    ax2[1, 0].scatter(R23[mass_bin_three], mass_ind_metal[mass_bin_three], facecolors = 'None', edgecolors = 'orange', label = 'Bin 3')
    ax2[1, 1].scatter(R23[mass_bin_four], mass_ind_metal[mass_bin_four], facecolors = 'None', edgecolors = 'magenta', label = 'Bin 4')
    ax2[2, 0].scatter(R23[mass_bin_five], mass_ind_metal[mass_bin_five], facecolors = 'None', edgecolors = 'green', label = 'Bin 5')
    ax2[2, 1].scatter(R23[mass_bin_six], mass_ind_metal[mass_bin_six], facecolors = 'None', edgecolors = 'blue', label = 'Bin 6')
    
    ax2[0, 1].scatter(R23[mass_bin_two], mass_ind_metal[mass_bin_two], marker='^', facecolors = 'None', edgecolors = 'purple', label = 'Bin 2')
    
    for ii in range(3):
        for jj in range(2):
            ax2[ii, jj].legend(loc = 'best')
    ax2[0, 0].set_ylabel('Metallicity')
    ax2[0, 0].set_xlabel('$R_{23}$')
    ax2[0, 0].set_title('$M_\star$ Bins: Metallicity vs. $R_{23}$')
    fig2.set_size_inches(8, 8)
    fig2.savefig(pdf_pages, format='pdf')
    
    
    #Mass bin Metallicity vs O32 (one plot)
    fig3, ax3 = plt.subplots()
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax3.scatter(O32[mass_bin_one], mass_ind_metal[mass_bin_one], facecolors = 'None', edgecolors = 'red', label = 'Bin 1')
    ax3.scatter(O32[mass_bin_three], mass_ind_metal[mass_bin_three], facecolors = 'None', edgecolors = 'orange', label = 'Bin 3')
    ax3.scatter(O32[mass_bin_four], mass_ind_metal[mass_bin_four], facecolors = 'None', edgecolors = 'magenta', label = 'Bin 4')
    ax3.scatter(O32[mass_bin_five], mass_ind_metal[mass_bin_five], facecolors = 'None', edgecolors = 'green', label = 'Bin 5')
    ax3.scatter(O32[mass_bin_six], mass_ind_metal[mass_bin_six], facecolors = 'None', edgecolors = 'blue', label = 'Bin 6')
    
    ax3.scatter(O32[mass_bin_two], mass_ind_metal[mass_bin_two], marker='^', facecolors = 'None', edgecolors = 'purple', label = 'Bin 2')
    
    ax3.legend(loc = 'best')
    ax3.set_ylabel('Metallicity')
    ax3.set_xlabel('$O_{32}$')
    ax3.set_title('$M_\star$ Bins: Metallicity vs. $O_{32}$')
    fig3.set_size_inches(8, 8)
    fig3.savefig(pdf_pages, format='pdf')
    
    
    #Mass bin Metallicity vs O32 (per bin)
    fig4, ax4 = plt.subplots(nrows = 3, ncols = 2, sharex = True, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax4[0, 0].scatter(O32[mass_bin_one], mass_ind_metal[mass_bin_one], facecolors = 'None', edgecolors = 'red', label = 'Bin 1')
    ax4[1, 0].scatter(O32[mass_bin_three], mass_ind_metal[mass_bin_three], facecolors = 'None', edgecolors = 'orange', label = 'Bin 3')
    ax4[1, 1].scatter(O32[mass_bin_four], mass_ind_metal[mass_bin_four], facecolors = 'None', edgecolors = 'magenta', label = 'Bin 4')
    ax4[2, 0].scatter(O32[mass_bin_five], mass_ind_metal[mass_bin_five], facecolors = 'None', edgecolors = 'green', label = 'Bin 5')
    ax4[2, 1].scatter(O32[mass_bin_six], mass_ind_metal[mass_bin_six], facecolors = 'None', edgecolors = 'blue', label = 'Bin 6')
    
    ax4[0, 1].scatter(O32[mass_bin_two], mass_ind_metal[mass_bin_two], marker='^', facecolors = 'None', edgecolors = 'purple', label = 'Bin 2')
    
    for ii in range(3):
        for jj in range(2):
            ax4[ii, jj].legend(loc = 'best')
    ax4[0, 0].set_ylabel('Metallicity')
    ax4[0, 0].set_xlabel('$O_{32}$')
    ax4[0, 0].set_title('$M_\star$ Bins: Metallicity vs. $O_{32}$')
    fig4.set_size_inches(8, 8)
    fig4.savefig(pdf_pages, format='pdf')
    
    
    
    
    #Mass-LHBeta bin Metallicity vs R23 (one plot)
    fig5, ax5 = plt.subplots()
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax5.scatter(R23[HB_bin_two], HB_ind_metal[HB_bin_two], facecolors = 'None', edgecolors = 'red', label = 'Bin 2')
    ax5.scatter(R23[HB_bin_four], HB_ind_metal[HB_bin_four], facecolors = 'None', edgecolors = 'orange', label = 'Bin 4')
    ax5.scatter(R23[HB_bin_six], HB_ind_metal[HB_bin_six], facecolors = 'None', edgecolors = 'magenta', label = 'Bin 6')
    ax5.scatter(R23[HB_bin_ten], HB_ind_metal[HB_bin_ten], facecolors = 'None', edgecolors = 'green', label = 'Bin 10')
    
    ax5.scatter(R23[HB_bin_three], HB_ind_metal[HB_bin_three], marker='^', facecolors = 'None', edgecolors = 'cyan', label = 'Bin 3')
    ax5.scatter(R23[HB_bin_seven], HB_ind_metal[HB_bin_seven], marker='^', facecolors = 'None', edgecolors = 'blue', label = 'Bin 7')
    ax5.scatter(R23[HB_bin_eight], HB_ind_metal[HB_bin_eight], marker='^', facecolors = 'None', edgecolors = 'purple', label = 'Bin 8')
    ax5.scatter(R23[HB_bin_twelve], HB_ind_metal[HB_bin_twelve], marker='^', facecolors = 'None', edgecolors = 'black', label = 'Bin 12')
    
    ax5.legend(loc = 'best')
    ax5.set_ylabel('Metallicity')
    ax5.set_xlabel('$R_{23}$')
    ax5.set_title('$M_\star$-LH$\\beta$ Bins: Metallicity vs. $R_{23}$')
    fig5.set_size_inches(8, 8)
    fig5.savefig(pdf_pages, format='pdf')
    
    
    #Mass-LHBeta bin Metallicity vs R23 (per bin)
    fig6, ax6 = plt.subplots(nrows = 4, ncols = 2, sharex = True, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)    
    ax6[0, 0].scatter(R23[HB_bin_two], HB_ind_metal[HB_bin_two], facecolors = 'None', edgecolors = 'red', label = 'Bin 2')
    ax6[1, 0].scatter(R23[HB_bin_four], HB_ind_metal[HB_bin_four], facecolors = 'None', edgecolors = 'orange', label = 'Bin 4')
    ax6[1, 1].scatter(R23[HB_bin_six], HB_ind_metal[HB_bin_six], facecolors = 'None', edgecolors = 'magenta', label = 'Bin 6')
    ax6[3, 0].scatter(R23[HB_bin_ten], HB_ind_metal[HB_bin_ten], facecolors = 'None', edgecolors = 'green', label = 'Bin 10')

    ax6[0, 1].scatter(R23[HB_bin_three], HB_ind_metal[HB_bin_three], marker='^', facecolors = 'None', edgecolors = 'cyan', label = 'Bin 3')
    ax6[2, 0].scatter(R23[HB_bin_seven], HB_ind_metal[HB_bin_seven], marker='^', facecolors = 'None', edgecolors = 'blue', label = 'Bin 7')
    ax6[2, 1].scatter(R23[HB_bin_eight], HB_ind_metal[HB_bin_eight], marker='^', facecolors = 'None', edgecolors = 'purple', label = 'Bin 8')
    ax6[3, 1].scatter(R23[HB_bin_twelve], HB_ind_metal[HB_bin_twelve], marker='^', facecolors = 'None', edgecolors = 'black', label = 'Bin 12')
    
    for ii in range(4):
        for jj in range(2):
            ax6[ii, jj].legend(loc = 'best')
    ax6[0, 0].set_ylabel('Metallicity')
    ax6[0, 0].set_xlabel('$R_{23}$')
    ax6[0, 0].set_title('$M_\star$-LH$\\beta$ Bins: Metallicity vs. $R_{23}$')
    fig6.set_size_inches(8, 8)
    fig6.savefig(pdf_pages, format='pdf')
    
    
    #Mass-LHBeta bin Metallicity vs O32 (one plot)
    fig7, ax7 = plt.subplots()
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)    
    ax7.scatter(O32[HB_bin_two], HB_ind_metal[HB_bin_two], facecolors = 'None', edgecolors = 'red', label = 'Bin 2')
    ax7.scatter(O32[HB_bin_four], HB_ind_metal[HB_bin_four], facecolors = 'None', edgecolors = 'orange', label = 'Bin 4')
    ax7.scatter(O32[HB_bin_six], HB_ind_metal[HB_bin_six], facecolors = 'None', edgecolors = 'magenta', label = 'Bin 6')
    ax7.scatter(O32[HB_bin_ten], HB_ind_metal[HB_bin_ten], facecolors = 'None', edgecolors = 'green', label = 'Bin 10')
    
    ax7.scatter(O32[HB_bin_three], HB_ind_metal[HB_bin_three], marker='^', facecolors = 'None', edgecolors = 'cyan', label = 'Bin 3')
    ax7.scatter(O32[HB_bin_seven], HB_ind_metal[HB_bin_seven], marker='^', facecolors = 'None', edgecolors = 'blue', label = 'Bin 7')
    ax7.scatter(O32[HB_bin_eight], HB_ind_metal[HB_bin_eight], marker='^', facecolors = 'None', edgecolors = 'purple', label = 'Bin 8')
    ax7.scatter(O32[HB_bin_twelve], HB_ind_metal[HB_bin_twelve], marker='^', facecolors = 'None', edgecolors = 'black', label = 'Bin 12')
    
    ax7.legend(loc = 'best')
    ax7.set_ylabel('Metallicity')
    ax7.set_xlabel('$O_{32}$')
    ax7.set_title('$M_\star$-LH$\\beta$ Bins: Metallicity vs. $O_{32}$')
    fig7.set_size_inches(8, 8)
    fig7.savefig(pdf_pages, format='pdf')
    
    
    #Mass-LHBeta bin Metallicity vs O32 (per bin)
    fig8, ax8 = plt.subplots(nrows = 4, ncols = 2, sharex = True, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)    
    ax8[0, 0].scatter(O32[HB_bin_two], HB_ind_metal[HB_bin_two], facecolors = 'None', edgecolors = 'red', label = 'Bin 2')
    ax8[1, 0].scatter(O32[HB_bin_four], HB_ind_metal[HB_bin_four], facecolors = 'None', edgecolors = 'orange', label = 'Bin 4')
    ax8[1, 1].scatter(O32[HB_bin_six], HB_ind_metal[HB_bin_six], facecolors = 'None', edgecolors = 'magenta', label = 'Bin 6')
    ax8[3, 0].scatter(O32[HB_bin_ten], HB_ind_metal[HB_bin_ten], facecolors = 'None', edgecolors = 'green', label = 'Bin 10')
    
    ax8[0, 1].scatter(O32[HB_bin_three], HB_ind_metal[HB_bin_three], marker='^', facecolors = 'None', edgecolors = 'cyan', label = 'Bin 3')
    ax8[2, 0].scatter(O32[HB_bin_seven], HB_ind_metal[HB_bin_seven], marker='^', facecolors = 'None', edgecolors = 'blue', label = 'Bin 7')
    ax8[2, 1].scatter(O32[HB_bin_eight], HB_ind_metal[HB_bin_eight], marker='^', facecolors = 'None', edgecolors = 'purple', label = 'Bin 8')
    ax8[3, 1].scatter(O32[HB_bin_twelve], HB_ind_metal[HB_bin_twelve], marker='^', facecolors = 'None', edgecolors = 'black', label = 'Bin 12')
    
    for ii in range(4):
        for jj in range(2):
            ax8[ii, jj].legend(loc = 'best')
    ax8[0, 0].set_ylabel('Metallicity')
    ax8[0, 0].set_xlabel('$O_{32}$')
    ax8[0, 0].set_title('$M_\star$-LH$\\beta$ Bins: Metallicity vs. $O_{32}$')
    fig8.set_size_inches(8, 8)
    fig8.savefig(pdf_pages, format='pdf')
    
    
    pdf_pages.savefig()
    
    pdf_pages.close()









