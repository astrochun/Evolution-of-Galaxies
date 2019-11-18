import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit 
      
    

def bin_derived_props_plots(metal_file, em_file, out_file):
    metal_table = asc.read(metal_file)
    em_table = asc.read(em_file)
    pdf_pages = PdfPages(out_file)
    
    com_O_log = metal_table['com_O_log'].data
    T_e = metal_table['Temperature'].data
    R23_composite = metal_table['R23_Composite'].data
    O32_composite = metal_table['O32_Composite'].data
    
    avg_mass = em_table['mass_avg'].data
    detection = em_table['Detection'].data
    OII = em_table['OII_3727_Flux_Observed'].data
    OIII4363 = em_table['Updated_OIII_4363_Flux_Observed'].data
    OIII4363_SN = em_table['OIII_4363_S/N'].data
    OIII5007 = em_table['OIII_5007_Flux_Observed'].data
    OIII5007_SN = em_table['OIII_5007_S/N'].data 
    HBeta = em_table['HBETA_Flux_Observed'].data
    
    #non_detect = np.where(detection != 1)[0]
    #detect = np.where(detection == 1)[0]
    detect = np.where((detection == 1) & (OIII4363_SN >= 3) & (OIII5007_SN > 100))[0]
    non_detect = np.where((detection == 0) & (OIII4363_SN < 3) & (OIII5007_SN > 100))[0]   #non-detection but may be reliable
        
    #Line ratios        
    OII_HBeta = OII / HBeta
    OIII_HBeta = (OIII5007 * (1 + 1/3.1)) / HBeta            
    OIII_ratio = OIII4363 / (OIII5007 * (1 + 1/3.1))            
    
    #Line Ratios vs Mass
    fig1, ax1 = plt.subplots(2, 2)
    ax1[0, 0].scatter(avg_mass, R23_composite, label = 'R_23')
    ax1[0, 0].scatter(avg_mass, O32_composite, label = 'O_32')
    ax1[0, 0].legend(loc = 'best')
    ax1[0, 0].set_xticklabels([])
    ax1[0, 1].scatter(avg_mass, np.log10(OII_HBeta), label = 'OII/HBeta')
    ax1[0, 1].legend(loc = 'best')
    ax1[0, 1].set_xticklabels([])
    ax1[1, 0].scatter(avg_mass, np.log10(OIII_HBeta), label = 'OIII/HBeta')
    ax1[1, 0].legend(loc = 'best')
    ax1[1, 0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1[1, 1].scatter(avg_mass, np.log10(OIII_ratio), label = '4363/(5007 * (1 + 1/3.1))')
    ax1[1, 1].legend(loc = 'best')
    ax1[1, 1].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    for t_ax in ax1:
        for tt in range(len(t_ax)):
            t_ax[tt].set_xlim(8.5,11.0)
    plt.subplots_adjust(left = 0.075, right = 0.97, bottom = 0.075, top = 0.99, wspace = 0.225, hspace = 0.05)
    pdf_pages.savefig()    
    
    #R23 vs Temperature
    fig2, ax2 = plt.subplots()
    ax2.scatter(np.log10(T_e[detect]), R23_composite[detect], marker = '.')
    ax2.scatter(np.log10(T_e[non_detect]), R23_composite[non_detect], marker = '^')
    ax2.set_xlabel('$T_e$ (K)')
    ax2.set_ylabel('$R_{23}$')
    ax2.set_title('$R_{23}$ vs. Temperatures')
    pdf_pages.savefig()
     
    #O32 vs Temperature
    fig3, ax3 = plt.subplots()
    ax3.scatter(np.log10(T_e[detect]), O32_composite[detect], marker = '.')
    ax3.scatter(np.log10(T_e[non_detect]), O32_composite[non_detect], marker = '^')
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
    ax6.scatter(avg_mass[non_detect], np.log10(T_e[non_detect]), marker = '^')
    ax6.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax6.set_ylabel('$T_e$ (K)')
    ax6.set_title('Temperatures vs. Avg Mass')
    pdf_pages.savefig()
    
    #Metallicity vs Mass
    fig7, ax7 = plt.subplots()
    ax7.scatter(avg_mass[detect], com_O_log[detect], marker = '.')
    #y = 8.798 - np.log10(1 + ((10**8.901)/(10**avg_mass))**0.640)
    #ax7.plot(avg_mass, y, color='g', linestyle = '-', marker = '*')
    ax7.scatter(avg_mass[non_detect], com_O_log[non_detect], marker = '^')
    ax7.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax7.set_ylabel('12+log(O/H) $T_e$')
    ax7.set_title('Composite Metallicity vs. Avg Mass')
    pdf_pages.savefig()
    
    pdf_pages.close()
    



def derived_properties_plots(fitspath):    
    metal_Te_file = fitspath + 'individual_derived_properties_metallicity.tbl'
    HB_valid_file = fitspath + 'hbeta_revised_75_112_113_300_600_1444_1444_massbin_validation.tbl'
    mass_valid_file = fitspath + 'revised_75_112_113_300_600_1444_1444_massbin_validation.tbl'
    HB_bin_metal_file = fitspath + 'hbeta_revised_75_112_113_300_600_1444_1444_updated_massbin_derived_properties_metallicity.tbl'
    mass_bin_metal_file = fitspath + 'revised_75_112_113_300_600_1444_1444_updated_massbin_derived_properties_metallicity.tbl'
    
    MT_ascii = asc.read(metal_Te_file)
    HB_valid_table = asc.read(HB_valid_file)
    mass_valid_table = asc.read(mass_valid_file)
    HB_bin_metal_table = asc.read(HB_bin_metal_file)
    mass_bin_metal_table = asc.read(mass_bin_metal_file)

    log_mass = MT_ascii['Log10(Mass)'].data
    LHbeta = MT_ascii['HBeta_Luminosity'].data
    mass_ind_metal = MT_ascii['Mass Bin com_O_log'].data
    HB_ind_metal = MT_ascii['LHBeta Bin com_O_log'].data
    R23 = MT_ascii['R23'].data
    O32 = MT_ascii['O32'].data
    indiv_detect = MT_ascii['Individual Detections'].data
    mass_bin_detect = MT_ascii['Mass Bin Detections'].data    #len = 4088  --> 1's and 0's populated over all individual
    HB_bin_detect = MT_ascii['LHBeta Bin Detections'].data    #len = 4088  --> 1's and 0's populated over all individual
    HB_bin_avg_mass = HB_valid_table['mass_avg'].data
    HB_bin_detect_valid = HB_valid_table['Detection'].data    #len = 14  --> 1's and 0's populated over all bins
    mass_bin_avg_mass = mass_valid_table['mass_avg'].data
    mass_bin_detect_valid = mass_valid_table['Detection'].data     #len = 7 --> 1's and 0's populated over all bins
    HB_bin_metal = HB_bin_metal_table['com_O_log'].data
    mass_bin_metal = mass_bin_metal_table['com_O_log'].data
    
    HB_detection = np.where((indiv_detect == 1.0) & (LHbeta > 0) & (HB_bin_detect == 1.0))[0]
    mass_detection = np.where((indiv_detect == 1.0) & (mass_bin_detect == 1.0))[0]
    HB_detect = np.where(HB_bin_detect_valid == 1.0)[0]
    HB_detect = HB_detect[:-1]
    HB_non_detect = np.where(HB_bin_detect_valid == 0.0)[0]
    mass_detect = np.where(mass_bin_detect_valid == 1.0)[0]
    mass_non_detect = np.where(mass_bin_detect_valid == 0.0)[0]

    pdf_pages = PdfPages(fitspath + 'individual_metal_plots.pdf')
    
    
    
    
    ##HBeta Luminosity vs Mass ColorMap=Metallicity##
    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot1 = ax1.scatter(log_mass[HB_detection], LHbeta[HB_detection], 0.8, c=HB_ind_metal[HB_detection], marker='*')
    cb = fig1.colorbar(plot1)
    cb.set_label('Metallicity')
    ax1.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1.set_ylabel('log(H$\\beta$ Luminosity)')
    ax1.set_title('$M_\star$-LH$\\beta$ Bins: H$\\beta$ Luminosity vs. Mass Colormap=Metallicity')
    fig1.set_size_inches(8, 8)
    fig1.savefig(pdf_pages, format='pdf')
    
    fig2, ax2 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot2 = ax2.scatter(log_mass[mass_detection], LHbeta[mass_detection], 0.8, c=mass_ind_metal[mass_detection], marker='*')
    cb = fig2.colorbar(plot2)
    cb.set_label('Metallicity')
    ax2.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax2.set_ylabel('log(H$\\beta$ Luminosity)')
    ax2.set_title('$M_\star$ Bins: H$\\beta$ Luminosity vs. Mass Colormap=Metallicity')
    fig2.set_size_inches(8, 8)
    fig2.savefig(pdf_pages, format='pdf')
    
    
    
    
    ##O32 vs R23 ColorMap=Metallicity##
    fig3, ax3 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot3 = ax3.scatter(R23[HB_detection], O32[HB_detection], 0.8, c=HB_ind_metal[HB_detection], marker='*')
    cb = fig3.colorbar(plot3)
    cb.set_label('Metallicity')
    ax3.set_xlabel('$R_{23}$')
    ax3.set_ylabel('$O_{32}$')
    ax3.set_title('$M_\star$-LH$\\beta$ Bins: $O_{32}$ vs. $R_{23}$ Colormap=Metallicity')
    fig3.set_size_inches(8, 8)
    fig3.savefig(pdf_pages, format='pdf')
    
    fig4, ax4 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot4 = ax4.scatter(R23[mass_detection], O32[mass_detection], 0.8, c=mass_ind_metal[mass_detection], marker='*')
    cb = fig4.colorbar(plot4)
    cb.set_label('Metallicity')
    ax4.set_xlabel('$R_{23}$')
    ax4.set_ylabel('$O_{32}$')
    ax4.set_title('$M_\star$ Bins: $O_{32}$ vs. $R_{23}$ Colormap=Metallicity')
    fig4.set_size_inches(8, 8)
    fig4.savefig(pdf_pages, format='pdf')
    
    
    
    
    ##Metallicity vs R23 and O32##
    fig5, ax5 = plt.subplots()
    ax5.scatter(R23[HB_detection], HB_ind_metal[HB_detection], facecolors = 'None', edgecolors = 'blue', label = '$R_{23}$')
    ax5.scatter(O32[HB_detection], HB_ind_metal[HB_detection], facecolors = 'None', edgecolors = 'red', label = '$O_{32}$')
    ax5.legend(loc = 'best')
    ax5.set_ylabel('Metallicity')
    ax5.set_title('$M_\star$-LH$\\beta$ Bins: Metallicity vs. $R_{23}$ and $O_{32}$')
    fig5.set_size_inches(8, 8)
    fig5.savefig(pdf_pages, format='pdf')
    
    fig6, ax6 = plt.subplots()
    ax6.scatter(R23[mass_detection], mass_ind_metal[mass_detection], facecolors = 'None', edgecolors = 'blue', label = '$R_{23}$')
    ax6.scatter(O32[mass_detection], mass_ind_metal[mass_detection], facecolors = 'None', edgecolors = 'red', label = '$O_{32}$')
    ax6.legend(loc = 'best')
    ax6.set_ylabel('Metallicity')
    ax6.set_title('$M_\star$ Bins: Metallicity vs. $R_{23}$ and $O_{32}$')
    fig6.set_size_inches(8, 8)
    fig6.savefig(pdf_pages, format='pdf')
    
    
    
    
    ##Metallcity vs Mass##
    fig7, ax7 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.98, wspace = 0.0)
    
    
    #Andrews&Martini fit
    mass = np.arange(7.5, 10, 0.05)
    y = 8.798 - np.log10(1 + ((10**8.901)/(10**mass))**0.640)
    ax7[0].plot(mass, y, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    ax7[1].plot(mass, y, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    
    
    #Individual detections
    ax7[0].scatter(log_mass[HB_detection], HB_ind_metal[HB_detection], s = 15, facecolors = 'None',
                   edgecolors = 'blue', label = 'Individual Detections')
    ax7[1].scatter(log_mass[mass_detection], mass_ind_metal[mass_detection], s = 15, facecolors = 'None',
                   edgecolors = 'blue', label = 'Individual Detections')
    
    
    #HBeta bin detections and non-detections
    ax7[0].scatter(HB_bin_avg_mass[HB_detect], HB_bin_metal[HB_detect], s = 25, color = 'red', label = 'Bin Detections')
    ax7[0].scatter(HB_bin_avg_mass[HB_non_detect], HB_bin_metal[HB_non_detect], s = 25, color = 'red', marker = '^', label = 'Bin Non-Detections')
    #Mass bin detections and non-detections
    ax7[1].scatter(mass_bin_avg_mass[mass_detect], mass_bin_metal[mass_detect], s = 25, color = 'red', label = 'Bin Detections')
    ax7[1].scatter(mass_bin_avg_mass[mass_non_detect], mass_bin_metal[mass_non_detect], s = 25, color = 'red', marker = '^', label = 'Bin Non-Detections')

    
    #Curve fit with restricted M_TO
    fail = False
    #p0 = [8.798, 8.901, 0.640]
    p0 = [8.798, 0.640]
    #para_bounds = ((8.0, 8.0, 0.0), (9.0, 9.5, 1.0))
    para_bounds = ((8.0, 0.0), (9.0, 1.0))
    try:
        o11, o21 = curve_fit(mass_metal_fit2, HB_bin_avg_mass[HB_detect], HB_bin_metal[HB_detect], p0 = p0, bounds = para_bounds)
        print(o11)
    except ValueError:
        print('fail HBeta Bins')
        fail = True
    try:
        o12, o22 = curve_fit(mass_metal_fit2, mass_bin_avg_mass[mass_detect], mass_bin_metal[mass_detect], p0 = p0, bounds = para_bounds)
        print(o12)
    except ValueError:
        print('fail Mass Bins')
        fail = True
        
    if not fail:
        ax7[0].plot(mass, mass_metal_fit2(mass, *o11), alpha = 0.5, label = 'Our Fit')
        ax7[1].plot(mass, mass_metal_fit2(mass, *o12), alpha = 0.5, label = 'Our Fit')
          
    ax7[0].legend(title = '$M_\star$-LH$\\beta$ Bins', fontsize = 5, loc = 'upper left')
    ax7[1].legend(title = '$M_\star$ Bins', fontsize = 5, loc = 'upper left')
    ax7[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax7[0].set_ylabel('12+log(O/H) $T_e$')
    ax7[1].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
   
        
    pdf_pages.savefig()
    
    pdf_pages.close()
    
    
    
def mass_metal_fit(mass, a, b, g):
    '''
    Andrews & Martini Mass-Metallicity Relation:
    8.798 - np.log10(1 + ((10**8.901)/(10**mass))**0.640)
    '''
    
    return a - np.log10(1 + ((10**b)/(10**mass))**g)   



def mass_metal_fit2(mass, a, g):
    '''
    Andrews & Martini Mass-Metallicity Relation:
    8.798 - np.log10(1 + ((10**8.901)/(10**mass))**0.640)
    '''
    
    return a - np.log10(1 + ((10**8.901)/(10**mass))**g) 



def find_outlier(fitspath):
    metal_Te_file = fitspath + 'individual_derived_properties_metallicity.tbl'
    
    MT_ascii = asc.read(metal_Te_file)
    
    ID = MT_ascii['Source_ID'].data
    log_mass = MT_ascii['Log10(Mass)'].data
    indiv_detect = MT_ascii['Individual Detections'].data
    LHbeta = MT_ascii['HBeta_Luminosity'].data
    
    outlier = np.where((log_mass < 7.0) & (indiv_detect == 1.0) & (LHbeta > 0))[0]
    
    print(log_mass[outlier]) 
    print(ID[outlier]) 



def lum_vs_mass(fitspath, binning_file):
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










