import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit 
      
    

def bin_derived_props_plots(metal_file, em_file, out_file):
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
    pdf_pages = PdfPages(out_file)
    
    com_O_log = metal_table['com_O_log'].data
    T_e = metal_table['Temperature'].data
    R23_composite = metal_table['R23_Composite'].data
    O32_composite = metal_table['O32_Composite'].data
    
    avg_mass = em_table['mass_avg'].data
    detection = em_table['Detection'].data
    OII = em_table['OII_3727_Flux_Observed'].data
    OIII4363 = em_table['Updated_OIII_4363_Flux_Observed'].data
    OIII5007 = em_table['OIII_5007_Flux_Observed'].data
    HBeta = em_table['HBETA_Flux_Observed'].data
    
    #non_detect = np.where(detection != 1)[0]
    #detect = np.where(detection == 1)[0]
    detect = np.where(detection == 1)[0]
    non_detect = np.where(detection == 0.5)[0]   #non-detection with reliable 5007
        
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
    #y = 8.798 - np.log10(1 + ((10**8.901)/(10**avg_mass))**0.640)
    #ax7.plot(avg_mass, y, color='g', linestyle = '-', marker = '*')
    ax7.scatter(avg_mass[non_detect], com_O_log[non_detect], marker = '^')
    ax7.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax7.set_ylabel('12+log(O/H) $T_e$')
    ax7.set_title('Composite Metallicity vs. Avg Mass')
    pdf_pages.savefig()
    
    pdf_pages.close()
    





def indiv_derived_props_plots(fitspath, metal_Te_file, mass_bin_file, HB_bin_file, mass_metal_file,
                              HB_metal_file, MTO, restrict_MTO = False):    
    '''
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
    
    MT_ascii = asc.read(metal_Te_file)
    mass_bin_tbl = asc.read(mass_bin_file)
    HB_bin_tbl = asc.read(HB_bin_file)
    mass_metal_tbl = asc.read(mass_metal_file)
    HB_metal_tbl = asc.read(HB_metal_file)


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



    pdf_pages = PdfPages(fitspath + 'individual_metal_plots' + MTO + '.pdf')
    
    
    
    
    ##HBeta Luminosity vs Mass ColorMap=Metallicity##
    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot1 = ax1.scatter(log_mass[mass_ind_detect], LHbeta[mass_ind_detect], 5.0,
                        c=mass_ind_metal[mass_ind_detect], marker='*')
    plot1 = ax1.scatter(log_mass[mass_ind_nondetect], LHbeta[mass_ind_nondetect], 5.0, facecolors = 'None',
                        c=mass_ind_metal[mass_ind_nondetect], marker='^')
    cb = fig1.colorbar(plot1)
    cb.set_label('Metallicity')
    ax1.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1.set_ylabel('log(H$\\beta$ Luminosity)')
    ax1.set_title('$M_\star$ Bins: H$\\beta$ Luminosity vs. Mass Colormap=Metallicity')
    fig1.set_size_inches(8, 8)
    fig1.savefig(pdf_pages, format='pdf')
    
    fig2, ax2 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot2 = ax2.scatter(log_mass[HB_ind_detect], LHbeta[HB_ind_detect], 5.0, c=HB_ind_metal[HB_ind_detect],
                        marker='*')
    plot2 = ax2.scatter(log_mass[HB_ind_nondetect], LHbeta[HB_ind_nondetect], 5.0, facecolors = 'None',
                        c=HB_ind_metal[HB_ind_nondetect], marker='^')
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
    plot3 = ax3.scatter(R23[mass_ind_detect], O32[mass_ind_detect], 5.0, c=mass_ind_metal[mass_ind_detect],
                        marker='*')
    plot3 = ax3.scatter(R23[mass_ind_nondetect], O32[mass_ind_nondetect], 5.0, facecolors = 'None',
                        c=mass_ind_metal[mass_ind_nondetect], marker='^')
    cb = fig3.colorbar(plot3)
    cb.set_label('Metallicity')
    ax3.set_xlabel('$R_{23}$')
    ax3.set_ylabel('$O_{32}$')
    ax3.set_title('$M_\star$ Bins: $O_{32}$ vs. $R_{23}$ Colormap=Metallicity')
    fig3.set_size_inches(8, 8)
    fig3.savefig(pdf_pages, format='pdf')
    
    fig4, ax4 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot4 = ax4.scatter(R23[HB_ind_detect], O32[HB_ind_detect], 5.0, c=HB_ind_metal[HB_ind_detect], marker='*')
    plot4 = ax4.scatter(R23[HB_ind_nondetect], O32[HB_ind_nondetect], 5.0, facecolors = 'None', 
                        c=HB_ind_metal[HB_ind_nondetect], marker='^')
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
    ax5[0].scatter(R23[mass_ind_detect], mass_ind_metal[mass_ind_detect], facecolors = 'None', edgecolors = 'blue',
                   label = '$R_{23}$')
    ax5[1].scatter(O32[mass_ind_detect], mass_ind_metal[mass_ind_detect], facecolors = 'None', edgecolors = 'red',
                   label = '$O_{32}$')
    ax5[0].scatter(R23[mass_ind_nondetect], mass_ind_metal[mass_ind_nondetect], marker='^', facecolors = 'None', edgecolors = 'blue',
                   label = '$R_{23}$', alpha = 0.5)
    ax5[1].scatter(O32[mass_ind_nondetect], mass_ind_metal[mass_ind_nondetect], marker='^', facecolors = 'None', edgecolors = 'red',
                   label = '$O_{32}$', alpha = 0.5)
    ax5[0].legend(loc = 'best')
    ax5[1].legend(loc = 'best')
    ax5[0].set_ylabel('Metallicity')
    ax5[0].set_title('$M_\star$ Bins: Metallicity vs. $R_{23}$')
    ax5[1].set_title('$M_\star$ Bins: Metallicity vs. $O_{32}$')
    fig5.set_size_inches(8, 8)
    fig5.savefig(pdf_pages, format='pdf')
    
    fig6, ax6 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax6[0].scatter(R23[HB_ind_detect], HB_ind_metal[HB_ind_detect], facecolors = 'None', edgecolors = 'blue',
                label = '$R_{23}$')
    ax6[1].scatter(O32[HB_ind_detect], HB_ind_metal[HB_ind_detect], facecolors = 'None', edgecolors = 'red',
                label = '$O_{32}$')
    ax6[0].scatter(R23[HB_ind_nondetect], HB_ind_metal[HB_ind_nondetect], marker='^', facecolors = 'None', edgecolors = 'blue',
                label = '$R_{23}$', alpha = 0.5)
    ax6[1].scatter(O32[HB_ind_nondetect], HB_ind_metal[HB_ind_nondetect], marker='^', facecolors = 'None', edgecolors = 'red',
                label = '$O_{32}$', alpha = 0.5)
    ax6[0].legend(loc = 'best')
    ax6[1].legend(loc = 'best')
    ax6[0].set_ylabel('Metallicity')
    ax6[0].set_title('$M_\star$-LH$\\beta$ Bins: Metallicity vs. $R_{23}$')
    ax6[1].set_title('$M_\star$-LH$\\beta$ Bins: Metallicity vs. $O_{32}$')
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
    
    
    #Individual detections and non-detections
    ax7[0].scatter(log_mass[mass_ind_detect], mass_ind_metal[mass_ind_detect], s = 15, facecolors = 'None',
                   edgecolors = 'blue', label = 'Individual Detections')
    ax7[1].scatter(log_mass[HB_ind_detect], HB_ind_metal[HB_ind_detect], s = 15, facecolors = 'None',
                   edgecolors = 'blue', label = 'Individual Detections')
    ax7[0].scatter(log_mass[mass_ind_nondetect], mass_ind_metal[mass_ind_nondetect], s = 15, marker='^',
                   facecolors = 'None', edgecolors = 'blue', label = 'Individual Non-Detections', alpha = 0.5)
    ax7[1].scatter(log_mass[HB_ind_nondetect], HB_ind_metal[HB_ind_nondetect], s = 15, marker='^',
                   facecolors = 'None', edgecolors = 'blue', label = 'Individual Non-Detections', alpha = 0.5)
    
    print('Number of mass bin individual sources plotted:', len(log_mass[mass_ind_detect]))
    print('Number of mass-LHbeta bin individual sources plotted:', len(log_mass[HB_ind_detect]))
    
    #Mass bin detections and non-detections
    ax7[0].scatter(mass_avg_mass[mass_detect], mass_metal[mass_detect], s = 25, color = 'red',
                   label = 'Bin Detections')
    ax7[0].scatter(mass_avg_mass[mass_nondetect], mass_metal[mass_nondetect], s = 25, color = 'red',
                   marker = '^', label = 'Bin Non-Detections', alpha = 0.5)
    #HBeta bin detections and non-detections
    ax7[1].scatter(HB_avg_mass[HB_detect], HB_metal[HB_detect], s = 25, color = 'red',
                   label = 'Bin Detections')
    ax7[1].scatter(HB_avg_mass[HB_nondetect], HB_metal[HB_nondetect], s = 25, color = 'red',
                   marker = '^', label = 'Bin Non-Detections', alpha = 0.5)

    
    ##Curve fit 
    #exclude last mass bin data point for curve-fitting purposes
    #mass_detect = mass_detect[:-1]
    
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
        o12, o22 = curve_fit(fit, HB_avg_mass[HB_detect], HB_metal[HB_detect], p0 = p0, bounds = para_bounds)
        print(o12)
    except ValueError:
        print('fail HBeta Bins')
        fail = True
        
    if not fail:
        if restrict_MTO == False:
            ax7[0].plot(mass, mass_metal_fit(mass, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
            ax7[1].plot(mass, mass_metal_fit(mass, *o12), alpha = 0.5, color = 'red', label = 'Our Fit')
        else:
            ax7[0].plot(mass, mass_metal_fit2(mass, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
            ax7[1].plot(mass, mass_metal_fit2(mass, *o12), alpha = 0.5, color = 'red', label = 'Our Fit')
          
    ax7[0].legend(title = '$M_\star$ Bins', fontsize = 5, loc = 'upper left')    
    ax7[1].legend(title = '$M_\star$-LH$\\beta$ Bins', fontsize = 5, loc = 'upper left')
    ax7[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax7[0].set_ylabel('12+log(O/H) $T_e$')
    ax7[1].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
   
        
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










