import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from Metallicity_Stack_Commons.column_names import indv_names0, temp_metal_names0, bin_mzevolve_names0
from Metallicity_Stack_Commons.column_names import bin_names0, filename_dict
from .relation_fitting import curve_fitting, mass_metal_fit, extract_error_bars, AM13




def indiv_derived_props_plots(fitspath, restrictMTO = False, revised = False, err_bars = False, hbeta_bin = False):    
    '''
    Purpose:
        This function creates plots for the individual measurements within each bin: HBeta Luminosity vs
        Stellar Mass (color Metallicity), O32 vs R23 (color Metallicity), Metallicity vs R23 and O32,
        Metallicity vs Stellar Mass.
           
    Params:
        fitspath --> string of the file path where data files are located.
        restrictMTO (OPTIONAL) --> if the mass turnover value should be held constant in the curve fit of
                                    Metallicity vs Mass, then restrictMTO = True.
        revised (OPTIONAL) --> if the revised data tables should be used, then revised = True.
        err_bars (OPTIONAL) --> if error bars for metallicity and temperature should be plotted, then 
                                err_bars = True.
        hbeta_bin (OPTIONAL) --> if the binning type is mass-LHbeta bins, then hbeta_bin = True.
        
    Returns:
        None
        
    Outputs:
        pdf_pages --> a pdf containing all of the plots for the individual spectra on separate pages.
    '''
    
    #Read in individual data tables
    indiv_derivedprops_tbl = asc.read(fitspath + filename_dict['indv_derived_prop'])
    indiv_props_tbl = asc.read(fitspath + filename_dict['indv_prop'])
    indiv_bininfo_tbl = asc.read(fitspath + filename_dict['indv_bin_info'])
    
    #Read in composite data tables
    bininfo_tbl = asc.read(fitspath + filename_dict['bin_info'])
    if revised == True:
        bin_derivedprops_tbl = asc.read(fitspath + filename_dict['bin_derived_prop_rev'])
        bin_valid_tbl = asc.read(fitspath + filename_dict['bin_valid_rev'])
    else:    
        bin_derivedprops_tbl = asc.read(fitspath + filename_dict['bin_derived_prop'])
        bin_valid_tbl = asc.read(fitspath + filename_dict['bin_valid'])


    #Read in individual data
    indiv_logM = indiv_props_tbl[indv_names0[3]].data
    indiv_logLHb = indiv_props_tbl[indv_names0[4]].data
    indiv_metal = indiv_derivedprops_tbl[temp_metal_names0[1]].data
    indiv_logR23 = indiv_derivedprops_tbl[indv_names0[1]].data
    indiv_logO32 = indiv_derivedprops_tbl[indv_names0[2]].data
    indiv_logTwoBeta = indiv_derivedprops_tbl[indv_names0[5]].data
    indiv_logThreeBeta = indiv_derivedprops_tbl[indv_names0[6]].data
    indiv_bin_detect = indiv_bininfo_tbl[bin_names0[2]].data
    
    #Read in composite data
    bin_logM = bininfo_tbl[bin_mzevolve_names0[2]].data      
    bin_metal = bin_derivedprops_tbl[temp_metal_names0[1]].data
    bin_detect_col = bin_valid_tbl[bin_names0[2]].data
    
    
    #Define detection and non-detection (with reliable 5007) arrays for bins
    bin_detect = np.where(bin_detect_col == 1.0)[0]
    bin_nondetect = np.where(bin_detect_col == 0.5)[0]
    
    #Define detection and non-detection (with reliable 5007) arrays for individual galaxies  
    indiv_detect, indiv_nondetect = get_indiv_detect(indiv_props_tbl, indiv_bin_detect, LHbeta_bins = hbeta_bin)

    
    #Define output file name
    if restrictMTO == True:
        pdf_pages = PdfPages(fitspath + filename_dict['indv_derived_prop'].replace('.tbl', 'constMTO.pdf'))
    else:
        pdf_pages = PdfPages(fitspath + filename_dict['indv_derived_prop'].replace('.tbl', '.pdf'))
    
    
    #HBeta Luminosity vs Mass ColorMap=Metallicity
    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot1 = ax1.scatter(indiv_logM[indiv_detect], indiv_logLHb[indiv_detect], 5.0, 
                        c=indiv_metal[indiv_detect], marker='*')
    plot1 = ax1.scatter(indiv_logM[indiv_nondetect], indiv_logLHb[indiv_nondetect], 5.0, 
                        facecolors = 'None', c=indiv_metal[indiv_nondetect], marker='^')
    cb = fig1.colorbar(plot1)
    cb.set_label('Metallicity')
    ax1.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1.set_ylabel('log(H$\\beta$ Luminosity)')
    ax1.set_title('H$\\beta$ Luminosity vs. Mass Colormap=Metallicity')
    fig1.set_size_inches(8, 8)
    fig1.savefig(pdf_pages, format='pdf')
    
    
    #log(O32) vs log(R23) ColorMap=Metallicity
    fig3, ax3 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot3 = ax3.scatter(indiv_logR23[indiv_detect], indiv_logO32[indiv_detect], 5.0, 
                        c=indiv_metal[indiv_detect], marker='*')
    plot3 = ax3.scatter(indiv_logR23[indiv_nondetect], indiv_logO32[indiv_nondetect], 5.0, 
                        facecolors = 'None', c=indiv_metal[indiv_nondetect], marker='^')
    cb = fig3.colorbar(plot3)
    cb.set_label('Metallicity')
    ax3.set_xlabel('log($R_{23}$)')
    ax3.set_ylabel('log($O_{32}$)')
    ax3.set_title('$O_{32}$ vs. $R_{23}$ Colormap=Metallicity')
    fig3.set_size_inches(8, 8)
    fig3.savefig(pdf_pages, format='pdf')    
    
    
    #Metallicity vs log(R23) and log(O32)
    fig5, ax5 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax5[0].scatter(indiv_logR23[indiv_detect], indiv_metal[indiv_detect], facecolors = 'None',
                   edgecolors = 'blue', label = 'log($R_{23}$)')
    ax5[1].scatter(indiv_logO32[indiv_detect], indiv_metal[indiv_detect], facecolors = 'None', 
                   edgecolors = 'red', label = 'log($O_{32}$)')
    ax5[0].scatter(indiv_logR23[indiv_nondetect], indiv_metal[indiv_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'blue', label = 'log($R_{23}$)', alpha = 0.5)
    ax5[1].scatter(indiv_logO32[indiv_nondetect], indiv_metal[indiv_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'red', label = 'log($O_{32}$)', alpha = 0.5)
    ax5[0].legend(loc = 'best')
    ax5[1].legend(loc = 'best')
    ax5[0].set_ylabel('Metallicity')
    ax5[0].set_title('Metallicity vs. $R_{23}$')
    ax5[1].set_title('Metallicity vs. $O_{32}$')
    fig5.set_size_inches(8, 8)
    fig5.savefig(pdf_pages, format='pdf')
    
    
    #log(R23) and log(O32) vs Mass
    fig7, ax7 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax7[0].scatter(indiv_logM[indiv_detect], indiv_logR23[indiv_detect], facecolors = 'None',
                   edgecolors = 'blue', label = 'log($R_{23}$)')    
    ax7[1].scatter(indiv_logM[indiv_detect], indiv_logO32[indiv_detect], facecolors = 'None', 
                   edgecolors = 'red', label = 'log($O_{32}$)')
    ax7[0].scatter(indiv_logM[indiv_nondetect], indiv_logR23[indiv_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'blue', label = 'log($R_{23}$)')
    ax7[1].scatter(indiv_logM[indiv_nondetect], indiv_logO32[indiv_nondetect], marker='^', 
                   facecolors = 'None', edgecolors = 'red', label = 'log($O_{32}$)')
    ax7[0].legend(loc = 'best')
    ax7[1].legend(loc = 'best')
    ax7[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax7[0].set_title('$R_{23}$ vs. Mass')
    ax7[1].set_title('$O_{32}$ vs. Mass')
    fig7.set_size_inches(8, 8)
    fig7.savefig(pdf_pages, format='pdf')
    
    
    #log(OII/HBeta) and log(OIII/HBeta) vs Mass
    fig9, ax9 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.95, wspace = 0.0)
    ax9[0].scatter(indiv_logM[indiv_detect], np.log10(indiv_logTwoBeta[indiv_detect]), 
                   facecolors = 'None', edgecolors = 'cyan', label = 'log($\\frac{OII}{H\\beta}$)')    
    ax9[1].scatter(indiv_logM[indiv_detect], np.log10(indiv_logThreeBeta[indiv_detect]), 
                   facecolors = 'None', edgecolors = 'orange', label = 'log($\\frac{OIII}{H\\beta}$)')
    ax9[0].scatter(indiv_logM[indiv_nondetect], np.log10(indiv_logTwoBeta[indiv_nondetect]), 
                   marker='^', facecolors = 'None', edgecolors = 'cyan', label = 'log($\\frac{OII}{H\\beta}$)')
    ax9[1].scatter(indiv_logM[indiv_nondetect], np.log10(indiv_logThreeBeta[indiv_nondetect]), 
                   marker='^', facecolors = 'None', edgecolors = 'orange', label = 'log($\\frac{OIII}{H\\beta}$)')
    ax9[0].legend(loc = 'best')
    ax9[1].legend(loc = 'best')
    ax9[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax9[0].set_title('$\\frac{OII}{H\\beta}$ vs. Mass')
    ax9[1].set_title('$\\frac{OIII}{H\\beta}$ vs. Mass')
    fig9.set_size_inches(8, 8)
    fig9.savefig(pdf_pages, format='pdf')
    
    
    #Metallcity vs Mass
    fig11, ax11 = plt.subplots()
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.14, top = 0.98, wspace = 0.0)
    
    #Andrews & Martini fit
    mass_range = np.arange(7.5, 10, 0.05)
    AM_relation = AM13(mass_range)
    ax11.plot(mass_range, AM_relation, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    
    #Plot individual detections and non-detections
    ax11.scatter(indiv_logM[indiv_detect], indiv_metal[indiv_detect], s = 15, facecolors = 'None',
                 edgecolors = 'blue', label = 'Individual Detections')
    ax11.scatter(indiv_logM[indiv_nondetect], indiv_metal[indiv_nondetect], s = 15, marker='^', 
                 facecolors = 'None', edgecolors = 'blue', label = 'Individual Non-Detections', alpha = 0.5)
    
    print('Number of individual sources plotted:', len(indiv_logM[indiv_detect]))
    
    #Plot bin detections and non-detections
    ax11.scatter(bin_logM[bin_detect], bin_metal[bin_detect], s = 25, color = 'red', label = 'Bin Detections')
    ax11.scatter(bin_logM[bin_nondetect], bin_metal[bin_nondetect], s = 25, color = 'red', marker = '^', 
                 label = 'Bin Non-Detections', alpha = 0.5)
    if err_bars == True:
        err_dict = extract_error_bars(fitspath)
        ax11.errorbar(bin_logM[bin_detect], bin_metal[bin_detect], yerr = err_dict['12+log(O/H)_lowhigh_error'], fmt = '.')
    
    #Fit bin detections and plot relation
    o11, o21, fail = curve_fitting(bin_logM[bin_detect], bin_metal[bin_detect], restrict_MTO = restrictMTO)
    if not fail:
        ax11.plot(mass_range, mass_metal_fit(mass_range, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
          
    ax11.legend(fontsize = 5, loc = 'upper left')    
    ax11.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax11.set_ylabel('12+log(O/H) $T_e$')    
    pdf_pages.savefig()
    
    pdf_pages.close()
    
    
    
    
    
def indiv_metal_mass_plot(Mbin_dict, MLHbbin_dict, restrictMTO = False, revised = False, err_bars = False):
    '''
    Purpose:
        This function creates a two-paneled Metallicity vs Mass plot containing individual measurements and 
        composite measurements. The left subplot is for mass bins, and the right subplot is for mass-LHbeta
        bins. It also plots a relation fit for the composite measurements.
           
    Params:
        Mbin_dict --> a dictionary containing the file path to the mass bin folder, composite and individual
                      masses, metallicities, detection arrays, non-detection arrays, and composite metallicity
                      errors.
        MLHbbin_dict --> a dictionary containing the file path to the mass bin folder, composite and 
                         individual masses, metallicities, detection arrays, non-detection arrays, and 
                         composite metallicity errors.
        restrictMTO (OPTIONAL) --> if the mass turnover value should be held constant in the curve fit of
                                    Metallicity vs Mass, then restrictMTO = True.
        revised (OPTIONAL) --> if the revised data tables should be used, then revised = True.
        err_bars (OPTIONAL) --> if error bars for metallicity and temperature should be plotted, then 
                                err_bars = True.
        
    Returns:
        None
        
    Outputs:
        Mbin_pdf_pages --> a pdf containing the two-paneled plot placed in the mass bin file path.
        MLHbbin_pdf_pages --> a pdf containing the two-paneled plot placed in the mass-LHbeta bin file path.
    '''
    
    #Define output file names
    if restrictMTO == True:
        Mbin_pdf_pages = PdfPages(Mbin_dict['path'] + 'combined_' + filename_dict['indv_derived_prop'].replace('.tbl', 'constMTO.pdf'))
        MLHbbin_pdf_pages = PdfPages(MLHbbin_dict['path'] + 'combined_' + filename_dict['indv_derived_prop'].replace('.tbl', 'constMTO.pdf'))
    else:
        Mbin_pdf_pages = PdfPages(Mbin_dict['path'] + 'combined_' + filename_dict['indv_derived_prop'].replace('.tbl', '.pdf'))
        MLHbbin_pdf_pages = PdfPages(MLHbbin_dict['path'] + 'combined_' + filename_dict['indv_derived_prop'].replace('.tbl', '.pdf'))
        
        
    #Metallcity vs Mass
    fig1, ax1 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.14, top = 0.98, wspace = 0.0)
    
    #Andrews & Martini fit
    mass_range = np.arange(7.5, 10, 0.05)
    AM_relation = AM13(mass_range)
    ax1[0].plot(mass_range, AM_relation, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    ax1[1].plot(mass_range, AM_relation, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    
    #Plot individual detections and non-detections
    ax1[0].scatter(Mbin_dict['indiv_logM'][Mbin_dict['indiv_detect']], Mbin_dict['indiv_metallicity'][Mbin_dict['indiv_detect']],
                   s = 15, facecolors = 'None', edgecolors = 'blue', label = 'Individual Detections')
    ax1[0].scatter(Mbin_dict['indiv_logM'][Mbin_dict['indiv_nondetect']], Mbin_dict['indiv_metallicity'][Mbin_dict['indiv_nondetect']], 
                   s = 15, marker='^', facecolors = 'None', edgecolors = 'blue', label = 'Individual Non-Detections', alpha = 0.5)
    ax1[1].scatter(MLHbbin_dict['indiv_logM'][MLHbbin_dict['indiv_detect']], MLHbbin_dict['indiv_metallicity'][MLHbbin_dict['indiv_detect']],
                   s = 15, facecolors = 'None', edgecolors = 'blue', label = 'Individual Detections')
    ax1[1].scatter(MLHbbin_dict['indiv_logM'][MLHbbin_dict['indiv_nondetect']], MLHbbin_dict['indiv_metallicity'][MLHbbin_dict['indiv_nondetect']],
                   s = 15, marker='^', facecolors = 'None', edgecolors = 'blue', label = 'Individual Non-Detections', alpha = 0.5)
    
    print('Number of mass bin individual sources plotted:', len(Mbin_dict['indiv_logM'][Mbin_dict['indiv_detect']]))
    print('Number of mass-LHbeta bin individual sources plotted:', len(MLHbbin_dict['indiv_logM'][MLHbbin_dict['indiv_detect']]))
    
    #Plot mass bin detections and non-detections
    ax1[0].scatter(Mbin_dict['composite_logM'][Mbin_dict['composite_detect']], Mbin_dict['composite_metallicity'][Mbin_dict['composite_detect']],
                   s = 25, color = 'red', label = 'Bin Detections')
    ax1[0].scatter(Mbin_dict['composite_logM'][Mbin_dict['composite_nondetect']], Mbin_dict['composite_metallicity'][Mbin_dict['composite_nondetect']],
                   s = 25, color = 'red', marker = '^', label = 'Bin Non-Detections', alpha = 0.5)
    
    #Plot mass-LHBeta bin detections and non-detections
    ax1[1].scatter(MLHbbin_dict['composite_logM'][MLHbbin_dict['composite_detect']], MLHbbin_dict['composite_metallicity'][MLHbbin_dict['composite_detect']], 
                   s = 25, color = 'red', label = 'Bin Detections')
    ax1[1].scatter(MLHbbin_dict['composite_logM'][MLHbbin_dict['composite_nondetect']], MLHbbin_dict['composite_metallicity'][MLHbbin_dict['composite_nondetect']], 
                   s = 25, color = 'red', marker = '^', label = 'Bin Non-Detections', alpha = 0.5)
    if err_bars == True:
        ax1[0].errorbar(Mbin_dict['composite_logM'][Mbin_dict['composite_detect']], Mbin_dict['composite_metallicity'][Mbin_dict['composite_detect']],
                        yerr = Mbin_dict['composite_metal_errors'], fmt = '.')
        ax1[1].errorbar(MLHbbin_dict['composite_logM'][MLHbbin_dict['composite_detect']], MLHbbin_dict['composite_metallicity'][MLHbbin_dict['composite_detect']],
                        yerr = MLHbbin_dict['composite_metal_errors'], fmt = '.')
     
    #Fit bin detections and plot relation    
    o11, o21, Mbin_fail = curve_fitting(Mbin_dict['composite_logM'][Mbin_dict['composite_detect']], 
                                        Mbin_dict['composite_metallicity'][Mbin_dict['composite_detect']], 
                                        restrict_MTO = restrictMTO)
    o12, o22, MLHbbin_fail = curve_fitting(MLHbbin_dict['composite_logM'][MLHbbin_dict['composite_detect']], 
                                           MLHbbin_dict['composite_metallicity'][MLHbbin_dict['composite_detect']], 
                                           restrict_MTO = restrictMTO)
    if Mbin_fail == False:
        ax1[0].plot(mass_range, mass_metal_fit(mass_range, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
    if MLHbbin_fail == False:
        ax1[1].plot(mass_range, mass_metal_fit(mass_range, *o12), alpha = 0.5, color = 'red', label = 'Our Fit')
          
    ax1[0].legend(title = '$M_\star$ Bins', fontsize = 5, loc = 'upper left')    
    ax1[1].legend(title = '$M_\star$-LH$\\beta$ Bins', fontsize = 5, loc = 'upper left')
    ax1[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1[0].set_ylabel('12+log(O/H) $T_e$')
    ax1[1].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    Mbin_pdf_pages.savefig()
    MLHbbin_pdf_pages.savefig()
    
    Mbin_pdf_pages.close()
    MLHbbin_pdf_pages.close()
    
    
    
    
def get_indiv_detect(indiv_props_tbl, bin_detect, LHbeta_bins = False):
    '''
    Purpose:
        This function creates index arrays of the individual detections and non-detections based on
        valid individual emission lines (OIII5007, OII, and HBETA) and bin detections and non-detections.
           
    Params:
        OIII5007 --> an array containing individual galaxies' OIII5007 flux values (length = # of galaxies).
        OII --> an array containing individual galaxies' OII flux values (length = # of galaxies).
        HBETA --> an array containing individual galaxies' HBETA flux values (length = # of galaxies).
        logLHb --> an array containing individual galaxies' HBeta Luminosity values (length = # of galaxies).
        LHbeta_bins (OPTIONAL) --> if the binning type is mass-LHbeta bins, then LHbeta_bins = True.
        
    Returns:
        combined_detect --> a numpy array of detection indices that pass all detection conditions.
        combined_nondetect --> a numpy array of non-detection indices that pass all non-detection conditions.
    '''
    
    #Read in individual emission lines and 
    logLHb = indiv_props_tbl[indv_names0[4]].data
    HBETA = indiv_props_tbl['HBETA_Flux_Observed'].data
    OII = indiv_props_tbl['OII_3727_Flux_Observed'].data
    OIII5007 = indiv_props_tbl['OIII_5007_Flux_Observed'].data
    
    #Get detection and non-detection indices
    OIII5007_idx = np.where((np.isfinite(OIII5007) == True) & (OIII5007 >= 1e-18) & (OIII5007 <= 1e-15))[0]
    OII_idx = np.where((np.isfinite(OII) == True) & (OII >= 1e-18) & (OII <= 1e-15))[0]
    HBETA_idx = np.where((np.isfinite(HBETA) == True) & (HBETA >= 1e-18) & (HBETA <= 1e-15))[0]
    detect_idx = np.where(bin_detect == 1.0)[0]
    non_detect_idx = np.where(bin_detect == 0.5)[0]
    
    combined_detect = set(OIII5007_idx) & set(OII_idx) & set(HBETA_idx) & set(detect_idx)
    combined_nondetect = set(OIII5007_idx) & set(OII_idx) & set(HBETA_idx) & set(non_detect_idx)
    if LHbeta_bins == True:
        LHb_idx = np.where(logLHb > 0)[0]
        combined_detect &= set(LHb_idx)
        combined_nondetect &= set(LHb_idx)
        
    print('combined_detect:', len(combined_detect))
    print('combined_nondetect:', len(combined_nondetect))
        
    return np.array(list(combined_detect)), np.array(list(combined_nondetect))
    
    
    
    
    
    
    
    
    
    
    
    