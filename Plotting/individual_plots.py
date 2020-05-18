import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from Metallicity_Stack_Commons.column_names import indv_names0, temp_metal_names0, bin_mzevolve_names0
from Metallicity_Stack_Commons.column_names import bin_names0, filename_dict, npz_filename_dict
from .relation_fitting import curve_fitting, mass_metal_fit, mass_metal_fit_constMTO, extract_error_bars




def indiv_derived_props_plots(fitspath, restrict_MTO = False, revised = False, err_bars = False, hbeta_bin = False):    
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
    
    indiv_derivedprops_tbl = asc.read(fitspath + filename_dict['indv_derived_prop'])
    indiv_props_tbl = asc.read(fitspath + filename_dict['indv_prop'])
    indiv_bininfo_tbl = asc.read(fitspath + filename_dict['indv_bin_info'])
    
    bin_valid_tbl = asc.read(fitspath + filename_dict['bin_valid'])
    bininfo_tbl = asc.read(fitspath + filename_dict['bin_info'])
    if revised == True:
        bin_derivedprops_tbl = asc.read(fitspath + filename_dict['bin_derived_prop_rev'])
    else:    
        bin_derivedprops_tbl = asc.read(fitspath + filename_dict['bin_derived_prop'])


    ##individual galaxy data, i.e. comes from MT_ascii
    indiv_logM = indiv_props_tbl[indv_names0[3]].data
    indiv_logLHb = indiv_props_tbl[indv_names0[4]].data
    
    indiv_HBETA = indiv_props_tbl['HBETA_Flux_Observed'].data
    indiv_OII = indiv_props_tbl['OII_3727_Flux_Observed'].data
    indiv_OIII5007 = indiv_props_tbl['OIII_5007_Flux_Observed'].data
    
    indiv_metal = indiv_derivedprops_tbl[temp_metal_names0[1]].data
    indiv_logR23 = indiv_derivedprops_tbl[indv_names0[1]].data
    indiv_logO32 = indiv_derivedprops_tbl[indv_names0[2]].data
    indiv_logTwoBeta = indiv_derivedprops_tbl[indv_names0[5]].data
    indiv_logThreeBeta = indiv_derivedprops_tbl[indv_names0[6]].data
    indiv_bin_detect = indiv_bininfo_tbl[bin_names0[2]].data
    
    
    ##bin data, i.e. comes from either mass or massLHb specific files
    bin_logM = bininfo_tbl[bin_mzevolve_names0[2]].data      
    bin_metal = bin_derivedprops_tbl[temp_metal_names0[1]].data
    bin_detect_col = bin_valid_tbl[bin_names0[2]].data
    
    
    ##detection determinations  
    #bins
    bin_detect = np.where(bin_detect_col == 1.0)[0]
    bin_nondetect = np.where(bin_detect_col == 0.5)[0]
    
    #individual    
    if hbeta_bin == True:
        indiv_detect = np.where((indiv_bin_detect == 1.0) & (np.isfinite(indiv_OIII5007) == True) & 
                                (indiv_OIII5007 >= 1e-18) & (indiv_OIII5007 <= 1e-15) & (np.isfinite(indiv_OII) == True) & 
                                (indiv_OII >= 1e-18) & (indiv_OII <= 1e-15) & (np.isfinite(indiv_HBETA) == True) & 
                                (indiv_HBETA >= 1e-18) & (indiv_HBETA <= 1e-15) & (indiv_logLHb > 0))[0]
        indiv_nondetect = np.where((indiv_bin_detect == 0.5) & (np.isfinite(indiv_OIII5007) == True) & 
                                   (indiv_OIII5007 >= 1e-18) & (indiv_OIII5007 <= 1e-15) & (np.isfinite(indiv_OII) == True) & 
                                   (indiv_OII >= 1e-18) & (indiv_OII <= 1e-15) & (np.isfinite(indiv_HBETA) == True) & 
                                   (indiv_HBETA >= 1e-18) & (indiv_HBETA <= 1e-15) & (indiv_logLHb > 0))[0]
    else:
        indiv_detect = np.where((indiv_bin_detect == 1.0) & (np.isfinite(indiv_OIII5007) == True) &
                                (indiv_OIII5007 >= 1e-18) & (indiv_OIII5007 <= 1e-15) & (np.isfinite(indiv_OII) == True) & 
                                (indiv_OII>= 1e-18) & (indiv_OII<= 1e-15) & (np.isfinite(indiv_HBETA) == True) & 
                                (indiv_HBETA >= 1e-18) & (indiv_HBETA <= 1e-15))[0]
        indiv_nondetect = np.where((indiv_bin_detect == 0.5) & (np.isfinite(indiv_OIII5007) == True) & 
                                   (indiv_OIII5007 >= 1e-18) & (indiv_OIII5007 <= 1e-15) & (np.isfinite(indiv_OII) == True) & 
                                   (indiv_OII>= 1e-18) & (indiv_OII<= 1e-15) & (np.isfinite(indiv_HBETA) == True) & 
                                   (indiv_HBETA >= 1e-18) & (indiv_HBETA <= 1e-15))[0]


    if restrict_MTO == True:
        pdf_pages = PdfPages(fitspath + filename_dict['indv_derived_prop'].replace('.tbl', 'constMTO.pdf'))
    else:
        pdf_pages = PdfPages(fitspath + filename_dict['indv_derived_prop'].replace('.tbl', '.pdf'))
    
    
    ##HBeta Luminosity vs Mass ColorMap=Metallicity##
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
    
    
    
    ##O32 vs R23 ColorMap=Metallicity##
    fig3, ax3 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot3 = ax3.scatter(indiv_logR23[indiv_detect], indiv_logO32[indiv_detect], 5.0, 
                        c=indiv_metal[indiv_detect], marker='*')
    plot3 = ax3.scatter(indiv_logR23[indiv_nondetect], indiv_logO32[indiv_nondetect], 5.0, 
                        facecolors = 'None', c=indiv_metal[indiv_nondetect], marker='^')
    cb = fig3.colorbar(plot3)
    cb.set_label('Metallicity')
    ax3.set_xlabel('log($R_{23}$)')
    ax3.set_ylabel('$O_{32}$')
    ax3.set_title('$O_{32}$ vs. $R_{23}$ Colormap=Metallicity')
    fig3.set_size_inches(8, 8)
    fig3.savefig(pdf_pages, format='pdf')    
    
    
    
    ##Metallicity vs R23 and O32##
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
    
    
    
    ##R23 and O32 vs Mass##
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
    
    
    
    ##OII/HBeta and OIII/HBeta vs Mass##
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
    
    
    
    ##Metallcity vs Mass##
    fig11, ax11 = plt.subplots()
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.14, top = 0.98, wspace = 0.0)
    
    
    #Andrews&Martini fit
    mass_range = np.arange(7.5, 10, 0.05)
    AM_relation = 8.798 - np.log10(1 + ((10**8.901)/(10**mass_range))**0.640)
    ax11.plot(mass_range, AM_relation, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    
    
    #Individual detections and non-detections
    ax11.scatter(indiv_logM[indiv_detect], indiv_metal[indiv_detect], s = 15, facecolors = 'None',
                 edgecolors = 'blue', label = 'Individual Detections')
    ax11.scatter(indiv_logM[indiv_nondetect], indiv_metal[indiv_nondetect], s = 15, marker='^', 
                 facecolors = 'None', edgecolors = 'blue', label = 'Individual Non-Detections', alpha = 0.5)
    
    print('Number of individual sources plotted:', len(indiv_logM[indiv_detect]))
    
    #Bin detections and non-detections
    ax11.scatter(bin_logM[bin_detect], bin_metal[bin_detect], s = 25, color = 'red', label = 'Bin Detections')
    ax11.scatter(bin_logM[bin_nondetect], bin_metal[bin_nondetect], s = 25, color = 'red', marker = '^', 
                 label = 'Bin Non-Detections', alpha = 0.5)
    if err_bars == True:
        err_dict = extract_error_bars(fitspath)
        ax11.errorbar(bin_logM[bin_detect], bin_metal[bin_detect], yerr = err_dict['12+log(O/H)_lowhigh_error'], fmt = '.')
    
    
    o11, o21, fail = curve_fitting(bin_logM[bin_detect], bin_metal[bin_detect], restrict_MTO = False)
        
    if not fail:
        if restrict_MTO == False:
            ax11.plot(mass_range, mass_metal_fit(mass_range, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
        else:
            ax11.plot(mass_range, mass_metal_fit_constMTO(mass_range, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
          
    ax11.legend(fontsize = 5, loc = 'upper left')    
    ax11.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax11.set_ylabel('12+log(O/H) $T_e$')
   
        
    pdf_pages.savefig()
    
    pdf_pages.close()
    
    
    
def indiv_metal_mass_plot(Mbin_dict, MLHbbin_dict, restrict_MTO = False, revised = False, err_bars = False):
    if restrict_MTO == True:
        Mbin_pdf_pages = PdfPages(Mbin_dict['path'] + 'combined_' + filename_dict['indv_derived_prop'].replace('.tbl', 'constMTO.pdf'))
        MLHbbin_pdf_pages = PdfPages(MLHbbin_dict['path'] + 'combined_' + filename_dict['indv_derived_prop'].replace('.tbl', 'constMTO.pdf'))
    else:
        Mbin_pdf_pages = PdfPages(Mbin_dict['path'] + 'combined_' + filename_dict['indv_derived_prop'].replace('.tbl', '.pdf'))
        MLHbbin_pdf_pages = PdfPages(MLHbbin_dict['path'] + 'combined_' + filename_dict['indv_derived_prop'].replace('.tbl', '.pdf'))
        
        
    ##Metallcity vs Mass##
    fig1, ax1 = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.14, top = 0.98, wspace = 0.0)
    
    
    #Andrews&Martini fit
    mass_range = np.arange(7.5, 10, 0.05)
    AM_relation = 8.798 - np.log10(1 + ((10**8.901)/(10**mass_range))**0.640)
    ax1[0].plot(mass_range, AM_relation, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    ax1[1].plot(mass_range, AM_relation, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    
    
    #Individual detections and non-detections
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
    
    #Mass bin detections and non-detections
    ax1[0].scatter(Mbin_dict['composite_logM'][Mbin_dict['composite_detect']], Mbin_dict['composite_metallicity'][Mbin_dict['composite_detect']],
                   s = 25, color = 'red', label = 'Bin Detections')
    ax1[0].scatter(Mbin_dict['composite_logM'][Mbin_dict['composite_nondetect']], Mbin_dict['composite_metallicity'][Mbin_dict['composite_nondetect']],
                   s = 25, color = 'red', marker = '^', label = 'Bin Non-Detections', alpha = 0.5)
    #HBeta bin detections and non-detections
    ax1[1].scatter(MLHbbin_dict['composite_logM'][MLHbbin_dict['composite_detect']], MLHbbin_dict['composite_metallicity'][MLHbbin_dict['composite_detect']], 
                   s = 25, color = 'red', label = 'Bin Detections')
    ax1[1].scatter(MLHbbin_dict['composite_logM'][MLHbbin_dict['composite_nondetect']], MLHbbin_dict['composite_metallicity'][MLHbbin_dict['composite_nondetect']], 
                   s = 25, color = 'red', marker = '^', label = 'Bin Non-Detections', alpha = 0.5)
    if err_bars == True:
        ax1[0].errorbar(Mbin_dict['composite_logM'][Mbin_dict['composite_detect']], Mbin_dict['composite_metallicity'][Mbin_dict['composite_detect']],
                        yerr = Mbin_dict['composite_metal_errors'], fmt = '.')
        ax1[1].errorbar(MLHbbin_dict['composite_logM'][MLHbbin_dict['composite_detect']], MLHbbin_dict['composite_metallicity'][MLHbbin_dict['composite_detect']],
                        yerr = MLHbbin_dict['composite_metal_errors'], fmt = '.')
        
    o11, o21, Mbin_fail = curve_fitting(Mbin_dict['composite_logM'][Mbin_dict['composite_detect']], 
                                        Mbin_dict['composite_metallicity'][Mbin_dict['composite_detect']], 
                                        restrict_MTO = False)
    o12, o22, MLHbbin_fail = curve_fitting(MLHbbin_dict['composite_logM'][MLHbbin_dict['composite_detect']], 
                                           MLHbbin_dict['composite_metallicity'][MLHbbin_dict['composite_detect']], 
                                           restrict_MTO = False)
        
    if Mbin_fail == False and MLHbbin_fail == False:
        if restrict_MTO == False:
            ax1[0].plot(mass_range, mass_metal_fit(mass_range, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
            ax1[1].plot(mass_range, mass_metal_fit(mass_range, *o12), alpha = 0.5, color = 'red', label = 'Our Fit')
        else:
            ax1[0].plot(mass_range, mass_metal_fit_constMTO(mass_range, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
            ax1[1].plot(mass_range, mass_metal_fit_constMTO(mass_range, *o12), alpha = 0.5, color = 'red', label = 'Our Fit')
          
    ax1[0].legend(title = '$M_\star$ Bins', fontsize = 5, loc = 'upper left')    
    ax1[1].legend(title = '$M_\star$-LH$\\beta$ Bins', fontsize = 5, loc = 'upper left')
    ax1[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1[0].set_ylabel('12+log(O/H) $T_e$')
    ax1[1].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
   
    Mbin_pdf_pages.savefig()
    MLHbbin_pdf_pages.savefig()
    
    Mbin_pdf_pages.close()
    MLHbbin_pdf_pages.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    