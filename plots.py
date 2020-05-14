import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit 
from Metallicity_Stack_Commons.column_names import indv_names0, temp_metal_names0, bin_mzevolve_names0
from Metallicity_Stack_Commons.column_names import bin_names0, filename_dict, npz_filename_dict



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
        
    if err_bars == True:
        derivedprops_errors = np.load(fitspath + npz_filename_dict['der_prop_errors'])
        
        metal_errors = derivedprops_errors['12+log(O/H)_lowhigh_error']
        
        metal_low_err = [metal_errors[ii][0] for ii in range(len(metal_errors))]
        metal_high_err = [metal_errors[ii][1] for ii in range(len(metal_errors))]
        
        metal_lowhigh_err = [metal_low_err, metal_high_err]


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
        ax11.errorbar(bin_logM[bin_detect], bin_metal[bin_detect], yerr = metal_lowhigh_err, fmt = '.')
    
    
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
        o11, o21 = curve_fit(fit, bin_logM[bin_detect], bin_metal[bin_detect], p0 = p0, bounds = para_bounds)
        print(o11)
    except ValueError:
        print('Failed curve fitting')
        fail = True
        
    if not fail:
        if restrict_MTO == False:
            ax11.plot(mass_range, mass_metal_fit(mass_range, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
        else:
            ax11.plot(mass_range, mass_metal_fit2(mass_range, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
          
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
            ax1[0].plot(mass_range, mass_metal_fit2(mass_range, *o11), alpha = 0.5, color = 'red', label = 'Our Fit')
            ax1[1].plot(mass_range, mass_metal_fit2(mass_range, *o12), alpha = 0.5, color = 'red', label = 'Our Fit')
          
    ax1[0].legend(title = '$M_\star$ Bins', fontsize = 5, loc = 'upper left')    
    ax1[1].legend(title = '$M_\star$-LH$\\beta$ Bins', fontsize = 5, loc = 'upper left')
    ax1[0].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1[0].set_ylabel('12+log(O/H) $T_e$')
    ax1[1].set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
   
    Mbin_pdf_pages.savefig()
    MLHbbin_pdf_pages.savefig()
    
    Mbin_pdf_pages.close()
    MLHbbin_pdf_pages.close()
    
    
    
def curve_fitting(x_array, y_array, restrict_MTO = False):
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
        o11, o21 = curve_fit(fit, x_array, y_array, p0 = p0, bounds = para_bounds)
        print(o11)
    except ValueError:
        print('Failed curve fitting!')
        fail = True
        
    return o11, o21, fail
    
    
    
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









