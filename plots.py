import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from scipy.optimize import curve_fit 
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





def derived_properties_plots():    
    metal_Te_file = fitspath2 + 'individual_derived_properties_metallicity.tbl'
    HB_valid_file = fitspath2 + 'hbeta_revised_75_112_113_300_600_1444_1444_massbin_validation.tbl'
    mass_valid_file = fitspath2 + 'revised_75_112_113_300_600_1444_1444_massbin_validation.tbl'
    HB_bin_metal_file = fitspath2 + 'hbeta_revised_75_112_113_300_600_1444_1444_updated_massbin_derived_properties_metallicity.tbl'
    mass_bin_metal_file = fitspath2 + 'revised_75_112_113_300_600_1444_1444_updated_massbin_derived_properties_metallicity.tbl'
    
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

    pdf_pages = PdfPages(fitspath2 + 'individual_metal_plots.pdf')
    
    
    
    
    ##HBeta Luminosity vs Mass ColorMap=Metallicity##
    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot1 = ax1.scatter(log_mass[HB_detection], LHbeta[HB_detection], 0.8, c=HB_ind_metal[HB_detection], marker='*')
    cb = fig1.colorbar(plot1)
    cb.set_label('Metallicity')
    ax1.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1.set_ylabel('log(H$\\beta$ Luminosity)')
    ax1.set_title('H$\\beta$ Luminosity Bins: H$\\beta$ Luminosity vs. Mass Colormap=Metallicity')
    fig1.set_size_inches(8, 8)
    fig1.savefig(pdf_pages, format='pdf')
    
    fig2, ax2 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot2 = ax2.scatter(log_mass[mass_detection], LHbeta[mass_detection], 0.8, c=mass_ind_metal[mass_detection], marker='*')
    cb = fig2.colorbar(plot2)
    cb.set_label('Metallicity')
    ax2.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax2.set_ylabel('log(H$\\beta$ Luminosity)')
    ax2.set_title('Mass Bins: H$\\beta$ Luminosity vs. Mass Colormap=Metallicity')
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
    ax3.set_title('H$\\beta$ Luminosity Bins: $O_{32}$ vs. $R_{23}$ Colormap=Metallicity')
    fig3.set_size_inches(8, 8)
    fig3.savefig(pdf_pages, format='pdf')
    
    fig4, ax4 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot4 = ax4.scatter(R23[mass_detection], O32[mass_detection], 0.8, c=mass_ind_metal[mass_detection], marker='*')
    cb = fig4.colorbar(plot4)
    cb.set_label('Metallicity')
    ax4.set_xlabel('$R_{23}$')
    ax4.set_ylabel('$O_{32}$')
    ax4.set_title('Mass Bins: $O_{32}$ vs. $R_{23}$ Colormap=Metallicity')
    fig4.set_size_inches(8, 8)
    fig4.savefig(pdf_pages, format='pdf')
    
    
    
    
    ##Metallicity vs R23 and O32##
    fig5, ax5 = plt.subplots()
    ax5.scatter(R23[HB_detection], HB_ind_metal[HB_detection], facecolors = 'None', edgecolors = 'blue', label = '$R_{23}$')
    ax5.scatter(O32[HB_detection], HB_ind_metal[HB_detection], facecolors = 'None', edgecolors = 'red', label = '$O_{32}$')
    ax5.legend(loc = 'best')
    ax5.set_ylabel('Metallicity')
    ax5.set_title('H$\\beta$ Luminosity Bins: Metallicity vs. $R_{23}$ and $O_{32}$')
    fig5.set_size_inches(8, 8)
    fig5.savefig(pdf_pages, format='pdf')
    
    fig6, ax6 = plt.subplots()
    ax6.scatter(R23[mass_detection], mass_ind_metal[mass_detection], facecolors = 'None', edgecolors = 'blue', label = '$R_{23}$')
    ax6.scatter(O32[mass_detection], mass_ind_metal[mass_detection], facecolors = 'None', edgecolors = 'red', label = '$O_{32}$')
    ax6.legend(loc = 'best')
    ax6.set_ylabel('Metallicity')
    ax6.set_title('Mass Bins: Metallicity vs. $R_{23}$ and $O_{32}$')
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
    #para_bounds = ((8.0, 0.0, 0.0), (9.0, 1.0, 1.0))
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
          
    ax7[0].legend(title = 'H$\\beta$ Lum Bins', fontsize = 5, loc = 'upper left')
    ax7[1].legend(title = 'Mass Bins', fontsize = 5, loc = 'upper left')
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



def find_outlier():
    metal_Te_file = fitspath2 + 'individual_derived_properties_metallicity.tbl'
    
    MT_ascii = asc.read(metal_Te_file)
    
    ID = MT_ascii['Source_ID'].data
    log_mass = MT_ascii['Log10(Mass)'].data
    indiv_detect = MT_ascii['Individual Detections'].data
    LHbeta = MT_ascii['HBeta_Luminosity'].data
    
    outlier = np.where((log_mass < 7.0) & (indiv_detect == 1.0) & (LHbeta > 0))[0]
    
    print(log_mass[outlier]) 
    print(ID[outlier]) 















