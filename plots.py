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
    valid_file = fitspath2 + 'hbeta_revised_75_112_113_300_600_1444_1444_massbin_validation.tbl'
    bin_metal_file = fitspath2 + 'hbeta_revised_75_112_113_300_600_1444_1444_updated_massbin_derived_properties_metallicity.tbl'
    
    MT_ascii = asc.read(metal_Te_file)
    valid_table = asc.read(valid_file)
    bin_metal_table = asc.read(bin_metal_file)

    log_mass = MT_ascii['Log10(Mass)'].data
    LHbeta = MT_ascii['HBeta_Luminosity'].data
    ind_metal = MT_ascii['com_O_log'].data
    R23 = MT_ascii['R23'].data
    O32 = MT_ascii['O32'].data
    indiv_detect = MT_ascii['Individual Detections'].data
    bin_avg_mass = valid_table['mass_avg'].data
    bin_detect = valid_table['Detection'].data
    bin_metal = bin_metal_table['com_O_log'].data
    
    detection = np.where((indiv_detect == 1.0) & (LHbeta > 0))[0]
    detect = np.where(bin_detect == 1.0)[0]
    detect = detect[:-1]
    non_detect = np.where(bin_detect == 0.0)[0]

    pdf_pages = PdfPages(fitspath2 + 'individual_metal_plots_ConstantMto.pdf')
    

    ##HBeta Luminosity vs Mass ColorMap=Metallicity##
    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot1 = ax1.scatter(log_mass[detection], LHbeta[detection], 0.8, c=ind_metal[detection], marker='*')
    cb = fig1.colorbar(plot1)
    cb.set_label('Metallicity')
    ax1.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax1.set_ylabel('log(H$\\beta$ Luminosity)')
    ax1.set_title('H$\\beta$ Luminosity vs. Mass Colormap=Metallicity')
    fig1.set_size_inches(8, 8)
    fig1.savefig(pdf_pages, format='pdf')
    
    
    ##O32 vs R23 ColorMap=Metallicity##
    fig2, ax2 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')
    plot2 = ax2.scatter(R23[detection], O32[detection], 0.8, c=ind_metal[detection], marker='*')
    cb = fig2.colorbar(plot2)
    cb.set_label('Metallicity')
    ax2.set_xlabel('$R_{23}$')
    ax2.set_ylabel('$O_{32}$')
    ax2.set_title('$O_{32}$ vs. $R_{23}$ Colormap=Metallicity')
    fig2.set_size_inches(8, 8)
    fig2.savefig(pdf_pages, format='pdf')
    
    
    ##Metallicity vs R23 and O32##
    fig3, ax3 = plt.subplots()
    ax3.scatter(R23[detection], ind_metal[detection], facecolors = 'None', edgecolors = 'blue', label = '$R_{23}$')
    ax3.scatter(O32[detection], ind_metal[detection], facecolors = 'None', edgecolors = 'red', label = '$O_{32}$')
    ax3.legend(loc = 'best')
    ax3.set_ylabel('Metallicity')
    ax3.set_title('Metallicity vs. $R_{23}$ and $O_{32}$')
    fig3.set_size_inches(8, 8)
    fig3.savefig(pdf_pages, format='pdf')
    
    
    ##Metallcity vs Mass##
    fig4, ax4 = plt.subplots()
    plt.subplots_adjust(left = 0.12, right = 0.98, bottom = 0.12, top = 0.98)
    
    #Andrews&Martini fit
    mass = np.arange(7.5, 10, 0.05)
    y = 8.798 - np.log10(1 + ((10**8.901)/(10**mass))**0.640)
    ax4.plot(mass, y, color='g', linestyle = '-', alpha = 0.5, label = 'Andrews & Martini (2013)')
    
    #Individual detections
    ax4.scatter(log_mass[detection], ind_metal[detection], s = 15, facecolors = 'None', edgecolors = 'blue',
                label = 'Individual Detections')
    
    #HBeta bin detections and non-detections subplot
    ax4.scatter(bin_avg_mass[detect], bin_metal[detect], s = 25, color = 'red', label = 'Bin Detections')
    ax4.scatter(bin_avg_mass[non_detect], bin_metal[non_detect], s = 25, color = 'red', marker = '^',
                label = 'Bin Non-Detections')
    
    #Mass bin detections and non-detections subplot
    
    
    ax4.legend(loc = 'best')
    ax4.set_xlabel('log($\\frac{M_\star}{M_{\odot}}$)')
    ax4.set_ylabel('12+log(O/H) $T_e$')
    
    fail = False
    #p0 = [8.798, 8.901, 0.640]
    p0 = [8.798, 0.640]
    para_bounds = ((8.0, 0.0), (9.0, 1.0))
    try:
        o1, o2 = curve_fit(mass_metal_fit2, bin_avg_mass[detect], bin_metal[detect], p0 = p0, bounds = para_bounds)
        print(o1)
    except ValueError:
        print('fail')
        fail = True
        
    if not fail:
        ax4.plot(mass, mass_metal_fit2(mass, *o1), alpha = 0.5)
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















