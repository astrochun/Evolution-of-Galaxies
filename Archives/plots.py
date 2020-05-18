import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages 



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









