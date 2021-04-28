import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages 
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from astropy.io import ascii as asc
from astropy.io import fits
import astropy.units as u
#from Metallicity_Stack_Commons.plotting.balmer import HbHgHd_fits


#TEMPORARY FILE TO RUN THESIS PLOTS (TO MAKE CHANGING THINGS EASIER & NOT WORRY ABOUT BREAKING CODE)
#WILL INTEGRATE THIS INTO MAIN CODES LATER

data_path = 'C:\\Users\\carol\\Google Drive\\MZEvolve\\massbin\\20210224\\75_112_113_300_600_1444_1443\\'

########################################################################################################
#Fig 4
def zoom_in_4363():
    fitspath = data_path
    bin_file = np.load(fitspath + "bin_info.npz")
    em_table = asc.read(data_path + 'bin_emission_line_fit.MC.tbl', format='fixed_width_two_line')
    flux = bin_file['flux']
    wavelength = bin_file['wavelength']
    OIII4363_SN = em_table['OIII_4363_S/N'].data
    N_stack = ['75', '112', '113', '300', '600', '1444', '1443']
    out_pdf = fitspath + 'comp_spec_zoom_in_4363_thesis.pdf'
    pdf_pages = PdfPages(out_pdf)
    
    flux_fits_file = fitspath + 'composite_spectra.fits'
    Spect_1D, header = fits.getdata(flux_fits_file, header=True)
    dispersion = header['CDELT1']
    wave = header['CRVAL1'] + dispersion*np.arange(header['NAXIS1'])

    # Plotting stacked spectra
    xlim = [4325, 4370]
    fig, ax = plt.subplots(2, 3, sharex=True, sharey=True)
    y = 0
    x = -1
    for i in range(0, 6):
        working_wave4363 = 4363.21
        #working_waveHG = 4340.46
        x_idx4363 = np.where((wave >= (working_wave4363 - 100)) & (wave <= (working_wave4363 + 100)))[0]
        x_idx_mask4363 = np.where((wave >= (working_wave4363 - 100)) & (wave <= (working_wave4363 + 100)))[0]
        #x_idxHG = np.where((wave >= (working_waveHG - 100)) & (wave <= (working_waveHG + 100)))[0]
        #x_idx_maskHG = np.where((wave >= (working_waveHG - 100)) & (wave <= (working_waveHG + 100)))[0]
        
        y0 = Spect_1D[i]
        y_norm = y0 / 1e-17
        o1_4363, med0_4363, max0_4363 = get_gaussian_fit(working_wave4363, wave, y0, y_norm, x_idx4363, 
                                                         x_idx_mask4363, 'Single', 5)
        #o1_HG, med0_HG, max0_HG = get_gaussian_fit(working_waveHG, wave, y0, y_norm, x_idxHG, x_idx_maskHG, 'Balmer', 5)
        gauss0_4363 = gauss(wave, *o1_4363)
        #gauss0_HG = double_gauss(wave, *o1_HG)
        
        wavelen_idx = np.where((wavelength[i] > xlim[0]) & (wavelength[i] < xlim[1]))[0]
        max0 = np.max((flux[i][wavelen_idx] / 1e-17))
        
        if x < 2:
            x += 1
        else:
            x = 0
            y = 1
        ax[y, x].axvline(x=4363.21, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #OIII4363
        ax[y, x].text(4364, 0.3, '$\\mathrm{[O}$ III]$\\lambda4363$', rotation=90, verticalalignment='center', size=8)
        ax[y, x].axvline(x=4340.46, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #HG
        ax[y, x].text(4342, 0.3, '$\\mathrm{H\\gamma}$', rotation=90, verticalalignment='center', size=8)
        
        ax[y, x].plot(wavelength[i], np.zeros(len(wavelength[i])), color='k', linestyle='dashed', linewidth=0.5)
        
        OIII4363range = np.where((wave >= 4358.21) & (wave <= 4368.21))[0]
        ax[y, x].plot(wave[OIII4363range], gauss0_4363[OIII4363range], color='r', linestyle='solid', linewidth=0.5)
        
        #HGrange = np.where((wave >= 4324) & (wave <= 4355))[0]
        #ax[y, x].plot(wave[HGrange], gauss0_HG[HGrange], color='r', linestyle='solid', linewidth=0.5)
        
        ax[y, x].plot(wave, Spect_1D[i] / 1e-17, color='b', linestyle='solid', linewidth=0.2)
        
        an_text = f"N = {int(N_stack[i])}"
        if i == 3 or i == 5:
            sn_text = "ND"
            ax[y, x].annotate(sn_text, [0.92, 0.96], xycoords='axes fraction', ha='center', va='top', size=8)
        else:
            sn_text = f"S/N = {np.round_(OIII4363_SN[i], decimals=2)}"
            ax[y, x].annotate(sn_text, [0.78, 0.96], xycoords='axes fraction', ha='center', va='top', size=8)

        
        ax[y, x].tick_params(direction='in')
        ax[y, x].set_xlim(xlim)
        ax[y, x].set_ylim(-0.05, max0*1.1)
        ax[y, x].annotate(an_text, [0.17, 0.96], xycoords='axes fraction', ha='center', va='top', size=8) #ax[y, x].annotate(an_text, [0.37, 0.96], xycoords='axes fraction', ha='center', va='top', size=8)
        
    
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.text(0.5, 0.04, 'Wavelength ($\\mathrm{\\AA}$)', ha='center', va='center')
    fig.text(0.06, 0.5, 'Intensity ($10^{-17} \\mathrm{erg}\\ \\mathrm{s^{-1}}\\ \\mathrm{cm^{-2}}\\ \\mathrm{\\AA^{-1}}$)',
             ha='center', va='center', rotation='vertical') 
    pdf_pages.savefig()
    
    bin_file.close()
    pdf_pages.close()
    
    
    
def get_gaussian_fit(working_wave, x0, y0, y_norm, x_idx, x_idx_mask,
                     line_type, s2):
    med0 = np.median(y_norm[x_idx_mask])
    max0 = np.max(y_norm[x_idx]) - med0

    fail = 0
    if line_type == 'Single': 
        p0 = [working_wave, 1.0, max0, med0]  # must have some reasonable values
        para_bound = ((working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0)), 
                      (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0)))
        try:
            o1, o2 = curve_fit(gauss, x0[x_idx], y_norm[x_idx],
                               p0=p0, sigma=None, bounds=para_bound)
        except ValueError:
            fail = 1
            
    if line_type == 'Balmer': 
        p0 = [working_wave, 1.0, max0, med0, s2, -0.05 * max0]  # must have some reasonable values
        para_bound = ((working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0), 0.0, -max0), 
                      (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0), 25.0, 0))
        try:
            o1, o2 = curve_fit(double_gauss, x0[x_idx], y_norm[x_idx],
                               p0=p0, sigma=None, bounds=para_bound)
        except ValueError:
            fail = 1

    if line_type == 'Oxy2':
        p0 = [working_wave, 1.0, 0.75 * max0, med0, 1.0, max0] # must have some reasonable values
        para_bound = ((working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0), 0.0, 0.0), 
                      (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0), 10.0, 100.0))
        try:
            o1, o2 = curve_fit(oxy2_gauss, x0[x_idx], y_norm[x_idx], p0=p0, sigma=None, 
                               bounds=para_bound)
        except ValueError:
            fail = 1
   
    if not fail:
        return o1, med0, max0
    else:
        return None, med0, max0
    
    
########################################################################################################
    
#Fig 3
def comp_spectra():
    fitspath = data_path
    bin_file = np.load(fitspath + "bin_info.npz")
    flux = bin_file['flux']
    wavelength = bin_file['wavelength']
    out_pdf = fitspath + 'composite_spectra_full1.pdf'
    pdf_pages = PdfPages(out_pdf)

    # Plotting stacked spectra
    xlim = [3700, 5100]
    fig, ax_arr = plt.subplots()
    ax_arr.axvline(x=3726.18, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #OII3727
    plt.text(3734, 3e-17, '$\\mathrm{[O}$ II]$\\lambda3727}$', rotation=90, verticalalignment='center', size=8)
    
    ax_arr.axvline(x=3868.74, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #NeIII3869
    plt.text(3828, 3e-17, '$\\mathrm{[Ne}$ III]$\\lambda3869}$', rotation=90, verticalalignment='center', size=8)
    
    ax_arr.axvline(x=3888.65, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #HeI
    ax_arr.axvline(x=3889.05, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #HZeta
    plt.text(3893, 3e-17, '$\\mathrm{He}$ I, $\\mathrm{H}\\zeta$', rotation=90, verticalalignment='center', size=8)
    
    ax_arr.axvline(x=3967.51, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #NeIII3968
    plt.text(3931, 3e-17, '$\\mathrm{[Ne}$ III]$\\lambda3968$', rotation=90, verticalalignment='center', size=8)
    
    ax_arr.axvline(x=3970.07, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #HEpsilon
    plt.text(3977, 3e-17, '$\\mathrm{H\\epsilon}$', rotation=90, verticalalignment='center', size=8)
    
    ax_arr.axvline(x=4101.73, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #HD
    plt.text(4109, 3e-17, '$\\mathrm{H\\delta}$', rotation=90, verticalalignment='center', size=8)
    
    ax_arr.axvline(x=4340.46, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #HG
    plt.text(4298, 3e-17, '$\\mathrm{H\\gamma}$', rotation=90, verticalalignment='center', size=8)
    
    ax_arr.axvline(x=4363.21, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #OIII4363
    plt.text(4370, 3e-17, '$\\mathrm{[O}$ III]$\\lambda4363$', rotation=90, verticalalignment='center', size=8)
    
    ax_arr.axvline(x=4861.32, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #HB
    plt.text(4814, 3e-17, '$\\mathrm{H\\beta}$', rotation=90, verticalalignment='center', size=8)
    
    ax_arr.axvline(x=4958.91, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #OIII4959
    plt.text(4912, 3e-17, '$\\mathrm{[O}$ III]$\\lambda4959$', rotation=90, verticalalignment='center', size=8)
    
    ax_arr.axvline(x=5006.84, color='k', linestyle='dashed', alpha=0.5, linewidth=0.5) #OIII5007
    plt.text(5014, 3e-17, '$\\mathrm{[O}$ III]$\\lambda5007$', rotation=90, verticalalignment='center', size=8)
    
    ax_arr.plot(wavelength[0], np.zeros(len(wavelength[0])), color='k', linestyle='dashed', linewidth=0.5)
    ax_arr.plot(wavelength[0], flux[0], color='b', linestyle='solid', linewidth=0.2)
    ax_arr.set_xlim(xlim)
    plt.xlabel('Wavelength ($\\mathrm{\\AA}$)')
    plt.ylabel('Intensity ($10^{-17} \\mathrm{erg}\\ \\mathrm{s^{-1}}\\ \\mathrm{cm^{-2}}\\ \\mathrm{\\AA^{-1}}$)')
        
    pdf_pages.savefig()
    
    bin_file.close()
    pdf_pages.close()
    
    
########################################################################################################   

#Fig 2
result = asc.read("C:\\Users\\carol\\Google Drive\\MZEvolve\\results_deeps_revised.tbl")
mass = result['best.stellar.m_star'].data
deep_fields = asc.read("C:\\Users\\carol\\Google Drive\\MZEvolve\\magfiles\\deep_fields.mag")

phot_cols = deep_fields.colnames[2:5]
phot_cols = [str0 for str0 in phot_cols if '_err' not in str0]

no_mass_idx = np.where((mass <= 0) | (np.isfinite(mass) == False))[0]
mag_no_mass = -2.5*np.log10(deep_fields['cfht_I'][no_mass_idx].data)+16.4


def interpolation(filename, band, mag_no_mass, no_mass_idx):
    npz_file = np.load(filename)
    grid, avg_mass, valid_mag, valid_mass, standard_dev, N_gals = npz_file['grid'], npz_file['average_mass'], npz_file['valid_mag'], npz_file['valid_mass'], npz_file['standard_deviation'], npz_file['N_gals']
    gals_in_bin = np.where(N_gals > 0)[0]
    interp_data = interp1d(grid[gals_in_bin], avg_mass[gals_in_bin], fill_value="extrapolate")
    exclude = np.where(valid_mag <= 24.1)[0]
    nonzero = np.where(avg_mass != 0)[0]
    nonzero = nonzero[:-2]
    out_pdf = data_path + "interpolation_thesis.pdf"
    pdf_pages = PdfPages(out_pdf)
    
    plt.figure(figsize = (15, 5))
    plt.scatter(valid_mag[exclude], valid_mass[exclude], s = 10, alpha = 0.5)
    plt.scatter(grid[nonzero] + 0.125, avg_mass[nonzero], s = 50, color = "black")
    plt.errorbar(grid[nonzero] + 0.125, avg_mass[nonzero], yerr = standard_dev[nonzero], color = "black", fmt = "none")
    plt.plot(grid[nonzero] + 0.125, avg_mass[nonzero], grid[nonzero] + 0.125, interp_data(grid[nonzero]), 'r-')
    plt.xlabel('I-Band Magnitude', fontsize=22)
    plt.ylabel('log($\\frac{M_\star}{M_{\odot}}$)', fontsize=22)
    plt.tick_params(labelsize=15)
    
    print(len(np.where(valid_mag <= 24.1)[0]))
    
    print(max(valid_mag))
    
    npz_file.close()
    
    #pdf_pages.savefig("C:\\Users\\carol\\Documents\\UA\\Research\\Thesis\\Images\\interpolation.pdf")
    pdf_pages.savefig()
    pdf_pages.close()


########################################################################################################
    
#Fig TBD
def mass_bin_cut_offs(filename, band):
    npz_file = np.load(filename)
    grid, avg_mass, valid_mag, valid_mass, standard_dev, N_gals = npz_file['grid'], npz_file['average_mass'], npz_file['valid_mag'], npz_file['valid_mass'], npz_file['standard_deviation'], npz_file['N_gals']
    gals_in_bin = np.where(N_gals > 0)[0]
    interp_data = interp1d(grid[gals_in_bin], avg_mass[gals_in_bin], fill_value="extrapolate")
    exclude = np.where(valid_mag <= 24.1)[0]
    #nonzero = np.where(avg_mass != 0)[0]
    #nonzero = nonzero[:-2]
    ##out_pdf = data_path + "interpolation.pdf"
    ##pdf_pages = PdfPages(out_pdf)
    
    plt.figure(figsize = (15, 5))
    plt.scatter(valid_mag[exclude], valid_mass[exclude], s = 10, alpha = 0.5)
    plt.annotate(band, [0.95, 0.95], xycoords = "axes fraction", ha = "right", va = "top")
    #plt.scatter(grid[nonzero] + 0.125, avg_mass[nonzero], s = 50, color = "black")
    #plt.errorbar(grid[nonzero] + 0.125, avg_mass[nonzero], yerr = standard_dev[nonzero], color = "black", fmt = "none")
    #plt.plot(grid[nonzero] + 0.125, avg_mass[nonzero], grid[nonzero] + 0.125, interp_data(grid[nonzero]), 'r-')
    plt.xlabel('I-Band Magnitude', fontsize=50)
    plt.ylabel('log($\\frac{M_\star}{M_{\odot}}$)', fontsize=50)
    plt.tick_params(labelsize=12)
    plt.subplots_adjust(wspace=0.1, hspace=0.08)
    
    print(len(valid_mass[exclude]))
    print(len(gals_in_bin))
    print(len(np.array(interp_data(mag_no_mass))))
    
    npz_file.close()
    
    ##pdf_pages.savefig()
    #pdf_pages.close()
    
    
########################################################################################################
    
#Fig 5
def plotting_gaussian_curves():
    pdf_pages = PdfPages(data_path + 'curve_fitting_thesis.pdf')
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.set_figheight(4)
    fig.set_figwidth(10)
    x = np.arange(1, 100)
    xbar = 50.0
    s = 15.0
    a = 20.0
    c = 0.0
    singlecurve = gauss(x, xbar, s, a, c)

    # Balmer Emission Lines
    x = np.arange(1, 100)
    xbar = 50.0
    s1 = 5.0 #15
    a1 = 20.0
    s2 = 15.0 #25
    a2 = -2.0
    c = 0.0

    doublecurve = double_gauss(x, xbar, s1, a1, c, s2, a2)

    positive = gauss(x, xbar, s1, a1, doublecurve[0])
    negative = gauss(x, xbar, s2, a2, c)

    # Oxygen Two Line
    x = np.arange(1, 100)
    xbar = 40.0
    s1 = 8.0
    a1 = 20.0
    s2 = 8.0
    a2 = 30.0
    c = 0.0

    oxycurve = oxy2_gauss(x, xbar, s1, a1, c, s2, a2)

    xbar3 = 40.0
    xbar4 = 63.5
    s3 = 8.0
    a3 = 20.0

    s4 = 8.0
    a4 = 30.0

    positive1 = gauss(x, xbar3, s3, a3, oxycurve[0])
    positive2 = gauss(x, xbar4, s4, a4, oxycurve[0])

    ax1.plot(x, singlecurve, 'k')
    ax2.plot(x, doublecurve, 'k')
    ax2.plot(x, positive, c='C0', linestyle='--')
    ax2.plot(x, negative, 'r', linestyle='--')
    ax3.plot(x, oxycurve, 'k')
    ax3.plot(x, positive1, c='C0', linestyle='--')
    ax3.plot(x, positive2, c='C0', linestyle='--')
    ax1.set_yticklabels([])
    ax2.set_yticklabels([])
    ax1.set_ylim(-3, 25.5)
    ax2.set_ylim(-3, 21.5)
    ax3.set_ylim(-3, 31.5)
    ax3.set_yticklabels([])
    ax1.set_title('Single Gaussian Fitting')
    ax2.set_title('Balmer Fitting')
    ax3.set_title('[OII]$\lambda3727$ Fitting')
    txt1 = '(a)'
    txt2 = '(b)'
    txt3 = '(c)'
    ax1.annotate(txt1, [0.95, 0.95], xycoords='axes fraction', va='top',
                 ha='right', fontsize='10')
    ax2.annotate(txt2, [0.95, 0.95], xycoords='axes fraction', va='top',
                 ha='right', fontsize='10')
    ax3.annotate(txt3, [0.95, 0.95], xycoords='axes fraction', va='top',
                 ha='right', fontsize='10')

    pdf_pages.savefig()
    pdf_pages.close()
    
    
def gauss(x: np.ndarray, xbar: float, s: float, a: float, c: float) \
        -> np.ndarray:
    return a * np.exp(-(x - xbar) ** 2 / (2 * s ** 2)) + c


def double_gauss(x: np.ndarray, xbar: float, s1: float, a1: float, c: float,
                 s2: float, a2: float) -> np.array:
    return a1 * np.exp(-(x - xbar) ** 2 / (2 * s1 ** 2)) + c + \
           a2 * np.exp(-(x - xbar) ** 2 / (2 * s2 ** 2))


def oxy2_gauss(x: np.ndarray, xbar: float, s1: float, a1: float, c: float,
               s2: float, a2: float) -> np.ndarray:
    con1 = 72.0/45.0
    return a1 * np.exp(-(x - xbar) ** 2 / (2 * s1 ** 2)) + c + \
           a2 * np.exp(-(x - (xbar * con1)) ** 2 / (2 * s2 ** 2))
    

########################################################################################################
           
#Fig 7    
def Z_vs_mass():
    # Read in data tables
    pdf_pages = PdfPages(data_path + 'MZ_plot_thesis.pdf')
    metal_table = asc.read(data_path + 'bin_derived_properties.MC.dustcorr.tbl', format='fixed_width_two_line')
    valid_table = asc.read(data_path + 'bin_validation.revised.tbl', format='fixed_width_two_line')
    bin_table = asc.read(data_path + 'bin_info.tbl', format='fixed_width_two_line')
    
    # Read in composite values
    com_O_log = metal_table['12+log(O/H)'].data
    logM_avg = bin_table['logM_avg'].data
    detection = valid_table['Detection'].data
    
    # Define detection and non-detection (with reliable 5007) arrays
    detect = np.where(detection == 1)[0]
    non_detect = np.where(detection == 0.5)[0]
    
    fig11, ax11 = plt.subplots()
    plt.subplots_adjust(top=0.99)
    
    # Andrews & Martini fit
    mass_range = np.arange(8.2, 9.9, 0.05)
    AM_relation = 8.798 - np.log10(1 + ((10**8.901)/(10**mass_range))**0.640)
    ax11.plot(mass_range, AM_relation, color='k', linestyle='dotted', label='A & M (2013)')
        
    # Plot bin detections and non-detections
    ax11.scatter(logM_avg[detect], com_O_log[detect], color='b', s=25, label='Detections')
    ax11.scatter(logM_avg[non_detect], com_O_log[non_detect], color='b', s=25, marker='^', 
                 label='Non-Detections')
    err_dict = extract_error_bars()
    ax11.errorbar(logM_avg[detect], com_O_log[detect], yerr=err_dict['12+log(O/H)_error'], 
                  color='b', fmt='.')
    
    # Fit bin detections and plot relation
    o11, o21, fail = curve_fitting(logM_avg[detect], com_O_log[detect], restrict_MTO=True)
    print("Fitting results [asm, gamma]:", o11)
    if not fail:
        ax11.plot(mass_range, mass_metal_fit(mass_range, *o11), alpha=0.5, color='b', 
                  linestyle='dashdot', label='Best Fit')
          
    ax11.legend(fontsize=8, loc='upper left')    
    ax11.set_xlabel('log($M_{\star}/M_{\odot}$)', fontsize=12)
    ax11.set_ylabel('12+log(O/H) $T_e$', fontsize=12)    
    pdf_pages.savefig()
    
    pdf_pages.close()
    
    

def extract_error_bars():    
    der_prop_err = np.load(data_path + 'derived_properties_errors.dustcorr.npz')
    der_prop_file = asc.read(data_path + 'bin_derived_properties.MC.dustcorr.tbl', format='fixed_width_two_line')
    valid_tbl = asc.read(data_path + 'bin_validation.revised.tbl', format='fixed_width_two_line')
        
    Te_err = der_prop_err['T_e_error']        
    metal_err = der_prop_err['12+log(O/H)_error']
    T_e = der_prop_file['T_e'].data
    detect = np.where(valid_tbl['Detection'].data == 1.0)[0]
        
    Te_low_err = -1*np.log10(1 - Te_err[:,0]/T_e[detect])
    Te_high_err = np.log10(1 + Te_err[:,1]/T_e[detect])
    
    err_dict = {'T_e_error': [Te_low_err, Te_high_err], 
                '12+log(O/H)_error': [metal_err[:,0], metal_err[:,1]]}
    
    der_prop_err.close()
    
    return err_dict



def curve_fitting(x_array, y_array, restrict_MTO=False):
    fail = False
    if not restrict_MTO:
        p0 = [8.798, 0.640, 8.901]
        para_bounds = ((8.0, 0.0, 8.0), (9.0, 1.0, 9.5))
    else:
        p0 = [8.798, 0.640]
        para_bounds = ((8.0, 0.0), (9.0, 1.0))
    try:
        o11, o21 = curve_fit(mass_metal_fit, x_array, y_array, p0=p0, bounds=para_bounds)
    except ValueError:
        print('Failed curve fitting!')
        fail = True
    print(o11)
    
    return o11, o21, fail


def mass_metal_fit(mass, a, g, b=8.901):    
    return a - np.log10(1 + ((10**b)/(10**mass))**g)


########################################################################################################
    
#Fig 6
def run_HbHgHd_plots():
    #HbHgHd_fits(data_path, out_pdf_prefix='HbHgHd_fits', use_revised=True)
    out_pdf_prefix = 'HbHgHd_fits'
    use_revised = True
    n_rows = 3
    n_cols = 3
    scalefact = 1e-17
    verbose = False
    
    comp_spec_file = data_path + 'composite_spectra.fits'
    stack2D, header = fits.getdata(comp_spec_file, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])

    if not use_revised:
        astropy_table_file = data_path + 'bin_emission_line_fit.tbl'
        out_pdf = data_path + out_pdf_prefix + '_thesis.pdf'
    else:
        astropy_table_file = data_path + 'bin_emission_line_fit.MC.tbl'
        out_pdf = data_path + out_pdf_prefix + '_thesis.MC.pdf'

    astropy_table = asc.read(astropy_table_file)
    ID = astropy_table['bin_ID'].data
    pdf_pages = PdfPages(out_pdf)

    for ii in range(len(ID)):
        if ii % n_rows == 0:
            fig, ax_arr = plt.subplots(nrows=n_rows, ncols=n_cols,
                                       squeeze=False)

        y0 = stack2D[ii]
        y_norm = y0/scalefact

        Hb_dict = extract_fit(astropy_table[ii], 'HBETA', balmer=True, verbose=False)
        Hg_dict = extract_fit(astropy_table[ii], 'HGAMMA', balmer=True, verbose=False)
        Hd_dict = extract_fit(astropy_table[ii], 'HDELTA', balmer=True, verbose=False)

        wave_beta  = 4861.32
        wave_gamma = 4340.46
        wave_delta = 4101.73

        # Beta
        Hb_fit_dict = fitting_result(wave, y_norm, wave_beta, **Hb_dict,
                                     use_revised=use_revised)

        # Gamma
        Hg_fit_dict = fitting_result(wave, y_norm, wave_gamma, **Hg_dict,
                                     use_revised=use_revised)

        # Delta
        Hd_fit_dict = fitting_result(wave, y_norm, wave_delta, **Hd_dict,
                                     use_revised=use_revised)

        # Calculate E(B-V)
        EBV_HgHb = compute_EBV(Hg_fit_dict['flux_gauss']/Hb_fit_dict['flux_gauss'],
                               source='HgHb', verbose=verbose)
        EBV_HdHb = compute_EBV(Hd_fit_dict['flux_gauss']/Hb_fit_dict['flux_gauss'],
                               source='HdHb', verbose=verbose)

        row = ii % n_rows

        # Label in upper left once
        ax_arr[row][0].annotate(f'ID: {ID[ii]}', [0.05, 0.95], va='top', ha='left',
                                xycoords='axes fraction', fontsize='8')

        # The below code could be refactored or simplified
        txt0 = "+$\\sigma$: " + str(np.round_(Hb_dict['line_fit'][1], decimals=3)) + \
               "; -$\\sigma$: " + str(np.round_(Hb_dict['line_fit_neg'][1], decimals=3)) + "\n"
        txt0 += "$F_G$: " + str(np.round_(Hb_fit_dict['flux_gauss'], decimals=3)) + \
                "; $F_S$: " + str(np.round_(Hb_fit_dict['flux_spec'], decimals=3)) + "\n"

        ax_arr[row][2].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][2].plot(wave, Hb_fit_dict['gauss'], 'm', linewidth=0.25, label='Beta Fit')
        ax_arr[row][2].set_xlim(4846.32, 4876.32)

        ax_arr[row][2].annotate(txt0, [0.95, 0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize='8')
        ax_arr[row][2].plot(wave[Hb_fit_dict['idx_sig']], Hb_fit_dict['residual'],
                            'r', linestyle='dashed', linewidth=0.2, label='Residuals')

        txt1 = "+$\\sigma$: " + str(np.round_(Hg_dict['line_fit'][1], decimals=3)) + \
               "; -$\\sigma$: " + str(np.round_(Hg_dict['line_fit_neg'][1], decimals=3)) + "\n"
        txt1 += "$F_G$: " + str(np.round_(Hg_fit_dict['flux_gauss'], decimals=3)) + \
                "; $F_S$: " + str(np.round_(Hg_fit_dict['flux_spec'], decimals=3)) + "\n"
        HgHb = Hg_fit_dict['flux_gauss']/Hb_fit_dict['flux_gauss']
        txt1 += "H$\\gamma$/H$\\beta$: " + str(np.round_(HgHb, decimals=2)) + "; E(B-V): " + \
                str(np.round_(EBV_HgHb, decimals=2))

        ax_arr[row][1].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][1].plot(wave, Hg_fit_dict['gauss'], 'm', linewidth=0.25,
                            label='Gamma Fit')
        ax_arr[row][1].set_xlim(4325.46, 4355.46)

        ax_arr[row][1].annotate(txt1, [0.95, 0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize='8')
        ax_arr[row][1].plot(wave[Hg_fit_dict['idx_sig']], Hg_fit_dict['residual'],
                            'r', linestyle='dashed', linewidth=0.2, label='Residuals')

        txt2 = "+$\\sigma$: " + str(np.round_(Hd_dict['line_fit'][1], decimals=3)) + \
               "; -$\\sigma$: " + str(np.round_(Hd_dict['line_fit_neg'][1], decimals=3)) + "\n"
        txt2 += "$F_G$: " + str(np.round_(Hd_fit_dict['flux_gauss'], decimals=3)) + \
                "; $F_S$: " + str(np.round_(Hd_fit_dict['flux_spec'], decimals=3)) + "\n"
        HdHb = Hd_fit_dict['flux_gauss']/Hb_fit_dict['flux_gauss']
        txt2 += "H$\\delta$/H$\\beta$: " + str(np.round_(HdHb, decimals=2)) + "; E(B-V): " + \
                str(np.round_(EBV_HdHb, decimals=2))

        ax_arr[row][0].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][0].plot(wave, Hd_fit_dict['gauss'], 'm', linewidth=0.25,
                            label='Delta Fit')
        ax_arr[row][0].set_xlim(4086.73, 4116.73)

        ax_arr[row][0].set_ylim(0, 1.6)
        
        ax_arr[row][0].annotate(txt2, [0.95, 0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize='8')
        ax_arr[row][0].plot(wave[Hd_fit_dict['idx_sig']], Hd_fit_dict['residual'],
                            'r', linestyle='dashed', linewidth=0.2, label='Residuals')
       
        ax_arr[row][0].set_yticklabels([0, 0.5, 1, 1.5])
        ax_arr[row][1].set_yticklabels([])
        ax_arr[row][2].set_yticklabels([])
        
        ax_arr[row][0].xaxis.set_major_locator(plt.MaxNLocator(4, prune='both'))
        ax_arr[row][1].xaxis.set_major_locator(plt.MaxNLocator(4, prune='both'))
        ax_arr[row][2].xaxis.set_major_locator(plt.MaxNLocator(4, prune='both'))
        
        if row != n_rows-1 and ii != stack2D.shape[0]-1:
            ax_arr[row][0].set_xticklabels([])
            ax_arr[row][1].set_xticklabels([])
            ax_arr[row][2].set_xticklabels([])
        else:
            ax_arr[row][0].set_xticklabels([4090, 4098, 4106, 4112])

        for col in range(n_cols):
            ax_arr[row][col].tick_params(direction='in')  # ticks on the inside

        if row == 1:
            ax_arr[row][0].set_ylabel('Intensity ($10^{-17} \\mathrm{erg}\\ \\mathrm{s^{-1}}\\ \\mathrm{cm^{-2}}\\ \\mathrm{\\AA^{-1}}$)',
                                      fontsize=12)

        if row == n_rows-1 or ii == stack2D.shape[0]-1:
            ax_arr[row][1].set_xlabel('Wavelength ($\\mathrm{\\AA}$)', fontsize=12)

        if ii % n_rows == n_rows-1:
            plt.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.99,
                                hspace=0.05, wspace=0.025)
            fig.savefig(pdf_pages, format='pdf')

    pdf_pages.close()
    
    
def fitting_result(wave: np.ndarray, y_norm: np.ndarray, lambda_cen: float,
                   line_fit: list, line_fit_neg: list, flux_gauss: float,
                   flux_spec: float, use_revised: bool = False) -> dict:
    dx = wave[2] - wave[1]

    fit_dict = dict()

    idx_sig = np.where(np.abs((wave - lambda_cen)) / line_fit[1] <= 2.5)[0]
    fit_dict['gauss'] = double_gauss(wave, *line_fit)
    fit_dict['negative'] = gauss(wave, *line_fit_neg)

    gauss0_diff = fit_dict['gauss'] - fit_dict['negative']
    y_norm_diff = y_norm[idx_sig] - fit_dict['negative'][idx_sig]

    # Residuals
    idx_sig_2 = np.where(np.abs((wave - lambda_cen)) / line_fit[1] <= 3.0)[0]
    fit_dict['residual'] = y_norm[idx_sig_2] - fit_dict['gauss'][idx_sig_2] + line_fit[3]

    fit_dict['idx_sig'] = idx_sig_2

    # Fluxes
    if not use_revised:
        fit_dict['flux_gauss'] = np.sum(gauss0_diff * dx)
        fit_dict['flux_spec']  = np.sum(y_norm_diff * dx)
    else:
        fit_dict['flux_gauss'] = flux_gauss
        fit_dict['flux_spec']  = flux_spec

    return fit_dict



def extract_fit(astropy_table, line_name, balmer=False, verbose=False):
    try:
        astropy_table[line_name + '_Center']
    except KeyError:
        return

    xbar = astropy_table[line_name + '_Center']
    sp   = astropy_table[line_name + '_Sigma']
    ap   = astropy_table[line_name + '_Norm']
    con  = astropy_table[line_name + '_Median']

    result_dict = dict()

    result_dict['line_fit'] = [xbar, sp, ap, con]

    result_dict['flux_gauss'] = astropy_table[line_name + '_Flux_Gaussian']
    result_dict['flux_spec']  = astropy_table[line_name + '_Flux_Observed']

    if balmer:
        sn = astropy_table[line_name + '_Abs_Sigma']
        an = astropy_table[line_name + '_Abs_Norm']

        result_dict['line_fit'] += [sn, an]

        result_dict['line_fit_neg'] = [xbar, sn, an, con]
    return result_dict



def compute_EBV(ratio, source='HgHb', zero_neg=True, verbose=False):
    if isinstance(ratio, list):
        raise TypeError("!!! Incorrect type for input [ratio].  Cannot be list !!!")

    if source not in ['HgHb', 'HdHb']:
        raise KeyError("!!! Incorrect [source]")
        
    line_name = ['OII_3727', 'HDELTA', 'HGAMMA', 'OIII_4363', 'HBETA', 'OIII_4958', 'OIII_5007']
    lambda0   = [3726.18, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]
    k_values = cardelli(lambda0 * u.Angstrom)
    k_dict = dict(zip(line_name, k_values))
        
    k_HBETA  = k_dict['HBETA']
    k_HGAMMA = k_dict['HGAMMA']
    k_HDELTA = k_dict['HDELTA']

    # For HgHb (default)
    ratio0 = 0.468
    k1 = k_HGAMMA

    if 'HdHb' in source:
        ratio0 = 0.259
        k1 = k_HDELTA

    EBV = -2.5 * np.log10(ratio/ratio0)/(k1 - k_HBETA)

    if zero_neg:
        if isinstance(EBV, float):
            if EBV < 0.0:
                EBV = 0.0
            return EBV
        else:
            if len(EBV.shape) == 1:
                neg_idx = np.where(EBV < 0.0)[0]
                if len(neg_idx) > 0:
                    EBV[neg_idx] = 0.0
                return EBV
            if len(EBV.shape) == 2:
                EBV_avg = np.average(EBV, axis=0)  # initial guess
                EBV_err, EBV_peak = compute_onesig_pdf(EBV, EBV_avg, usepeak=True)
                neg_idx = np.where(EBV_peak < 0.0)[0]
                if len(neg_idx) > 0:
                    EBV[neg_idx, :] -= EBV_peak[neg_idx].reshape((len(neg_idx), 1))
                    EBV_peak[neg_idx] = 0.0
                return EBV, EBV_peak
    else:
        return EBV
    
    
def cardelli(lambda0, R=3.1): #, extrapolate=False):
    t_lam = lambda0.to(u.nm).value

    ## Handles individual values, x
    if type(t_lam) == 'list':
        t_lam = np.array(t_lam)
    else:
        if isinstance(t_lam, (np.ndarray, np.generic)) == False:
            t_lam = np.array([t_lam])

    x = 1.0/(t_lam/1000.0) #in micron^-1

    k = np.zeros(np.size(t_lam), dtype=np.float64)

    mark = np.where((x <= 1.10) & (x >= 0.30))[0]
    if len(mark) > 0: k[mark] = ir_func_mw(x[mark], R)

    mark = np.where((x <= 3.30) & (x > 1.10))[0]
    if len(mark) > 0: k[mark] = opt_func_mw(x[mark], R)

    mark = np.where((x <= 8.00) & (x > 3.3))[0]
    if len(mark) > 0: k[mark] = uv_func_mw(x[mark], R)

    mark = np.where((x <= 10.00) & (x > 8.0))[0]
    if len(mark) > 0: k[mark] = fuv_func_mw(x[mark], R)

    k = k * R
    if np.size(x) == 1: k = k[0]
    return k 


def compute_onesig_pdf(arr0, x_val, usepeak=False, silent=True, verbose=False):
    import numpy as np

    if silent == False: print('### Begin compute_onesig_pdf | ' + systime())

    len0 = arr0.shape[0]  # arr0.shape[1] # Mod on 29/06/2016

    err = np.zeros((len0, 2))  # np.zeros((2,len0)) # Mod on 29/06/2016
    xpeak = np.zeros(len0)

    conf = 0.68269  # 1-sigma

    for ii in range(len0):
        test = arr0[ii]  # arr0[:,ii] # Mod on 29/06/2016
        good = np.where(np.isfinite(test) == True)[0]
        if len(good) > 0:
            v_low = np.percentile(test[good], 15.8655)
            v_high = np.percentile(test[good], 84.1345)

            xpeak[ii] = np.percentile(test[good], 50.0)
            if usepeak == False:
                t_ref = x_val[ii]
            else:
                t_ref = xpeak[ii]

            err[ii, 0] = t_ref - v_low
            err[ii, 1] = v_high - t_ref
    if silent == False: print('### End compute_onesig_pdf | ' + systime())
    return err, xpeak


def systime():
    import time
    return time.strftime("%d_%b_%Y_%H:%M:%S", time.localtime())


class TimerClass:
    from datetime import datetime as dt

    def __init__(self):
        self.start = 0
        self.stop = 0
        self.delta = 0
        self.format = ""

    def _start(self):
        self.start = self.dt.now()

    def _stop(self):
        self.stop = self.dt.now()
        self.delta = self.stop - self.start
        sec = self.delta.seconds
        HH = sec // 3600
        MM = (sec // 60) - (HH * 60)
        SS = sec - (HH * 3600) - (MM * 60)
        self.format = "%i hours  %i minutes  %i seconds" % (HH, MM, SS)
        
        
def uv_func_mw(x, R):
  
    Fa = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
    Fb = 0.2130*(x-5.9)**2 - 0.1207*(x-5.9)**3
    
    mark1 = np.where(x < 5.9)[0]
    if len(mark1) > 0:
        Fa[mark1] = 0.0
        Fb[mark1] = 0.0
        
    a = 1.752 - 0.316 * x - 0.104/((x-4.67)**2 + 0.341) + Fa
    b = -3.090 + 1.825 * x + 1.206/((x-4.62)**2+0.263) + Fb

    return a + b/R


def fuv_func_mw(x, R):

    a = -1.073 - 0.628 * (x-8) + 0.137*(x-8)**2 - 0.070*(x-8)**3
    b = 13.670 + 4.257*(x-8) - 0.420*(x-8)**2 + 0.374*(x-8)**3
    return a + b/R


def ir_func_mw(x, R):

    a = 0.574 * x**1.61
    b = -0.527 * x**1.61
    return a + b/R


def opt_func_mw(x, R):
    y = x - 1.82

    a = 1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4 + \
        0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7
    b = 1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4 - \
        0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7
    return a + b/R


########################################################################################################


def Te_vs_R():
    # Read in data tables
    pdf_pages = PdfPages(data_path + 'TevsR_plot_thesis.pdf')
    der_props_table = asc.read(data_path + 'bin_derived_properties.MC.dustcorr.tbl', format='fixed_width_two_line')
    valid_table = asc.read(data_path + 'bin_validation.revised.tbl', format='fixed_width_two_line')
    
    # Read in composite values
    logTe = np.log10(der_props_table['T_e'].data)
    logR = np.log10(der_props_table['R_composite'].data)
    detection = valid_table['Detection'].data
    
    # Define detection and non-detection (with reliable 5007) arrays
    detect = np.where(detection == 1)[0]
    non_detect = np.where(detection == 0.5)[0]
    
    fig1, ax1 = plt.subplots()
    plt.subplots_adjust(top=0.99)
    
    logR_range = np.logspace(-2.4, -1.5, 10)
    T_e = 13025 * (-logR_range - 0.92506) ** (-1 * 0.98062)
    #sorted_logR, sorted_logT_e = zip(*sorted(zip(logR_range, np.log10(T_e))))
        
    # Plot bin detections and non-detections
    ax1.scatter(logR[detect], logTe[detect], color='b', s=25, label='Detections')
    ax1.scatter(logR[non_detect], logTe[non_detect], color='b', s=25, marker=(3, 0, 135), label='Non-Detections')
    ax1.plot(logR_range, np.log10(T_e), color='k', linestyle='dashed', label='Eq. 3')
    err_dict = extract_error_bars()
    ax1.errorbar(logR[detect], logTe[detect], yerr=err_dict['T_e_error'], color='b', fmt='.')
          
    ax1.legend(fontsize=8, loc='upper left')    
    ax1.set_xlabel('log($R$)', fontsize=12)
    ax1.set_ylabel('log($T_e$)', fontsize=12)    
    pdf_pages.savefig()
    
    pdf_pages.close()
    
    
#comp_spectra()  
#zoom_in_4363()    
#interpolation("C:\\Users\\carol\\Google Drive\\MZEvolve\\revised_mag_mass_cfht_I.npz", phot_cols[2], mag_no_mass, no_mass_idx)
#mass_bin_cut_offs("C:\\Users\\carol\\Google Drive\\MZEvolve\\revised_mag_mass_cfht_I.npz", phot_cols[2])
#Z_vs_mass()
#plotting_gaussian_curves()
#run_HbHgHd_plots()
#Te_vs_R()
    
    
    
    
    
    
    
    
    
    
    