import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from pylab import subplots_adjust
from scipy.optimize import curve_fit 
from collections import OrderedDict

from Metallicity_Stack_Commons.analysis.fitting import movingaverage_box1D, gauss, double_gauss, oxy2_gauss
from Metallicity_Stack_Commons.analysis.fitting import rms_func, con1
from Metallicity_Stack_Commons import scalefact
from Metallicity_Stack_Commons.column_names import gauss_lines_names0, filename_dict, bin_names0
from .. import create_empty_dict



def get_gaussian_fit(working_wave, x0, y0, y_norm, x_idx, x_idx_mask, line_type, s2):
    med0 = np.median(y_norm[x_idx_mask])
    max0 = np.max(y_norm[x_idx]) - med0

    fail = 0
    if line_type == 'Single': 
        p0 = [working_wave, 1.0, max0, med0] # must have some reasonable values
        para_bound = ((working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0)), 
                      (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0)))
        try:
            o1, o2 = curve_fit(gauss, x0[x_idx], y_norm[x_idx], p0=p0, sigma=None, bounds=para_bound)
        except ValueError:
            print('fail')
            fail = 1
            
    if line_type == 'Balmer': 
        p0 = [working_wave, 1.0, max0, med0, s2, -0.05 * max0] # must have some reasonable values
        para_bound = ((working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0), 0.0, -max0), 
                      (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0), 25.0, 0))
        try:
            o1, o2 = curve_fit(double_gauss, x0[x_idx], y_norm[x_idx], p0=p0, sigma=None, bounds=para_bound)
        except ValueError:
            print('fail')
            fail = 1

    if line_type == 'Oxy2':
        p0 = [working_wave, 1.0, 0.75 * max0, med0, 1.0, max0] # must have some reasonable values
        para_bound = ((working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0), 0.0, 0.0), 
                      (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0), 10.0, 100.0))
        try:
            o1, o2 = curve_fit(oxy2_gauss, x0[x_idx], y_norm[x_idx], p0=p0, sigma=None, 
                               bounds=para_bound)
        except ValueError:
            print('fail')
            fail = 1
   
    if not fail:
        return o1, med0, max0
    else:
        return None, med0, max0



def zoom_gauss_plot(line_dict, col_names, pdf_pages, N, wave, Spect_1D, dispersion, s2, lambda0, 
                    working_wave, line_type='', outpdf='', line_name='', hbeta_bin=False):
    '''
    Purpose:
        This function fits each emission line with a Gaussian curve. It also calculates and saves in a table
        the Gaussian flux, observed flux, sigma, median, norm, RMS, and S/N for each provided emission line.
        
    Usage:
        emission_line_fit.zoom_gauss_plot(pdf_pages, N, wave, Spect_1D, dispersion, s2, lambda0, working_wave,
                                          line_type='', outpdf='', line_name='', hbeta_bin=False)
        
    Params:
        pdf_pages --> the pdf that the emission line fit plots are saved to.
        N --> the number of galaxies in each bin.
        Spect_1D --> an array of spectral data.
        lambda0 --> a list of wavelengths (in Angstroms) that all the emission lines are at.
        working_wave --> the wavelength (in Angstroms) that the emission line is at.
        line_type (OPTIONAL) --> a string of the emission line type: Oxy2, Balmer, or Single. Default is ''.
        outpdf (OPTIONAL) --> a string naming the pdf containing the emission line fits. Default is ''.
        line_name (OPTIONAL) --> a string of the emission line name. Default is ''.
        hbeta_bin (OPTIONAL) --> True if the binning used was mass-Hbeta luminosity binning. False if the 
            binning used was mass binning (default).
        
    Returns:
        tab0 --> the ascii table containing each emission line's Gaussian flux, observed flux, sigma, median,
            norm, RMS, and S/N values.
        
    Outputs:
        None
    '''
    
    if hbeta_bin:
        nrows = 4
        ncols = 4
    else:
        nrows = 2
        ncols = 4
    
    mask_flag = np.zeros(len(wave))
    for t_lam in lambda0:
        em_idx = np.where((wave >= (t_lam-5)) & (wave <= (t_lam+5)))[0]
        if len(em_idx) > 0: mask_flag[em_idx] = 1

    x_idx = np.where((wave >= (working_wave - 100)) & (wave <= (working_wave + 100)))[0] 
    x_idx_mask = np.where((wave >= (working_wave - 100)) & (wave <= (working_wave + 100)) &
                          (mask_flag == 0))[0]
    
    x0 = wave   
    
    
    fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey='row', squeeze=False) 
    fig.text(0.5, 0.98, line_name)
    
    ref_ymax = 0
    ref_ymed = -0.25
    ref_ymax_list = []
    ref_ymed_list = []
    count = 1

    for rr in range(Spect_1D.shape[0]):       
        y0 = Spect_1D[rr]
        y_norm = y0 / scalefact

        row = np.int(np.floor(rr / ncols))
        col = np.int(np.floor(rr % ncols))
        
        x1 = working_wave - 100
        x2 = working_wave + 100

        y_smooth = movingaverage_box1D(Spect_1D[rr] / scalefact, 2, boundary='extend')

        o1, med0, max0 = get_gaussian_fit(working_wave, x0, y0, y_norm, x_idx, x_idx_mask, line_type, s2)
        
    
        # Determines y limits 
        if med0 + max0 > ref_ymax:
            ref_ymax = med0 + max0
        if med0 < ref_ymed:
            ref_ymed = med0
            
        if (col == ncols - 1) or ((row == nrows - 1) and (rr == Spect_1D.shape[0] - 1)):
            ref_ymax_list.append(ref_ymax)
            ref_ymed_list.append(ref_ymed)
            ref_ymax = 0
            ref_ymed = -0.25
        
        
        # Calculating Flux: Signal Line Fit
        if type(o1) != type(None):
            dx = x0[2] - x0[1]
            if line_type == 'Single':
                x_sigsnip = np.where((np.abs((x0 - working_wave)) / o1[1]) <= 2.5)[0]
                gauss0 = gauss(x0,*o1)
            
            if line_type == 'Balmer':
                x_sigsnip = np.where(np.abs((x0 - working_wave)) / o1[1] <= 2.5)[0] 
                gauss0 = double_gauss(x0, *o1)
                o1_neg = [o1[0], o1[4], o1[5], o1[3]]
                neg0 = gauss(x0, *o1_neg)
                gauss0_diff = gauss0 - neg0
                y_norm_diff = y_norm[x_sigsnip] - neg0[x_sigsnip]

            if line_type == 'Oxy2':
                x_sigsnip = np.where(((x0 - working_wave) / o1[1] >= -2.5) & 
                                     ((x0 - working_wave * con1) / o1[4] <= 2.5))[0]
                gauss0 = oxy2_gauss(x0, *o1)
            
            if line_type == 'Single' or line_type == 'Oxy2':
                flux_g = np.sum((gauss0 - o1[3]) * dx)             # flux from gaussian distribution 
                flux_s = np.sum((y_norm[x_sigsnip] - o1[3]) * dx)  # flux from snipping method

            if line_type == 'Balmer':
                flux_g = np.sum(gauss0_diff * dx)
                flux_s = np.sum(y_norm_diff * dx)

            # Calculating RMS
            ini_sig1, RMS_pix = rms_func(wave, dispersion, working_wave, y0, o1[1], mask_flag)
            
            # Filling In Arrays
            line_dict[col_names[0]][rr] = flux_g                #Flux Gaussian
            line_dict[col_names[1]][rr] = flux_s                #Flux Observed
            line_dict[col_names[2]][rr] = (flux_s / ini_sig1)   #S/N
            line_dict[col_names[3]][rr] = o1[0]                 #Center
            line_dict[col_names[4]][rr] = max0                  #Norm
            line_dict[col_names[5]][rr] = o1[3]                 #Median
            line_dict[col_names[6]][rr] = o1[1]                 #Sigma
            line_dict[col_names[7]][rr] = ini_sig1              #RMS
            if line_type == 'Balmer':
                line_dict[col_names[8]][rr] = o1[5]             #Abs Norm
                line_dict[col_names[9]][rr] = o1[4]             #Abs Sigma

            resid = y_norm[x_sigsnip] - gauss0[x_sigsnip] + o1[3]  

            # Plotting
            t_ax = ax_arr[row, col]
            t_ax.plot(wave, y_norm, 'k', linewidth=0.6, label='Emission')
            t_ax.set_xlim([x1 + 45,x2 - 45])
            t_ax.plot(x0, gauss0, 'b--', linewidth=0.5, label='Gauss Fit')
            t_ax.plot(x0[x_sigsnip], resid, 'r', linestyle='dashed', linewidth=0.2, label='Residuals')
    
            
            if line_type == 'Balmer' or line_type == 'Oxy2':
                str1 = 'Sigma:%.2f' % (o1[1])
                str2 = 'Sigma:%.2f' % (o1[4])
                str3 = 'S/N:%.2f' % (np.round_((flux_s / ini_sig1), decimals=2))
                str4 = 'FluxO:%.2f' % (flux_s) 
                str5 = 'FluxG:%.2f' % (flux_g) 
                str6 = 'NGal:%.2f' % (N[rr])
                str7 = 'MBin:%.2f' % (rr + 1)
                if hbeta_bin:
                    t_ax.annotate('HbBin:%.2f' % (rr + 1), [0.45, 0.84], xycoords='axes fraction', 
                                  va='top', ha='right', fontsize='8')
                    str7 = 'MBin:%.2f' % (count)
                    if rr % 2 != 0:
                        count += 1
                t_ax.annotate(str1, [0.97, 0.98], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
                t_ax.annotate(str2, [0.97, 0.91], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
                t_ax.annotate(str3, [0.97, 0.84], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
                t_ax.annotate(str4, [0.97, 0.77], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
                t_ax.annotate(str5, [0.97, 0.70], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
                t_ax.annotate(str6, [0.45, 0.98], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
                t_ax.annotate(str7, [0.45, 0.91], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
            if line_type == 'Single':
                str1 = 'Sigma:%.2f' % (o1[1])
                str2 = 'S/N:%.2f' % (np.round_((flux_s / ini_sig1), decimals=2))
                str3 = 'FluxO:%.2f' % (flux_s)
                str4 = 'FluxG:%.2f' % (flux_g)
                str5 = 'NGal:%.2f' % (N[rr])
                str6 = 'MBin:%.2f' % (rr + 1)
                if hbeta_bin:
                    t_ax.annotate('HbBin:%.2f' % (rr + 1), [0.45, 0.84], xycoords='axes fraction', 
                                  va='top', ha='right', fontsize='8')
                    str6 = 'MBin:%.2f' % (count)
                    if rr % 2 != 0:
                        count += 1
                t_ax.annotate(str1, [0.97, 0.98], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
                t_ax.annotate(str2, [0.97, 0.91], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
                t_ax.annotate(str3, [0.97, 0.84], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
                t_ax.annotate(str4, [0.97, 0.77], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
                t_ax.annotate(str5, [0.45, 0.98], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
                t_ax.annotate(str6, [0.45, 0.91], xycoords='axes fraction', va='top', ha='right', 
                              fontsize='8')
            
            for x in lambda0:
                t_ax.axvline(x=x, linewidth=0.3, color='k')

            if rr == Spect_1D.shape[0] - 1 and rr % (nrows * ncols) != nrows * ncols - 1:
                for jj in range(col + 1, ncols):
                    ax_arr[row, jj].axis('on')
            for kk in range(row + 1, nrows):
                for zz in range(ncols):
                    ax_arr[kk, zz].axis('on')

                
        if (rr % (nrows * ncols) == nrows * ncols - 1) or rr == Spect_1D.shape[0] - 1: 
            subplots_adjust(left=0.1, right=0.98, bottom=0.06, top=0.97, hspace=0.05)
            row_count = 0
            for ii in range(0, nrows):
                for jj in range(0, ncols):
                    ax_arr[ii, jj].set_ylim([ref_ymed_list[row_count], ref_ymax_list[row_count] * 1.05])
                    if jj == ncols - 1:
                        row_count += 1
            
            fig.set_size_inches(8, 8)
            fig.savefig(pdf_pages, format='pdf')
     
    print('Done!')
    return line_dict

   
    

def zm_general(fitspath, Spect_1D, dispersion, wave, lambda0, line_type, line_name, s, a,
               c, s1, a1, s2, a2, hbeta_bin=False):
    '''
    Purpose:
        This function calls the emission line fitting and plotting functions, producing a pdf of the all
        the fits and an ascii table containing measured data about each line.
        
    Usage:
        emission_line_fit.zm_general(fitspath, Spect_1D, dispersion, wave, lambda0,
                                     line_type, line_name, s, a, c, s1, a1, s2, a2, hbeta_bin=False)
        
    Params:
        fitspath --> a string of the file path where the output files will be placed.
        lambda0 --> a list of wavelengths (in Angstroms) that all the emission lines are at.
        line_type --> a list of strings denoting the emission line type: Oxy2, Balmer, or Single.
        line_name --> a list of strings denoting the emission line name.
        hbeta_bin (OPTIONAL) --> True if the binning used was mass-Hbeta luminosity binning. False if the 
            binning used was mass binning (default).        
        
    Returns:
        None
        
    Outputs:
        table_stack --> an ascii table containing each emission line's Gaussian flux, observed flux, sigma,
            median, norm, RMS, and S/N values.
        pdf_pages --> a pdf of all the emission line fits.
    '''
    
    outpdf = fitspath + filename_dict['bin_fit'].replace('.tbl', '.pdf')
    pdf_pages = PdfPages(outpdf)
    table0 = asc.read(fitspath + filename_dict['bin_info'], format='fixed_width_two_line')
    out_ascii = fitspath + filename_dict['bin_fit']
    ID = table0[bin_names0[0]].data
    N = table0[bin_names0[1]].data
    
    em_dict = OrderedDict([(bin_names0[0], ID)])
    for ii in range(len(lambda0)):
        curr_line_cols = [line_col for line_col in gauss_lines_names0 if line_col.startswith(line_name[ii])]
        line_dict = create_empty_dict(curr_line_cols, arr_size=len(N), ordered=True)
        one_line_dict = zoom_gauss_plot(line_dict, curr_line_cols, pdf_pages, N, wave, Spect_1D, dispersion, s2, 
                                   lambda0, lambda0[ii], line_type=line_type[ii], line_name=line_name[ii], 
                                   hbeta_bin=hbeta_bin)
        em_dict.update(one_line_dict)
    asc.write(em_dict, out_ascii, format='fixed_width_two_line', overwrite=True)

    
    pdf_pages.close()
