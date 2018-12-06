import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from os.path import exists
from pylab import subplots_adjust
from scipy.optimize import curve_fit 
from astropy.convolution import Box1DKernel, convolve
from getpass import getuser


if getuser() == 'carol':
    fitspath = "C:\\Users\\carol\\Google Drive\\MZEvolve\\"
    fitspath2 = fitspath
else:
    fitspath = "../DEEP2/" 
    fitspath2 = "../"

flux = fitspath + 'flux_400_bin_4089.fits'                                      

'''stack2D, header = fits.getdata(flux, header = True)
wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])'''

xcoor = []

RestframeMaster = fitspath + 'Master_Grid.fits'
s = 1.0
a = 1.0
c = 2.0
s1 = 1.3
a1 = 1.5
s2 = 5
a2 = 1.8

lambda0 = [3726.18, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]
line_type = ['Oxy2', 'Balmer', 'Balmer', 'Single', 'Balmer', 'Single', 'Single']
line_name = ['OII_3727', 'HDELTA', 'HGAMMA', 'OIII_4363', 'HBETA', 'OIII_4958', 'OIII_5007'] 


def movingaverage_box1D(values, width, boundary = 'fill', fill_value = 0.0):
    box_kernel = Box1DKernel(width)
    smooth = convolve(values, box_kernel, boundary = boundary, fill_value = fill_value)
    return smooth


def gauss(x, xbar, s, a, c):
    return a * np.exp(-(x - xbar)**2 / s**2) + c


def double_gauss(x, xbar, s1, a1, c, s2, a2): 
    return a1 * np.exp(-(x - xbar)**2 / (2 * s1**2)) + c + a2 * np.exp(-(x - xbar)**2 / (2 * s2**2))


def oxy2_gauss(x, xbar, s1, a1, c, s2, a2):
    con1 = 3728.91/3726.16
    return a1 * np.exp(-(x - xbar)**2 / (2 * s1**2)) + c + a2 * np.exp(-(x - (xbar * con1))**2 / (2 * s2**2)) 


def get_gaussian_fit(working_wave, x0, y0, y_norm, x_idx, x_idx_mask, line_type):
    med0 = np.median(y_norm[x_idx_mask])
    max0 = np.max(y_norm[x_idx]) - med0

    fail = 0
    if line_type == 'Single': 
        p0 = [working_wave, 1.0, max0, med0] #must have some reasonable values
        para_bound = ((working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0)), (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0)))
        try:
            o1, o2 = curve_fit(gauss, x0[x_idx], y_norm[x_idx], p0 = p0, sigma = None, bounds = para_bound)
        except ValueError:
            print('fail')
            fail = 1
            
    if line_type == 'Balmer': 
        p0 = [working_wave, 1.0, max0, med0, s2, -0.05 * max0] #must have some reasonable values
        para_bound = (working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0), 0.0, -max0), (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0), 25.0, 0)
        try:
            o1, o2 = curve_fit(double_gauss, x0[x_idx], y_norm[x_idx], p0 = p0, sigma = None, bounds = para_bound)
        except ValueError:
            print('fail')
            fail = 1

    if line_type == 'Oxy2':
        p0 = [working_wave, 1.0, 0.75 * max0, med0, 1.0, max0] #must have some reasonable values
        para_bound = (working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0), 0.0, 0.0), (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0), 10.0, 100.0)
        try:
            o1, o2 = curve_fit(oxy2_gauss, x0[x_idx], y_norm[x_idx], p0 = p0, sigma = None, bounds = para_bound)
        except ValueError:
            print('fail')
            fail = 1
   
    if not fail:
        return o1, med0, max0
    else:
        return None, med0, max0


#calculating rms
def rms_func(wave, dispersion, lambda_in, y0, sigma_array, scalefact):
    x_idx = np.where((np.abs(wave - lambda_in) <= 100))[0]

    sigma = np.std(y0[x_idx])
    pix = 5 * sigma_array / dispersion 
    s_pix = np.sqrt(pix)

    ini_sig = s_pix * sigma * dispersion
    RMS_pix = sigma * dispersion / scalefact
   
    return ini_sig / scalefact, RMS_pix


#for each individual stack
#electron temperature and the R23 and O32 values 
#Plotting Zoomed 
def zoom_gauss_plot(dataset, working_wave, line_type = '', outpdf = '', line_name = ''):
    Spect_1D, header = fits.getdata(flux, header=True)
    dispersion = header['CDELT1']
    wave = header['CRVAL1'] + dispersion*np.arange(header['NAXIS1'])
    if outpdf == '':
        outpdf = fitspath + 'mass_bin_400.pdf'
    
    pdf_pages = PdfPages(outpdf)
    nrows = 2
    ncols = 5
    
    mask_flag = np.zeros(len(wave))
    for t_lam in lambda0:
        em_idx = np.where(wave >= (t_lam-5) & (wave <= (t_lam+5)))[0]
        if len(em_idx) > 0: mask_flag[em_idx] = 1

    x_idx = np.where((wave >= (working_wave - 100)) & (wave <= (working_wave + 100)))[0] 
    x_idx_mask = np.where((wave >= (working_wave - 100)) & (wave <= (working_wave + 100)) &
                          (mask_flag == 0))[0]
    x0 = wave   
    scalefact = 1e-17
    
    
    #mask_flag[x_idx_mask] = [1 for ii in range(len(wave))]
    
    #Initializing Arrays
    flux_g_array = np.zeros(Spect_1D.shape[0])
    flux_s_array = np.zeros(Spect_1D.shape[0])
    sigma_array = np.zeros(Spect_1D.shape[0])
    median_array = np.zeros(Spect_1D.shape[0])
    norm_array = np.zeros(Spect_1D.shape[0])
    RMS_array = np.zeros(Spect_1D.shape[0])
    SN_array = np.zeros(Spect_1D.shape[0])
    N_gal_array = np.zeros(Spect_1D.shape[0])
    R_23_array = np.zeros(Spect_1D.shape[0])
    O_32_array = np.zeros(Spect_1D.shape[0])

    '''#Initialize Error Propogation
    flux_g_err = np.zeros(Spect_1D.shape[0])
    flux_s_err = np.zeros(Spect_1D.shape[0])
    sigma_err = np.zeros(Spect_1D.shape[0])
    median_err = np.zeros(Spect_1D.shape[0])
    norm_err = np.zeros(Spect_1D.shape[0])
    RMS_err = np.zeros(Spect_1D.shape[0])
    SN_err = np.zeros(Spect_1D.shape[0])'''

    for rr in range(Spect_1D.shape[0]):
        y0 = Spect_1D[rr]
        y_norm = y0 / scalefact

        row = np.int(np.floor(rr / ncols))
        col = np.int(np.floor(rr % ncols))
        
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey='row', squeeze=False) 
       
        t_ax = ax_arr[row, col]
        
        x1 = working_wave - 100
        x2 = working_wave + 100

        y_smooth = movingaverage_box1D(Spect_1D[rr] / scalefact, 2, boundary = 'extend')

        o1, med0, max0 = get_gaussian_fit(working_wave, x0, y0, y_norm, x_idx, x_idx_mask, line_type)
     
        #Calculating Flux: Signal Line Fit
        if type(o1) != type(None):
            dx = x0[2] - x0[1]  #0.16866575
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
                con1 = 3728.91/3726.16
                x_sigsnip = np.where(((x0 - working_wave) / o1[1] >= -2.5) & ((x0 - working_wave * con1) / o1[4] <= 2.5))[0]
                gauss0 = oxy2_gauss(x0, *o1)
            
            if line_type == 'Single' or line_type == 'Oxy2':
                flux_g = np.sum((gauss0 - o1[3]) * dx)    #flux from gaussian distribution 
                flux_s = np.sum((y_norm[x_sigsnip] - o1[3]) * dx)  #flux from snipping method (spectral flux)where snip off sigma >2.5

            if line_type == 'Balmer':
                flux_g = np.sum(gauss0_diff * dx)
                flux_s = np.sum(y_norm_diff * dx)

            #Calculating RMS
            ini_sig1, RMS_pix = rms_func(wave, dispersion, working_wave, y0, o1[1], scalefact)
            
            #Filling In Arrays
            flux_g_array[rr] = flux_g 
            flux_s_array[rr] = flux_s
            sigma_array[rr] = o1[1]
            median_array[rr] = o1[3]
            norm_array[rr] = max0
            RMS_array[rr] = ini_sig1
            SN_array[rr] = (flux_s / ini_sig1)

            resid = y_norm[x_sigsnip] - gauss0[x_sigsnip] + o1[3]  
            
            '''#Error Propogation 
            flux_g_err[rr] = error_prop_chuncodes(flux_g,1)
            flux_s_err[rr] = error_prop_chuncodes(flux_s,1)
            sigma_err[rr] =  error_prop_chuncodes(o1[1],1)
            median_err[rr] = error_prop_chuncodes(o1[3],1)
            norm_err[rr] =   error_prop_chuncodes(max0,1)
            RMS_err[rr] =    error_prop_chuncodes(ini_sig1,1)
            SN_err[rr] =     error_prop_chuncodes(flux_s/ini_sig1,1)'''

        #Plotting
            t_ax.plot(wave, y_norm, 'k', linewidth = 0.6, label = 'Emission')
        
            t_ax.set_xlim([x1 + 45,x2 - 45])

            t_ax.plot(x0, gauss0, 'b--', linewidth = 0.5, label = 'Gauss Fit')
            t_ax.plot(x0[x_sigsnip], resid, 'r', linestyle = 'dashed', linewidth = 0.2, label = 'Residuals')
            t_ax.legend(bbox_to_anchor = (0.25,0.1), borderaxespad = 0, ncol = 2, fontsize = 3)
            #t_ax.annotate(np.nanmax(gauss0[x_idx]), [0.95,0.95], xycoords = 'axes fraction', va = 'top', ha = 'right', fontsize = '5')
            if line_type == 'Balmer' or line_type == 'Oxy2': 
                t_ax.annotate(np.round_(o1[0], decimals=2), [0.95,0.95], xycoords = 'axes fraction', va = 'top', ha = 'right', fontsize = '5')
                t_ax.annotate(np.round_(o1[1], decimals=2), [0.95,0.90], xycoords = 'axes fraction', va = 'top', ha = 'right', fontsize = '5')
                t_ax.annotate(np.round_(o1[2], decimals=2), [0.95,0.85], xycoords = 'axes fraction', va = 'top', ha = 'right', fontsize = '5')
                t_ax.annotate(np.round_(o1[3], decimals=2), [0.95,0.80], xycoords = 'axes fraction', va = 'top', ha = 'right', fontsize = '5')
                t_ax.annotate(np.round_(o1[4], decimals=2), [0.95,0.75], xycoords = 'axes fraction', va = 'top', ha = 'right', fontsize = '5')
                t_ax.annotate(np.round_(o1[5], decimals=2), [0.95,0.70], xycoords = 'axes fraction', va = 'top', ha = 'right', fontsize = '5')
            if line_type == 'Single':
                t_ax.annotate(np.round_(o1[0], decimals=2), [0.95,0.95], xycoords = 'axes fraction', va = 'top', ha = 'right', fontsize = '5')
                t_ax.annotate(np.round_(o1[1], decimals=2), [0.95,0.90], xycoords = 'axes fraction', va = 'top', ha = 'right', fontsize = '5')
                t_ax.annotate(np.round_(o1[2], decimals=2), [0.95,0.85], xycoords = 'axes fraction', va = 'top', ha = 'right', fontsize = '5')
                t_ax.annotate(np.round_(o1[3], decimals=2), [0.95,0.80], xycoords = 'axes fraction', va = 'top', ha = 'right', fontsize = '5')
            
            for x in lambda0:
                t_ax.axvline(x = x, linewidth = 0.3, color = 'k')

            '''if col != 0:
                #t_ax.set_ylabel('Spect_1D')
                t_ax.set_yticklabels([])  #sets y-tick labels '''

            if rr == Spect_1D.shape[0] - 1 and rr % (nrows * ncols) != nrows * ncols - 1:
                for jj in range(col + 1, ncols):
                    ax_arr[row, jj].axis('on')
            for kk in range(row + 1, nrows):
                for zz in range(ncols):
                    ax_arr[kk, zz].axis('on')
                
        if (rr % (nrows * ncols) == nrows * ncols - 1) or rr == Spect_1D.shape[0] - 1: 
            subplots_adjust(left = 0.1, right = 0.98, bottom = 0.06, top = 0.97, hspace = 0.05)
            
            fig.set_size_inches(8, 8)
            fig.savefig(pdf_pages, format ='pdf')
        
    #endfor
    '''
    #Error Propogation 
    flux_g_err = error_prop_chuncodes(flux_g_array,1)
    flux_s_err = error_prop_chuncodes(flux_s_array,1)
    sigma_err =  error_prop_chuncodes(sigma_array,1)
    median_err = error_prop_chuncodes(median_array,1)
    norm_err =   error_prop_chuncodes(norm_array,1)
    RMS_err =    error_prop_chuncodes(RMS_array,1)
    SN_err =     error_prop_chuncodes(SN_array,1)'''
     
    #Writing Ascii Tables and Fits Tables
    out_ascii = fitspath + dataset + '_gaussian_' + str(np.int(working_wave)) + '.tbl'  
    out_fits = fitspath + dataset + '_gaussian_' + line_name + '.fits'
    n = ('Flux_Gaussian', 'Flux_Observed', 'Sigma', 'Median', 'Norm', 'RMS', 'S/N')
    n = tuple([line_name + '_' + val for val in n])
    tab0 = Table([flux_g_array, flux_s_array, sigma_array, median_array, norm_array, RMS_array, SN_array], names = n)
    asc.write(tab0, out_ascii, format = 'fixed_width_two_line')
         
    out_ascii_single = fitspath + dataset + '_Average_R23_O32_Values.tbl'
    if not exists(out_ascii_single):
        n2 = ('R_23_Average', 'O_32_Average', 'N_Galaxies', 'RMS')
        tab1 = Table([R_23_array, O_32_array, N_gal_array, RMS_array], names = n2)
        asc.write(tab1, out_ascii_single, format = 'fixed_width_two_line') 
       
        '''
        print 'Starting fits'
        #header['Comment'] = 'This fits file contains the calculated and observed fluxs for emission line stated in file name'
        c1 = fits.Column(name='Flux_Gaussian'+line_name, format= 'K', array=flux_g_array)
        c2 = fits.Column(name='Flux_Observed'+line_name, format= 'K', array=flux_s_array)
        c3 = fits.Column(name='Sigma'+line_name,  format= 'K', array=sigma_array)
        c4 = fits.Column(name='Median'+line_name,  format= 'K', array=median_array)
        c5 = fits.Column(name='Norm'+line_name, format= 'K',  array=norm_array)
        c6 = fits.Column(name='RMS'+line_name, format= 'K',  array=ini_sig1)

        #fits_cols = fits.ColDefs([c1,c2,c3,c4,c5])
        f_tab= fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6])
        
        f_tab.writeto(fitspath+'/Mar19_outputs/Flux_Outputs'+line_name+'.fits', overwrite= True)'''
     
    pdf_pages.close()
    print('Done!')



def zm_general(dataset):
    for ii in range(len(lambda0)):
        if line_type[ii] == 'Single':
            outpdf = fitspath + 'mass_bin_400_' + line_name[ii] + '.pdf'   
            print(outpdf)
            zoom_gauss_plot(dataset, lambda0[ii], line_type = line_type[ii], outpdf = outpdf, line_name = line_name[ii])
            
        if line_type[ii] == 'Balmer': 
            outpdf = fitspath + 'mass_bin_400_' + line_name[ii] + '.pdf'
            print(outpdf)
            zoom_gauss_plot(dataset, lambda0[ii], line_type = line_type[ii], outpdf = outpdf, line_name = line_name[ii])
            
        if line_type[ii] == 'Oxy2': 
            outpdf = fitspath + 'mass_bin_400_' + line_name[ii] + '.pdf'
            print(outpdf)
            zoom_gauss_plot(dataset, lambda0[ii], line_type = line_type[ii], outpdf = outpdf, line_name = line_name[ii])



 