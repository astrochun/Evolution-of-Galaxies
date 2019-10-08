from astropy.io import fits, ascii
from collections import Counter
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from getpass import getuser
from os.path import exists
import os



if getuser() == 'carol':
    path = "C:\\Users\\carol\\Google Drive\\MZEvolve\\"
    path2 = path + "massbin\\"
else:
    path = "../DEEP2/" 
    path2 = "../"
    

def exclude_outliers(objno):
    flag = np.zeros(len(objno), dtype=int)
    bad_data = np.array(['32007727', '32101412', '42006031a', '32035286', '14023705'])
    for ii in range(len(bad_data)):
        idx = [xx for xx in range(len(objno)) if bad_data[ii] == objno[xx]]
        flag[idx] = 1
    
    return flag


def stack_spectra(fname, mname='', plot = False, indices = 'placeholder'):
    hdu = fits.open(fname)
    image = hdu[0].data
    header = hdu[0].header
    wavelen = header["CRVAL1"] + np.array(header["CDELT1"])*range(header['NAXIS1'])
    
    if mname != '':
        hdu = fits.open(mname)
        mask_image = hdu[0].data
        mask_header = hdu[0].header
        
        image = np.ma.array(image, mask=mask_image)
    
    if (indices == 'placeholder'):
        flux = np.mean(image, axis = 0)
    else:
        flux = np.mean(image[indices], axis = 0)
    if plot == True:
        plt.plot(wavelen, flux)
        plt.xlabel("Wavelength")
        plt.ylabel("Flux")
        plt.title("Stacked Spectra of %s" % fname)
    
    hdu.close()
    
    return image, flux, wavelen


def interpolate_data(interp_file):
    npz_file = np.load(interp_file)
    interp_data, mag_no_mass, no_mass_idx = npz_file['interp_data'][0], npz_file['mag_no_mass'], npz_file['no_mass_idx']
    interp_mass = interp_data(mag_no_mass)
    
    plt.scatter(mag_no_mass, interp_mass)
    
    return interp_mass, no_mass_idx
    


def binning(temp_x, objno, bin_pts_input, interp_file, bin_pts_fname, mname = '', fitspath0 = '',
            spectra_plot = False, filename = False, adaptive = False, hbeta_bin = False, lum = []):
    """
    temp_x = quantity to be divided into bins --> unsorted
    bin_pts_input = number of points in each bin
    x_sort = masses sorted by mass value
    ind_sort = indices corresponding to sorted mass values in same order as x_sort
    logx_sort = log of masses sorted by mass value
    bin_start = lowest mass value for each bin 
    """
    
    interp_mass, no_mass_idx = interpolate_data(interp_file)    
    interp_mass = 10**(interp_mass)
    temp_x[no_mass_idx] = interp_mass
    
    temp_flag = exclude_outliers(objno)
    
    #Remove bad data and initially excluded data
    valid_ind = np.where((np.isfinite(temp_x)==True) & (temp_x>0) & (temp_x < 1e13) & (temp_flag == 0))[0]
    
    x_sort = np.sort(temp_x[valid_ind])
    ind_sort = np.argsort(temp_x[valid_ind])
    objno = objno[valid_ind][ind_sort]
   
    logx = np.log10(temp_x)
    logx_sort = np.log10(x_sort)
    y = range(len(logx_sort))
    
            
    out_file = fitspath0 + bin_pts_fname + '.npz'
        
    if exists(out_file):
        print('File exists')
        rinput = input('Do you want to delete file? Yes or no ')
        if rinput.lower() == 'yes':
            os.remove(out_file)
            print(fitspath0 + bin_pts_fname + '.npz deleted.')
        else:
            idx_file = np.load(out_file) 
            bin_ind = idx_file['bin_ind'] 
            bin_start = idx_file['bin_start']
            bin_edge = idx_file['bin_edge']
            bin_redge = idx_file['bin_redge']
            flux = idx_file['flux']
            wavelength = idx_file['wavelength']
            mass_avg = idx_file['mass_avg']
            bin_ID = idx_file['bin_ID']
            count = len(bin_ID)
            N = idx_file['N']
            idx_file.close()
            
    if not exists(out_file):
        bin_ind = []
        start = 0
        bin_start = logx_sort[start]
        bin_edge = []
        bin_redge = []
        flux = []
        wavelength = []
        mass_avg = []
        bin_ID = []
        N = []
        lowest_hbeta = []
        highest_hbeta = []
        count = 0
        count2 = 0
            
        while (bin_start < logx_sort[-1]):
            if adaptive == True:
                bin_pts = bin_pts_input[count]
            else:
                bin_pts = bin_pts_input
                    
            stop = start + bin_pts
            if ((stop + bin_pts) > len(logx_sort)):
                stop = len(logx_sort) - 1
            count += 1
            bin_stop = logx_sort[stop]
            bin_edge.append(bin_start)
                
            if hbeta_bin == False:
                bin_ind.append(ind_sort[start:stop])
                if filename != False:
                    _, flx, wave = stack_spectra(filename, mname, indices = ind_sort[start:stop])
                    flux.append(flx)
                    wavelength.append(wave)
                N.append(len(ind_sort[start:stop]))
                mass_avg.append(np.mean(logx_sort[start:stop]))
                bin_ID.append(count)
            else:
                lum_sort = lum[valid_ind][ind_sort]
                valid_hbeta = np.where(lum_sort[start:stop] < 44)[0]
                median0 = np.median(lum_sort[start:stop][valid_hbeta])
                invalid_hbeta = np.where((lum_sort[start:stop] > 44) | (np.isfinite(lum_sort[start:stop]) == False))[0]
                lum[valid_ind[ind_sort[start:stop][invalid_hbeta]]] = -1
                temp_lower = np.where(lum <= median0)[0]
                temp_upper = np.where(lum > median0)[0]
                lower_idx = np.array(list(set(valid_ind[ind_sort][start:stop]) & set(temp_lower)))
                upper_idx = np.array(list(set(valid_ind[ind_sort][start:stop]) & set(temp_upper)))
                    
                non_neg = np.where(lum[lower_idx] != -1)[0]
                non_neg = lower_idx[non_neg]
                lowest_hbeta.append(np.min(lum[non_neg]))
                highest_hbeta.append(np.max(lum[non_neg]))
                lowest_hbeta.append(np.min(lum[upper_idx]))
                highest_hbeta.append(np.max(lum[upper_idx]))
                    
                bin_redge.append(logx_sort[start + len(lower_idx) - 1])
                bin_edge.append(logx_sort[start + len(lower_idx)])
                    
                bin_ind.append(lower_idx)
                bin_ind.append(upper_idx)
                if filename != False:
                    _, lower_flx, lower_wave = stack_spectra(filename, mname, indices = lower_idx)
                    _, upper_flx, upper_wave = stack_spectra(filename, mname, indices = upper_idx)
                    flux += [lower_flx] + [upper_flx]
                    wavelength += [lower_wave] + [upper_wave]
                    N += [len(lower_idx), len(upper_idx)]
                mass_avg += [np.mean(logx[lower_idx])] + [np.mean(logx[upper_idx])]
                bin_ID.append(count2)
                count2 += 1
                bin_ID.append(count2)
                count2 += 1
                    
            start, bin_start = stop, bin_stop
            bin_redge.append(bin_stop)
            
        '''
        The mass array saved here contains the 4140 masses (not the log(mass)) with the no mass indices replaced
        with interpolated values; it does NOT exclude sources excluded for binning (indices relative to
        original table). 
        The luminosity array saved here contains the 4140 luminosities (log(lum)) with the >44 and nans replaced
        to -1.
        '''
        np.savez(out_file, mass = temp_x, lum = lum, bin_ind = bin_ind, bin_start = bin_start, bin_edge = bin_edge, 
                 bin_redge = bin_redge, flux = flux, wavelength = wavelength, mass_avg = mass_avg,
                 bin_ID = bin_ID, N = N, lowest_hbeta = lowest_hbeta, highest_hbeta = highest_hbeta)
        
        if adaptive == False:
            out_ascii = fitspath0 + str(bin_pts_input) + '_binning.tbl' 
        else:
            out_ascii = fitspath0 + bin_pts_fname + '_binning.tbl'
        if hbeta_bin == False:
            lowest_hbeta = np.zeros(len(bin_pts_input))
            highest_hbeta = np.zeros(len(bin_pts_input))
        n = ('ID', 'mass_min', 'mass_max', 'mass_avg', 'Number of Galaxies', 'Lowest Hbeta', 'Highest Hbeta')
        table_stack = Table([bin_ID, bin_edge, bin_redge, mass_avg, N, lowest_hbeta, highest_hbeta], names = n)
        ascii.write(table_stack, out_ascii, format = 'fixed_width_two_line', overwrite = True)
        
    
    xlim = [4250,4450]
    if (spectra_plot == True):
        for i in range(count):
            plt.subplot(np.ceil(count/2.0), 2, i+1)
            if hbeta_bin == False:
                plt.plot(wavelength[i], flux[i], color = 'b', linestyle = 'solid')
                wavelen_idx = np.where((wavelength[i] > xlim[0]) & (wavelength[i] < xlim[1]))[0]
                max0 = np.max(flux[i][wavelen_idx])   
            else:
                plt.plot(wavelength[i*2], flux[i*2], color = 'b', linestyle = 'solid')
                wavelen_idx = np.where((wavelength[i*2] > xlim[0]) & (wavelength[i*2] < xlim[1]))[0]
                max0 = np.max(flux[i*2][wavelen_idx])
                plt.plot(wavelength[i*2+1], flux[i*2+1], color = 'orange', linestyle = 'solid')
                max0 = np.max(flux[i*2+1][wavelen_idx])
            plt.ylim(-0.05e-17, max0*1.1)
            plt.xlim(xlim)
            plt.axvline(x=5007, color='k', linestyle = 'dashed', alpha=0.5)
            plt.axvline(x=4363, color='r', linestyle = 'dashed', alpha=0.5)
            plt.suptitle(str(i))
            if hbeta_bin == False:
                an_text = 'N = '+str(N[i])
            else:
                an_text = 'N = [%i,%i]' %(N[i*2], N[i*2+1]) 
            plt.annotate(an_text, [0.05,0.95], xycoords='axes fraction', ha='left', va='top')
    else:
        plt.bar(bin_edge, N, align = 'edge', width = bin_redge)
        
    
                
    
    return bin_edge, flux




def gen_hist(fname, hist, nbins = 10): 

#hist options: mass, chi2

    result = Table(ascii.read(fname))

    if hist == 'mass':
        mass = result['best.stellar.m_star']
        ind = np.where((np.isfinite(mass) == True) & (mass>0))[0]
        plt.hist(np.log10(mass[ind]), bins = nbins)
        plt.title('Distribution of stellar mass')
        
    if hist == 'chi2':
        chi2 = result['best.reduced_chi_square']
        ind = np.where(np.isfinite(chi2) == True)[0]
        plt.hist(np.log10(chi2[ind]))




"""
def duplicates(column):
    col = list(column)
    
    counts = Counter(col)
    ascii.write(counts.items(), 'count.csv')
    for s,num in counts.items():
        #print s, num
        if num > 1:
            for suffix in list(string.ascii_lowercase)[0:num]:
                ind = col.index(s)
                col[ind] = str(s) + suffix
                #pdb.set_trace()
                print s, 'converted to', col[ind] 
    return np.array(col)
"""


def duplicates(column):
    col = [str(a) for a in column.data]
    counts  = Counter(col)
    N       = np.array(list(counts.values()))
    idx_dup = np.where(N > 1)[0]
    keys    = np.array(list(counts.keys()))
    for ii in range(len(idx_dup)):
        s = keys[idx_dup[ii]]
        print(ii, s)
        t_dup = [cc for cc in range(len(col)) if s in col[cc]]
        suffix0 = list(string.ascii_lowercase[0:len(t_dup)])
        for ind,suffix in zip(t_dup,suffix0):
            col[ind] = s + suffix
            print(ii, s, 'converted to', col[ind])
    return np.array(col)