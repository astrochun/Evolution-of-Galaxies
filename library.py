from astropy.io import fits, ascii
from astropy.table import Table, hstack
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import glob, pdb, string
from getpass import getuser
from os.path import exists



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
    #temp1 = fname.split('/')
    #temp2 = temp1[1].split('.')
    #field = temp2[0]
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
    


def binning(temp_x, objno, bin_pts_input, interp_file, bin_pts_fname, mname = '', bin_array_file = '',
            spectra_plot = False, filename = False, adaptive = False, hbeta_bin = False, lum = []):
    """
    temp_x = quantity to be divided into bins [must NOT be sorted]
    bin_pts_input = Number of points in each bin
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
            
    lum_sort = lum[valid_ind][ind_sort]
    
    if exists(bin_array_file + bin_pts_fname + '.npz'):
        idx_file = np.load(bin_array_file + bin_pts_fname + '.npz') 
        bin_ind = idx_file['bin_ind'] 
        bin_start = idx_file['bin_start']
        bin_edge = idx_file['bin_edge']
        bin_redge = idx_file['bin_redge']
        distribution = idx_file['distribution']
        flux = idx_file['flux']
        wavelength = idx_file['wavelength']
        mass_avg = idx_file['mass_avg']
        bin_ID = idx_file['bin_ID']
        count = len(bin_ID)
        N = idx_file['N']
        
    else:
        bin_ind = []
        start = 0
        bin_start = logx_sort[start]
        bin_edge = []
        bin_redge = []
        distribution = []
        flux = []
        wavelength = []
        mass_avg = []
        bin_ID = []
        N = []
        count = 0
        
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
            dist = len(y[start:stop])
            distribution.append(dist)
            bin_edge.append(bin_start)
            
            bin_ind.append([])
            if hbeta_bin == False:
                bin_ind[-1].append(ind_sort[start:stop])
                if filename != False:
                    _, flx, wave = stack_spectra(filename, mname, indices = ind_sort[start:stop])
                    flux.append([flx])
                    wavelength.append([wave])
                N.append(len(ind_sort[start:stop]))
                mass_avg.append(np.mean(logx_sort[start:stop]))
            else:
                valid_hbeta = np.where(lum_sort[start:stop] < 44)[0]
                median0 = np.median(lum_sort[start:stop][valid_hbeta])
                invalid_hbeta = np.where((lum_sort[start:stop] > 44) | (np.isfinite(lum_sort[start:stop]) == False))[0]
                lum[valid_ind[ind_sort[start:stop][invalid_hbeta]]] = -1
                temp_lower = np.where(lum <= median0)[0]
                temp_upper = np.where(lum > median0)[0]
                lower_idx = list(set(valid_ind[ind_sort][start:stop]) & set(temp_lower))
                upper_idx = list(set(valid_ind[ind_sort][start:stop]) & set(temp_upper))
                bin_ind[-1].append(lower_idx)
                bin_ind[-1].append(upper_idx)
                if filename != False:
                    _, lower_flx, lower_wave = stack_spectra(filename, mname, indices = lower_idx)
                    _, upper_flx, upper_wave = stack_spectra(filename, mname, indices = upper_idx)
                    flux.append([lower_flx, upper_flx])
                    wavelength.append([lower_wave, upper_wave])
                N.append([len(lower_idx), len(upper_idx)])
                mass_avg.append([np.mean(logx[lower_idx]), np.mean(logx[upper_idx])])
            bin_ID.append(count)
            
            start, bin_start = stop, bin_stop
            bin_redge.append(bin_stop)
    
    
    np.savez(bin_array_file + bin_pts_fname + '.npz', bin_ind = bin_ind, bin_start = bin_start, bin_edge = bin_edge,
             bin_redge = bin_redge, distribution = distribution, flux = flux, wavelength = wavelength,
             mass_avg = mass_avg, bin_ID = bin_ID, N = N)
    
    if adaptive == False:
        out_ascii = path2 + str(bin_pts_input) + '_massbin.tbl' 
    else:
        out_ascii = path2 + bin_pts_fname + '_massbin.tbl'
    n = ('ID', 'mass_min', 'mass_max', 'mass_avg', 'Number of Galaxies')
    table_stack = Table([bin_ID, bin_edge, bin_redge, mass_avg, N], names = n)
    ascii.write(table_stack, out_ascii, format = 'fixed_width_two_line', overwrite = True)
        
    if (spectra_plot == True):
        for i in range(count):
            plt.subplot(np.ceil(count/2.0), 2, i+1)
            #plt.plot(wavelength[i], flux[i])
            plt.plot(wavelength[i][0], flux[i][0], color = 'b', linestyle = 'solid')
            if hbeta_bin == True:
                plt.plot(wavelength[i][1], flux[i][1], color = 'orange', linestyle = 'solid')
            plt.ylim(-0.05e-17, 0.5e-17)
            plt.xlim(4250,4450)
            plt.axvline(x=5007, color='k', linestyle = 'dashed', alpha=0.5)
            plt.axvline(x=4363, color='r', linestyle = 'dashed', alpha=0.5)
            plt.suptitle(str(i))
            plt.annotate('N = '+str(N[i]), [0.05,0.95], xycoords='axes fraction',
                         ha='left', va='top')
    else:
        plt.bar(bin_edge, distribution, align = 'edge', width = bin_redge)
                
    
    return distribution, bin_edge, flux




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