from astropy.io import fits, ascii
from astropy.table import Table, hstack
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import glob, pdb, string
from getpass import getuser



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
            spectra_plot = False, filename = False, adaptive = False):
    """
    temp_x = quantity to be divided into bins [must NOT be sorted]
    bin_pts_input = Number of points in each bin
    """
    
    interp_mass, no_mass_idx = interpolate_data(interp_file)    
    interp_mass = 10**(interp_mass)
    temp_x[no_mass_idx] = interp_mass
    
    x_sort = np.sort(temp_x)
    ind_sort = np.argsort(temp_x)
    objno = objno[ind_sort]

    temp_flag = exclude_outliers(objno)
    

    #Remove bad data and initially excluded data
    valid_ind = np.where((np.isfinite(x_sort)==True) & (x_sort>0) & (x_sort < 1e13) & (temp_flag == 0))[0]


    '''excluded = [ii for ii in range(len(x_sort)) if ii not in valid_ind]
    print("total =", len(excluded))
    #con1 = is not finite, con2 = less than (a) or equal to (b) zero, con3 = greater than (a) or equal to (b) 1e13
    con1 = np.where(np.isfinite(x_sort)==False)[0]
    con2a = np.where(x_sort < 0)[0]
    con2b = np.where(x_sort == 0)[0]
    con3a = np.where(x_sort > 1e13)[0]
    con3b = np.where(x_sort == 1e13)[0]
    print("condition 1 =", len(con1))
    print("condition 2a =", len(con2a))
    print("condition 2b =", len(con2b))
    print("condition 3a =", len(con3a))
    print("condition 3b =", len(con3b))''' 
     
    x = np.log10(x_sort[valid_ind])
    ind = ind_sort[valid_ind]
    y = range(len(x))

    
    start = 0
    bin_start = x[start]
    bin_edge = []
    bin_redge = []
    distribution = []
    flux = []
    wavelength = []
    mass_avg = []
    bin_ID = []
    N = []
    count = 0

            
    while (bin_start < x[-1]):
        if adaptive == True:
            bin_pts = bin_pts_input[count]
        else:
            bin_pts = bin_pts_input
            
        stop = start + bin_pts
        if ((stop + bin_pts) > len(x)):
            stop = len(x) - 1
        count += 1
        bin_stop = x[stop]
        dist = len(y[start:stop])
        distribution.append(dist)
        bin_edge.append(bin_start)
        if filename != False:
            _, flx, wave = stack_spectra(filename, mname, indices = ind[start:stop])
            flux.append(flx)
            wavelength.append(wave)
        N.append(len(ind[start:stop]))
        mass_avg.append(np.mean(x[start:stop]))
        bin_ID.append(count)
        if bin_array_file != '':
            if adaptive == False:
                np.savez(bin_array_file + '_' + str(count) + '.npz', indices=ind[start:stop])
            else:
                np.savez(bin_array_file + bin_pts_fname + '_' + str(count) + '.npz', indices=ind[start:stop])
        
        start, bin_start = stop, bin_stop
        bin_redge.append(bin_stop)
    
    if adaptive == False:
        out_ascii = path2 + str(bin_pts) + '_massbin.tbl' 
    else:
        out_ascii = path2 + bin_pts_fname + '_massbin.tbl'
    n = ('ID', 'mass_min', 'mass_max', 'mass_avg', 'Number of Galaxies')
    table_stack = Table([bin_ID, bin_edge, bin_redge, mass_avg, N], names = n)
    ascii.write(table_stack, out_ascii, format = 'fixed_width_two_line', overwrite = True)
        
    if (spectra_plot == True):
        for i in range(count):
            plt.subplot(np.ceil(count/2.0), 2, i+1)
            plt.plot(wavelength[i], flux[i])
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