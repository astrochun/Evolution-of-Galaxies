from astropy.io import fits, ascii
from astropy.table import Table
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import glob, pdb, string





def stack_spectra(fname, plot = False, indices = 'placeholder'):
    hdu = fits.open(fname)
    temp1 = fname.split('/')
    temp2 = temp1[1].split('.')
    field = temp2[0]
    image = hdu[0].data
    header = hdu[0].header
    wavelen = header["CRVAL1"] + np.array(header["CDELT1"])*range(header['NAXIS1'])
    if (indices == 'placeholder'):
        flux = np.mean(image, axis = 0)
    else:
        flux = np.mean(image[indices], axis = 0)
    if plot == True:
        plt.plot(wavelen, flux)
        plt.xlabel("Wavelength")
        plt.ylabel("Flux")
        plt.title("Stacked Spectra of %s" % field)
    
    hdu.close()
    
    return image, flux, wavelen








def binning(temp_x, bin_pts, spectra_plot = False, filename = False):
    """
    temp_x = quantity to be divided into bins [must NOT be sorted]
    bin_pts = Number of points in each bin
    """
    
    temp_x1 = np.sort(temp_x)
    temp_ind = np.argsort(temp_x)
    temp_ind1 = np.where((np.isfinite(temp_x1)==True) & (temp_x1>0) & (temp_x1 < 1e13))[0]
    x = np.log10(temp_x1[temp_ind1])
    ind = temp_ind[temp_ind1]
    y = range(len(x))
    
    start = 0
    bin_start = x[start]
    bin_edge = []
    bin_size = []
    distribution = []
    flux = []
    wavelength = []
    count = 0
    while (bin_start < x[-1]):
        stop = start + bin_pts
        if (stop > len(x)):
            stop = len(x) - 1
        count += 1
        bin_stop = x[stop]
        dist = len(y[start:stop])
        distribution.append(dist)
        bin_edge.append(bin_start)
        bin_size.append((bin_stop - bin_start))
        if filename != False:
            _, flx, wave = stack_spectra(filename, indices = ind[start:stop])
            flux.append(flx)
            wavelength.append(wave)
        start, bin_start = stop, bin_stop
    
    if (spectra_plot == True):
        for i in range(count):
            plt.subplot(count/2, 2, i+1)
            plt.plot(wavelength[i], flux[i])
            plt.ylim(-0.2e-17, 2.0e-17)
            plt.axvline(x=5007, color='k', linestyle = 'dashed')
            plt.axvline(x=4363, color='r', linestyle = 'dashed')
    else:
        plt.bar(bin_edge, distribution, align = 'edge', width = bin_size)
        
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
    N       = np.array(counts.values())
    idx_dup = np.where(N > 1)[0]
    keys    = np.array(counts.keys())
    for ii in xrange(len(idx_dup)):
        s = keys[idx_dup[ii]]
        print(ii, s)
        t_dup = [cc for cc in xrange(len(col)) if s in col[cc]]
        suffix0 = list(string.ascii_lowercase[0:len(t_dup)])
        for ind,suffix in zip(t_dup,suffix0):
            col[ind] = s + suffix
            print(ii, s, 'converted to', col[ind])
    return np.array(col)