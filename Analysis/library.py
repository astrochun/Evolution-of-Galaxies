from astropy.io import fits, ascii
from collections import Counter
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from os.path import exists
import os

from Metallicity_Stack_Commons import exclude_outliers
from Metallicity_Stack_Commons.column_names import filename_dict, bin_names0, bin_mzevolve_names0
from Metallicity_Stack_Commons.column_names import merge_column_names, remove_from_list

    

def stack_spectra(fname, mname='', plot=False, indices='placeholder'):
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
        flux = np.mean(image, axis=0)
    else:
        flux = np.mean(image[indices], axis=0)
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
    


def binning(temp_x, objno, bin_pts_input, interp_file, mname='', fitspath0='',
            spectra_plot=False, filename=False, adaptive=False, hbeta_bin=False, lum=[]):
    '''
    Purpose:
        This function excludes sources with invalid mass, sorts the sources based on mass, and stacks
        spectra according to the provided number of sources in each bin. This function can bin the data
        according to mass or subdivide the mass bins according to HBeta luminosity values.
        
    Params:
        *See section below for any parameters not listed here*
        interp_file --> an npz file containing the interpolated mass results for sources with no mass.
        mname (OPTIONAL) --> the MastermaskArray.fits file.
        fitspath0 (OPTIONAL) --> a string of the file path where the output files will be placed.
        spectra_plot (OPTIONAL) --> True if a plot of stacked spectra is wanted. False otherwise (default).
        filename (OPTIONAL) --> the Master_Grid.fits file.
        adaptive (OPTIONAL) --> True if bin sizes vary. False if constant bin size (default).
        hbeta_bin (OPTIONAL) --> True if HBeta luminosity sub-binning is wanted. False otherwise (default). 
        lum (OPTIONAL) --> if HBeta luminosity binning is wanted, lum is a list of luminosity values for
            all sources. The default is an empty list.
        
    Returns:
        bin_edge --> left edge of each bin.
        flux --> array of flux values for each bin.
        
    Outputs:
        out_file --> an npz file containing bin information such as the indices of sources relative to 
            temp_x that are in each bin, the first and last sources in each bin, average values, etc.
        table_stack --> an ascii table with bin IDs and minimum, maximum, and average mass/HBeta luminosity
            values for each bin.
        composite_spectra --> a plot of the stacked spectra for each bin.
    '''
    
    """    
    Parameters/Variables:
        
    temp_x = array of stellar masses to be divided into bins --> unsorted
    objno = array of ID numbers --> unsorted
    bin_pts_input = number of points in each bin
    lum = array of HBeta luminosities
    
    interp_mass = array of interpolated mass values
    no_mass_idx = indices in temp_x that have negative or nan masses
    temp_flag = array of 0's and 1's that is the size of objno where 0 = valid object and 1 = invalid object
    
    valid_ind = indices relative to temp_x of valid masses
    x_sort = only valid masses sorted by mass value
    argsort_valid_x = indices relative to array containing only the valid masses (< 4140) 
    valid_ind_sort = indices in temp_x corresponding to sorted mass values put in same order as x_sort
    logx = log of entire mass array
    logx_sort = log of only valid masses sorted by mass value
    """
    
    # Replace indices of no mass with interpolated mass
    interp_mass, no_mass_idx = interpolate_data(interp_file)    
    interp_mass = 10**(interp_mass)
    temp_x[no_mass_idx] = interp_mass
    
    temp_flag = exclude_outliers(objno)    
    
    # Remove invalid data and initially excluded data
    valid_ind = np.where((np.isfinite(temp_x)==True) & (temp_x > 0) & (temp_x < 1e13) & (temp_flag == 0))[0]
    
    # sort valid masses by mass and ID array by sorted mass
    x_sort = np.sort(temp_x[valid_ind])
    argsort_valid_x = np.argsort(temp_x[valid_ind])
    valid_ind_sort = valid_ind[argsort_valid_x]  # values in valid_ind_sort are indices relative to full 4140
    #objno = objno[valid_ind_sort]               # IDs of valid masses in order of sorted valid masses
    
    # take log of all masses and of valid masses
    logx = np.log10(temp_x)
    logx_sort = np.log10(x_sort)
    
           
    out_file = fitspath0 + filename_dict['bin_info'].replace('.tbl', '.npz')
    
    # Check if binning npz file already exists    
    if exists(out_file):
        print('File exists')
        rinput = input('Do you want to delete file? Yes or no ')
        if rinput.lower() == 'yes':
            os.remove(out_file)
            print(fitspath0 + 'binning.npz deleted.')
        else:
            idx_file = np.load(out_file) 
            bin_edge = idx_file['bin_edge']
            bin_redge = idx_file['bin_redge']
            flux = idx_file['flux']
            wavelength = idx_file['wavelength']
            bin_ID = idx_file['bin_ID']
            count = len(bin_ID)
            N = idx_file['N']
            idx_file.close()
    
    '''
    Parameters/Variables (cont.):
        
    bin_ind = list of lists where each inner list contains the indices relative to temp_x that are in each bin
    start = number that starts each bin --> i.e. if bin_pts = 75, 112, 113, 300, 600, 1444, 1444 then
                                            start = 0, 75, 112, 113, 300, 600, 1444
    stop = number that ends each bin (exclusive) --> i.e. if bin_pts = 75, 112, 113, 300, 600, 1444, 1444 then
                                                     stop = 75, 112, 113, 300, 600, 1444, 1444
    bin_start = first mass value for each bin
    bin_stop = last mass value in each bin
    
    flux = 2D list of flux values from the stack_spectra function where each inner list is one bin
    wavelength = 2D list of wavelength values from the stack_spectra function where each inner list is one bin
    
    bin_edge = list of mass values that appear first (left edge) in each bin
    bin_redge = list of mass values that appear last (right edge) in each bin
    lowest_hbeta = list of Hbeta luminosity values that appear first (left edge) in each mass-LHbeta bin 
    highest_hbeta = list of Hbeta luminosity values that appear last (right edge) in each mass-LHbeta bin
    mass_avg = list of average masses from each bin
    bin_ID = list of bin numbers --> goes from 1 to number of bins
    N = list of number of values in each bin (i.e. same numbers as passed in bin_pts)
    count = number of mass bins
    count2 = number of mass-LHbeta bins 
    '''
    
    # If binning npz file doesn't exist, do binning           
    if not exists(out_file):
        bin_ind = []
        start = 0
        bin_start = logx_sort[start]
        bin_edge = []
        bin_redge = []
        flux = []
        wavelength = []
        mass_avg = []
        mass_median = []
        bin_ID = []
        N = []
        lowest_hbeta = []
        highest_hbeta = []
        lum_avg = []
        lhb_median = []
        count = 0
        count2 = 0
            
        while (bin_start < logx_sort[-1]):
            if adaptive:
                bin_pts = bin_pts_input[count]
            else:
                bin_pts = bin_pts_input
                    
            stop = start + bin_pts
            if ((stop + bin_pts) > len(logx_sort)):
                stop = len(logx_sort) - 1
            count += 1
            bin_stop = logx_sort[stop]
            bin_edge.append(bin_start)
            
            # for mass bins
            if hbeta_bin == False:
                bin_ind.append(valid_ind_sort[start:stop])
                
                # stack spectra
                if filename != False:
                    _, flx, wave = stack_spectra(filename, mname, indices=valid_ind_sort[start:stop])
                    flux.append(flx)
                    wavelength.append(wave)
                N.append(len(valid_ind_sort[start:stop]))
                mass_avg.append(np.mean(logx_sort[start:stop]))
                mass_median.append(np.median(logx_sort[start:stop]))
                bin_ID.append(count)    
                    
            # for mass-LHbeta bins
            else:
                '''
                Parameters/Variables (cont.):
                
                lum_sort = Hbeta luminosity values of the valid masses in mass-sorted order
                valid_hbeta = valid Hbeta luminosity indices that also have valid masses for one bin
                              --> sorted by valid mass
                median0 = median Hbeta luminosity value of the mass bin
                invalid_hbeta = invalid Hbeta luminosity indices
                temp_lower = indices of Hbeta luminosity values less than or equal to the bin median value
                temp_upper = indices of Hbeta luminosity values greater than the bin median value
                lower_idx = indices of samples with valid mass values and Hbeta luminosity values less than
                            or equal to the bin median for a given bin
                upper_idx = indices of samples with valid mass values and Hbeta luminosity values greater than
                            the bin median for a given bin
                '''
                
                lum_sort = lum[valid_ind_sort]
                valid_hbeta = np.where(lum_sort[start:stop] < 44)[0]
                median0 = np.median(lum_sort[start:stop][valid_hbeta])
                invalid_hbeta = np.where((lum_sort[start:stop] > 44) | 
                                         (np.isfinite(lum_sort[start:stop]) == False))[0]
                lum[valid_ind_sort[start:stop][invalid_hbeta]] = -1
                temp_lower = np.where(lum <= median0)[0]
                temp_upper = np.where(lum > median0)[0]
                lower_idx = np.array(list(set(valid_ind_sort[start:stop]) & set(temp_lower)))
                upper_idx = np.array(list(set(valid_ind_sort[start:stop]) & set(temp_upper)))
                    
                non_neg = np.where(lum[lower_idx] != -1)[0]
                non_neg = lower_idx[non_neg]
                lowest_hbeta.append(np.min(lum[non_neg]))
                highest_hbeta.append(np.max(lum[non_neg]))
                lowest_hbeta.append(np.min(lum[upper_idx]))
                highest_hbeta.append(np.max(lum[upper_idx]))
                lum_avg += [np.mean(lum[non_neg])] + [np.mean(lum[upper_idx])]
                lhb_median += [np.median(lum[non_neg])] + [np.median(lum[upper_idx])]
                bin_redge.append(logx_sort[start + len(lower_idx) - 1])
                bin_edge.append(logx_sort[start + len(lower_idx)])
                bin_ind.append(lower_idx)
                bin_ind.append(upper_idx)
                
                # stack spectra
                if filename != False:
                    _, lower_flx, lower_wave = stack_spectra(filename, mname, indices=lower_idx)
                    _, upper_flx, upper_wave = stack_spectra(filename, mname, indices=upper_idx)
                    flux += [lower_flx] + [upper_flx]
                    wavelength += [lower_wave] + [upper_wave]
                    N += [len(lower_idx), len(upper_idx)]
                mass_avg += [np.mean(logx[lower_idx])] + [np.mean(logx[upper_idx])]
                mass_median += [np.median(logx[lower_idx])] + [np.median(logx[upper_idx])]
                count2 += 1
                bin_ID.append(count2)
                count2 += 1
                bin_ID.append(count2)
                    
            start, bin_start = stop, bin_stop
            bin_redge.append(bin_stop)
            
        '''
        Note:
            
        The mass array saved here contains the 4140 masses (log(mass)) with the no mass indices replaced
        with interpolated values; it does NOT exclude sources excluded for binning (indices relative to
        original table). 
        The luminosity array saved here contains the 4140 luminosities (log(lum)) with the >44 and nans 
        replaced with -1.
        '''
        np.savez(out_file, mass=np.log10(temp_x), lum=lum, bin_ind=bin_ind, bin_start=bin_start, 
                 logM_min=bin_edge, logM_max=bin_redge, flux=flux, wavelength=wavelength, 
                 logM_avg=mass_avg, logM_median=mass_median, bin_ID=bin_ID, N_stack=N, 
                 logLHb_min=lowest_hbeta, logLHb_max=highest_hbeta, logLHb_avg=lum_avg, 
                 logLHb_median=lhb_median)
        
        if not adaptive:
            out_ascii = fitspath0 + str(bin_pts_input) + '_binning.tbl' 
        else:
            out_ascii = fitspath0 + filename_dict['bin_info']
        if not hbeta_bin:
            lowest_hbeta = np.zeros(len(bin_pts_input))
            highest_hbeta = np.zeros(len(bin_pts_input))
            lum_avg = np.zeros(len(bin_pts_input))
            lhb_median = np.zeros(len(bin_pts_input))
        column_names0 = remove_from_list(bin_names0, [bin_names0[-1]])
        all_column_names = tuple(merge_column_names(column_names0, bin_mzevolve_names0))
        table_stack = Table([bin_ID, N, bin_edge, bin_redge, mass_avg, mass_median, lowest_hbeta, 
                             highest_hbeta, lum_avg, lhb_median], names=all_column_names)        
        ascii.write(table_stack, out_ascii, format='fixed_width_two_line', overwrite=True)
        
    
    # Plotting stacked spectra
    xlim = [4250,4450]
    if (spectra_plot == True):
        for i in range(count):
            plt.subplot(np.ceil(count/2.0), 2, i+1)
            if not hbeta_bin:
                plt.plot(wavelength[i], flux[i], color='b', linestyle='solid')
                wavelen_idx = np.where((wavelength[i] > xlim[0]) & (wavelength[i] < xlim[1]))[0]
                max0 = np.max(flux[i][wavelen_idx])   
            else:
                plt.plot(wavelength[i*2], flux[i*2], color='b', linestyle='solid')
                wavelen_idx = np.where((wavelength[i*2] > xlim[0]) & (wavelength[i*2] < xlim[1]))[0]
                max0 = np.max(flux[i*2][wavelen_idx])
                plt.plot(wavelength[i*2+1], flux[i*2+1], color='orange', linestyle='solid')
                max0 = np.max(flux[i*2+1][wavelen_idx])
            plt.ylim(-0.05e-17, max0*1.1)
            plt.xlim(xlim)
            plt.axvline(x=5007, color='k', linestyle='dashed', alpha=0.5)
            plt.axvline(x=4363, color='r', linestyle='dashed', alpha=0.5)
            plt.suptitle(str(i))
            if hbeta_bin == False:
                an_text = 'N = '+str(N[i])
            else:
                an_text = 'N = [%i,%i]' %(N[i*2], N[i*2+1]) 
            plt.annotate(an_text, [0.05,0.95], xycoords='axes fraction', ha='left', va='top')
    else:
        plt.bar(bin_edge, N, align='edge', width=bin_redge)
        
    
    return bin_edge, flux




def gen_hist(fname, hist, nbins=10): 

#hist options: mass, chi2

    result = Table(ascii.read(fname))

    if hist == 'mass':
        mass = result['best.stellar.m_star']
        ind = np.where((np.isfinite(mass) == True) & (mass>0))[0]
        plt.hist(np.log10(mass[ind]), bins=nbins)
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