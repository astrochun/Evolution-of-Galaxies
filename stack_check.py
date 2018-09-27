from astropy.io import fits, ascii
import numpy as np
import matplotlib.pyplot as plt

from pylab import subplots_adjust

from matplotlib.backends.backend_pdf import PdfPages

nrows = 10
ncols = 1
def stack_check(fname, bin_array_file, out_pdf, mname=''):    
    hdu = fits.open(fname)
    image = hdu[0].data
    header = hdu[0].header
    wavelength = header["CRVAL1"] + np.array(header["CDELT1"])*range(header['NAXIS1'])
    
    npz_zero = np.load(bin_array_file)
    indices = npz_zero["indices"]
    
    if mname != '':
        hdu = fits.open(mname)
        mask_image = hdu[0].data
        mask_header = hdu[0].header
        
        zero_x, zero_y = np.where(image == 0)
        mask_image[zero_x, zero_y] = 1
        
        image = np.ma.array(image, mask=mask_image)
    
    sigma_array = np.zeros(len(indices))
    
    pp = PdfPages(out_pdf)
    for ii in range(len(indices)):
        row = np.int(np.floor((ii/ncols) % nrows))
        if ii % (nrows*ncols) == 0:
            fig, ax = plt.subplots(nrows=10, ncols=1)
        
        flux = image[indices[ii]]
        #if ii == 0: print(flux.shape, row)
        ax[row].axvline(x=5006.8, color="black", linestyle="dashed", linewidth=0.5)
        ax[row].axvline(x=3726.16, color="black", linestyle="dashed", linewidth=0.5)
        ax[row].axvline(x=3728.91, color="black", linestyle="dashed", linewidth=0.5)
        ax[row].axvline(x=4861.32, color="black", linestyle="dashed", linewidth=0.5)
        ax[row].plot(wavelength, flux/1e-17)
        
        ax[row].set_xlabel("Wavelength")
        ax[row].set_ylabel("Flux")
        ax[row].set_xlim(3700, 5300)
        #ax[row].set_title("Stacked Spectra of %s" % fname)
        
        sigma_array[ii] = np.ma.std(flux/1e-17)
        ax[row].annotate('sigma = %.3f' % sigma_array[ii], [0.95,0.95], xycoords='axes fraction',
                         ha='right', va='top')

        
        if (ii % (nrows*ncols) == nrows*ncols-1) or (ii == len(indices)):
            subplots_adjust(left=0.05, bottom=0.05, top=0.99, right=0.99, hspace=0.03)
            fig.set_size_inches(8,11)
            fig.savefig(pp, format='pdf')
            
    pp.close()
    hdu.close()
    
    fig, ax = plt.subplots()
    ax.hist(sigma_array, bins=50)
    fig.savefig(out_pdf.replace(".pdf", "_hist.pdf"))