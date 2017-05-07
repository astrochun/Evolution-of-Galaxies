from astropy.io import fits, ascii
from astropy.table import Table, vstack, Column
import numpy as np
import library as lib


# Loading zcat files
zcat1 = Table(fits.getdata('../DEEP2/DEEP2_Field1_zcat_ext.fits'))
zcat2 = Table(fits.getdata('../DEEP2/DEEP2_Field2_zcat_ext.fits'))
zcat3 = Table(fits.getdata('../DEEP2/DEEP2_Field3_zcat_ext.fits'))
zcat4 = Table(fits.getdata('../DEEP2/DEEP2_Field4_zcat_ext.fits'))


# Loading all_line files
all_line1 = Table(fits.getdata('../DEEP2/DEEP2_Field1_all_line_fit.cat.fits'))
all_line2 = Table(fits.getdata('../DEEP2/DEEP2_Field2_all_line_fit.cat.fits'))
all_line3 = Table(fits.getdata('../DEEP2/DEEP2_Field3_all_line_fit.cat.fits'))
all_line4 = Table(fits.getdata('../DEEP2/DEEP2_Field4_all_line_fit.cat.fits'))


# Creating 1 table of each for easier cross-matching
zcat = vstack([zcat1, zcat2, zcat3, zcat4])
all_line = vstack([all_line1, all_line2, all_line3, all_line4])


# Replacing RA_DEEP column (the 2nd column) by redshift values
# The column will be renamed along with other columns
zcat['RA_DEEP'] = all_line['ZBEST']


# Deleting unnecessary columns 
# Renaming the remaining columns
colnames = zcat.colnames
colnames = colnames[0:2] + colnames[5:21] + colnames[28:29]
zcat.keep_columns(colnames)
new_names = ['id', 'redshift', 'cfht_B', 'cfht_R', 'cfht_I', 'cfht_B_err', 'cfht_R_err', 'cfht_I_err', 'CF_u', 'CF_g', 'CF_r', 'CF_i', 'CF_z', 'CF_u_err', 'CF_g_err', 'CF_r_err', 'CF_i_err', 'CF_z_err', 'source', 'u_prime', 'g_prime', 'r_prime','i_prime', 'z_prime', 'u_prime_err', 'g_prime_err', 'r_prime_err', 'i_prime_err', 'z_prime_err']
for nm, new_nm in zip(colnames, new_names):
    zcat[nm].name = new_nm


#Replacing duplicate id's
ID = lib.duplicates(zcat['id'])
zcat.remove_column('id')
zcat.add_column(Column(ID, name='id'), 0) 


# Differentiating between sloan & cfht ugriz
sdss_ind = [ii for ii in xrange(len(zcat)) if 'SDSS' in zcat['source'][ii]]
for a, b in zip(new_names[19:29], new_names[8:18]):
    # a = SDSS ugriz and their errors
    # b = CFHT ugriz and their errors
    zcat[a] = [-9999.0]*len(zcat)            #creating new column for SDSS data
    zcat[a][sdss_ind] = zcat[b][sdss_ind]    #at sdss_ind, the value is sloan's
    zcat[b][sdss_ind] = -9999.0       #so at sdss_ind there is no data for CFHT
zcat.remove_column('source')


# Converting flux to mJy
colnames = zcat.colnames
fluxes = [ii for ii in xrange(len(colnames)) if 'err' not in colnames[ii]][2:]
errs = [ii for ii in xrange(len(colnames)) if 'err' in colnames[ii]]
for ff, ee in zip(fluxes, errs):
	w_data = np.where(zcat[colnames[ff]] >= 0)[0]                      #avoiding -9999.0
	m_col = zcat[colnames[ff]][w_data]
	del_m = zcat[colnames[ee]][w_data]
	
	flux = 10**(-0.4*(m_col-23.9))                                     #arbitrary unit to microJy
	pos = 10**(-0.4*(m_col-delm-23.9)) - flux
	neg = flux - 10**(-0.4*(m_col+delm-23.9))
	err = np.mean(pos,neg)
	zcat[colnames[ff]][w_data] = flux/1000                                   #microJy to mJy
    
def sadalk():
    

print 'Table conversion completed'

ascii.write(zcat, 'stacked_deep_fields.mag', format = 'commented_header')
print 'wrote original file'
ascii.write(zcat[0:30], 'stacked_deep_fields_sh.mag', format = 'commented_header')
print 'wrote short file'
