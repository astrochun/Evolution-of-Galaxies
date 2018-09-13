from astropy.io import fits, ascii
from astropy.table import Table, Column, vstack
import astropy.units as u
import astropy.coordinates as crd
import numpy as np

from getpass import getuser
if getuser() == 'carol':
    path = "C:\\Users\\carol\\Google Drive\\MZEvolve\\"
    path2 = path
else:
    path = "../DEEP2/"
    path2 = "../"    


#loading catalogs
aegis = Table(ascii.read(path + 'aegis-n2.deblend.v5.1.cat'))
zcat = Table(fits.getdata(path + 'DEEP2_Field1_zcat_ext.fits'))
all_line = Table(fits.getdata(path + 'DEEP2_Field1_all_line_fit.cat.fits'))


#cross matching coordinates
aegis_crd = crd.SkyCoord(ra = aegis['ra']*u.degree, dec = aegis['dec']*u.degree)
zcat_crd = crd.SkyCoord(ra = zcat['RA_DEEP']*u.degree, dec = zcat['DEC_DEEP']*u.degree)
indexa, indexz, d2d, d3d = zcat_crd.search_around_sky(aegis_crd, 1*u.arcsec)
print('aeigs length', len(indexa), 'zcat length', len(indexz))


#making a table with objno (zcat), id (aegis) & redshift (all_line)
field1 = aegis[indexa]
obj = Column(zcat[indexz]['OBJNO'])
redshift = Column(all_line[indexz]['ZBEST'])
field1.add_column(obj, index = 0)
field1.add_column(redshift, index = 2)
ascii.write(field1, path2 + 'magfiles/aegis_objno_id.mag', format = 'commented_header')


#deleting unnecessary columns
field1.remove_columns(['id', 'x', 'y', 'ra', 'dec', 'Kaper', 'eKaper', 'wCH4', 'wCH3', 'wCH2', 'wCH1', 'wKs', 'wK', 'wH', 'wH2', 'wH1', 'wJ', 'wJ3', 'wJ2', 'wJ1', 'wz', 'wI', 'wR', 'wG', 'wU', 'wNUV', 'wFUV', 'ftot24um_uJy', 'f24um_uJy', 'e24um_uJy', 'w24um', 'wmin', 'wmin_irac', 'z_spec', 'star_flag', 'ap_col', 'ap_tot', 'totcor', 'K_ellip', 'K_theta_J2000', 'K_R50', 'K_class_star', 'K_flags', 'UH2_flags', 'Near_Star', 'CH1_contam', 'CH2_contam', 'CH3_contam', 'CH4_contam', 'NUV_contam', 'FUV_contam',  'contam_flag', 'nchild', 'id_parent', 'use'])


#renaming columns for pcigale
field1['ZBEST'].name = 'redshift'
field1['OBJNO'].name = 'id'
ind = np.arange(2,len(field1.colnames),2)               #only taking band names
names = np.array(field1.colnames)[ind]
new_names = ['IRAC4', 'IRAC3', 'IRAC2', 'IRAC1', 'WIRCam_cfh8302_Ks', 'NEWFIRM_k', 'WIRCam_cfh8201_H', 'NEWFIRM_h2', 'NEWFIRM_h1', 'WIRCam_cfh8101_J', 'NEWFIRM_j3', 'NEWFIRM_j2', 'NEWFIRM_j1', 'CF_z', 'CF_i', 'CF_r', 'CF_g', 'CF_u', 'NUV', 'FUV']
for a,b in zip(names, new_names):
    if a != b:
        field1[a].name = b
    field1['e'+a].name = b+'_err'


#converting flux to mJy
print(field1.colnames)
for col in new_names:
    w_data = np.where(field1[col] >= 0)[0]                      #avoiding -9999
    new_col = field1[col][w_data]*(10**-0.44)       #arbitrary units to microJy
    field1[col][w_data] = new_col/1000                          #microJy to mJy
    #same for error columns
    w_data_err = np.where(field1[col+'_err'] >= 0)[0]
    new_err = field1[col+'_err'][w_data_err]*(10**-0.44)
    field1[col+'_err'][w_data_err] = new_err/1000


ascii.write(field1, path2 + 'magfiles/aegis5.1_deep2_crossmatch.mag', format = 'commented_header')
print('completed writing the mag file. Length:', len(field1))

