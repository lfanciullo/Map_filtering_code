from astropy.io import fits
import numpy as np
import os
from starlink import kappa
from starlink import convert

map1d = fits.open('ifake_4as_final.fits')
map1v = fits.open('ivarfake_4as_final.fits')
map2d = fits.open('qfake_4as_final.fits')
map2v = fits.open('qvarfake_4as_final.fits')
map3d = fits.open('ufake_4as_final.fits')
map3v = fits.open('uvarfake_4as_final.fits')
iref = fits.open('iext_serp.fits')
qref = fits.open('qext_serp.fits')
uref = fits.open('uext_serp.fits')

n1 = map1d[0].data.shape[0]
n2 = map1d[0].data.shape[1]
n1r = iref[0].data.shape[1]
n2r = iref[0].data.shape[2]

if n1r > n1 and n2r > n2:
	delt_c = n2r-n2
	delt_r = n1r-n1
	half_delt_c = np.floor(0.5*(delt_c))
	half_delt_r = np.floor(0.5*(delt_r))
	ind1_1 = int(half_delt_r-1)
	ind1_2 = int((delt_r - half_delt_r + 1)*(-1))
	ind2_1 = int(half_delt_c-1)
	ind2_2 = int((delt_c - half_delt_c + 1)*(-1))
	iblank = np.empty((n1r,n2r))*np.nan
	ivarblank = np.empty((n1r,n2r))*np.nan
	qblank = np.empty((n1r,n2r))*np.nan
	qvarblank = np.empty((n1r,n2r))*np.nan
	ublank = np.empty((n1r,n2r))*np.nan
	uvarblank = np.empty((n1r,n2r))*np.nan
	iblank[ind1_1:ind1_2,ind2_1:ind2_2] = map1d[0].data
	ivarblank[ind1_1:ind1_2,ind2_1:ind2_2] = map1v[0].data
	qblank[ind1_1:ind1_2,ind2_1:ind2_2] = map2d[0].data
	qvarblank[ind1_1:ind1_2,ind2_1:ind2_2] = map2v[0].data
	ublank[ind1_1:ind1_2,ind2_1:ind2_2] = map3d[0].data
	uvarblank[ind1_1:ind1_2,ind2_1:ind2_2] = map3v[0].data
	iref[0].data[0,:,:] = iblank
	iref[1].data[0,:,:] = ivarblank
	qref[0].data[0,:,:] = qblank
	qref[1].data[0,:,:] = qvarblank
	uref[0].data[0,:,:] = ublank
	uref[1].data[0,:,:] = uvarblank
elif n1 > n1r and n2 > n2r:
	delt_c = n2-n2r
	delt_r = n1-n1r
	half_delt_c = np.floor(0.5*(delt_c))
	half_delt_r = np.floor(0.5*(delt_r))
	ind1_1 = int(half_delt_r-1)
	ind1_2 = int((delt_r - half_delt_r + 1)*(-1))
	ind2_1 = int(half_delt_c-1)
	ind2_2 = int((delt_c - half_delt_c + 1)*(-1))
	iref[0].data[0,:,:] = map1d[0].data[ind1_1:ind1_2,ind2_1:ind2_2]
	iref[1].data[0,:,:] = map1v[0].data[ind1_1:ind1_2,ind2_1:ind2_2]
	qref[0].data[0,:,:] = map2d[0].data[ind1_1:ind1_2,ind2_1:ind2_2]
	qref[1].data[0,:,:] = map2v[0].data[ind1_1:ind1_2,ind2_1:ind2_2]
	uref[0].data[0,:,:] = map3d[0].data[ind1_1:ind1_2,ind2_1:ind2_2]
	uref[1].data[0,:,:] = map3v[0].data[ind1_1:ind1_2,ind2_1:ind2_2]


iref.writeto('fakei_trans_pW_sc_4as_al.fits')
qref.writeto('fakeq_trans_pW_sc_4as_al.fits')
uref.writeto('fakeu_trans_pW_sc_4as_al.fits')

