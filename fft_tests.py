
import numpy as np
import copy
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from aplpy import FITSFigure


### READ MAPS (filtered + unfiltered) ###
folder_maps = 'Map_files/'
folder_image_results = 'Figures_other/'

# Original (simulated) data
mapfile = folder_maps + 'Synthobs_filament_nonoise_converted_i+var_Jy-arcs2.fits'       # Filament, I (no noise), unfiltered
#mapfile = folder_maps + 'Synthobs_core_nonoise_converted_i+var_Jy-arcs2.fits'           # Core, I (no noise), unfiltered
#mapfile = folder_maps + 'Synthobs_filament+cores_nonoise_converted_i+var_Jy-arcs2.fits' # Clumps+filament, I (no noise), unfiltered
synthmap = fits.open(mapfile)
synthmap_data = copy.deepcopy(synthmap[0].data)
wcs = WCS(synthmap[0].header)

# Data filtered with POL-2 pipeline
# To add: core
mapfile_filt = folder_maps + 'Synth_filament_filtered_iext_Jysqa.fits' # Filament, I (no noise), filtered
#mapfile_filt = folder_maps + 'Synth_filament+cores_filtered_iext_filtered_model.fits' # Clumps+filament, I (no noise), filtered
synthmap_filt = fits.open(mapfile_filt)
synthmap_filt_data = copy.deepcopy(np.squeeze(synthmap_filt[0].data))
wcs_filt = WCS(synthmap_filt[0].header).dropaxis(2)  # Dropping shallow 3rd axis


### PREPARE DATA ###

# NOTA: Need to edit WCS for proper coordinates in plot
nside_init = 273                         # Size of original (unfiltered) FITS image
nside_init_filt = 267                    # Size of filtered FITS image
nside = 151                              # Size of region selected for analysis
#image = synthmap_data[0:nside,0:nside]
image = synthmap_data[int((nside_init-nside)/2):int(nside+(nside_init-nside)/2),  # Select square region
                      int((nside_init-nside)/2):int(nside+(nside_init-nside)/2)]
image_filt = synthmap_filt_data[int((nside_init_filt-nside)/2):int(nside+(nside_init_filt-nside)/2),  # Select square region
                                int((nside_init_filt-nside)/2):int(nside+(nside_init_filt-nside)/2)]
#image[np.where(~np.isfinite(image))] = 0. # Get rid of NaNs
# Some relevant numbers:
#   Pixel size: 3.4 arcsec (synthmap[0].header['CDELT1']*3600)
#   Map size: 273 * 3.4 ~ 928 arcsec
#   180" scale ~ map size / 5


### Fourier transforms ###
fourier_image = np.fft.fft2(image)                      # FFT
fourier_image_shifted = np.fft.fftshift(fourier_image)  # Shift --> low spatial frequencies now at the center of the image
fourier_image_filt = np.fft.fft2(image_filt)
fourier_image_filt_shifted = np.fft.fftshift(fourier_image_filt)

#image_filt(where(~np.isfinite(image_filt)) = np.nanmin(image_filt))  # Getting rid of NaNs

# Update this
reconstructed_image = np.fft.ifft2(fourier_image)                         # Sanity check
reconstructed_image_filt = np.fft.ifft2(fourier_image_filt)
max_deviation = np.max(np.angle(reconstructed_image[where(image != 0)]))
max_deviation_filt = np.max(np.angle(reconstructed_image_filt[where(image_filt != 0)]))
print('Max deviation in inverse FFT (original image) = ', max_deviation)
print('Max deviation in inverse FFT (original image) = ', max_deviation_filt)

x, y = np.meshgrid(np.arange(nside), np.arange(nside))     # Create array of "distance" from center (--> frequencies)
dist = np.sqrt((x-int(nside/2))**2 + (y-int(nside/2))**2)

# Filters
filtsize = 3. #5.
floor = 1e-6
filt_boxcar = 1 - (dist < filtsize).astype(int) * (1-floor)
filt_gauss = 1 - exp(-(dist/filtsize)**2) * (1-floor)
filt_quartic = 1 - exp(-(dist/filtsize)**4) * (1-floor)
# Visual test of filters
#range_min = 69; range_max = 82
#plt.subplot()
#plt.imshow(dist[range_min:range_max, range_min:range_max], cmap = 'viridis', vmin = 0., vmax = filtsize)
#cbar = plt.colorbar()
#cbar.set_label('Distance (pixels)', size = 'large')
#plt.contour(filt_boxcar[range_min:range_max, range_min:range_max], [.1, .5, .9], colors = 'k')

fourier_image_shifted *= filt_boxcar                     # Manually eliminating large scales
image_FFTfiltered = np.fft.ifft2(fourier_image_shifted)  # Inverse Fourier transform --> FFT-filtered image

#kfreq = np.fft.fftfreq(nside) * nside


### PLOTS ###

# Plot original image
#wcs = WCS(synthmap[0].header)
plt.subplot(projection = wcs)
plt.imshow(image, cmap = 'rainbow', norm = colors.LogNorm(vmin = 1e-2, vmax = 2.))
#plt.imshow(image/np.nanmax(image), cmap = 'rainbow', norm = colors.LogNorm(vmin = 1e-2, vmax = 1.)) # Normalized version
plt.xlabel('RA (J2000)', fontsize = 'large')
plt.ylabel('Dec (J2000)', fontsize = 'large')
cbar = plt.colorbar()
cbar.set_label('Intensity', size = 'large')
#plt.savefig(folder_image_results + 'Image_orig.png')

# Filtered / unfiltered ratio
plt.subplot()
plt.title('$I$ ratio for filament (filtered/unfiltered)', size = 'x-large')
plt.imshow(image_filt/image, cmap = 'viridis')#, norm = colors.LogNorm(vmin = 1e-2, vmax = 2.))
cbar = plt.colorbar()
cbar.set_label('Ratio', size = 'large')
#plt.savefig(folder_image_results + 'Image_ratio.png')


# Plot Fourier amplitude
fourier_amplitudes = np.abs(fourier_image_shifted)**2
fourier_amplitudes_filt = np.abs(fourier_image_filt_shifted)**2
plt.subplot()
plt.imshow(fourier_amplitudes, cmap = 'viridis', norm = colors.LogNorm(vmin = 1e-3, vmax = 1e5))
#plt.imshow(fourier_amplitudes_filt/fourier_amplitudes, cmap = 'viridis', norm = colors.LogNorm(vmin = 1e-6, vmax = 1)) # Ratio
#plt.imshow(fourier_amplitudes/np.nanmax(fourier_amplitudes), cmap = 'viridis', norm = colors.LogNorm(vmin = 1e-6, vmax = 1)) # Normalized
plt.title('FFT for filament ($I$, unfiltered)', size = 'x-large')
cbar = plt.colorbar()
cbar.set_label('Fourier amplitude', size = 'large')
#cbar.set_label('Fourier amplitude ratio: filtered/unfiltered', size = 'large') # Ratio
#cbar.set_label('Fourier amplitude (normalized)', size = 'large') # Normalized
#plt.savefig(folder_image_results + 'Image_FFT-amplitudes_shifted.png')

# Plot inverse Fourier transform
data2plot = \
    np.abs(image_FFTfiltered) #/ image
wcs = WCS(synthmap[0].header)
plt.subplot(projection = wcs)
plt.imshow(data2plot, cmap = 'rainbow', norm = colors.LogNorm(vmin = 1e-2, vmax = 2.))
#plt.imshow(np.abs(reconstructed_image)/np.nanmax(np.abs(reconstructed_image)), cmap = 'rainbow',
#           norm = colors.LogNorm(vmin = 1e-2, vmax = 1.))
plt.xlabel('RA (J2000)', fontsize = 'large')
plt.ylabel('Dec (J2000)', fontsize = 'large')
cbar = plt.colorbar()
cbar.set_label('Intensity', size = 'large')
#plt.savefig(folder_image_results + 'Image_iFFT.png')

'''
# APLPY version -- colorbar is not flexible :(
fig1 = FITSFigure(data)
fig1.show_colorscale(cmap='rainbow') #, vmin = -2.5e-2, vmax = 2.5e-1)
#fig1.set_title('Original')
#fig1.axis_labels.set_font(size='xx-large')
#fig1.tick_labels.set_font(size='x-large')
fig1.add_colorbar()
fig1.colorbar.set_font(size='x-large')
#fig1.colorbar.set_axis_label_font(size='xx-large')
#fig1.colorbar.set_axis_label_text('...')
'''
