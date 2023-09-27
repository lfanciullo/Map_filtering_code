
import numpy as np
import copy
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from aplpy import FITSFigure

folder_maps = 'Map_files/'
folder_image_results = 'Figures_other/'

mapfile = folder_maps + \
    #'Synthobs_core_nonoise_converted_i+var_Jy-arcs2.fits' # Core, I (no noise), unfiltered
    'Synthobs_filament_nonoise_converted_i+var_Jy-arcs2.fits' # Filament, I (no noise), unfiltered


synthmap = fits.open(mapfile)
synthmap_data = copy.deepcopy(synthmap[0].data)

data = synthmap_data[0:273,0:273] # Select square region
data[np.where(~np.isfinite(data))] = 0.
# 

fourier_image = np.fft.fftn(data)
fourier_image_shifted= np.fft.fftshift(fourier_image)
fourier_amplitudes = np.abs(fourier_image_shifted)**2


# Plot original image
wcs = WCS(synthmap[0].header)
plt.subplot(projection = wcs)
plt.imshow(data/np.nanmax(data), cmap = 'rainbow', norm = colors.LogNorm(vmin = 1e-3, vmax = 1.))
plt.xlabel('RA (J2000)', fontsize = 'large')
plt.ylabel('Dec (J2000)', fontsize = 'large')
cbar = plt.colorbar()
cbar.set_label('Intensity (normalized)', size = 'large')
#plt.savefig('Image_orig.png')


# Plot Fourier amplitude
plt.subplot()
plt.imshow(fourier_amplitudes/np.nanmax(fourier_amplitudes), cmap = 'viridis', norm = colors.LogNorm(vmin = 1e-6, vmax = 1))
cbar = plt.colorbar()
cbar.set_label('Fourier amplitude (normalized)', size = 'large')
plt.savefig(folder_image_results + 'Image_FFT-amplitudes_shifted.png')

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
