# Map_filtering_code

This page hosts the results of our efforts to filter the large scales from SOFIA maps and compare them to POL-2.

The code for comparing filtered and non-filtered maps (including Fourier transforms) is fft_tests.py

The code for the filtering itself can be found in autofilter_script_Janik/auto_script.txt

## Available folders

### autofilter_script_Janik

Contains a semi-automated script for filtering IQU maps using Starlink, as well as the FITS and SDF files for the Serpens region used as "empty field".

### Map_files

FITS files of maps (real or simulated). 

Contains the following:
* I, Q, U for "toy" synthetic observations, to be filtered
    * Three structure types: filament, radially symmetric core, elongated cores + filament
    * Noiseless and noisy versions of each
    * Units of Jy/arcsec^2, realistic fluxes for HAWC+ band D (154 um)
* Filtered versions of the I, Q, U toy maps (noiseless version only at the moment)
    * Note: the method used for the filament (by Kate) is slightly different from the one for the other maps (by Janik), most notably the scaling used.

### Figures

Figures for maps, plots comparing before/after filtering, and other forms of data visualization.

Current content:
* Cropped I maps of the three synthetic observations, before and after filtering;
* Maps of the before/after filtering ratio for I;
* FFT of the filtered and unfiltered I maps;
* Comparison of the radial profiles for the FFT-transformed I maps. 

#### Subfolder: Synthetic_observation_maps_original

Quiver plots for the three synthetic structures (filament, radially symmetric core, elongated cores + filament), noiseless.

#### Subfolder: Old_pictures

Older and no longer relevant figures are moved here to avoid clutter. 

Current content:
* "Filtering simulations" obtained by masking FFT-transformed map (main issue: severe ringing).

### Reports_and_presentations

Presentations on the intermediate analysis results, used for internal reports or to present at DR meetings. 
