# Map_filtering_code
Code to filter the large scales from SOFIA maps and compare them to POL-2.

## Available folders
### Map_files
Contains the FITS files of maps (real or simulated). 

Current content:
* I, Q, U for "toy" synthetic observations
    * Three structure types: filament, radially symmetric core, elongated cores + filament
    * Noiseless and noisy versions of each
    * Units of Jy/arcsec^2, realistic fluxes for HAWC+ band D (154 um)

### Map_figures

Current content:
* Quiver plots for the three synthetic structures (filament, radially symmetric core, elongated cores + filament), noiseless

### Figures_other
Other forms of data visualization

Current content:
* Plots of masks for FFT filtering
* FFT-filtered maps for one synthetic structure (filament)
