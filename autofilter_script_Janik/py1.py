from astropy.io import fits
import numpy as np
import os
from starlink import kappa
from starlink import convert
import sys

map1 = fits.open(str(sys.argv[1]))
map2 = fits.open(str(sys.argv[2]))
map3 = fits.open(str(sys.argv[3]))
scalei = float(sys.argv[4])
scaleq = float(sys.argv[5])
scaleu = float(sys.argv[6])

fits.writeto('Synthobs_i_sc.fits',data=map1[0].data/(2.795*scalei),header=map1[0].header)
fits.writeto('Synthobs_ivar_sc.fits',data=map1[1].data/(2.795*2.795*scalei),header=map1[1].header)
fits.writeto('Synthobs_q_sc.fits',data=map2[0].data/(2.795*scaleq),header=map2[0].header)
fits.writeto('Synthobs_qvar_sc.fits',data=map2[1].data/(2.795*2.795*scaleq),header=map2[1].header)
fits.writeto('Synthobs_u_sc.fits',data=map3[0].data/(2.795*scaleu),header=map3[0].header)
fits.writeto('Synthobs_uvar_sc.fits',data=map3[1].data/(2.795*2.795*scaleu),header=map3[1].header)

