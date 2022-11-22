import numpy as np
from scipy.io import FortranFile
import xarray as xr

'''
This script converts the netcdf file to fortran binary.
interpolate_sam.py prepares the SAM output for this script.
I think the same interpoloate_sam.py, or .ipynb at the time, is
used to produce input to CLUBB for the input_* experiments.
-03/01/2019
'''

sam = xr.open_dataset('atex_r1_i.nc')

t1, t2 = 0, 480

rtm = sam.QTO.T[0,0,:,t1:t2]
thlm = sam.THETAL.T[0,0,:,t1:t2]
wp2 = sam.W2.T[0,0,:,t1:t2]
wp3 = sam.W3.T[0,0,:,t1:t2]
wprtp = sam.QTOW.T[0,0,:,t1:t2]
wpthlp = sam.THELW.T[0,0,:,t1:t2]
rtp2 = sam.QT2.T[0,0,:,t1:t2]*1.0e-6
thlp2 = sam.THEL2.T[0,0,:,t1:t2]
rtpthlp = sam.QTOTHEL.T[0,0,:,t1:t2]*1.0e-3
rtp3 = sam.QTO3.T[0,0,:,t1:t2]*1.0e-9
thlp3 = sam.THEL3.T[0,0,:,t1:t2]
pressure = sam.PRES.T[0,0,:,t1:t2]*1.0e2

f = FortranFile('input_atex_r1_skx.bin', 'w')
for var in [rtm, thlm, wp2, wp3, wprtp, wpthlp,
            rtp2, thlp2, rtpthlp, rtp3, thlp3, pressure]:
# for var in [rtm, thlm, wp2, wp3, wprtp, wpthlp,
            # rtp2, thlp2, rtpthlp, pressure]:
    f.write_record(var.values.T.astype(np.float64))
f.close()
