import numpy as np
from scipy.io import FortranFile
import xarray as xr

"""
This script converts the WRFSTAT netcdf file generated by Tak Yamaguchi's WRF-LES package to fortran binary.
These runs were made for the 2014 paper (Xiao, Gustafson and Wang, 2014, in JGR-Atmos).
The main work is to convert from LIWSE (liquid-ice water static energy) to THETAL (liquid water potential temperature).
"""

# constants
lv = 2.501e6 # J/kg
g = 9.81 # m/s**2
cp = 1005.7 # J/K/kg
p0 = 1.0e5 # Pa
kappa = 2./7. #

wrfstat = xr.open_dataset("/global/homes/x/xiao169/asrRE2/old_comptran/cam_wdm6_100_100h192_15v_0.2s/stat/comporef_wrf.0.5s.nc")

t1, t2 = 0, 1081
pressure = wrfstat.QT.T[:, t1:t2]*0.0 + wrfstat.p * 1.0e2
# print(pressure[:,2])
exner = (p0/pressure)**kappa
rho = wrfstat.RHO.T[:,t1:t2]
rtm = wrfstat.QT.T[:, t1:t2] * 1.0e-3
thlm = wrfstat.TPL.T[:, t1:t2]
wp2 = wrfstat.W2.T[:, t1:t2]
wp3 = wrfstat.W3.T[:, t1:t2]
wprtp = wrfstat.QTFLUX.T[:,t1:t2] # w/m^2
wprtp = wprtp/rho/lv # kg/kg m/s
wpthlp = wrfstat.TPLFLUX.T[:,t1:t2] # w/m^2
wpthlp = wpthlp/rho/cp # K m/s
rtp2 = wrfstat.QT2.T[:, t1:t2]*1.0e-6 # (g/kg)**2 to (kg/kg)**2
thlp2 = wrfstat.LIWSE2.T[:, t1:t2] # K**2
thlp2 = thlp2*(exner**2) # K**2
rtpthlp = wrfstat.LIWSEQT.T[:, t1:t2] # K**2
rtpthlp = rtpthlp*(exner*cp/lv) # K kg/kg
rtp3 = wrfstat.QT3.T[:, t1:t2] * 1.0e-9
thlp3 = wrfstat.LIWSE3.T[:, t1:t2]
thlp3 = thlp3*(exner**3) # K**3

f = FortranFile(
    # "input_comptran.bin",
    "input_comptran_skx.bin",
    "w",
)
for var in [
    rtm,
    thlm,
    wp2,
    wp3,
    wprtp,
    wpthlp,
    rtp2,
    thlp2,
    rtpthlp,
    rtp3,
    thlp3,
    pressure,
]:
    f.write_record(var.values.T.astype(np.float64))
f.close()
