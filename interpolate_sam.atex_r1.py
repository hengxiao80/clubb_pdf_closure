import numpy as np
from netCDF4 import Dataset as ncfile
import subprocess

# get a copy of the SAM output to change
subprocess.call(['cp', '/global/homes/x/xiao169/cscratch/asrRE1/sam_runs/atex_r1/OUT_STAT/ATEX_r1.nc', 'atex_r1_i.nc'])
ds = ncfile('atex_r1_i.nc', 'r+')
# variables on zi-levels
varszi = ['W2', 'W3', 'THELW', 'THELWS', 'QTOW', 'QTOWS',
          'WWQC', 'WQQ', 'WTT', 'WQT', 'WWQ', 'WWT']
# There are two ways to do it,
# 1. interpolate all the variables that are output on zi-levels to z-levels
for var in varszi:
    data = ds.variables[var]
    temp = np.zeros((data.shape[0], data.shape[1]+1, 1, 1))
    temp[:,:-1,:,:] = data[:,:,:,:]
    data[:] = (temp[:,1:,:,:] + temp[:,:-1,:,:])*0.5
# add SGS contributions to get the total
ds.variables['THELW'][:] = ds.variables['THELW'][:] \
        + ds.variables['THELWS'][:]
ds.variables['W2'][:] = ds.variables['W2'][:] \
        + ds.variables['TKES'][:]*(2./3.)
ds.variables['QTOW'][:] = (ds.variables['QTOW'][:] \
        + ds.variables['QTOWS'][:])
# ds.variables['QTOW'][:] = (ds.variables['QTOFLX'][:] \
        # + ds.variables['QTOFLXS'][:])*1.0e-3
# We also need to process THETAL and QTO in SAM output
# rgas = 287.0
# cp = 1004.0
# lcond = 2.5104e+06
# pres = ds.variables['PRES'][:]
# qci = ds.variables['QCI'][:]*1.0e-3
# prespot = (1000./pres)**(rgas/cp)
# THETAL in the SAM output is actually liquid/ice potential temperature
# THEL on the other hand, is actual liquid potential temperature
# BUT only THELW, THELWS, THEL2 and THEL3 are output, not THEL itself.
# TL is also output, which is the actual prognostic variable in SAM.
# It is liquid/ice water static energy, so it has qci in it also.
# BUT since doicemicro is turned off in this run. So qci = 0.
# ds.variables['THETAL'][:] = ds.variables['THETAL'][:] - qci*prespot
ds.variables['QTO'][:] = ds.variables['QTO'][:]*1.0e-3
ds.variables['QTOW'][:] = ds.variables['QTOW'][:]*1.0e-3
ds.close()