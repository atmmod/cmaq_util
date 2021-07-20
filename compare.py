### Python script for comparing CMAQ and CMAQ-hyd results
###! ~/miniconda3/bin/python
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import pprint
import re

# File path directory for original and hyd run (aero + vdiff only)
fname_orig = '/glade/work/edliu/models/CMAQ_v5.3.2/data/output_CCTM_v532_intel_Bench_2016_12SE1_cpdcp/CCTM_CONC_v532_intel_Bench_2016_12SE1_20160701.nc'
fname_hyd = '/glade/work/edliu/models/CMAQ_v5.3.2_hyd/data/output_CCTM_v532_intel_Bench_2016_12SE1_hyd_cpdcp/CCTM_CONC_v532_intel_Bench_2016_12SE1_20160701.nc'


# Import the netcdf4 files as Datasets
f_orig = Dataset(fname_orig)
f_hyd = Dataset(fname_hyd)

attr_dict = f_orig.__dict__  # show all the attributes in the cgrid file
varstr = attr_dict['VAR-LIST']  # select the variable list in the map, this is a string
varstr = re.sub("\s+", ",", varstr.strip())  # replace whitespaces with comma
varlst = varstr.split(',')

# The hyperdual perturbation
pert_hyperdual = 1e-3

def separate_var(conc, runtime=6):
    '''Separate concentration to x, dx1, dx2, dx1x2 information
    from the hyperdual netCDF file'''
    
#   For 6-hour run
    if runtime == 6: 
    	x_conc = conc[:7][:][:][:]
    	dx1_conc = conc[7:13][:][:][:]
    	dx2_conc = conc[13:19][:][:][:]
    	dx1x2_conc = conc[19:][:][:][:]
    
    
#   For 1-hour Run
    elif runtime == 1: 
    	x_conc = conc[1][:][:][:]
    	dx1_conc = conc[2][:][:][:]
    	dx2_conc = conc[3][:][:][:]
    	dx1x2_conc = conc[4][:][:][:]    
    	
    return x_conc, dx1_conc, dx2_conc, dx1x2_conc

diff_absolute = {}
diff_relative = {}
sens = {}

for name in varlst:
    orig_conc = f_orig.variables[name][:][:][:][:]  # (time + 1), L, C, R
    hyd_conc, dx1, dx2, dx1x2 = separate_var(np.squeeze(f_hyd.variables[name][:][:][:][:]), 6)
    diff_absolute[name] = np.average(np.abs(orig_conc[-1,:,:,:] - hyd_conc[-1,:,:,:]))
    diff_relative[name] = np.average(np.abs((orig_conc[-1,:,:,:] - hyd_conc[-1,:,:,:]) / orig_conc[-1,:,:,:]))
    if name == 'VLVOO1':
    	diff_matrix =np.abs((orig_conc[-1,:,:,:] - hyd_conc[-1,:,:,:]) / orig_conc[-1,:,:,:])
#     	for i in range(diff_matrix.shape[0]):
#     		for j in range(diff_matrix.shape[1]):
#     			for k in range(diff_matrix.shape[2]):
#     				if diff_matrix[i,j,k] > 1e8: 
#     					print(i)
#     					print(j)
#     					print(k)
#     					break
    
    	print('ORIG_CONC_VLVOO1')
    	print(orig_conc[-1, 34, 32, 0])
    	print('HYD_CONC_VLVOO1')
    	print(hyd_conc[-1, 34, 32, 0])
    
# 
pprint.pprint('Average Absolute Difference')
pprint.pprint(diff_absolute)
print('\n')
# 
pprint.pprint('Average Relative Difference')
pprint.pprint(diff_relative)
print('\n')
# 
# pprint.pprint('Sensitivity from hyd')
# pprint.pprint(sens)
# print('\n')
# 
