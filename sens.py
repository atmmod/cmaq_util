import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from netCDF4 import Dataset
import re
import pprint

# This script verify the perturbation done to real mode CMAQ
# and does the sensitivity comparison between real and hyperdual
# versions of CMAQ

def get_namelist(f):
	'''
	Get the species namelist from a netCDF file
	'''
	
	attr_dict = f.__dict__  # show all the attributes in the cgrid file
	varstr = attr_dict['VAR-LIST']  # select the variable list in the map, this is a string
	varstr = re.sub("\s+", ",", varstr.strip())  # replace whitespaces with comma
	varlst = varstr.split(',')
	
	return varlst
	
	
def calc_pert(f_orig, f_pert, pert_code):
	'''Calculate the perturbation amount for the 
	original CMAQ
	'''
	
	# (2, 35, 80, 100)
	orig_conc = f_orig.variables[pert_code][:][:][:][:]
	pert_conc = f_pert.variables[pert_code][:][:][:][:]	
	
	# Perturbation happens at timestep = 0, layer = 0
	pert = np.sum(pert_conc[0, :, :, :], axis=0) - np.sum(orig_conc[0, :, :, :], axis=0)
	
	# For                                                                    1-layer run:  
# 	pert = pert_conc[0, 0, :, :] - orig_conc[0, 0, :, :]

	return pert
	
	
def real_sens(orig_dir, pert_dir, pert_code):
	''' 
	This function checks whether the perturbation is done correctly
	and outputs a sensitivity map and a sensitivity matrix for each species  
	'''
	
	# import the two netCDF files
	f_orig = Dataset(orig_dir)
	f_pert = Dataset(pert_dir)
	
	# Get the species namelist from the file, (224, 1)
	namelist = get_namelist(f_orig)
	
	# Get the perturbation, (80, 100)
	perturbation = calc_pert(f_orig, f_pert, pert_code)
	
	sens_map = {}
	aecj_conc = f_orig.variables['AECJ'][:][:][:][:]
	
	# (221, 35, 80, 100) matrix
	# Initialize matrices for first and second order derivatives
	sens_matrix = np.zeros((len(namelist), aecj_conc.shape[2], aecj_conc.shape[3]))

	# Start looping over species
	for i, name in enumerate(namelist):
		orig_conc = f_orig.variables[name][:][:][:][:]
		pert_conc = f_pert.variables[name][:][:][:][:]
		
		# For each layer ... 
#		for l in range(orig_conc.shape[1]):
		sens_matrix[i, :, :] = (np.sum(pert_conc[-1, :, :, :], axis=0) - np.sum(orig_conc[-1, :, :, :], axis=0)) / perturbation[:, :]
		
	 	# Sens matrix: (221, 35, 80, 100)
	 	# Average the sensitivity over all layers, columns, and rows
		sens_map[name] = np.average(sens_matrix[i, :, :])
	
	return sens_matrix, sens_map


pert_hyperdual = 35 * 1e-3


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



def hyd_sens(hyd_dir, pert, runtime):
	'''
	This function calculates the hyperdual CMAQ sensitivity
	and outputs a sensitivity map for each species
	'''
	
	f_hyd = Dataset(hyd_dir)
	namelist = get_namelist(f_hyd)
	
	sens_map = {}
	hyd_aecj_conc, _, _, _ = separate_var(f_hyd.variables['AECJ'][:][:][:][:], runtime)
	# 221 x 35 x 80 x 100
	
	# For 6-hour run
	hyd_sens_matrix = np.zeros((len(namelist), hyd_aecj_conc.shape[2], hyd_aecj_conc.shape[3]))
	
	# For 1-hour run: 
# 	hyd_sens_matrix = np.zeros((len(namelist), hyd_aecj_conc.shape[1], hyd_aecj_conc.shape[2]))
	
	for i, name in enumerate(namelist): 
		hyd_conc, dx1, dx2, dx1x2 = separate_var(f_hyd.variables[name][:][:][:][:], runtime)
		# 221x80x100 result: 
# 		for l in range(hyd_conc.shape[0]):
		
		# For 6-hour run'
		hyd_sens_matrix[i, :, :] = np.sum(dx1[-1, :, :, :], axis=0) / pert_hyperdual
		
		# For 1-hour run
# 		hyd_sens_matrix[i, :, :] = np.sum(dx1[:, :, :], axis=0) / pert_hyperdual
		sens_map[name] = np.average(hyd_sens_matrix[i, :, :])        	    
	
	return hyd_sens_matrix, sens_map


# with open('cmaq_sens_vdiffonly_cd.txt', 'w') as orig_sens_file: 
# 	for k, v in cmaq_sens.items():
# 		orig_sens_file.write('%s:%s\n' % (k, v))
	
# with open('cmaq_hyd_sens_aero_wo_orgaer.txt', 'w') as hyd_sens_file: 
# 	for k, v in cmaq_hyd_sens.items():
# 		hyd_sens_file.write('%s:%s\n' % (k, v))


def output_comparison_cell(orig_directory, sens_matrix, sens_matrix_hyd):
	'''Compare the sensitivity coefficients for each cell, 
	The output will be in 221x80x100 format
	sens_matrix: 221x80x100
	sens_matrix_hyd: 221x80x100
	output: 221x80x100
	'''
	# Get the namelist for all species
	f_orig = Dataset(orig_directory)
	namelist = get_namelist(f_orig)
	
	result = np.zeros((len(namelist), sens_matrix.shape[1], sens_matrix.shape[2]))
	
	for i, name in enumerate(namelist):
		result[i, :, :] = (sens_matrix[i, :, :] - sens_matrix_hyd[i, :, :]) / sens_matrix_hyd[i, :, :]
	
	return result



def output_comparison(sens_map, sens_map_hyd):
	'''Output statistics comparing the averaged sensitivity coefficients
	for the FD and HYD sensitivities'''
	
	# Create a diff map
	diff_map = {}
	
	for key, val in sens_map.items(): 
	    
	    # set a threshold for the sensitivity coefficients
		if abs(val) > 1e-5:
			print(key)
		
		# Technically, the mean normalized bias is calculated
			diff_map[key] = (val - sens_map_hyd[key]) / sens_map_hyd[key]
	
	return diff_map


def get_latlon():
	latlonfilename = '/glade/work/edliu/models/latlon.1day'
	latlonfile = Dataset(latlonfilename, mode='r', open=True)

	gridLat =  np.squeeze(latlonfile.variables['LAT'][:][:][:][:])
	gridLon =  np.squeeze(latlonfile.variables['LON'][:][:][:][:])

	last1ind,last2ind = gridLon.shape

	loncrn = [gridLon[0,0],gridLon[last1ind-1,last2ind-1]]
	latcrn = [gridLat[0,0],gridLat[last1ind-1,last2ind-1]]
	lonmid = gridLon[last1ind//2,last2ind//2]
	latmid = gridLat[last1ind//2-1,last2ind//2-1]
	
	return loncrn, latcrn, lonmid, latmid

# Make the species that satisfy the conditions a list of names
# species_list = list(diff_map.keys())
# print(species_list)

def find_species_index(orig_dir, species_list):
	'''Find the indices of species in the species list'''
	f_orig = Dataset(orig_dir)
	namelist = get_namelist(f_orig)
	
	numlist = []
	for i, name in enumerate(namelist):
		if name in species_list:
			numlist += [i]
	return numlist
	
	
def draw_pdf(orig_directory, species_list, diff_matrix):
	
	numlist = find_species_index(orig_directory, species_list)

	# Get the latlon information from the function above
	loncrn, latcrn, lonmid, latmid = get_latlon()
	
	for i in range(diff_matrix.shape[0]):
		if i in numlist: 
			with PdfPages('/glade/work/edliu/models/analysis'+ species_list[numlist.index(i)]+'.pdf') as pdf:
			
				# Formulate the map
				m = Basemap(projection='lcc', llcrnrlon=loncrn[0], llcrnrlat=latcrn[0],\
				urcrnrlon=loncrn[1], urcrnrlat=latcrn[1], lon_0=lonmid,\
				lat_0=latmid, resolution='l')	
            
            # draw coastlines, state and country boundaries, edge of map.
				m.drawcoastlines()
				m.drawstates()
				m.drawcountries()
				ny = diff_matrix.shape[1]
				nx = diff_matrix.shape[2]
				lons, lats = m.makegrid(nx, ny)
				x, y = m(lons, lats)
				
				upperbound = max(np.abs(np.nanmax(diff_matrix[i,:,:])), np.abs(np.nanmin(diff_matrix[i,:,:])))
				print(f'The upperbound for {species_list[numlist.index(i)]} is: {upperbound:1.1f}')
				lowerbound = -upperbound
				clevs = np.arange(lowerbound,upperbound+(upperbound - lowerbound)/20,(upperbound - lowerbound)/20)
				cmap=cm.seismic

				cs = m.contourf(x,y,diff_matrix[i,:,:],clevs,cmap=cmap,extend='both') 
				cbar = m.colorbar(cs,location='right',pad="0.00%")
				cbar.formatter.set_powerlimits((0, 0))
				cbar.update_ticks()
				cbar.ax.tick_params(labelsize=14)
				cbar.set_label('rel_diff', fontsize=14)

				plt.title('SENS '+species_list[numlist.index(i)]+' (real - hyd) / hyd')
				plt.tight_layout()
				pdf.savefig()
				plt.close()


def main():
	'''The main function of the file. Call the functions needed here'''
	orig_directory = '/glade/work/edliu/models/CMAQ_v5.3.2/data/output_CCTM_v532_intel_Bench_2016_12SE1_cpdcp_dec2p5/CCTM_CONC_v532_intel_Bench_2016_12SE1_20160701.nc'
	pert_directory = '/glade/work/edliu/models/CMAQ_v5.3.2/data/output_CCTM_v532_intel_Bench_2016_12SE1_cpdcp_inc2p5/CCTM_CONC_v532_intel_Bench_2016_12SE1_20160701.nc'
	hyd_directory = '/glade/work/edliu/models/CMAQ_v5.3.2_hyd/data/output_CCTM_v532_intel_Bench_2016_12SE1_hyd_cpdcp/CCTM_CONC_v532_intel_Bench_2016_12SE1_20160701.nc'	
	
	# Calculate the sensitivity matrices and maps
	cmaq_sens_matrix, cmaq_sens_map = real_sens(orig_directory, pert_directory, 'AECJ')
	cmaq_hyd_sens_matrix, cmaq_hyd_sens_map = hyd_sens(hyd_directory, pert_hyperdual, 6)
	
	# Calculate the differences between sensitivity matrices and maps
	diff_map = output_comparison(cmaq_sens_map, cmaq_hyd_sens_map)
	diff_matrix = output_comparison_cell(orig_directory, cmaq_sens_matrix, cmaq_hyd_sens_matrix)
	
			
	print('Sens diff:')
	pprint.pprint(diff_map)
	print('FD sens:')
	pprint.pprint(cmaq_sens_map)
	print('HYD sens:')
	pprint.pprint(cmaq_hyd_sens_map)
	
	
# 	for i in range(cmaq_sens_matrix.shape[1]):
# 		for j in range(cmaq_sens_matrix.shape[2]):
# 			v = (cmaq_sens_matrix[159, i, j] - cmaq_hyd_sens_matrix[159, i, j]) / cmaq_hyd_sens_matrix[159, i, j]
# 			if abs(v) > 1000:
# 				print('CD')
# 				print(cmaq_sens_matrix[159, i, j])
# 				print('HYD')
# 				print(cmaq_hyd_sens_matrix[159, i, j])
# 				print(i)
# 				print(j)

# AECJ index is 139
# ANH4J index is 115
# ASQTJ index is 137
# AH2OJ index is 159
# ACLJ index is 165

	
#	species_list = list(diff_map.keys())
# 	draw_pdf(orig_directory, species_list, diff_matrix)
	
	
main()















	
	
	
	