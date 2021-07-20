#!/home/cl3293/anaconda3/bin/python

### Script for plotting and evaluating the ACONC in CMAQ v5.3 beta copied from Dan's directory

### CHANGE above: specify local python installation

import matplotlib as mpl
mpl.use('Agg')
from netCDF4 import Dataset
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from math import fsum, exp
from random import uniform
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import matplotlib.ticker
from scipy import stats

#236 variables, below are species related to NH3
lst_var = ['NH3', 'HNO3', 'HCL', 'ASO4J', 'ASO4I', 'ANH4J', 'ANH4I', 'ANO3J', 'ANO3I', 'ANAJ', 'ANAI', 'ACLJ', 'ACLI']

# Open the DS, *_off is the official result
ACONCname = "/glade/work/edliu/models/CMAQ_TEST/data/output_CCTM_v532_intel_Bench_2016_12SE1/CCTM_ACONC_v532_intel_Bench_2016_12SE1_20160701.nc"
ACONCname_off = "/glade/work/edliu/models/inputs/SEv5.3.2.BENCH/output_CCTM_v532_ISAM_gcc_Bench_2016_12SE1_opt/CCTM_ACONC_v532_ISAM_gcc_Bench_2016_12SE1_20160701.nc"
ACONCds = Dataset(ACONCname, "r")
ACONCds_off = Dataset(ACONCname_off, "r")

latlonfilename = '/glade/work/edliu/models/CMAQ_TEST/CCTM/scripts/latlon.1day'
latlonfile = Dataset(latlonfilename, mode='r', open=True)

gridLat =  np.squeeze(latlonfile.variables['LAT'][:][:][:][:])
gridLon =  np.squeeze(latlonfile.variables['LON'][:][:][:][:])

last1ind,last2ind = gridLon.shape

loncrn = [gridLon[0,0],gridLon[last1ind-1,last2ind-1]]
latcrn = [gridLat[0,0],gridLat[last1ind-1,last2ind-1]]
lonmid = gridLon[last1ind//2,last2ind//2]
latmid = gridLat[last1ind//2-1,last2ind//2-1]

for var in lst_var:
    ACONCds_var = np.squeeze(ACONCds.variables[var])             # 24 * 80 * 100
    ACONCds_off_var = np.squeeze(ACONCds_off.variables[var])     # 24 * 80 * 100
    ACONds_diff = ACONCds_var - ACONCds_off_var                  # 24 * 80 * 100

    #ACONC scatter plot
    unit = 'ppmV'
    pdf_file = '/glade/work/edliu/models/analysis/ed_vs_official' + var + '.pdf'
    with PdfPages(pdf_file) as pdf:
        nArray = ACONCds_var.flatten() #x
        fArray = ACONCds_off_var.flatten() #y
    
        # Calculate Statistical Parameters
        slope, intercept, r_value, p_value, std_err = stats.linregress(nArray, fArray)
        plt.rcParams["font.size"] = 14
    
        cs = plt.scatter(nArray,fArray, cmap=cm.Spectral_r)
    
        rangemin = min(np.min(nArray),np.min(fArray))
        rangemax = max(np.max(nArray),np.max(fArray))
        v = [rangemin, rangemax, rangemin, rangemax]
        plt.axis(v)
    
        plt.plot([rangemin, rangemax], [rangemin, rangemax])
        ax = plt.gca()
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.set_ylabel(r'ACONC official (' + unit + ')')
        ax.set_xlabel(r'ACONC el662 (' + unit + ')')
        ax.minorticks_on()
        ax.yaxis.grid(True, which='major', color='black', linestyle='--', linewidth='0.4')
        ax.yaxis.grid(True, which='minor', color='gray', linestyle=':', linewidth='0.3',)
        ax.xaxis.grid(True, which='major', color='black', linestyle='--', linewidth='0.4')
        ax.xaxis.grid(True, which='minor', color='gray', linestyle=':', linewidth='0.3')
    
        ax.annotate('R$^2$ = {:0.2F}\nSlope = {:0.2F}\nIntercept = {:0.2E}\nBlue Line = 1:1'.format(r_value**2., slope, intercept), xy=(0.05, 0.7), fontsize=14, fontweight='bold', xycoords='axes fraction',                \
                        horizontalalignment='left', verticalalignment='bottom')
    
        plt.tight_layout()
    
        # add title
        plt.title( r'ACONC evaluation for '+var)
        pdf.savefig()
        plt.close()
        plt.hist2d(nArray, fArray, bins =(10,10), cmap=plt.cm.BuPu)
        plt.colorbar()
        pdf.savefig()
        plt.close()
        
   #ACONC
    with PdfPages('/glade/work/edliu/models/analysis/ed.'+var+'map.pdf') as pdf:
    
        # Hour range is the timeframe in the ingested CMAQ values to be plotted and compared
        ###CHANGE range for that date
        j = 0
        for i in range(0, 24):
            ###CHANGE date
            rng = pd.date_range('7/1/2016', periods=24, freq='H')
    
            # Put the monitoring data on top of the concentration from CMAQ
            # Establish basemap, this is copied from the GRIDDESC file
            m = Basemap(projection='lcc', llcrnrlon=loncrn[0], llcrnrlat=latcrn[0],\
                        urcrnrlon=loncrn[1], urcrnrlat=latcrn[1], lon_0=lonmid,\
                        lat_0=latmid, resolution='l')
    
            # draw coastlines, state and country boundaries, edge of map.
            m.drawcoastlines()
            m.drawstates()
            m.drawcountries()
    
            # Try plotting hourly ozone from CMAQ
            ny = ACONCds_var.shape[1]
            nx = ACONCds_var.shape[2]
            lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
            x, y = m(lons, lats) # compute map proj coordinates.
    
            # draw filled contours.
            clevs = np.arange(np.min(ACONCds_var[0:24, :, :]), np.max(ACONCds_var[0:24, :, :])+(np.max(ACONCds_var[0:24, :, :])-np.min(ACONCds_var[0:24, :, :]))/100, (np.max(ACONCds_var[0:24, :, :])-np.min(ACONCds_var[0:24, :, :]))/100)
            print("max", np.max(ACONCds_var[0:24, :, :]))
            print("min", np.min(ACONCds_var[0:24, :, :]))
            # shift the hours by a fixed amount to the start of the CMAQ day of interest
            cs = m.contourf(x, y, ACONCds_var[i, :, :], clevs, cmap=cm.Spectral_r, vmin=np.min(ACONCds_var), vmax=np.max(ACONCds_var))
    
            # add colorbar.
            cbar = m.colorbar(cs, location='right', pad="0.01%")
            cbar.set_label('(ppmV)')
            plt.tight_layout()
    
            # add title
            #CHANGE TITLE
            plt.title('Hourly '+var+' for '+str(rng[j]))
    
            # iterate over dates
            j = j + 1
    
            #plt.show()
    
            pdf.savefig()
            plt.close()
   #
    #ACONC_diff
    with PdfPages('/glade/work/edliu/models/analysis/ACONC_'+var+'.diff.ed.vs.official_map.pdf') as pdf:
        j = 0
        ###CHANGE range for that date
        for i in range(0,24):
            ###CHANGE date
            rng = pd.date_range('7/1/2016',periods=24,freq='H')

            # Establish basemap
            m = Basemap(projection='lcc', llcrnrlon=loncrn[0], llcrnrlat=latcrn[0],\
                        urcrnrlon=loncrn[1], urcrnrlat=latcrn[1], lon_0=lonmid,\
                        lat_0=latmid, resolution='l')

            # draw coastlines, state and country boundaries, edge of map.
            m.drawcoastlines()
            m.drawstates()
            m.drawcountries()

            ny = ACONds_diff.shape[1]; nx = ACONds_diff.shape[2]
            lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
            x, y = m(lons, lats) # compute map proj coordinates.

            # draw filled contours.
            upperbound = max(np.abs(np.nanmax(ACONds_diff[i,:,:])), np.abs(np.nanmin(ACONds_diff[i,:,:])))
            print("upperbound", upperbound)
            lowerbound = -upperbound
            clevs = np.arange(lowerbound,upperbound+(upperbound - lowerbound)/20,(upperbound - lowerbound)/20)
            cmap=cm.seismic

            # shift the hours by a fixed amount to the start of the CMAQ day of interest
            cs = m.contourf(x,y,ACONds_diff[i,:,:],clevs,cmap=cmap,extend='both')        # add colorbar.
            cbar = m.colorbar(cs, location='right',pad="0.00%")
            #scientific notation
            cbar.formatter.set_powerlimits((0, 0))
            cbar.update_ticks()
            # change tick font size
            cbar.ax.tick_params(labelsize=14)
            # add label
            cbar.set_label('(ppmV)', fontsize=14)


            # add title
            plt.title('ACONC '+var+'ed vs. official')

            # iterate over dates
            j = j+1

            plt.tight_layout()
            pdf.savefig()
            plt.close()