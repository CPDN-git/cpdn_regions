###########################################################################
# Program:  blank_map.py
# Author:   Sarah Sparrow
# Date:     21/09/2013
# Purpose:  To produce a blank global map for use
###########################################################################
import sys
import os
import numpy as np
import math
import datetime
import fnmatch
import matplotlib
import glob
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

from scipy.io import netcdf

ddir='data/'
region_files=['eu25_region.nc']
region_cols=['SpringGreen']

def get_rot_global_coords(region_file):
    f=netcdf.netcdf_file(ddir+region_file,'r')
    glat=f.variables['global_latitude0']
    glon=f.variables['global_longitude0']
    f.close()
    return glat,glon   

def get_global_coords(region_file):
    f=netcdf.netcdf_file(ddir+region_file,'r')
    lat=f.variables['latitude0']
    lon=f.variables['longitude0']
    f.close()
    return lat,lon
 
def map_plot():
    # Plot a map of the seasonal ocean climatology
    # Set the plot font size
    font = {'family' : 'sans-serif',
            'size'   : 14}
    matplotlib.rc('font', **font)

    # Produce the map plot
    fig = plt.figure()
    m = Basemap(llcrnrlon=-15.5,llcrnrlat=10.,urcrnrlon=118.566,urcrnrlat=75.352,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=38.,lon_0=25.)
#    m = Basemap(projection='ortho',lat_0=35,lon_0=-100,resolution='c')
    #m = Basemap(projection='robin',lon_0=35,resolution='c')
    m.bluemarble()
    m.drawmapboundary()
    for ir,region_file in enumerate(region_files):
	print region_file
	if region_file=="afr50_region.nc" or region_file=="dub25_region.nc" or region_file=="cafr25_region.nc":
	    lat,lon=get_global_coords(region_file)
            for i in range(0,4):
               if i==0:
                   longitude=np.repeat(lon[0],lat.shape)
		   latitude=lat[:]
               elif i==1:
                   longitude=lon[:]
                   latitude=np.repeat(lat[0],lon.shape)
               elif i==2:
                   longitude=np.repeat(lon[-1],lat.shape)
                   latitude=lat[:]
               elif i==3:
                   longitude=lon[:]
                   latitude=np.repeat(lat[-1],lon.shape)
               x, y = m(longitude, latitude)
               plt.scatter(x,y,lw=1,s=2, c=region_cols[ir],edgecolor=region_cols[ir],zorder=1)	
            for il in range(1,8):
               for i in range(0,4):
                if i==0:
                   longitude=np.repeat(lon[0+il],lat.shape)
                   latitude=lat[:]
                elif i==1:
                   longitude=lon[:]
                   latitude=np.repeat(lat[0+il],lon.shape)
                elif i==2:
                   longitude=np.repeat(lon[-1-il],lat.shape)
                   latitude=lat[:]
                elif i==3:
                   longitude=lon[:]
                   latitude=np.repeat(lat[-1-il],lon.shape)
                x, y = m(longitude, latitude)
                plt.scatter(x,y,lw=1,s=2, c=region_cols[ir],edgecolor=region_cols[ir],linewidths=0.0,alpha=0.1,zorder=1)
	else:
            glat,glon=get_rot_global_coords(region_file)
            for i in range(0,4):
    	       if i==0:
	           longitude=glon[0,:]
	           latitude=glat[0,:]
	       elif i==1:
	           longitude=glon[:,0]
                   latitude=glat[:,0]
               elif i==2:
	           longitude=glon[-1,:]
                   latitude=glat[-1,:]
               elif i==3:
                   longitude=glon[:,-1]
                   latitude=glat[:,-1]
           #Convert latitude and longitude to coordinates X and Y
	       x, y = m(longitude, latitude)
	       plt.plot(x,y,lw=1, c=region_cols[ir],zorder=2)
	    il=7
	    for i in range(0,4):
               if i==0:
	           longitude=glon[0,:]
		   latitude=glat[0,:]
	           longitude1=glon[0+il,:]
                   latitude1=glat[0+il,:] 
               elif i==1:
	           longitude=glon[:,0]
		   latitude=glat[:,0]
                   longitude1=glon[:,0+il]
                   latitude1=glat[:,0+il]
               elif i==2:
		   longitude=glon[-1,:]
		   latitude=glat[-1,:]
		   longitude1=glon[-1-il,:]
                   latitude1=glat[-1-il,:]
               elif i==3:
	           longitude=glon[:,-1]
		   latitude=glat[:,-1]
                   longitude1=glon[:,-1-il]
                   latitude1=glat[:,-1-il]
                #Convert latitude and longitude to coordinates X and Y
               print longitude,longitude1
	       x, y = m(longitude, latitude)
	       x1, y1 = m(longitude1, latitude1)
 	       plt.fill_between(x1,y,y1,facecolor=region_cols[ir],alpha=0.5,zorder=1)
	       plt.fill_betweenx(y1,x,x1,facecolor=region_cols[ir],alpha=0.5,zorder=1)

    plt.tight_layout()
    #plt.title("weather@home Regions")
    fig.savefig("region_plot_eu25.png")

#Main controling function
def main():
    map_plot()
    print 'Finished!'

#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()
