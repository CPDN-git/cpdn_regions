#############
# Author: Sihan Li
# Date: 13/09/2017
# Region map with spong zone
##############
import sys
import os
import numpy as np
import math
import datetime
import fnmatch
import matplotlib
from matplotlib.patches import Polygon
import glob
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from scipy.io import netcdf

ddir='data/'
region_files=['eas50_region.nc']
region_cols=['Red']

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

def get_rotated_pole(region_file,field16):
	# get the rotated pole longitude / latitude (for calculating weights)
	try:
		grid_map_name = getattr(nc_in_var,"grid_mapping")
		grid_map_var = nc_in_file.variables[grid_map_name]	
		plon = getattr(grid_map_var,"grid_north_pole_longitude")
		plat = getattr(grid_map_var,"grid_north_pole_latitude")
	except:
		plon = 0.0
		plat = 90.0
	return plon, plat

def draw_box( lat,lon, m, col):
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
               plt.plot(x,y,lw=1, c=col,zorder=1)

def fill_box(lat,lon, m, col):

    lons=np.repeat(lon[0],lat.shape)
    lats=lat[:]

    print len(lats),len(lons)

    lons=np.append(lons,lon[:])
    lats=np.append(lats,np.repeat(lat[-1],lon.shape))

    print len(lats),len(lons)

    lons=np.append(lons,np.repeat(lon[-1],lat.shape))
    lats=np.append(lats,lat[::-1])

    print len(lats),len(lons)

    lons=np.append(lons,lon[::-1])
    lats=np.append(lats,np.repeat(lat[0],lon.shape))

    print len(lats),len(lons)

    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, facecolor=col, alpha=0.4 )
    plt.gca().add_patch(poly)


def map_plot():
    # Plot a map of the seasonal ocean climatology
    # Set the plot font size
	font = {'family' : 'sans-serif',
            'size'   : 14}
	matplotlib.rc('font', **font)
    
fig = plt.figure()
	
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
            polelat,polelon=get_rotated_pole(region_file,'field16')
            ny,nx=glat.shape
            print glon[-1,0]
            print glat[-1,0]
            print glat[0,-1]
            print glon[0,-1]
            
            #m = Basemap(llcrnrlon=74,llcrnrlat=-28,urcrnrlon=190,urcrnrlat=57,\
            #			rsphere=(6378137.00,6356752.3142),\
            #            resolution='l',area_thresh=1000.,projection='lcc',\
            #            lat_1=20.,lon_0=125.)
            m = Basemap(llcrnrlon=glon[-1,0]-10,llcrnrlat=glat[-1,0]-20,urcrnrlon=glon[0,-1]+30,urcrnrlat=glat[0,-1]+10,\
            			rsphere=(6378137.00,6356752.3142),\
                        resolution='l',area_thresh=1000.,projection='lcc',\
                        #lat_1=-10.,lon_0=135)
                        lat_1=(glat[-1,0]+glat[0,-1])/2,lon_0=(glon[-1,0]+glon[0,-1])/2)
            m.bluemarble()
            m.drawmapboundary()
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
               x, y = m(longitude, latitude)
               x1, y1 = m(longitude1, latitude1)

               # plt.scatter(x,y,lw=0,s=0.9, c=region_cols[ir],edgecolor=None,linewidths=0.0,alpha=0.3,zorder=1)
               plt.fill_between(x1,y,y1,facecolor=region_cols[ir],alpha=0.5,zorder=1)
               plt.fill_betweenx(y1,x,x1,facecolor=region_cols[ir],alpha=0.5,zorder=1)


	    #Add in the China boxes
	    box2017=[25.,37.5,106.,122.]
	    lats = np.linspace(box2017[0], box2017[1])
	    lons = np.linspace(box2017[2], box2017[3])
	    draw_box( lats,lons, m, "Gold" )
	    fill_box( lats,lons, m, "Gold" )
#	    CC_box=[25,35,110,117]
#	    EC_box=[25,35,117,122]
	  
#	    lats = np.linspace(CC_box[0], CC_box[1])
#	    lons = np.linspace(CC_box[2], CC_box[3])
#	    draw_box( lats,lons, m, "Gold" )
   
#	    lats = np.linspace(EC_box[0], EC_box[1])
#        lons = np.linspace(EC_box[2], EC_box[3])
#        draw_box( lats,lons, m, "SpringGreen" )

plt.tight_layout()
fig.savefig("region_plot_eas50_v3.png")

#Main controling function
def main():
    map_plot()
    print 'Finished!'

#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()
