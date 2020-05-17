#!/usr/bin/env python2.7
#
###########################################################################
# Program:  region_map.py
# Author:   Sarah Sparrow
# Date:     03/01/2017
# Purpose:  To plot weather@home regions
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
import csv

from scipy.io import netcdf


# Dictionary detailing region names, colour and whether to display the outline as dotted
region_dict={'afr50':('SaddleBrown',False),'nawa25':('HotPink',False),'anz50':('Gold',False),'eas50':('Red',False),'eu25':('SpringGreen',False),'eu50r':('RoyalBlue',False),'cam50':('DarkOrange',False),'cam25':('DeepSkyBlue',False),'pnw25':('ForestGreen',False),'sas50':('BlueViolet',False),'wus25':('YellowGreen',False),'sam50':('LemonChiffon',False),'cafr25':('SandyBrown',False),'nam50':('Plum',False),'sam25':('Olive',False),'safr50':('Teal',False),'cari25':('DeepPink',False)}

def get_rot_global_coords(region_file):
    f=netcdf.netcdf_file('data/'+region_file,'r')
    glat=f.variables['global_latitude0']
    glon=f.variables['global_longitude0']
    f.close()
    return glat,glon   

def get_global_coords(region_file):
    f=netcdf.netcdf_file('data/'+region_file,'r')
    lat=f.variables['latitude0']
    lon=f.variables['longitude0']
    f.close()
    return lat,lon
 
def map_plot(plot_regions,plot_collabLines,plot_collab,plot_volunteers,idx):
    # Plot a map of the seasonal ocean climatology
    # Set the plot font size
    font = {'family' : 'sans-serif',
            'size'   : 14}
    matplotlib.rc('font', **font)

    # Produce the map plot
    fig = plt.figure()
    m = Basemap(projection='robin',lon_0=35,resolution='c')
    m.bluemarble()
    m.drawmapboundary()
    if plot_regions:
    	for region,region_style in region_dict.iteritems():
    		region_file=region+"_region.nc"
		region_col=region_style[0]
		region_line=region_style[1]
		print region_file
		if region=="afr50" or region=="nawa25" or region=="cafr25" or region=="safr50":
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
	       			if region_line==True:
	           			plt.scatter(x[::5],y[::5],lw=1,s=2, c=region_col,edgecolor=region_col,zorder=2,alpha=0.3)
               			else:
	           			plt.scatter(x,y,lw=1,s=2, c=region_col,edgecolor=region_col,zorder=2,alpha=0.3)	
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
	       			if region_line==True:
					plt.scatter(x[::5],y[::5],lw=1,s=2, c=region_col,edgecolor=region_col,zorder=2,alpha=0.3)
	       			else:
    	       				plt.scatter(x,y,lw=1,s=2, c=region_col,edgecolor=region_col,zorder=2,alpha=0.3)
    if plot_volunteers:    

    	with open('locations.csv') as csvfile:
  		readCSV = csv.reader(csvfile, delimiter=',')
    		clons = []
    		clats = []
    		for row in readCSV:
  			clon = float(row[1])
        		clat = float(row[0])

        		clons.append(clon)
        		clats.append(clat)

    	x,y = m(clons, clats)
    	m.plot(x,y,'o',c='Yellow',markersize=1,zorder=1)
    


    if plot_collab:
    	clons1=[360.-1.2577,144.9631,174.7762,360.-123.2620,36.8219,38.7578,73.8567,18.4241,360.-99.1332,118.7969,129.3435,90.4125,8.5417,360.-2.5879,360.-3.5339,360.-3.1883,13.0645,360-46.6333]
    	clats1=[51.752,-37.8136,-41.2865,44.5646,-1.2921, 8.9806,18.5204,-33.9249,19.4326,32.0603,36.0190,23.8103,47.3769,51.4545, 50.7184,55.9533,52.3906,-23.5505]

    	x,y = m(clons1, clats1)
    	m.plot(x,y,'wo',markersize=5)
    	if plot_collabLines:
        	for i in range(0,len(clons1)):
                	m.drawgreatcircle(360.-1.2577, 51.752, clons1[i], clats1[i], del_s=100.0,c='White',lw=1)


#    clons1=[360.-1.2577,144.9631,174.7762,360.-123.2620,36.8219,38.7578,73.8567,18.4241,360.-99.1332,118.7969,129.3435,90.4125,8.5417,360.-2.5879,360.-3.5339,360.-3.1883,13.0645]
#    clats1=[51.752,-37.8136,-41.2865,44.5646,-1.2921, 8.9806,18.5204,-33.9249,19.4326,32.0603,36.0190,23.8103,47.3769,51.4545, 50.7184,55.9533,52.3906]
#
#    x,y = m(clons1, clats1)
#    m.plot(x,y,'wo',markersize=5,zorder=2)
#    if plot_collabLines:
#    	m.drawgreatcircle(360.-1.2577, 51.752, 144.9631, -37.8136, del_s=100.0,c='White',lw=1) #Melbourne
#    	m.drawgreatcircle(360.-1.2577, 51.752, 174.7762, -41.2865, del_s=100.0,c='White',lw=1) #Wellington
#    	m.drawgreatcircle(360.-1.2577, 51.752, 360.-123.2620, 44.5646, del_s=100.0,c='White',lw=1) # Corvallis
#    	m.drawgreatcircle(360.-1.2577, 51.752, 36.8219, -1.2921, del_s=100.0,c='White',lw=1) # Nairobi
#    	m.drawgreatcircle(360.-1.2577, 51.752, 38.7578, 8.9806, del_s=100.0,c='White',lw=1) #Addis Ababa
#    	m.drawgreatcircle(360.-1.2577, 51.752, 73.8567, 18.5204, del_s=100.0,c='White',lw=1) #Pune
#    	m.drawgreatcircle(360.-1.2577, 51.752, 18.4241, -33.9249, del_s=100.0,c='White',lw=1) #Cape Town
#    	m.drawgreatcircle(360.-1.2577, 51.752, 360.-99.1332, 19.4326, del_s=100.0,c='White',lw=1) #Mexico City
#    	m.drawgreatcircle(360.-1.2577, 51.752, 118.7969, 32.0603, del_s=100.0,c='White',lw=1) # Nanjing
#    	m.drawgreatcircle(360.-1.2577, 51.752, 129.3435, 36.0190, del_s=100.0,c='White',lw=1) #Pohang
#    	m.drawgreatcircle(360.-1.2577, 51.752, 90.4125, 23.8103, del_s=100.0,c='White',lw=1) #Dhaka
#    	#m.drawgreatcircle(360.-1.2577, 51.752, 5.1810, 52.1093, del_s=100.0,c='White',lw=1) #De Bilt (KNMI)
#    	m.drawgreatcircle(360.-1.2577, 51.752, 8.5417, 47.3769, del_s=100.0,c='White',lw=1) #Zurich
#    	#m.drawgreatcircle(360.-1.2577, 51.752, 10.7522, 59.9139, del_s=100.0,c='White',lw=1) #Oslo
#    	m.drawgreatcircle(360.-1.2577, 51.752, 360.-2.5879, 51.4545, del_s=100.0,c='White',lw=1) #Bristol
#    	m.drawgreatcircle(360.-1.2577, 51.752, 360.-3.5339, 50.7184, del_s=100.0,c='White',lw=1) #Exeter
#    	m.drawgreatcircle(360.-1.2577, 51.752, 360.-3.1883, 55.9533, del_s=100.0,c='White',lw=1) #Edinburgh
#    	m.drawgreatcircle(360.-1.2577, 51.752, 13.0645, 52.3906, del_s=100.0,c='White',lw=1) #Potsdam

    plt.tight_layout()
    #plt.title("weather@home Regions")
    print "Saving region_network_plot_"+idx+".png"
    fig.savefig("region_network_plot_"+idx+".png",dpi=80)

#Main controling function
def main():
    plot_collab=True
    plot_collabLines=True
    plot_regions=False
    plot_volunteers=False
    map_plot(plot_regions,plot_collabLines,plot_collab,plot_volunteers,"01")

    plot_collab=True
    plot_collabLines=False
    plot_regions=True
    plot_volunteers=False
    map_plot(plot_regions,plot_collabLines,plot_collab,plot_volunteers,"02")

    plot_collab=False
    plot_collabLines=False
    plot_regions=True
    plot_volunteers=True
    map_plot(plot_regions,plot_collabLines,plot_collab,plot_volunteers,"03")

    plot_collab=False
    plot_collabLines=False
    plot_regions=False
    plot_volunteers=True
    map_plot(plot_regions,plot_collabLines,plot_collab,plot_volunteers,"04")

    plot_collab=True
    plot_collabLines=True
    plot_regions=False
    plot_volunteers=True
    map_plot(plot_regions,plot_collabLines,plot_collab,plot_volunteers,"05")
    print 'Finished!'

#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()
