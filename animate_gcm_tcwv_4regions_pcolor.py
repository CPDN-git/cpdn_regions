# script to animate the 4 region figure (with boundary, colorbar, and dates etc.)
# Author : Sihan Li
# Some time Sep 2017
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animat
import glob
from mpl_toolkits.basemap import *
import matplotlib.pyplot as plt
from scipy.io import netcdf
from netCDF4 import *
from datetime import date, timedelta
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import datetime, time
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset as NetCDFFile, date2index, num2date
import matplotlib.animation as animation
ddir='batch_d811_nam50/atmos/' 
ddir1= 'batch_d811_nam50/region/' 
def get_rot_global_coords(in_file):
    f=netcdf.netcdf_file(ddir1+in_file,'r')
    glat=f.variables['global_latitude2'][:]
    glon=f.variables['global_longitude2'][:]
    f.close()
    return glat,glon  
 
def get_rot_global_coords1(in_file):
    f=netcdf.netcdf_file(ddir1+in_file,'r')
    glat1=f.variables['global_latitude3'][:]
    glon1=f.variables['global_longitude3'][:]
    f.close()
    return glat1,glon1
def get_global_coords(in_file):
    f=netcdf.netcdf_file(ddir+in_file,'r')
    lat=f.variables['latitude3'][:]
    lon=f.variables['longitude3'][:]
    f.close()
    return lat,lon
def get_global_wind_coords(in_file):
    f=netcdf.netcdf_file(ddir+in_file,'r')
    lat1=f.variables['latitude2'][:]
    lon1=f.variables['longitude2'][:]
    f.close()
    return lat1,lon1
def plot_global_winds(nt,fig,Varin,lon,lat,rglon1,rglat1,rVarin,srglon1,srglat1,sVarin,anzrglon1,anzrglat1,anzVarin,sasrglon1,sasrglat1,sasVarin,datelabel,fname_out=False):
	print nt
	fig = plt.figure(figsize=(8,4.5))
	ax = fig.add_axes([0.1,0.1,0.8,0.8])
	pos = ax.get_position()
	l, b, w, h = pos.bounds
	map = Basemap(projection='cyl',llcrnrlat= -90.,urcrnrlat= 90.,\
              resolution='c',  llcrnrlon=0.,urcrnrlon=360.)
	map.drawcoastlines()
	map.drawcountries()
	map.drawcoastlines(color = '0.15')
	map.drawparallels(np.arange( -90., 91.,30.),labels=[1,0,0,0],fontsize=10,linewidth=0.001, zorder = 0)
	map.drawmeridians(np.arange(0.,360.,30.),labels=[0,0,0,1],fontsize=10,linewidth=0.001, zorder = 0)
	clevs = np.arange(0,80,5)
	x1,y1=map(*np.meshgrid(lon,lat))
	CS2 = map.pcolor(x1,y1,Varin,cmap=plt.cm.PuBuGn, vmin=0, vmax=80)
	cax = plt.axes([l+w, b+0.01, 0.01, h-0.03]) # setup colorbar axes
	cbar=plt.colorbar(CS2,drawedges=True, cax=cax) # draw colorbar
	cax.text(0.0,-0.05,'mm')
	cbar.set_ticks([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80])
	cbar.set_ticklabels(['0','','10','','20','','30','','40','','50','','60','','70','','80'])
	cbar.update_ticks()
	plt.axes(ax) # reset current axes
	for i in range(0,4):
    		if i==0:
	    		longitude=rglon1[0,:]
	    		latitude=rglat1[0,:]
#		   #print glon[0,:]
    		elif i==1:
	    		longitude=rglon1[:,0]
            		latitude=rglat1[:,0]
		   #print glon[:,0]
    		elif i==2:
	    		longitude=rglon1[-1,:]
            		latitude=rglat1[-1,:]
		   #print glon[-1,:]
    		elif i==3:
            		longitude=rglon1[:,-1]
            		latitude=rglat1[:,-1]
		   #print glon[:,-1]
           #Convert latitude and longitude to coordinates X and Y
    		x, y = map(longitude, latitude)
    		plt.plot(x,y,lw=3, c='Plum',zorder=2)
    #sam50
    	for i in range(0,4):
		if i==0:
	           longitude=srglon1[0,:]
	           latitude=srglat1[0,:]
	    	elif i==1:
	           longitude=srglon1[:,0]
                   latitude=srglat1[:,0]
        	elif i==2:
	           longitude=srglon1[-1,:]
                   latitude=srglat1[-1,:]
        	elif i==3:
                   longitude=srglon1[:,-1]
                   latitude=srglat1[:,-1]
		
           #Convert latitude and longitude to coordinates X and Y
	    	x, y = map(longitude, latitude)
    		plt.plot(x,y,lw=3, c='LemonChiffon',zorder=3)
    #anz50
    	for i in range(0,4):
    		if i==0:
	           longitude=anzrglon1[0,:]
	           latitude=anzrglat1[0,:]
	    	elif i==1:
	           longitude=anzrglon1[:,0]
                   latitude=anzrglat1[:,0]
        	elif i==2:
	           longitude=anzrglon1[-1,:]
                   latitude=anzrglat1[-1,:]
        	elif i==3:
                   longitude=anzrglon1[:,-1]
                   latitude=anzrglat1[:,-1]
		
           #Convert latitude and longitude to coordinates X and Y
	    	x, y = map(longitude, latitude)
    		plt.plot(x,y,lw=3, c='Gold',zorder=4)
    #sas50
    	for i in range(0,4):
    		if i==0:
	           longitude=sasrglon1[0,:]
	           latitude=sasrglat1[0,:]
	    	elif i==1:
	           longitude=sasrglon1[:,0]
                   latitude=sasrglat1[:,0]
        	elif i==2:
	           longitude=sasrglon1[-1,:]
                   latitude=sasrglat1[-1,:]
        	elif i==3:
                   longitude=sasrglon1[:,-1]
                   latitude=sasrglat1[:,-1]
		
           #Convert latitude and longitude to coordinates X and Y
	    	x, y = map(longitude, latitude)
    		plt.plot(x,y,lw=3, c='BlueViolet',zorder=5)
	map.drawcoastlines(color = '0.15')
	clevs = np.arange(0,80,5)
	x2,y2 = map(rglon1[:],rglat1[:])

	CS2 = map.pcolor(x2,y2,rVarin,cmap=plt.cm.Greens, vmin=0, vmax=80)
	x3,y3 = map(srglon1[:],srglat1[:])
	CS2 = map.pcolor(x3,y3,sVarin,cmap=plt.cm.Greens,vmin=0,vmax=80)
	x4,y4 = map(anzrglon1[:],anzrglat1[:])
	CS2 = map.pcolor(x4,y4,anzVarin,cmap=plt.cm.Greens,vmin=0,vmax=80)
	x5,y5 = map(sasrglon1[:],sasrglat1[:])
	CS2 = map.pcolor(x5,y5,sasVarin,cmap=plt.cm.Greens,vmin=0,vmax=80)
	map.drawcoastlines(color = '0.15',linewidth=2)
	txt = plt.title('Total Column Water Vapor ' + str(datelabel), fontsize=14, fontweight='bold' )
	cax = plt.axes([l+w+0.05, b+0.01, 0.01, h-0.03]) # setup colorbar axes
	cbar=plt.colorbar(CS2,drawedges=True, cax=cax) # draw colorbar
	cax.text(0.0,-0.05,'mm')
	cbar.set_ticks([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80])
	cbar.set_ticklabels(['0','','10','','20','','30','','40','','50','','60','','70','','80'])
	cbar.update_ticks()
	plt.axes(ax) # reset current axes
	if fname_out:
    		fig.savefig(fname_out,dpi=100)


def main():
	ddir='batch_d811_nam50/atmos/'  
	in_file1='item3225_daily_mean_a000_2014-12_2015-12.nc'
	in_file2='item3226_daily_mean_a000_2014-12_2015-12.nc'
	in_file3='item10_daily_mean_a000_2014-12_2015-12.nc'
	in_file4='item1_daily_mean_a000_2014-12_2015-12.nc'
	lat1,lon1=get_global_wind_coords(in_file1)
	f=netcdf.netcdf_file(ddir+in_file1,'r')
	u=f.variables['item3225_daily_mean'][:]
	f.close()
	f=netcdf.netcdf_file(ddir+in_file2,'r')
	v=f.variables['item3226_daily_mean'][:]
	f.close()
	lat,lon=get_global_coords(in_file3)
	nc = Dataset(ddir+in_file3, mode='r')
	shum = nc.variables['item10_daily_mean'][:]
	lat = nc.variables['latitude3']
	lon = nc.variables['longitude3']
	nc = Dataset(ddir+in_file4, mode='r')
	surfp = nc.variables['item1_daily_mean'][:]
	tcwv=shum*surfp/9.81
	U=u[:,0,:,:]
	V=v[:,0,:,:]
	U=np.squeeze(U)
	V=np.squeeze(V)
	ddir1='batch_d811_nam50/region/'  
	in_file1='nam50_item3225_daily_mean_a000_2014-12_2015-12.nc'
	in_file2='nam50_item3226_daily_mean_a000_2014-12_2015-12.nc'
	in_file3='nam50_item10_daily_mean_a000_2014-12_2015-12.nc'
	in_file4='nam50_item1_daily_mean_a000_2014-12_2015-12.nc'
	rglat,rglon=get_rot_global_coords(in_file1)
	nc = Dataset(ddir1+in_file1, mode='r')
	ru = nc.variables['item3225_daily_mean']
	rlat = nc.variables['latitude2']
	rlon = nc.variables['longitude2']
	nc = Dataset(ddir1+in_file2, mode='r')
	rv = nc.variables['item3226_daily_mean']
	rglat1,rglon1=get_rot_global_coords1(in_file3)
	nc = Dataset(ddir1+in_file3, mode='r')
	rshum = nc.variables['item10_daily_mean'][:]
	rlat1 = nc.variables['latitude3']
	rlon1 = nc.variables['longitude3']
	nc = Dataset(ddir1+in_file4, mode='r')
	rsurfp = nc.variables['item1_daily_mean'][:]
	rtcwv=rshum*rsurfp/9.81
	ddir1='batch_d811_nam50/region/'  
	in_file1='sam50_item3225_daily_mean_a000_2014-12_2015-12.nc'
	in_file2='sam50_item3226_daily_mean_a000_2014-12_2015-12.nc'
	in_file3='sam50_item10_daily_mean_a000_2014-12_2015-12.nc'
	in_file4='sam50_item1_daily_mean_a000_2014-12_2015-12.nc'
	srglat,srglon=get_rot_global_coords(in_file1)
	nc = Dataset(ddir1+in_file1, mode='r')
	sru = nc.variables['item3225_daily_mean']
	srlat = nc.variables['latitude2']
	srlon = nc.variables['longitude2']
	nc = Dataset(ddir1+in_file2, mode='r')
	srv = nc.variables['item3226_daily_mean']
	srglat1,srglon1=get_rot_global_coords1(in_file3)
	nc = Dataset(ddir1+in_file3, mode='r')
	srshum = nc.variables['item10_daily_mean'][:]
	srlat1 = nc.variables['latitude3']
	srlon1 = nc.variables['longitude3']
	nc = Dataset(ddir1+in_file4, mode='r')
	srsurfp = nc.variables['item1_daily_mean'][:]
	srtcwv=srshum*srsurfp/9.81
	# read in anz50 data
	ddir1='batch_d811_nam50/region/'  
	in_file1='anz50_item3225_daily_mean_a000_2014-12_2015-12.nc'
	in_file2='anz50_item3226_daily_mean_a000_2014-12_2015-12.nc'
	in_file3='anz50_item10_daily_mean_a000_2014-12_2015-12.nc'
	in_file4='anz50_item1_daily_mean_a000_2014-12_2015-12.nc'
	anzrglat,anzrglon=get_rot_global_coords(in_file1)
	nc = Dataset(ddir1+in_file1, mode='r')
	anzru = nc.variables['item3225_daily_mean']
	anzrlat = nc.variables['latitude2']
	anzrlon = nc.variables['longitude2']
	nc = Dataset(ddir1+in_file2, mode='r')
	anzrv = nc.variables['item3226_daily_mean']
	anzrglat1,anzrglon1=get_rot_global_coords1(in_file3)
	nc = Dataset(ddir1+in_file3, mode='r')
	anzrshum = nc.variables['item10_daily_mean'][:]
	anzrlat1 = nc.variables['latitude3']
	anzrlon1 = nc.variables['longitude3']
	nc = Dataset(ddir1+in_file4, mode='r')
	anzrsurfp = nc.variables['item1_daily_mean'][:]
	anzrtcwv=anzrshum*anzrsurfp/9.81
	# read in sas50 data
	ddir1='batch_d811_nam50/region/'  
	in_file1='sas50_item3225_daily_mean_a000_2014-12_2015-12.nc'
	in_file2='sas50_item3226_daily_mean_a000_2014-12_2015-12.nc'
	in_file3='sas50_item10_daily_mean_a000_2014-12_2015-12.nc'
	in_file4='sas50_item1_daily_mean_a000_2014-12_2015-12.nc'
	sasrglat,sasrglon=get_rot_global_coords(in_file1)
	nc = Dataset(ddir1+in_file1, mode='r')
	sasru = nc.variables['item3225_daily_mean']
	sasrlat = nc.variables['latitude2']
	sasrlon = nc.variables['longitude2']
	nc = Dataset(ddir1+in_file2, mode='r')
	sasrv = nc.variables['item3226_daily_mean']
	sasrglat1,sasrglon1=get_rot_global_coords1(in_file3)
	nc = Dataset(ddir1+in_file3, mode='r')
	sasrshum = nc.variables['item10_daily_mean'][:]
	sasrlat1 = nc.variables['latitude3']
	sasrlon1 = nc.variables['longitude3']
	nc = Dataset(ddir1+in_file4, mode='r')
	sasrsurfp = nc.variables['item1_daily_mean'][:]
	sasrtcwv=sasrshum*sasrsurfp/9.81
	# read in sas50 data
	ddir1='batch_d811_nam50/region/'  
	in_file1='sas50_item3225_daily_mean_a000_2014-12_2015-12.nc'
	in_file2='sas50_item3226_daily_mean_a000_2014-12_2015-12.nc'
	in_file3='sas50_item10_daily_mean_a000_2014-12_2015-12.nc'
	in_file4='sas50_item1_daily_mean_a000_2014-12_2015-12.nc'
	sasrglat,sasrglon=get_rot_global_coords(in_file1)
	nc = Dataset(ddir1+in_file1, mode='r')
	sasru = nc.variables['item3225_daily_mean']
	sasrlat = nc.variables['latitude2']
	sasrlon = nc.variables['longitude2']
	nc = Dataset(ddir1+in_file2, mode='r')
	sasrv = nc.variables['item3226_daily_mean']
	sasrglat1,sasrglon1=get_rot_global_coords1(in_file3)
	nc = Dataset(ddir1+in_file3, mode='r')
	sasrshum = nc.variables['item10_daily_mean'][:]
	sasrlat1 = nc.variables['latitude3']
	sasrlon1 = nc.variables['longitude3']
	nc = Dataset(ddir1+in_file4, mode='r')
	sasrsurfp = nc.variables['item1_daily_mean'][:]
	sasrtcwv=sasrshum*sasrsurfp/9.81
	
	d1 = date(2014, 12, 01)
	d2 = date(2015, 12, 30)
	delta = d2 - d1         # timedelta
	print delta.days
	datelist=[]
	for i in range(delta.days + 1):
    #print(d1 + timedelta(days=i))
     		dates=d1 + timedelta(days=i)
        	#print dates
        	if dates.day==31:
        		print 'haha' 
        	else:
        		datelist.append(dates)
        		if dates==date(2015, 2, 28):
        			print dates
        			datelist.append('2015-02-29')
        			datelist.append('2015-02-30')
    			
	for nt in range(0,90):
		if nt<10:
    	    		fout='batch_d811_nam50/atmos/atmos_4regions_tcwv_pcolor_00'+str(nt)
        	elif nt<100:
            		fout='batch_d811_nam50/atmos/atmos_4regions_tcwv_pcolor_0'+str(nt)
        	else:
            		fout='batch_d811_nam50/atmos/atmos_4regions_tcwv_pcolor_'+str(nt)
		f_ext='.png'
		Varin=np.squeeze(tcwv[nt,0,:,:])
		rVarin=np.squeeze(rtcwv[nt,0,:,:])
		sVarin=np.squeeze(srtcwv[nt,0,:,:])
		anzVarin=np.squeeze(anzrtcwv[nt,0,:,:])
		sasVarin=np.squeeze(sasrtcwv[nt,0,:,:])
		datelabel=datelist[nt]
		print datelabel
		plot_global_winds(nt,plt.figure(),Varin,lon,lat,rglon1,rglat1,rVarin,srglon1,srglat1,sVarin,anzrglon1,anzrglat1,anzVarin,sasrglon1,sasrglat1,sasVarin,datelabel,fname_out=fout+f_ext)
		print "saved as",fout+f_ext
		
	print 'Finished!'
	
#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()
