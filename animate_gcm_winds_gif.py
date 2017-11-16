# script to animate global 10-m winds, plain white vectors
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
def get_rot_global_coords(in_file):
    f=netcdf.netcdf_file(ddir+in_file,'r')
    glat=f.variables['global_latitude0']
    glon=f.variables['global_longitude0']
    f.close()
    return glat,glon   

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
    

def plot_global_winds(nt,fig,Uin,Vin,lon1,lat1,datelabel,fname_out=False):
	print nt
	fig = plt.figure(figsize=(8,4.5))
	ax = fig.add_axes([0.1,0.1,0.8,0.8])
	#nframe = 0; n1 = 390
	#nt = 0
	pos = ax.get_position()
	l, b, w, h = pos.bounds
	map = Basemap(projection='cyl',llcrnrlat= -90.,urcrnrlat= 90.,\
              resolution='c',  llcrnrlon=-180.,urcrnrlon=180.)
	map.bluemarble()
	map.drawcoastlines()
	map.drawcountries()
	map.drawcoastlines(color = '0.15')
	map.drawparallels(np.arange( -90., 91.,30.),labels=[1,0,0,0],fontsize=10,linewidth=0.001, zorder = 0)
	map.drawmeridians(np.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=10,linewidth=0.001, zorder = 0)
	clevs = np.arange(0,50,5)
	lonout, uu = map.shiftdata(lon1,Uin,lon_0=0)
	lonout, vv = map.shiftdata(lon1,Vin,lon_0=0)
	speed=np.sqrt(uu*uu+vv*vv)
	yy = np.arange(0, Y4.shape[0], 4)
	xx = np.arange(0, X4.shape[1], 4)
	points = np.meshgrid(yy, xx)
	vectplot=map.quiver(X4[points],Y4[points],uu[points],vv[points],color='lightgray')
	plt.quiverkey(vectplot, 0.1, -0.1, 10, '10 m/s', fontproperties={'weight': 'bold'}, labelpos='W',color='0.15')  #-- position and reference label
	txt = plt.title('10-m Wind ' + str(datelabel), fontsize=14, fontweight='bold' )

	if fname_out:
    		fig.savefig(fname_out,dpi=100)
def main():
	ddir='batch_d811_nam50/atmos/'  
	in_file1='item3225_daily_mean/item3225_daily_mean_a000_2014-12_2015-12.nc'
	in_file2='item3226_daily_mean/item3226_daily_mean_a000_2014-12_2015-12.nc'
	in_file3='item10_daily_mean/item10_daily_mean_a000_2014-12_2015-12.nc'
	in_file4='item1_daily_mean/item1_daily_mean_a000_2014-12_2015-12.nc'
	lat1,lon1=get_global_wind_coords(in_file1)
	f=netcdf.netcdf_file(ddir+in_file1,'r')
	u=f.variables['item3225_daily_mean'][:]
	f.close()
	f=netcdf.netcdf_file(ddir+in_file2,'r')
	v=f.variables['item3226_daily_mean'][:]
	f.close()
	U=u[:,0,:,:]
	V=v[:,0,:,:]
	U=np.squeeze(U)
	V=np.squeeze(V)
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
        			
	for nt in range(330,360):
		if nt<10:
    	    		fout='batch_d811_nam50/atmos/gcm_winds_white_00'+str(nt)
        	elif nt<100:
            		fout='batch_d811_nam50/atmos/gcm_winds_white_0'+str(nt)
        	else:
            		fout='batch_d811_nam50/atmos/gcm_winds_white_'+str(nt)
		f_ext='.png'
		Uin=np.squeeze(U[nt,:,:])
		Vin=np.squeeze(V[nt,:,:])
		datelabel=datelist[nt]
		print datelabel
		plot_global_winds(nt,plt.figure(),Uin,Vin,lon1,lat1,datelabel,fname_out=fout+f_ext)
		print "saved as",fout+f_ext
		
	print 'Finished!'
	
#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()
