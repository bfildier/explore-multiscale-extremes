import sys,os,glob
import psutil

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
from pprint import pprint

from matplotlib.colors import LogNorm
# from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import datetime as dt

import re
import gc
import warnings

thismodule = sys.modules[__name__]

#-- to redirect print output to standard output

class Unbuffered(object):

    def __init__(self, stream):
        self.stream = stream
    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
    def writelines(self, datas):
        self.stream.writelines(datas)
        self.stream.flush()
    def __getattr__(self, attr):
        return getattr(self.stream, attr)

sys.stdout = Unbuffered(sys.stdout)



# Global variables
from settings import *

# Own functions
from fcns_load_DYAMOND_SAM import *


def importData(i_t,varid,df):
    
    if varid == 'prec':     # Load DYAMOND precip data
        
        Prec_DYAMOND = loadPrec(i_t,df)
        
        return Prec_DYAMOND
    
    elif varid == 'PW':
        
        root_DYAMOND = df.iloc[i_t]['path_dyamond'] + '.%s.2D.nc'
        file_DYAMOND = root_DYAMOND%varid
        var_DYAMOND = xr.open_dataarray(os.path.join(DIR_DYAMOND,file_DYAMOND))[0]
        
        return var_DYAMOND

    elif varid == 'mcs':    # Load TOOCAN data
        
        path_TOOCAN = '/'+df.iloc[i_t]['img_seg_path']
        img_TOOCAN = xr.open_dataarray(path_TOOCAN)
        
        return img_TOOCAN

def getTimeStr(i_t):

    timestamp = df.path_dyamond[i_t].split('_')[-1]
    delta = dt.timedelta(seconds=int(int(timestamp)*7.5))
    date_t = dt.datetime(2016,8,1) + delta
    time_str = dt.datetime.strftime(date_t,"%h %d %Y, %H:%M")

    return time_str

def getCoords2D(dataset,slice_lon,slice_lat):
    
    # get correct coordinate names in dataset
    for prefix in 'lat','lon':
        r = re.compile("%s.*"%prefix)
        coord = list(filter(r.match,list(dataset.coords.dims)))[0]
        setattr(thismodule,'%s_coord'%prefix,coord)
    
    # extract coordinates
    lat_1D = dataset[lat_coord].sel({lat_coord:slice_lat})
    lon_1D = dataset[lon_coord].sel({lon_coord:slice_lon})

    # compute 2D meshgrid of coordinates
    lonarray,latarray = np.meshgrid(lon_1D,lat_1D)
    
    return lonarray,latarray

def showColorBar(fig,ax,im,varid):
    
    if varid == 'prec':
        ylabel = '30-mn Precipitation (mm)'

    elif varid == 'PW':
        ylabel = 'Column humidity (mm)'
    
    x,y,w,h = ax.get_position().bounds
    dx = w/60
    cax = plt.axes([x+w+2*dx,y,dx,h])
    cbar = fig.colorbar(im, cax=cax, orientation='vertical')
    cbar.ax.set_ylabel(ylabel)

def make_snapshot(i_t0,Lx_fig,Ly_fig,filename,show_mcs=False,title=None,norm=None):
    
    # load data
    var_DYAMOND = importData(i_t0,varid,df)
    
    if show_mcs:
        im_TOOCAN = importData(i_t0,'mcs',df)
    
    # initialize figure
    fig = plt.figure(figsize=(Lx_fig,Ly_fig))
    
    # remove margins~ish
    fig.subplots_adjust(left=0.02)
    
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))

    ims = []

    for slice_lon in slice(lon_lim[0],360),slice(0,lon_lim[1]):

        # coords
        lonarray_dyamond,latarray_dyamond = getCoords2D(var_DYAMOND,slice_lon,slice_lat)            
        # data
        Z = var_DYAMOND.sel(lon=slice_lon,lat=slice_lat)
        
        # show
        im = ax.pcolormesh(lonarray_dyamond,latarray_dyamond,Z,transform=ccrs.PlateCarree(),alpha=0.9,cmap=cmap,norm=norm)
        
        #- MCSs
        
        if show_mcs:

            # coords
            lonarray_toocan,latarray_toocan = getCoords2D(img_TOOCAN,slice_lon,slice_lat)            
            # data
            IMG_SEG = img_TOOCAN.sel(longitude=slice_lon,latitude=slice_lat)[0]%10    # display cyclic colors with period 10 (same as in corresponding cmap)
            # show
            im_MCS = ax.pcolormesh(lonarray_toocan,latarray_toocan,IMG_SEG,transform=ccrs.PlateCarree(),cmap=cmap_mcs,alpha=1)

            im_store = [im,im_MCS]
            
        else:
            
            im_store = im
        
        # store image placeholder for later updating
        ims.append(im_store)

    # delete data and remove from memory
    del var_DYAMOND
    del Z
    if show_mcs:
        del IMG_TOOCAN
        del IMG_SEG
    gc.collect()
    
    # cosmetics
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(title)
    showColorBar(fig,ax,im,varid)
        
    ax.set_extent([lon_lim[0]-360,lon_lim[1], *lat_lim],crs=ccrs.PlateCarree(central_longitude=0))
    ax.coastlines('110m')
    ax.gridlines()
    
    plt.savefig(os.path.join(DIR_FIG,filename),dpi=140)#,bbox_inches='tight')
    
    return fig, ax, ims



if __name__ == "__main__":
    
    # time step between 832 and 1917
    i_t = 914
    varid = 'PW'

    ## image
    # range
    clim = clim_specs[varid]
    # colors
    norm = norm_specs[varid]
    cmap = cmap_specs[varid]
    # range
    lon_lim = (280,100)
    lat_lim = (-10,30)
    slice_lat = slice(*lat_lim)
    
    ## compute figure size
    dlon = np.diff(lon_lim)[0] % 360
    dlat = np.diff(lat_lim)[0]
    Lx_fig = 15
    Lx_cbar = 1.5
    Ly_title = 1
    Ly_fig = (Lx_fig-Lx_cbar)/dlon*dlat + Ly_title
    
    # correspondence table
    df = pd.read_csv(os.path.join(DIR_DATA,'relation_2_table_UTC_dyamond_segmentation.csv'))
    df.sort_values(by='UTC',ignore_index=True,inplace=True)
    
    # output
    filename = '%s_DYAMOND_SAM_%d.png'%(varid,i_t)
    
    title_root = 'DYAMOND-Summer SAM-4km, %s'
    t_str = getTimeStr(i_t)
    title = title_root%t_str
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        make_snapshot(i_t,Lx_fig,Ly_fig,filename,title=title,norm=norm)
    