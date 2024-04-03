



import sys,os,glob
import psutil

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
from pprint import pprint

from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
import matplotlib.animation as animation
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


#-- INITIALIZATION

def importData(i_t):
    
    # paths
    root_DYAMOND = df.iloc[i_t]['path_dyamond'] + '.%s.2D.nc'
    file_PW_DYAMOND = root_DYAMOND%'PW'
    file_Precac_DYAMOND = root_DYAMOND%'Precac'
    path_TOOCAN = '/'+df.iloc[i_t]['img_seg_path']

    # Load DYAMOND data
    PW_DYAMOND = xr.open_dataarray(os.path.join(DIR_DYAMOND,file_PW_DYAMOND))
    
    # Load TOOCAN data
    img_TOOCAN = xr.open_dataarray(path_TOOCAN)
    
    return PW_DYAMOND, img_TOOCAN

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

def showColorBar(fig,ax,im):
    
    x,y,w,h = ax.get_position().bounds
    dx = w/60
    cax = plt.axes([x+w+2*dx,y,dx,h])
    cbar = fig.colorbar(im, cax=cax, orientation='vertical')
    cbar.ax.set_ylabel('Column humidity (mm)')

def initFigure(i_t0,Lx_fig,Ly_fig,title=None,norm=None):
    
    # load data
    PW_DYAMOND, img_TOOCAN = importData(i_t0)
    
    # initialize figure
    fig = plt.figure(figsize=(Lx_fig,Ly_fig))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))

    ims = []

    for slice_lon in slice(lon_lim[0],360),slice(0,lon_lim[1]):

        #- background

        # coords
        lonarray_dyamond,latarray_dyamond = getCoords2D(PW_DYAMOND,slice_lon,slice_lat)            
        # data
        Z = PW_DYAMOND.sel(lon=slice_lon,lat=slice_lat)[0]
        # show
        im = ax.pcolormesh(lonarray_dyamond,latarray_dyamond,Z,transform=ccrs.PlateCarree(),alpha=0.9,cmap=cmap)
        im.set_clim(*clim)

        #- MCSs
        
        # coords
        lonarray_toocan,latarray_toocan = getCoords2D(img_TOOCAN,slice_lon,slice_lat)            
        # data
        IMG_SEG = img_TOOCAN.sel(longitude=slice_lon,latitude=slice_lat)[0]%10    
        # show
        im_MCS = ax.pcolormesh(lonarray_toocan,latarray_toocan,IMG_SEG,transform=ccrs.PlateCarree(),cmap=cmap_mcs,alpha=1)

        # store image placeholders for later updating
        ims.append([im,im_MCS])

    # delete data and remove from memory
    del PW_DYAMOND
    del img_TOOCAN
    del Z
    del IMG_SEG
    gc.collect()
    
    # cosmetics
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(title)
    showColorBar(fig,ax,im)
        
    ax.set_extent([lon_lim[0]-360,lon_lim[1], *lat_lim],crs=ccrs.PlateCarree(central_longitude=0))
    ax.coastlines('110m')
    ax.gridlines()
    
    return fig, ax, ims



#-- MOVIE

def makeMovie(movie_path,j_start,j_end,Lx_fig,Ly_fig):
    
    Nt = len(df)
    
    # -- initialize figure
    
    # first time step
    i_t0 = j_row_t0
    
    title_root = 'DYAMOND-Summer SAM-4km, %s'
    # t_str = "%3.2f"%PW_DYAMOND.time.data[0]
    t_str = getTimeStr(i_t0)
    title = title_root%t_str
    # title = ''
    
    # initialize
    fig, ax, ims = initFigure(i_t0,Lx_fig,Ly_fig,title=title)

    # remove margins~ish
    fig.subplots_adjust(left=0.02)
    
    # -- define movie loop
    
    print()
    print('i_t:  MB used:')
    def updateImage(i_t):
        
        print(i_t,end=' ')
        
        # load data at i_t
        PW_DYAMOND, img_TOOCAN = importData(i_t)
        
        t_str = getTimeStr(i_t)
        
        for slice_lon,i_ims in zip([slice(lon_lim[0],360),slice(0,lon_lim[1])],list(np.arange(2))):
        
            Z = PW_DYAMOND.sel(lon=slice_lon,lat=slice_lat)[0]
            IMG_SEG = img_TOOCAN.sel(longitude=slice_lon,latitude=slice_lat)[0]%10    

            # update images
            ims[i_ims][0].set_array(np.ravel(Z[:-1,:-1].data))
            ims[i_ims][1].set_array(np.ravel(IMG_SEG[:-1,:-1].data))
        
        ax.set_title(title_root%t_str)
        
        MB_used = psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2
        print(', %7.1f'%MB_used)

        # delete data and remove from memory
        del PW_DYAMOND
        del img_TOOCAN
        del Z
        del IMG_SEG
        gc.collect()
        
        return [ims]

    print()
    anim = animation.FuncAnimation(fig, updateImage,
                                   frames=range(j_start,j_end), interval=interval, blit=False)
    
    writer = animation.writers['ffmpeg'](fps=frame_rate)

    anim.save(movie_path,writer=writer,dpi=150)#,savefig_kwargs={'bbox_inches':'tight'})
    
    
    
    
if __name__ == "__main__":
    
    #-- Paths
    
    DIR_DYAMOND = '/bdd/DYAMOND/SAM-4km/OUT_2D'
    DIR_DATA = '../input'
    
    #-- PARAMETERS

    ## animation
    N_frames_per_day = 48
    N_days_per_movie_second = 0.25
    frame_rate = N_frames_per_day * N_days_per_movie_second
    interval = int(1000/frame_rate)

    ## image
    cmap = plt.cm.RdBu
    # cmap_mcs = plt.cm.get_cmap('rainbow', 10)
    cmap_mcs = plt.cm.get_cmap('Accent', 10)
    clim = (10,70)
    # lon_lim = (280,100)
    lon_lim = (100,280)
    lat_lim = (-10,30)
    slice_lat = slice(*lat_lim)
    
    ## compute figure size
    dlon = np.diff(lon_lim)[0] % 360
    dlat = np.diff(lat_lim)[0]
    Lx_fig = 15
    Lx_cbar = 0
    Ly_title = 0.5
    Ly_fig = (Lx_fig-Lx_cbar)/dlon*dlat + Ly_title
    
    ## time bounds
    i_start = 832
    # i_start = 860    
    i_end = 1917
    # i_end = 836
    
    #-- check time bounds
    
    # correspondence table
    df = pd.read_csv(os.path.join(DIR_DATA,'relation_2_table_UTC_dyamond_segmentation.csv'))
    df.sort_values(by='UTC',ignore_index=True,inplace=True)
    
    # Find first available index
    t0_SAM = 200400
    root_SAM_t0 = 'DYAMOND_9216x4608x74_7.5s_4km_4608_0000%d'%t0_SAM
    j_row_t0 = df.index[df.path_dyamond == root_SAM_t0][0]
    print('first available index:',j_row_t0)
    
    if i_start < j_row_t0:
        err_str = "start index is too small (minimum is %d)"%j_row_t0
        raise ValueError(srr_str)
    
    # Find last available index
    t1_SAM = 460800
    root_SAM_t1 = 'DYAMOND_9216x4608x74_7.5s_4km_4608_0000%d'%t1_SAM
    j_row_t1 = df.index[df.path_dyamond == root_SAM_t1][0]
    print('last available index:',j_row_t1)
    
    if i_end > j_row_t1:
        err_str = "end index is too large (maximum is %d)"%j_row_t1
        raise ValueError(err_str)
    
    #-- Make movie
    
    moviedir = '../movies/'
    movie_name = 'toocan_DYAMOND_SAM_%dE%dE'%(lon_lim)
    movie_path = os.path.join(moviedir, '%s.mp4'%(movie_name))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        makeMovie(movie_path,i_start,i_end,Lx_fig,Ly_fig)
