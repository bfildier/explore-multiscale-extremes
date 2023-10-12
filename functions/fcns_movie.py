# REVISE THIS

import matplotlib.pyplot as plt
import xarray as xr
import datetime as dt
import os,glob,re
import numpy as np

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
lon_lim = (280,100)
lat_lim = (-10,30)
slice_lat = slice(*lat_lim)

# compute figure size
dlon = np.diff(lon_lim)[0] % 360
dlat = np.diff(lat_lim)[0]
Lx_fig = 15
Lx_cbar = 1.5
Ly_title = 1
Ly_fig = (Lx_fig-Lx_cbar)/dlon*dlat + Ly_title

# lat_dyamond = PW_DYAMOND.lat.sel(lat=slice_lat)
# lon_dyamond = PW_DYAMOND.lon.sel(lon=slice_lon)

# lonarray_dyamond,latarray_dyamond = np.meshgrid(lon_dyamond,lat_dyamond)


#-- INITIALIZATION

def importData(i_t):
    
    # paths
    root_DYAMOND = df.iloc[i_t]['path_dyamond'] + '.%s.2D.nc'
    file_PW_DYAMOND = root_DYAMOND%'PW'
    path_TOOCAN = '/'+df.iloc[i_t]['img_seg_path']
    path_PW_DYAMOND = os.path.join(DIR_DYAMOND,file_PW_DYAMOND)
    
    if len(glob.glob(path_PW_DYAMOND)) == 0:
        print(file_PW_DYAMOND,'does not exist')

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

def showSnapshot(ax,slice_lon,slice_lat,data=None,segmask=None):
    
    ims = []

    for slice_lon in slice(lon_lim[0],360),slice(0,lon_lim[1]):

        #- background

        # coords
        lonarray_dyamond,latarray_dyamond = getCoords2D(data,slice_lon,slice_lat)            
        # data
        Z = data.sel(lon=slice_lon,lat=slice_lat)[0]
        # show
        # im = ax.pcolormesh(np.ravel(lonarray_dyamond),np.ravel(latarray_dyamond),np.ravel(Z),transform=ccrs.PlateCarree(),alpha=0.9,cmap=cmap)
        im = ax.pcolormesh(lonarray_dyamond,latarray_dyamond,Z,transform=ccrs.PlateCarree(),alpha=0.9,cmap=cmap)
        im.set_clim(*clim)

        #- MCSs
        
        # coords
        lonarray_toocan,latarray_toocan = getCoords2D(segmask,slice_lon,slice_lat)            
        # data
        IMG_SEG = segmask.sel(longitude=slice_lon,latitude=slice_lat)[0]%10    
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
    
    return ims

def initFigure(i_t0,Lx_fig=15,Ly_fig=4,title=None,norm=None):
    
    # load data
    PW_DYAMOND, img_TOOCAN = importData(i_t0)
    
    # initialize figure
    fig = plt.figure(figsize=(Lx_fig,Ly_fig))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))

    ims = showSnapshot(ax,slice_lon,slice_lat,data=PW_DYAMOND,segmask=img_TOOCAN)
    
    return fig, ax, ims