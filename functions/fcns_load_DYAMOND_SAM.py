import xarray as xr
import gc,os
import pandas as pd
import re
import numpy as np
import warnings

from settings import *

## current script object
fcnload_module = sys.modules[__name__]


def loadRelTable(which='DYAMOND_SEG'):
    
    # relation table DYAMOND-SAM -- TOOCAN segmentation masks
    if which == 'DYAMOND_SEG':
        df = pd.read_csv(os.path.join(DIR_DATA,'relation_2_table_UTC_dyamond_segmentation.csv'))
        df.sort_values(by='UTC',ignore_index=True,inplace=True)

    return df

def loadPrecac(i_t,df=None):
    
    if df is not None:

        root = df.iloc[i_t]['path_dyamond']
        file_precac = root+'.Precac.2D.nc'
        # load
        precac = xr.open_dataarray(os.path.join(DIR_DYAMOND,file_precac)).load()[0]
    
    else:
        
        pass
    
    return precac

def loadPrec(i_t,df):
    
    # Load DYAMOND-SAM Precac
    precac_prev = loadPrecac(i_t-1,df)
    precac_current = loadPrecac(i_t,df)
    
    # get 30mn precipitation from difference
    prec = precac_current - precac_prev
    
    # free up memory
    del precac_prev
    del precac_current
    gc.collect()
    
    return prec

def loadTOOCANSeg(i_t,df,toocan_version='v2.08',keep_path=True):

    if toocan_version == 'v2.07':
            
        full_path = '/'+df.iloc[i_t]['img_seg_path']
        print(full_path)
            
        if DIR_TOOCANSEG_DYAMOND is None or keep_path:
            path_TOOCAN = full_path
        else :
            filename = os.path.basename(full_path)
            date = os.path.basename(os.path.dirname(full_path))
            path_TOOCAN = os.path.join(DIR_TOOCANSEG_DYAMOND,date,filename)
            
        # Load TOOCAN data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            img_TOOCAN = xr.open_dataarray(path_TOOCAN)
        
    elif toocan_version == 'v2.08':
        
        assert toocan_version in DIR_TOOCANSEG_DYAMOND, 'wrong TOOCAN version'
        date_key = '-0'.join(df['img_seg_path'][i_t].split('_')[-1].split('-'))
        year = date_key[:4]
        month = date_key[4:6]
        day = date_key[6:8]
        filename = 'mcs_mask_TOOCAN_SAM_'+date_key
        path_TOOCAN = os.path.join(DIR_TOOCANSEG_DYAMOND,'%s_%s_%s'%(year,month,day),filename)
        
        # Load TOOCAN data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            img_TOOCAN = xr.open_dataset(path_TOOCAN).cloud_mask
    
    return img_TOOCAN

def loadVar(i_t,df,varid):
    
    # get filename
    root_DYAMOND = df.iloc[i_t]['path_dyamond']
    file_DYAMOND = root_DYAMOND+'.%s.2D.nc'%varid
    
    # Load DYAMOND data
    var_DYAMOND = xr.open_dataarray(os.path.join(DIR_DYAMOND,file_DYAMOND))
    
    return var_DYAMOND

def getCoords2D(dataset,slice_lon,slice_lat):
    
    # get correct coordinate names in dataset
    for prefix in 'lat','lon':
        r = re.compile("%s.*"%prefix)
        coord = list(filter(r.match,list(dataset.coords.dims)))[0]
        setattr(fcnload_module,'%s_coord'%prefix,coord)
    
    # extract coordinates
    lat_1D = dataset[lat_coord].sel({lat_coord:slice_lat})
    lon_1D = dataset[lon_coord].sel({lon_coord:slice_lon})

    # compute 2D meshgrid of coordinates
    lonarray,latarray = np.meshgrid(lon_1D,lat_1D)
    
    return lonarray,latarray

def loadAllMCSs(path_to_dir,fun_load_TOOCAN):
    """load TOOCAN objects"""
    
    # where stored
    paths_TOOCAN = glob.glob(os.path.join(path_to_dir,'*.gz'))
    N_paths = len(paths_TOOCAN)
    # load
    toocan = []
    for i_p in range(N_paths):
        path = paths_TOOCAN[i_p]
        print('load %s'%path)
        toocan.extend(fun_load_TOOCAN(path))
    
    return toocan

def regionNameFromCoord(box):
    
    for coord,coordname in zip(box,['lonmin','lonmax','latmin','latmax']):
        
        if coord < 0 : 

            if re.compile('lon.*').match(coordname): setattr(fcnload_module,coordname,"%sW"%abs(coord))
            if re.compile('lat.*').match(coordname): setattr(fcnload_module,coordname,"%sS"%abs(coord))
            
        else:
            
            if re.compile('lon.*').match(coordname): setattr(fcnload_module,coordname,"%sE"%coord)
            if re.compile('lat.*').match(coordname): setattr(fcnload_module,coordname,"%sN"%coord)
    
    name = "%s_%s_%s_%s"%(lonmin,lonmax,latmin,latmax)
    
    return name

def regionCoordsFromName(name):
    
    # split into strings
    coord_str = name.split('_')
    # get absolute value
    coord_abs = [int(coord_str[i][:-1]) for i in range(len(coord_str))]
    # get direction
    coord_dir = [coord_str[i][-1:] for i in range(len(coord_str))]
    # adjust sign
    coord = [(-1)**(coord_dir[i] in ['S','W'])*coord_abs[i] for i in range(len(coord_str))]
    
    return coord

