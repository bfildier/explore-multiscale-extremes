import xarray as xr
import gc,os
import pandas as pd
import re

from settings import *

## current script object
thismodule = sys.modules[__name__]

def loadPrecac(i_t,df):
    
    root = df.iloc[i_t]['path_dyamond']
    file_precac = root+'.Precac.2D.nc'
    # load
    precac = xr.open_dataarray(os.path.join(DIR_DYAMOND,file_precac)).load()[0]
    
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

def loadRelTable(which='DYAMOND_SEG'):
    
    # relation table DYAMOND-SAM -- TOOCAN segmentation masks
    if which == 'DYAMOND_SEG':
        df = pd.read_csv(os.path.join(DIR_DATA,'relation_2_table_UTC_dyamond_segmentation.csv'))
        df.sort_values(by='UTC',ignore_index=True,inplace=True)

    return df

def loadTOOCANSeg(i_t,df):
    
    path_TOOCAN = '/'+df.iloc[i_t]['img_seg_path']
    # Load TOOCAN data
    img_TOOCAN = xr.open_dataarray(path_TOOCAN)
    
    return img_TOOCAN

def loadVar(i_t,df,varid):
    
    # get filename
    root_DYAMOND = df.iloc[i_t]['path_dyamond']
    file_DYAMOND = root_DYAMOND+'.%s.2D.nc'%varid
    
    # Load DYAMOND data
    var_DYAMOND = xr.open_dataarray(os.path.join(DIR_DYAMOND,file_DYAMOND))
    
    return var_DYAMOND

def regionNameFromCoord(box):
    
    for coord,coordname in zip(box,['lonmin','lonmax','latmin','latmax']):
        
        if coord < 0 : 

            if re.compile('lon.*').match(coordname): setattr(thismodule,coordname,"%sW"%abs(coord))
            if re.compile('lat.*').match(coordname): setattr(thismodule,coordname,"%sS"%abs(coord))
            
        else:
            
            if re.compile('lon.*').match(coordname): setattr(thismodule,coordname,"%sE"%coord)
            if re.compile('lat.*').match(coordname): setattr(thismodule,coordname,"%sN"%coord)
    
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
    