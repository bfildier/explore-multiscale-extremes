import xarray as xr
import gc,os
import pandas as pd

from settings import *

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