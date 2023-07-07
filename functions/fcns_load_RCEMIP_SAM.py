import xarray as xr
import gc,os
import pandas as pd
import re

from settings import *

## current script object
thismodule = sys.modules[__name__]

def loadVarRCEMIP(varid,SST):
    
    fileroot = '%s_rcemip_large_2048x128x74_3km_12s_%sK_64.2Dcom.nc'
    filename = fileroot%(varid,SST)
    var = xr.open_dataset(os.path.join(DIR_RCEMIP,"%sK"%SST,filename))[varid]
    
    return var

def loadTOOCANSEG_RCEMIP(SST):
    
    fileroot = 'TOOCAN_2.07_SAM_RCE_large%s_2D_irtb.nc'
    filename = fileroot%SST
    toocanseg = xr.open_dataset(os.path.join(DIR_TOOCANSEG_RCEMIP,filename))['MCS_label']
    
    return toocanseg
    