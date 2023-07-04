"""script computeVarDistributionEachT

load any variable at one time slice from DYAMOND and compute distribution of extremes on inverse logarithmic bins.

B. Fildier July 2023
"""

##--- Python modules

import sys,os,glob

import numpy as np
import xarray as xr
import pandas as pd
from pprint import pprint
import warnings

import datetime as dt
import re
import gc
import warnings
import argparse
import pickle

##--- Own settings and modules

from settings import *

## current script object
thismodule = sys.modules[__name__]
workdir = os.getcwd()
# Add own module library to path
moduledir, fcndir = defineDir(workdir,verbose=False)

# load own libraries
from conditionalstats import *
from plot1D import *

# own DYAMOND functions
from fcns_load_DYAMOND_SAM import *

def getTOOCANmask(i_t,df,lon_slice,lat_slice):
    
    # TOOCAN segmentation mask
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        toocan_seg = loadTOOCANSeg(i_t,df)[0].sel(latitude=lat_slice).sel(longitude=lon_slice)

    # boolean mask
    toocan_mask = toocan_seg > 0
    
    return toocan_mask

def getVarThresholdMask(i_t,df,lon_slice,lat_slice,mask):

    # assume mask string format is $varid_$thres
    varid_thres,thres = tuple(mask.split('_'))

    # get varid values
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        var_thres = loadVar(i_t,df,varid_thres)[0].sel(lat=lat_slice).sel(lon=lon_slice)

    # define boolean mask
    var_mask = var_thres > float(thres)
    
    return var_mask
    
def computeAtSlice(i_t,varid,lon_slice,lat_slice,region,mask):

    print('compute distribution @ i_t = %d ...'%i_t, end=' ')

    # relation table DYAMOND-SAM - TOOCAN segmentation
    reltab_dyam_seg = loadRelTable('DYAMOND_SEG')

    #-- calculation
    var = loadVar(i_t,reltab_dyam_seg,varid).sel(lat=lat_slice).sel(lon=lon_slice)
    # data
    var_1D = var.data.flatten()
    
    # select data in subregion
    if mask == 'all':
        
        var_1D_masked = var_1D

    elif mask == 'TOOCAN':
      
        toocan_mask = getTOOCANmask(i_t,reltab_dyam_seg,lon_slice,lat_slice).data.flatten()
        
        # valid pr values
        var_1D_masked = var_1D[toocan_mask]
        
    elif re.compile('TOOCAN_.*').match(mask):
        
        # get TOOCAN mask
        toocan_mask = getTOOCANmask(i_t,reltab_dyam_seg,lon_slice,lat_slice).data.flatten()
        # get var mask
        var_mask = getVarThresholdMask(i_t,reltab_dyam_seg,lon_slice,lat_slice,mask.replace("TOOCAN_","",1)).data.flatten()
        # merge masks (boolean intersection)
        full_mask = np.logical_and(toocan_mask,var_mask)
        
        # valid pr values
        var_1D_masked = var_1D[full_mask]
    
    else:
        
        var_mask = getVarThresholdMask(i_t,reltab_dyam_seg,lon_slice,lat_slice,mask).data.flatten()
        
        # valid pr values
        var_1D_masked = var_1D[var_mask]
        
    # init
    dist_pr_t = Distribution(name='%s, DYAMOND-SAM %s, i_t=%d'%(varid,region,i_t),nbins=50,bintype='invlogQ')
    # compute distribution
    dist_pr_t.computeDistribution(sample=var_1D_masked)
    # compute mean value
    dist_pr_t.computeMean(sample=var_1D_masked)

    # save on disk
    save_dir = os.path.join(DIR_OUT,region,varid,mask,'time_slices')
    os.makedirs(save_dir,exist_ok=True)
    pickle.dump(dist_pr_t,open(os.path.join(save_dir,'dist_%s_t_%d.pickle'%(varid,i_t)),'wb'))

    print('saved')

##--- Main script

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Compute input var distribution at time slice index provided.')
    parser.add_argument('-v','--varid', type=str, nargs=1, default='PW',
                    help='variable name of interest')
    parser.add_argument('-i','--index', type=int, nargs=1, default=[832],
                    help='index in relation table DYAMOND-SAM:segmentationTOOCAN')
    parser.add_argument('-n','--ninds', type=int, nargs=1, default=[1086],
                    help='number of indices in a row to analyze')
    parser.add_argument('-r','--region',type=str,default='tropics',
                    help='region of analysis')
    parser.add_argument('-m','--mask',type=str,default='all',
                    help='mask of subregion to analyze')
    
    args = parser.parse_args()
    i_t0 = args.index[0]
    Ni = args.ninds[0]
    varid = args.varid[0]

    # geographical parameters
    if args.region == 'tropics':
        lon_slice = slice(None,None)
        lat_slice = slice(-30,30)
    else:
        reg_coords = regionCoordsFromName(args.region)
        lon_slice = slice(reg_coords[0]%360,reg_coords[1]%360)
        lat_slice = slice(reg_coords[2],reg_coords[3])
    
    # compute
    for i_t in range(i_t0,i_t0+Ni):
        
        computeAtSlice(i_t,varid,lon_slice,lat_slice,args.region,args.mask)
        
    sys.exit(0)
    