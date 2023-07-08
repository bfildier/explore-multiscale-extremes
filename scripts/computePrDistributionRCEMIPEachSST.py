"""script computePrDistributionRCEMIPEachTime

load prec at one time slice from RCEMIP and compute distribution of extremes on inverse logarithmic bins.

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


workdir = os.getcwd()
# Add own module library to path
moduledir, fcndir = defineDir(workdir,verbose=False)

# load own libraries
from conditionalstats import *
# from plot1D import *

# own DYAMOND functions
from fcns_load_RCEMIP_SAM import *

## current script object
thismodule = sys.modules[__name__]
print('thismodule :',thismodule)

def getTOOCANmask(toocan_seg):

    # boolean mask
    toocan_mask = toocan_seg.data.flatten() > 0
    
    return toocan_mask

def getVarThresholdMask(var,thres):

    # # assume mask string format is $varid_$thres
    # varid,thres = tuple(mask.split('_'))

    # define boolean mask
    var_mask = var.data.flatten() > float(thres)
    
    return var_mask
    
def computeDistribution(name,var,mask,bintype):
    
    # initialize
    dist = Distribution(name=name,nbins=50,bintype=bintype)
    # data
    data_1D = var.data.flatten()
    # mask data
    data_1D_masked = data_1D[mask]
    # compute distribution
    dist.computeDistribution(sample=data_1D_masked)
    # compute mean value
    dist.computeMean(sample=data_1D_masked)
    
    return dist
    
def computeAllStats(SST):

    # print('compute all distributions @ i_t = %d ...'%i_t, end=' ')
    
    mmd_to_mmhr = 1/24

    #-- load all variables
    
    Prec = loadVarRCEMIP('Prec',SST)*mmd_to_mmhr
    PW = loadVarRCEMIP('PW',SST)
    TOOCANSEG = loadTOOCANSEG_RCEMIP(SST)

    #-- generate masks
    
    mask_PW_40 = getVarThresholdMask(PW,40)
    mask_TOOCAN = getTOOCANmask(TOOCANSEG)
    mask_PW_40_TOOCAN = np.logical_and(mask_PW_40,mask_TOOCAN)
    
    #-- prepare save directory
    save_dir = os.path.join(DIR_OUT,'RCEMIP',"%sK"%SST)
    os.makedirs(save_dir,exist_ok=True)
    
    #-- calculations
    
    varids = 'Prec',     'Prec',   'Prec',        'PW',    'PW',     'PW',     'PW'
    maskids = 'all',     'PW_40',  'PW_40_TOOCAN','all',    'all',   'TOOCAN', 'PW_40'
    bintypes = 'invlogQ','invlogQ','invlogQ',     'linear','invlogQ','invlogQ','invlogQ'
    
    # compute fraction in mask
    print('compute fraction in masks')
    for maskid in maskids:
        
        # load & compute
        if maskid == 'all':
            mask = slice(None)
            frac = 1.
        else:
            # mask = getattr(thismodule,'mask_%s'%maskid)
            mask = locals()['mask_%s'%maskid]
            frac = np.sum(np.array(mask,dtype=int))/mask.size
        # save
        pickle.dump(frac,open(os.path.join(save_dir,'frac_%s.pickle'%(maskid)),'wb'))
    
    # compute distributions
    for varid,maskid,bintype in zip(varids,maskids,bintypes):
        
        print('compute distribution for varid %s, mask %s and bintype %s'%(varid,maskid,bintype))
        
        name = '%s, RCEMIP-SAM %sK, mask=%s,bins=%s'%(varid,SST,mask,bintype)
        # var = getattr(thismodule,varid)
        var = locals()[varid]
        if maskid == 'all':
            mask = slice(None)
        else:
            # mask = getattr(thismodule,'mask_%s'%maskid)
            mask = locals()['mask_%s'%maskid].flatten()
    
        # compute
        dist = computeDistribution(name,var,mask,bintype)
        
        # save on disk
        pickle.dump(dist,open(os.path.join(save_dir,'dist_%s_%s_%s.pickle'%(varid,maskid,bintype)),'wb'))

    print('saved')

##--- Main script

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Compute Prec and PW distributions over the full dataset.')
    parser.add_argument('--sst',type=int,default='300',
                    help='Sea Surface Temperature (K)')
    # parser.add_argument('-m','--mask',type=str,default='all',
    #                 help='mask of subregion to analyze')
    
    args = parser.parse_args()
    
    print('SST = %dK'%args.sst)
    print()
    
    computeAllStats(args.sst)
        
    sys.exit(0)
    