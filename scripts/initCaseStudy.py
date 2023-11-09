"""script initCaseStudy

initialize CaseStudy Object and merge precipitation distribution as DistributionChunked object.

B. Fildier October 2023
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

# case 
from casestudy import *
# to access segmentation files and simulation outputs
from fcns_load_DYAMOND_SAM import *
# to access TOOCAN objects


##-- load chunked distribution

def defineCaseStudy(varid,region,mask,name,relation_table):

    cs = CaseStudy(name=name,
                   region=region,
                   rel_tab_dyam_seg=relation_table)

    cs.setSimulationSpecs(i_t_min = 832,
                          i_t_max = 1917)

    # where sliced distributions are stored
    cs.setDirectories(varid,mask)

    # load all distributions
    cs.loadDistSliced(varid,mask)

    # compute mean precip time series
    cs.computeMean(varid='Prec',mask='all')

    # find "missing" time indices
    cs.findTimeIndToIgnore()

    # store times
    cs.storeTimes()

    # compute full distributions
    cs.combineDistributions(varid,mask)
    
    return cs


def getDistFromCaseStudy(cs,varid,mask):
    
    varid_str = cs.getVaridStr(varid)
    dist_name = 'dist_%s_%s'%(varid_str,mask)
    dist_var = getattr(cs,dist_name)
    
    return dist_var

def savePrDistribution(cs,varid,mask):
    
    varid_str = cs.getVaridStr(varid)
    # get distribution
    dist = getDistFromCaseStudy(cs,varid,mask)
    # output dir
    dir_out = getattr(cs,'dir_dist_%s_%s_sliced'%(varid_str,mask))
    # save
    pickle.dump(dist,open(os.path.join(dir_out,'dist_%s.pickle'%varid_str),'wb'))
    
def saveCaseStudy(cs,varid,mask):
    
    varid_str = cs.getVaridStr(varid)
    # output dir
    dir_out = getattr(cs,'dir_dist_%s_%s_sliced'%(varid_str,mask))
    # save
    pickle.dump(cs,open(os.path.join(dir_out,'case_study.pickle'),'wb'))
    

##--- Main script

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Define case study and merge rain distribution from all time slices.')
    parser.add_argument('-r','--region',type=str,default='tropics',
                    help='region of analysis')
    parser.add_argument('--name', type=str,default='DYAMOND-SAM',
                    help='case name, typically <project>-<model>')
    parser.add_argument('-m','--mask',type=str,default='all',
                    help='mask of subregion to analyze')

    args = parser.parse_args()
    
    ##--- arguments ---##

    region = args.region
    mask = args.mask
    name = args.name
    
    # variables
    varid = 'Prec'
    
    ##--- Prepare ---##
    
    # load relation table
    print('- load relation table')
    relation_table = loadRelTable('DYAMOND_SEG')
    
    ##--- main ---##
    
    # init case study
    print('- create case study')
    cs = defineCaseStudy(varid,region,mask,name,relation_table)
    
    # save distribution to disk
    print('- save merged distribution to disk')
    savePrDistribution(cs,varid,mask)

    # save case study to disk
    print('- save case study to disk')
    saveCaseStudy(cs,varid,mask)
    
    print('Done. √')
    sys.exit(0)