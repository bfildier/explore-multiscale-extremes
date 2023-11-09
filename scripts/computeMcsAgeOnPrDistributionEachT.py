"""script computeMcsAgeOnPrDistributionEachT

load prec at one time slice from DYAMOND and a TOOCAN segmentation mask
to compute MCS age characteristics in extreme bins of the rain distribution.

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
from load_TOOCAN_DYAMOND_modif_BF import *
# mapping function
from lifecycle_mapping import *



def loadToocanTimeTable(timetable_file='TOOCAN_time_table.csv',overwrite = False):
    
    # set location on disk
    timetable_path = os.path.join(DIR_DATA,timetable_file)

    # check if table exists on disk
    table_exists = len(glob.glob(timetable_path)) > 0 

    if table_exists and not overwrite:

        #-- read from disk
        print('reading from %s'%timetable_path)
        toocan_timetable = pd.read_csv(timetable_path)

        # execution takes:
        # CPU times: user 98 ms, sys: 23.7 ms, total: 122 ms
        # Wall time: 159 ms

    else:

        print('n_MCS :',len(toocan))

        #-- initialize
        toocan_timetable = pd.DataFrame(columns=['label','i_t_min','i_t_max','duration'], index=np.arange(np.nanmax(labels_toocan)))

        #-- fill
        for i_MCS in range(len(toocan)):

            if i_MCS%10000 == 0:
                print(i_MCS,end='..')

            MCS = toocan[i_MCS]

            # label
            label = MCS.label

            # birth
            i_t_min = np.where(relation_table.UTC == MCS.Utime_Init)[0][0]

            # death
            i_t_max = np.where(relation_table.UTC == MCS.Utime_End)[0][0]

            # duration (equal (i_t_max-i_t_min+1)*0.5 hrs)
            duration = MCS.duration

            # save
            toocan_timetable.loc[label] = pd.Series({'label':label,'i_t_min':i_t_min,'i_t_max':i_t_max,'duration':duration})

        #--- Save to disk 
        if overwrite or not table_exists:
            print('saving to %s'%timetable_path)
            toocan_timetable.to_csv(timetable_path)
        
        # loop execution takes:
        # CPU times: user 2min 47s, sys: 4.53 s, total: 2min 51s
        # Wall time: 2min 47s
        
    return toocan_timetable


##-- load chunked distribution

# def defineCaseStudy(varid,region,mask,name,relation_table):
    
#         # TODO: DECIDE IF REDEFINE EACH TIME OR LOAD FROM DISK

#     cs = CaseStudy(name=name,
#                    region=region,
#                    rel_tab_dyam_seg=relation_table)

#     cs.setSimulationSpecs(i_t_min = 832,
#                           i_t_max = 1917)

#     # where sliced distributions are stored
#     cs.setDirectories(varid,mask)

#     # load all distributions
#     cs.loadDistSliced(varid,mask)

#     # compute mean precip time series
#     cs.computeMean(varid='Prec',mask='all')

#     # find "missing" time indices
#     cs.findTimeIndToIgnore()

#     # store times
#     cs.storeTimes()

#     # compute full distributions
#     cs.combineDistributions(varid,mask)
    
#     return cs

def loadCaseStudy(varid,region,mask):
    
    varid_str = varid
    if varid == 'Prec':
        varid_str = 'pr'
    # output dir
    dir_in = os.path.join(DIR_OUT,region,varid,mask)
    # save
    cs = pickle.load(open(os.path.join(dir_in,'case_study.pickle'),'rb'))
    # output
    return cs

def getDistFromCaseStudy(cs,varid,mask):
    
    varid_str = cs.getVaridStr(varid)
    dist_name = 'dist_%s_%s'%(varid_str,mask)
    dist_var = getattr(cs,dist_name)
    
    return dist_var

def writeDistInCaseStudy(cs,dist_var,varid,mask):
    
    varid_str = cs.getVaridStr(varid)
    dist_name = 'dist_%s_%s'%(varid_str,mask)
    setattr(cs,dist_name,dist_var)

def loadMaskXY(maskname):
    
    if maskname == 'allXY':
        
        return True
    
    elif maskname in ['land','ocean']:
        
        landfile = 'DYAMOND_9216x4608x74_7.5s_4km_4608_0000200160.LANDMASK.2D.nc'
        landpath = os.path.join(DIR_DYAMOND,landfile)
        landmask = xr.open_dataarray(landpath)[0]
        
        if maskname == 'land':
            mask = landmask > 0.5
        elif maskname == 'ocean':
            mask = landmask < 0.5
            
        return mask

def computeAtSlice(i_t,cs,lcm,dist_var,mask_xy,mask_T,diag,metric,save_dir):
    """Loads TOOCAN segmentation mask and varid at time slice i_t and computes conditional Age on varid Distribution.
    
    Arguments:
    - cs: CaseStudy object.
    - lcm: LifeCycleMapping object.
    - dist_var: reference DistributionChunked onto which to compute the diagnostic.
    - mask_xy: name of spatial mask on which to restrict calculation of diagnostic (e.g. 'land','ocean', etc.)
    - mask_T: name of mask on TOOCAN MCS lifetime on which to restrict calculation of diagnostic (e.g. allT, min5hr, etc.)
    - diag: diagnostic in argument of current script. Must be linear.
    - metric: age metric in argument of current script.
    - save_dir: path to output directory
    """

    print('compute age composite @ i_t = %d ...'%i_t, end=' ')

    # load segmentation mask for that date
    segmask = loadTOOCANSeg(i_t,lcm.relation_table)[0].sel(latitude=lat_slice).values

    # load prec data for that date
    sample = loadPrec(i_t,lcm.relation_table).sel(lat=lat_slice).values
    
    # load mask
    xymask = loadMaskXY(mask_xy)
    if xymask.__class__ is not bool:
        xymask = xymask.sel(lat=lat_slice).values.flatten()

    if i_t in cs.times_to_ignore:     # skip if index to ignore (corrupted data)
        
        diag_all_bins_t,N_all_bins_t,N_valid_all_bins_t = None, None, None

    else:    # compute if valid data
        
        diag_all_bins_t,N_all_bins_t,N_valid_all_bins_t = lcm.compositeMcsAgeOnDist(i_t,
                                                                                    segmask,
                                                                                    sample,
                                                                                    dist_var,
                                                                                    xymask,
                                                                                    mask_T,
                                                                                    diag,
                                                                                    metric)

    # save on disk
    pickle.dump(diag_all_bins_t,open(os.path.join(save_dir,'%s_t_%d.pickle'%(diag_name,i_t)),'wb'))
    pickle.dump(N_all_bins_t,open(os.path.join(save_dir,'N_all_bins_t_%d.pickle'%(i_t)),'wb'))
    pickle.dump(N_valid_all_bins_t,open(os.path.join(save_dir,'N_valid_all_bins_t_%d.pickle'%(i_t)),'wb'))

    

##--- Main script

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Compute rain distribution at time slice index provided.')
    parser.add_argument('-i','--index', type=int, nargs=1, default=[832],
                    help='index in relation table DYAMOND-SAM:segmentationTOOCAN')
    parser.add_argument('-n','--ninds', type=int, nargs=1, default=[1086],
                    help='number of indices in a row to analyze')
    parser.add_argument('-r','--region',type=str,default='tropics',
                    help='region of analysis')
    parser.add_argument('--metric',type=str,default='norm_age',
                    help='age metric to compute, e.g. age, norm_age')
    parser.add_argument('-d','--diag',type=str,default='mean',
                    help='diagnostic, e.g. mean, max, var, min')
    parser.add_argument('--name', type=str,default='DYAMOND-SAM',
                    help='case name, typically <project>-<model>')
    parser.add_argument('-m','--mask_ref',type=str,default='all',
                    help='mask over which the reference rain distribution is computed in subregion of analysis')
    parser.add_argument('--mask_xy',type=str,default='allxy',
                    help='name of spatial mask on which to restrict calculation of diagnostic (e.g. land,ocean, etc.)')
    parser.add_argument('--mask_T',type=str,default='allT',
                    help='name of mask on TOOCAN MCS lifetime on which to restrict calculation of diagnostic (e.g. allT, min5hr, etc.)')
    parser.add_argument('--load_toocan',type=bool,default=False,
                    help='whether or not loading toocan data is necessary')
    
    args = parser.parse_args()
    
    ##--- arguments ---##
    
    i_t0 = args.index[0]
    Ni = args.ninds[0]
    
    # what to compute
    metric = args.metric
    diag = args.diag
    region = args.region
    mask_ref = args.mask_ref
    mask_xy = args.mask_xy
    mask_T = args.mask_T
    name = args.name
    
    # others
    varid = 'Prec'
    load_toocan = args.load_toocan
    print('load_toocan : ',load_toocan)
    print('class:',load_toocan.__class__)

    
    ##--- Prepare ---##
    
    # load relation table
    print('- load relation table')
    relation_table = loadRelTable('DYAMOND_SEG')
    
    # MCSs
    print('- load TOOCAN')
    if load_toocan:
        toocan = loadAllMCSs(DIR_TOOCAN_DYAMOND,load_TOOCAN_DYAMOND)
    else:
        toocan = None
    
    # TOOCAN time table
    print('- load TOOCAN time table')
    toocan_timetable = loadToocanTimeTable()
    
    # create case & get distribution chunked 
    print('- create casestudy & load distribution')
    # cs = defineCaseStudy(varid,region,mask,name,relation_table)
    cs = loadCaseStudy(varid,region,mask_ref)
    dist_var = getDistFromCaseStudy(cs,varid,mask_ref)
    
    # create MCS mapping object
    print('- create MCS mapping')
    lcm = LifeCycleMapping(relation_table,toocan_timetable,toocan)

    # geography
    lat_slice = coord_slices['lat'][region]

    # define saving info
    diag_name = '%s_%s_%s_%s'%(diag,metric,mask_xy,mask_T)
    save_dir = os.path.join(DIR_OUT,region,diag_name,mask_ref,'time_slices')
    os.makedirs(save_dir,exist_ok=True)
    
    
    ##--- Compute ---##
    
    print('- compute MCS age diagnostic at each time')
    
    # 1. at each time slice
    
    diag_all_bins = [] # store diagnostics vs. bin at each t
    N_all_bins = [] # store sample size vs. bin each t
    N_valid_all_bins = [] # store sample size vs. bin each t
    
    for i_t in range(i_t0,i_t0+Ni):
                        
        computeAtSlice(i_t,cs,lcm,dist_var,mask_xy,mask_T,diag,metric,save_dir)
        
        
    
#     # ADAPT
#         # open file
#         file = tarfile.open(os.path.join(getattr(self,'dir_dist_%s_%s_sliced'%(varid_str,mask)),'time_slices.tar.gz'))
#         # extract file to temporary folder
#         file.extractall(DIR_TEMPDATA)
        
#         # remove files from temporary folder
#         shutil.rmtree(glob.glob(os.path.join(DIR_TEMPDATA,'*'))[0])


    print('Finished. Check. √')
        
    sys.exit(0)
    
    
    
