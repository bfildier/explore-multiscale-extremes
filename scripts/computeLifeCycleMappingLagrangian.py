"""script computeLifeCycleMappingLagrangian.py

load precipitation and toocan, and map MCS properties and rain maximum onto MCS lifecycle

! Works for version TOOCAN v2.07 and v2.08 (chosen manually). For newest versions, update the following:
- TOOCAN time table
- relation_table

B. Fildier February 2024 -- April 2024
"""

##--- Python modules

import sys,os,glob

import numpy as np
import xarray as xr
import pandas as pd
from pprint import pprint
import warnings
import pickle
# from tqdm import tqdm

import datetime as dt
import re
import gc
import warnings
import argparse
# time
import time
start_time = time.time()

import multiprocessing
from multiprocessing import Pool

##--- Own settings and modules

from settings import *

## current script object
thismodule = sys.modules[__name__]
workdir = os.getcwd()
# Add own module library to path
moduledir, fcndir = defineDir(workdir,verbose=False)

def defineDir(workdir,verbose=True):

    moduledir = os.path.join(os.path.dirname(workdir),'modules')
    fcndir = os.path.join(os.path.dirname(workdir),'functions')
    sys.path.insert(0,moduledir)
    sys.path.insert(0,fcndir)
    
    if verbose:
        for includedir in [moduledir,fcndir]:
            print("Own modules available:", [os.path.splitext(os.path.basename(x))[0]
                                             for x in glob.glob(os.path.join(includedir,'*.py'))])
    
    return moduledir, fcndir

#-- load own libraries

# to access segmentation files and simulation outputs
from fcns_load_DYAMOND_SAM import *
# to access TOOCAN objects
# from load_TOOCAN_DYAMOND_modif_BF import *
from load_TOOCAN_v208_DYAMOND import * # !! for version v2.08 !
# mapping function
from lifecycle_mapping import *

# TOOCAN specs

toocan_version = 'v2.08'
toocan_version_str = toocan_version.replace('.','')

if toocan_version == 'v2.07':
    label_key = 'label'
    duration_key = 'duration'
    area_key = 'surfkm2_172Wm2'
    utime_init_key = 'Utime_Init'
    utime_end_key = 'Utime_End'
elif toocan_version == 'v2.08':
    label_key = 'DCS_number'
    duration_key = 'INT_duration'
    area_key = 'LC_surfkm2_235K'
    utime_init_key = 'INT_UTC_timeInit'
    utime_end_key = 'INT_UTC_timeEnd'

##-- functions for output setup

def loadTOOCANTimeTable():
    
    # set location on disk and check if exists
    timetable_file = 'TOOCAN_time_table%s.csv'%toocan_version_str
    timetable_path = os.path.join(DIR_DATA,timetable_file)

    # check if table exists on disk
    table_exists = len(glob.glob(timetable_path)) > 0 

    # choose to overwrite
    overwrite = False

    # Load or compute time table

    if table_exists and not overwrite:

        #-- read from disk
        print('reading from %s'%timetable_path)
        toocan_timetable = pd.read_csv(timetable_path)

    else:

         toocan_timetable = makeTOOCANTimeTable()
            
    return toocan_timetable

def makeTOOCANTimeTable():
    
    #-- initialize
    toocan_timetable = pd.DataFrame(columns=['label','i_t_min','i_t_max','duration'], index=np.arange(np.nanmax(labels_toocan)))
    
    #-- fill
    for i_MCS in tqdm(range(len(toocan))):

        # if i_MCS%10000 == 0:
        #     print(i_MCS,end='..')

        MCS = toocan[i_MCS]

        # label
        label = getattr(MCS,label_key)

        # birth
        i_t_min = utime2i_t(getattr(MCS,utime_init_key),toocan_version=toocan_version)

        # death
        i_t_max = utime2i_t(getattr(MCS,utime_end_key),toocan_version=toocan_version)

        # duration (equal (i_t_max-i_t_min+1)*0.5 hrs)
        duration = getattr(MCS,duration_key)

        # save
        toocan_timetable.loc[label] = pd.Series({'label':label,'i_t_min':i_t_min,'i_t_max':i_t_max,'duration':duration})
    
    #--- Save to disk 
    if overwrite or not table_exists:
        print('saving to %s'%timetable_path)
        toocan_timetable.to_csv(timetable_path)
        
    return toocan_timetable

def initLifeCycleMappingObject():
    
    #-- prerequisites to LCM object
    
    # load relation table
    relation_table = loadRelTable('DYAMOND_SEG')
    
    # load TOOCAN (.gz version)
    toocan = loadAllMCSs(DIR_TOOCAN_DYAMOND,load_TOOCAN_DYAMOND)
    
    # list of TOOCAN labels, for quicker mapping on toocan list
    labels_toocan = [getattr(toocan[i],label_key) for i in range(len(toocan))]
    
    # load TOOCAN timetable
    toocan_timetable = loadTOOCANTimeTable()
    
    #-- initialize LCM object
    lcm = LifeCycleMapping(relation_table,toocan_timetable,toocan,timestep)
    
    # compute mapping label to index for faster fetching of MCS in toocan list by label
    lcm.defineMappingLabelToIndex()
    
    # compute unique labels
    lcm.findValidLabels()
    
    return lcm


def initOutputArray(label_toocan,age):
    """Initialize output xarray 'ds' containing relevant fields.
    Transfer relevant TOOCAN data.
    
    Returns:
    """

    N_labels = len(label_toocan)
    N_ages = len(age)
    
    # Initialize dataset with pre-computed MCS characteristics

    ds = xr.Dataset(
        data_vars=dict(
            valid_flag=(["MCS_label"],np.full((N_labels),np.nan)),
            area=(["MCS_label", "age"], np.full((N_labels,N_ages),np.nan)),
            time_birth=(["MCS_label"], np.full((N_labels),np.nan,dtype=np.datetime64)),
            time_death=(["MCS_label"], np.full((N_labels),np.nan,dtype=np.datetime64)),
            duration=(["MCS_label"],  np.full((N_labels),np.nan)),
            prec_max=(["MCS_label", "age"],  np.full((N_labels,N_ages),np.nan)),
            longitude_max=(["MCS_label", "age"],  np.full((N_labels,N_ages),np.nan)),
            latitude_max=(["MCS_label", "age"],  np.full((N_labels,N_ages),np.nan)),
            land_mask=(["MCS_label", "age"], np.full((N_labels,N_ages),np.nan)),
        ),
        coords=dict(
            MCS_label=label_toocan, # MCS label, removing duplicates
            age=age,#'hours from time_birth'), # time from MCS birth
        ),
        attrs=dict(description="Lagrangian MCS data"),
    )

    # attributes
    ds['valid_flag'].attrs['description'] = 'True if MCS label is valid, False if duplicate'
    ds['area'].attrs['description'] = 'area in km2 from 172W/m2 threshold in TOOCAN'
    # ds['time_birth'].attrs['units'] = "hours since 2000-01-01"
    ds['MCS_label'].attrs['description'] = 'TOOCAN labels for non-duplicate MCS labels only'
    ds['MCS_label'].attrs['units'] = 'None'
    ds['age'].attrs['description'] = 'age coordinate for anvil lifecycle'
    ds['age'].attrs['units'] = 'hours from time of birth'

    return ds

def fillWithTrackingData(ds,label_toocan,lcm):
    """Fill output array with relevant lifecycle data retrieved during tracking.
    
    Arguments:
    - mask_labels_valid: mask of valid MCS labels (not duplicates)
    - 
    """

    # mask for valid labels
    ds['valid_flag'].values = np.in1d(label_toocan,lcm.labels_valid)

    #-- then, only where labels are valid:
    # birth times
    ds['time_birth'][ds['valid_flag']] = np.array([utime2Datetime(getattr(lcm.toocan[lcm.toocan_index_of_label[lcm.labels_valid[i]]],utime_init_key),toocan_version) for i in range(len(lcm.labels_valid))])

    # death times
    ds['time_death'][ds['valid_flag']] = np.array([utime2Datetime(getattr(lcm.toocan[lcm.toocan_index_of_label[lcm.labels_valid[i]]],utime_end_key),toocan_version) for i in range(len(lcm.labels_valid))])

    # durations
    ds['duration'][ds['valid_flag']] = np.array([getattr(lcm.toocan[lcm.toocan_index_of_label[lcm.labels_valid[i]]],duration_key) for i in range(len(lcm.labels_valid))])

    # areas
    ds['area'][ds['valid_flag']] = formatLifecycleData(area_key,lcm,N_ages=ds['area'].shape[1])
    
    
def utime2Datetime(Utime,toocan_version):
    """Convert Utime (as in TOOCAN objects) into datetime.datetime objects"""
    
    if toocan_version == 'v2.07':
        
        days = int(str(Utime).split('.')[0])
        n_steps = int(str(Utime).split('.')[1])
    
        date = np.datetime64(dt.datetime(1970,1,1) + dt.timedelta(seconds=days*86400+n_steps*30*60))
        
    elif toocan_version == 'v2.08':
        
        date = np.datetime64(dt.datetime(1970,1,1) + dt.timedelta(seconds=Utime))
    
    return date

def datetime2Utime(date,toocan_version):
    """Convert datetime.datetime object into Utime (as in TOOCAN objects)"""
    
    delta_t = pd.to_datetime(date) - dt.datetime(1970,1,1)
    n_days = delta_t.days
    n_steps = int(delta_t.seconds/60/30)

    if toocan_version == 'v2.07': 
        
        utime = n_days+n_steps/100
        
    elif toocan_version == 'v2.08':
        
        utime = n_days*24*3600 + n_steps*30*60
    
    return utime

def utime2i_t(utime,toocan_version): # NOW REDUNDANT with LCM object attribute
    """Convert utime value (i.e. TOOCAN object attribute) into time index since beginning"""
    
    utime_v207_0 = 17014.03
    
    if toocan_version == 'v2.07':
        
        delta_utime = np.round(utime-utime_v207_0,2)
        n_days = int(str(delta_utime).split('.')[0])
        n_steps = int(str(delta_utime).split('.')[1])
    
        i_t = n_days*48 + n_steps
        
    elif toocan_version == 'v2.08':
        
        # reference time
        utime_v208_0 = datetime2Utime(utime2Datetime(utime_v207_0,
                                                 toocan_version='v2.07'),
                                  toocan_version='v2.08')
        
        # calculate
        delta_utime = utime-utime_v208_0
        i_t = int(delta_utime/30/60)
        
    else:
        print('TOOCAN version not recognized')
        
    return i_t

def formatLifecycleData(TOOCAN_attr,lcm,N_ages):
    """Convert variables stored along lifecycle in TOOCAN and format them as matrix to feed in output file"""
    
    # get values
    var_list = [getattr(lcm.toocan[lcm.toocan_index_of_label[lcm.labels_valid[i]]].clusters,TOOCAN_attr) for i in range(len(lcm.labels_valid))]
    # get length
    v_length = max(len(v) for v in var_list)
    # pad the different sizes of lists with nans
    var = np.vstack([v + [np.nan] * (max(N_ages,v_length) - len(v)) for v in var_list])
    
    return var
    
    
##-- Functions for time iteration
    
def loadDataAtSlice(i_t,lcm,lon_slice,lat_slice):
    """Load relevant fields : TOOCAN segmentation mask, precipitation values"""
    
    # date
    date = lcm.relation_table.loc[i_t].str_code

    # load segmentation mask for that date
    segmask = np.squeeze(loadTOOCANSeg(i_t,lcm.relation_table)).sel(latitude=lat_slice).data

    # load precipitation data
    prec_data = loadPrec(i_t,lcm.relation_table).sel(lat=lat_slice)
    prec = prec_data.data

    # coordinates
    longitude_2D,latitude_2D = getCoords2D(prec_data,lon_slice,lat_slice)
    del prec_data

    # compute all MCS ages (2D)
    MCS_age = lcm.computeAgesFromSegMask(i_t,segmask,metric='age')
    
    # age indices (2D)
    i_age_t = np.asarray(MCS_age/lcm.timestep,dtype=int)
    del MCS_age
    
    return date, segmask, prec, longitude_2D, latitude_2D, i_age_t
    
# def computeMaxPrecAtSlice(segmask, prec, longitude_2D, latitude_2D, i_age_t):
#     """Compute the maximum precipitation and its coordinates within each MCS label.
    
#     Returns:
#     - df_prec_max_label_ordered: pandas dataframe, with labels of MCS at t, prec_max, coordinates
#     """
    
#     ##-- extract maximum precipitation in mask at time t (spatial maximum), 
#     #    using only matrix manipulations and reordering, for better efficiency

#     # Create a DataFrame
#     df = pd.DataFrame({'label': segmask.flatten(),
#             'i_age': i_age_t.flatten(),
#             'prec': prec.flatten(),
#             'longitude': longitude_2D.flatten(),
#             'latitude': latitude_2D.flatten(),
#             })

#     #- remove nan labels
#     df_clean = df.dropna()

#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")

#         # replace all negative precip by zeros
#         df_clean['prec'][df_clean['prec'] < 0 ] = 0

#     # order array wrt. precip
#     df_ordered_prec = df_clean.sort_values(by=['label','prec'])

#     # find indices of last occurrence of each label
#     df_ordered_prec['label_diff'] = np.concatenate([np.diff(df_ordered_prec['label']),[1]])

#     # retain only the rows with max precip values (last in ordered values)
#     df_prec_max_label = df_ordered_prec[df_ordered_prec['label_diff'] > 0]

#     # reorder array wrt. label (so that, when extracting indices of labels, corresponding prec values are not swapped)
#     df_prec_max_label_ordered = df_prec_max_label.sort_values(by='label')
    
#     # remove rows with wrong ages ; these labels do not correctly appear in TOOCAN object
#     k_issue = np.where(df_prec_max_label_ordered['i_age'] < 0)[0]
#     df_prec_max_label_ordered.drop([df_prec_max_label_ordered.index[k] for k in k_issue],inplace = True)

#     # fix type changes, floats --> int
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")

#         # make some values integers
#         df_prec_max_label_ordered['label'] = df_prec_max_label_ordered['label'].astype(int)
#         df_prec_max_label_ordered['i_age'] = df_prec_max_label_ordered['i_age'].astype(int)
    
#     return df_prec_max_label_ordered

def computeMaxVarAtSlice(var,varid,segmask, longitude_2D, latitude_2D, i_age_t, i_t): 
    # Works similarly as following function, but is slower. Does it use less memory?
    """Compute the maximum value and its coordinates within each MCS label.
    
    Arguments:
    - var: variable,
    - varid: variable name.
    
    Returns:
    - df_var_max: dataframe of label-wise maximum variable values and coordinates
    """
    
    # Create a DataFrame
    df = pd.DataFrame({'label': segmask.flatten(),
            'i_age': i_age_t.flatten(),
            varid: var.flatten(),
            'longitude': longitude_2D.flatten(),
            'latitude': latitude_2D.flatten(),
            })
    
    #- remove nan labels
    df_clean = df.dropna()
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # replace all negative precip by zeros
        df_clean[varid][df_clean[varid] < 0 ] = 0
    
    #- iterate over labels and store max value
    all_labels = np.unique(df_clean.label)
    all_labels = np.array(all_labels[~np.isnan(all_labels)],dtype=int)

    df_list = []

    # for label in tqdm(all_labels):
    for label in all_labels:

        # (sub)frame at label
        df_ilab = df_clean.loc[df_clean.label==label]
        # extract max var
        var_max = np.max(df_ilab[varid])
        # save line of max
        df_list.append(df_ilab.loc[df_ilab[varid]==var_max])
        
    df_var_max = pd.concat(df_list)
    
    # remove duplicate labels
    df_var_max = df_var_max.drop_duplicates(keep='first',subset='label')
    
    # add current time index as attribute
    df_var_max.attrs['i_t'] = i_t
    
    return df_var_max
    
    
def assignVarMaxToOutput(ds,varnames,lcm,df):
    
    #- compute indices of ordered valid labels (where to assign prec_max) in labels_valid
    
    print('store df data at time index %d'%(df.attrs['i_t']))
    
    # a = lcm.labels_toocan
    a = lcm.labels_valid
    b = df['label'].values
    
    # show labels that are in df_max but not in lcm.labels_toocan
    ind_intruders = np.where(~np.isin(b,a))[0]
    label_intruders = b[ind_intruders] # just out of curiosity
    if len(label_intruders) > 0:
        print('intruder labels (present in df but not in TOOCAN object):',label_intruders)

    # ignore intruder labels (not found in b)
    ind_common = np.where(np.isin(b,a))[0]
    b = b[ind_common]
    df = df.iloc[ind_common]

    # get indices of labels at t
    ind_lab_t = np.argwhere(np.isin(a, b))[:,0]

    # get indices of ages at t
    ind_age_t = df['i_age'].values
    
    # create coordinates
    coords = np.vstack([ind_lab_t,ind_age_t])
    assert len(a) == ds['prec_max'].shape[0], "number of labels differ"
    flat_indices = np.ravel_multi_index(coords, ds['prec_max'].shape)

    #- assign in output xarray
    for key in varnames:
        
        try:
            
            # assign
            np.put(ds['%s_max'%key].values,flat_indices,df[key])

        except IndexError:
            
            print('IndexError ! ')
            print('occured during assignment assignVarMaxToOutput, for rel_table time index %d'%(df.attrs['i_t']))
            print('for key %s'%(key))
            print('labels :')
            print(a[ind_age_t])
            print('ind_lab_t :')
            print(ind_lab_t)
            print('ind_age_t:')
            print(ind_age_t)
            
            continue
    
    
def iteration(i_t):
    """load data, compute and return results at time step i_t. """
    
    print("time step #%d ; %s seconds "% (i_t,time.time() - start_time))
    
    # load data at time step t
    date, segmask, prec, longitude_2D, latitude_2D, i_age_t = loadDataAtSlice(i_t,lcm,lon_slice,lat_slice)

    # compute variables at maximum precipitation in each MCS segmentation mask
    # df_max = computeMaxPrecAtSlice(segmask, prec, longitude_2D, latitude_2D, i_age_t)
    df_max = computeMaxVarAtSlice(prec,'prec',segmask, longitude_2D, latitude_2D, i_age_t, i_t)
    
    del prec, segmask, longitude_2D, latitude_2D, i_age_t

    return df_max
    
    
##--- Main script

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Compute maximum precipitation in anvil lifecycle, and save as netcdf file.')
    parser.add_argument('-i','--index', type=int, nargs=1, default=[1],
                    help='index in relation table DYAMOND-SAM:segmentationTOOCAN')
    parser.add_argument('-n','--ninds', type=int, nargs=1, default=[1916],
                    help='number of indices in a row to analyze')
    parser.add_argument('-r','--region',type=str,default='tropics',
                    help='region of analysis')
    parser.add_argument('--n_ages',type=int,default=96,
                        help='number of time steps in age coordinate')
    parser.add_argument('--dt',type=float,default=0.5,
                        help='model time step (hours)')
    parser.add_argument('--n_proc',type=int,default=1,
                        help='number of simultaneous processes for parallelization')

    args = parser.parse_args()

    #-- arguments
    i_t0 = args.index[0]
    Ni = args.ninds[0]
    region = args.region
    N_ages = args.n_ages
    timestep = args.dt
    n_proc = args.n_proc
    
    # output directory
    out_dir = os.path.join(DIR_OUT,region,'life_cycle_mapping')
    os.makedirs(out_dir,exist_ok=True)

    # geographical parameters
    if region == 'tropics':
        lon_slice = slice(0,360)
        lat_slice = slice(-30,30)
        
    print("--- %s seconds ---" % (time.time() - start_time))
    print('> create LifeCycleMapping object')
    
    
    #- define LifeCycleMappingObject
    
    lcm_obj_name = 'lcm_%s.pickle'%toocan_version_str
    lcm_path = os.path.join(out_dir,lcm_obj_name)

    if len(glob.glob(lcm_path)) > 0: # then already exists
        
        # load
        lcm = pickle.load(open(lcm_path,'rb'))
        
    else:
        
        # initialize
        lcm = initLifeCycleMappingObject()

        # save
        pickle.dump(lcm,open(lcm_path,'wb'))

        
    #- define coordinates
    
    # define label coordinate
    # label_toocan = lcm.labels_toocan
    label_toocan = lcm.labels_valid  # !!! ignore duplicate labels
    # define age
    age = np.arange(N_ages)*lcm.timestep # in hours
    
    print("--- %s seconds ---" % (time.time() - start_time))
    print('> initialize output array')
    
    #- Initalize output
    ds = initOutputArray(label_toocan,age)
    
    #- Fill output array with toocan characteristics
    fillWithTrackingData(ds,label_toocan,lcm)
    
    print("--- %s seconds ---" % (time.time() - start_time))
    print('> start time loop')    
    
    #-- time iteration
    
    i_t_all = range(i_t0,i_t0+Ni)
    
    #- without multiprocessing
    # for i_t in i_t_all:
        # iteration(ds,i_t)
    
    #- with multiprocessing   
    # These two steps (loading and computing) are the most time consuming, hence the (small) parallelization
    
    # get the current start method
    method = multiprocessing.get_start_method()
    print('multiprocessing default start method is',method)
    
    # Create a multiprocessing Pool with n_proc processes
    pool = Pool(processes=n_proc)
    
    # compute
    df_all = pool.map(iteration,i_t_all)
    
    # Close the pool
    pool.close()
    pool.join()

    # assign
            
    for df in df_all:
            
        # vars to assign at location of MCS-wise max precipitation
        varnames = ['prec','longitude','latitude']
        
        # assign
        assignVarMaxToOutput(ds,varnames,lcm,df)
        
        # free up space
        gc.collect()
            
    print("Loop stopped ; %d seconds "% (time.time() - start_time))
    
    #-- save output xarray
    outfile_root = os.path.join(out_dir,'lifecycle_diag_DYAMOND_Summer_SAM_%s'%toocan_version_str)
    N_file = len(glob.glob(outfile_root+'_*.nc'))
    outfile_name = outfile_root+'_%d.nc'%N_file
    print('Save output to %s'%outfile_name)
    ds.to_netcdf(outfile_name)

    print("--- %s seconds ---" % (time.time() - start_time))
    print('Script ended successfully.')
        
    sys.exit(0)
        
        
        