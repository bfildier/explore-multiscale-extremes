import pickle
import os, sys, glob
import tarfile
import shutil

import numpy as np
import datetime as dt

# load own modules
from conditionalstats_chunked import *
# load project functions
from settings import *


class CaseStudy():
    
    """Documentation for class CaseStudy
    """
        
    ##-- Class constructor
    
    def __init__(self,name,region,mask='all',model='SAM',rel_tab_dyam_seg=None):
        
        self.name = name
        self.region = region
        self.model = model
        self.rel_tab_dyam_seg = rel_tab_dyam_seg
        
    def __repr__(self):
        """Creates a printable version of the Distribution object. Only prints the 
        attribute value when its string fits is small enough."""

        out = '< CaseStudy object:\n'
        # print keys
        for k in self.__dict__.keys():
            out = out+' . %s: '%k
            if k in ['dist_chunks','chunks_to_ignore']:
                # show type
                out = out+'%s\n'%str(getattr(self,k).__class__)
            else:
                if len(str(getattr(self,k))) < 80:
                    # show value
                    out = out+'%s\n'%str(getattr(self,k))
                else:
                    # show type
                    out = out+'%s\n'%getattr(self,k).__class__
        out = out+' >'

        return out
        
    def setSimulationSpecs(self,i_t_min,i_t_max,lat_slice=None,lon_slice=None):
        """Set time bounds and space bounds for analysis. Could be specified in a yaml file instead"""
        
        self.i_t_min = i_t_min
        self.i_t_max = i_t_max
        self.range_t = range(i_t_min,i_t_max+1)
        
        # latitude
        if lat_slice is None:
            coord_slices['lat'][self.region] # ! adapt here for other regions than tropics
        else:
            self.lat_slice = lat_slice

        # longitude
        if lon_slice is None:
            coord_slices['lon'][self.region] # ! adapt here for other regions than tropics
        else:
            self.lon_slice = lon_slice
        
    def getVaridStr(self,varid):
        
        varid_str = varid
        if varid == 'Prec':
            varid_str = 'pr'
            
        return varid_str
        
    def setDirectories(self,varid,mask):
        
        dirname = os.path.join(DIR_OUT,self.region,varid,mask)
        
        if varid == 'Prec':
            
            setattr(self,'dir_dist_pr_%s_sliced'%mask,dirname)
        
        else:
            
            setattr(self,'dir_dist_%s_%s_sliced'%(varid,mask),dirname)
    
    def storeTimes(self):
        
        # starting date
        date_ref = dt.datetime(year=2016,month=8,day=1)
        
        def getDateCurrent(i_t):
    
            # print(i_t,end='..')
        
            # number of seconds elapsed
            dt_dyam = int(self.rel_tab_dyam_seg.iloc[i_t]['path_dyamond'].split('_')[-1])
            # corresponding delta
            delta_t = dt.timedelta(seconds=int(dt_dyam*7.5))
            # current date
            date_current = date_ref + delta_t

            return date_current
        
        # time dimension
        self.time = np.array([getDateCurrent(i_t) for i_t in self.range_t])
        
    def getDistIndexFromTimeIndex(self,i_t):
        
        return i_t - self.i_t_min
    
    def getTimeIndexFromDistIndex(self,i_d):
        
        return i_d + self.i_t_min
    
    def loadDistSliced(self,varid,mask):
        
        varid_str = self.getVaridStr(varid)
        
        # open file
        file = tarfile.open(os.path.join(getattr(self,'dir_dist_%s_%s_sliced'%(varid_str,mask)),'time_slices.tar.gz'))
        # extract file to temporary folder
        file.extractall(DIR_TEMPDATA)
        
        # print(glob.glob(DIR_TEMPDATA+'/*'))
        
        # initialize
        dict_dist = {}

        # do operations on files
        for i_t in self.range_t:
            
            # load
            if len(glob.glob(os.path.join(DIR_TEMPDATA,'time_slices'))) > 0:
                dict_dist[i_t] = pickle.load(open(os.path.join(DIR_TEMPDATA,'time_slices','dist_%s_t_%d.pickle'%(varid_str,i_t)),'rb'))
            elif len(glob.glob(os.path.join(DIR_TEMPDATA,mask,'time_slices'))) > 0:
                dict_dist[i_t] = pickle.load(open(os.path.join(DIR_TEMPDATA,mask,'time_slices','dist_%s_t_%d.pickle'%(varid_str,i_t)),'rb'))
            else:
                dict_dist[i_t] = pickle.load(open(os.path.join(DIR_TEMPDATA,self.region,varid,mask,'time_slices','dist_%s_t_%d.pickle'%(varid_str,i_t)),'rb'))
                   
        # remove files from temporary folder
        shutil.rmtree(glob.glob(os.path.join(DIR_TEMPDATA,'*'))[0])
        
        # store dictionary
        setattr(self,'dict_dist_%s_%s_sliced'%(varid_str,mask),dict_dist)
    
    def loadDistPrSliced(self,mask):
        
        self.loadDistSliced(varid='Prec',mask=mask)
        
    def loadDistMerged(self,varid,mask):
        
        varid_str = self.getVaridStr(varid)
        dir_load = 'dir_dist_%s_%s_sliced'%(varid_str,mask)
        
        # load
        dist = pickle.load(open(os.path.join(dir_load,'dist_%s.pickle'%varid_str),'rb'))
        
        # store for later use
        dist_name = 'dist_%s_%s'%(varid_str,mask)
        setattr(self,dist_name,dist)
        
        
    def computeMean(self,varid,mask='all'):
        
        varid_str = self.getVaridStr(varid)
            
        setattr(self,'%s_%s_mean'%(varid_str,mask),np.array([getattr(self,'dict_dist_%s_%s_sliced'%(varid_str,mask))[i_t].mean for i_t in self.range_t]))
        
    def findTimeIndToIgnore(self):
        # Get indices where precip data is wrong
        #
        # choose to call k_t indices in range [0,1085] and i_t indices in range [832,1917]
        #
    
        k_t_skip = np.where(self.pr_all_mean < 0.01)[0] # find where precip is too small
        k_t_skip = np.hstack([k_t_skip-1,k_t_skip]) # also flag the previous index
        k_t_skip.sort() # sort
        
        # store
        self.times_to_ignore = k_t_skip + self.i_t_min
        
    def combineDistributions(self,varid,mask):
        
        varid_str = self.getVaridStr(varid)
        
        list_dist = list(getattr(self,'dict_dist_%s_%s_sliced'%(varid_str,mask)).values())
        inds_to_ignore = list(self.times_to_ignore - self.i_t_min)

        #-- initialize global distribution
        dist = DistributionChunked(name='%s, %s, %s %s, all times'%(varid_str,mask,self.name,self.region),
                                          dist_chunks=list_dist,
                                          chunks_to_ignore=inds_to_ignore,
                                          bintype='invlogQ',nd=8)

        # compute ranks
        dist.getInvLogRanks()
        # compute all percentiles with new algorithm
        dist.computeDistribution()
        
        # save
        dist_name = 'dist_%s_%s'%(varid_str,mask)
        setattr(self,dist_name,dist)
        
    
    def loadMcsAgeDiagnostic(self):
        
        pass

    
    def combineMcsAgeDiagnostic(self):
        
        pass
    
#    WRITE HERE function that loads MCS age diagnostics and save them in each distribution at time t
#    merge all times to get full mapping of MCS age onto distribution and save it in cs.dist_XX (DistributionChunked)
#
#    then test it in explore_lifecycle_mapping.ipynb
    