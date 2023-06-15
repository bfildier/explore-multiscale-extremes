import pickle
import os
import numpy as np

# load own modules
from conditionalstats_chunked import *


class CaseStudy():
    
    """Documentation for class CaseStudy

    """
        
    ##-- Class constructor
    
    def __init__(self,name,region,rel_tab_dyam_seg=None):
        
        self.name = name
        self.region = region
        self.rel_tab_dyam_seg = rel_tab_dyam_seg
        
    def __repr__(self):
        """Creates a printable version of the Distribution object. Only prints the 
        attribute value when its string fits is small enough."""

        out = '< DistributionChunked object:\n'
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
        
        self.i_t_min = i_t_min
        self.i_t_max = i_t_max
        self.range_t = range(i_t_min,i_t_max+1)
        self.lat_slice = lat_slice
        self.lon_slice = lon_slice
        
    def setDirectories(self,dir_dist_sliced):
        
        self.dir_dist_sliced = dir_dist_sliced
    
    def storeTimes(self):
        
        # starting date
        date_ref = dt.datetime(year=2016,month=8,day=1)
        
        def getDateCurrent(i_t):
    
            # number of seconds elapsed
            dt_dyam = int(self.reltab_dyam_seg.iloc[i_t]['path_dyamond'].split('_')[-1])
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
        
    def loadDistPrSliced(self):
        
        self.dict_dist_pr_sliced = {}

        for i_t in self.range_t:

            # load
            dist_pr_t = pickle.load(open(os.path.join(self.dir_dist_sliced,'dist_pr_t_%d.pickle'%i_t),'rb'))
            # store
            self.dict_dist_pr_sliced[i_t] = dist_pr_t
    
    def computePrMean(self):
        
        self.pr_mean = np.array([self.dict_dist_pr_sliced[i_t].mean for i_t in self.range_t])
        
    def findTimeIndToIgnore(self):
        # Get indices where precip data is wrong
        #
        # choose to call k_t indices in range [0,1085] and i_t indices in range [832,1917]
        #
    
        k_t_skip = np.where(self.pr_mean < 0.01)[0] # find where precip is too small
        k_t_skip = np.hstack([k_t_skip-1,k_t_skip]) # also flag the previous index
        k_t_skip.sort() # sort
        
        # store
        self.times_to_ignore = k_t_skip + self.i_t_min
        
    def combinePrDistributions(self):
        
        list_dist = list(self.dict_dist_pr_sliced.values())
        inds_to_ignore = list(self.times_to_ignore - self.i_t_min)

        #-- initiate global distribution
        self.dist_pr = DistributionChunked(name='pr, %s %s, all times'%(self.name,self.region),
                                          dist_chunks=list_dist,
                                          chunks_to_ignore=inds_to_ignore,
                                          bintype='invlogQ',nd=8)

        # compute ranks
        self.dist_pr.getInvLogRanks()
        # compute all percentiles with new algorithm
        self.dist_pr.computeDistribution()
    
            
    