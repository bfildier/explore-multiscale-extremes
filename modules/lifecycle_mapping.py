"""@package lifecyclemapping
Documentation for module mapping.

Class LifeCycleMapping defines methods to establish correspondence between a statistical distribution from
class Distribution, and (Lagrangian) MCS objects.
"""
import numpy as np
from math import log10,ceil,floor,exp
import time
import sys

from datetime import datetime as dt
from datetime import timedelta

class LifeCycleMapping():
    
    ###--- Class constructor
    
    def __init__(self,dist,relation_table,toocan):

        """Constructor for class LifeCycleMapping.
        Arguments:
        - dist: object of class Distribution.
        - relation_table: data frame mapping simulation files with TOOCAN segmentation masks.
        - toocan: list of TOOCAN objects.
        """
        
        self.dist = dist
        self.relation_table = relation_table
        self.toocan = toocan
        self.N_MCS = len(self.toocan)
        
        # store list of toocan attributes to speed up later calculation
        self.labels_toocan = [self.toocan[i].label for i in range(self.N_MCS)]
        
        # store list of toocan labels that appear for several MCSs
        self.findDuplicateLabels()
        self.findValidLabels()
        
        # # store list of indices accessible from labels --- TAKES TOO LONG
        # self.defineMappingLabelToIndex()
        
        
    def findDuplicateLabels(self):
        """Save TOOCAN labels appearing several times"""
        
        labels_unique, labels_counts = np.unique(self.labels_toocan,return_counts=True)
        ind = np.where(labels_counts > 1)[0]
        # store labels
        self.labels_duplicate = [labels_unique[i] for i in ind]
        
    def findValidLabels(self):
        """Save TOOCAN labels appearing only once"""
        
        # count label occurrence in TOOCAN
        labels_unique, labels_counts = np.unique(self.labels_toocan,return_counts=True)
        
        # valid labels are those appearing only once
        ind_unique = np.where(labels_counts == 1)[0]
        self.labels_valid = [labels_unique[i] for i in ind_unique]
        
    def defineMappingLabelToIndex(self):
        """Define array of indices arr such that for any label l, arr[l] is the index of l in toocan list.
        All other indices are nans; and duplicate labels in toocan are ignored.
        """
        
        # init
        self.toocan_index_of_label = np.full(np.nanmax(self.labels_toocan)+1,np.nan)
        
        # fill 
        for i_toocan in range(len(self.labels_toocan)):
            
            label = self.labels_toocan[i_toocan]
            
            if label in self.labels_valid:
                
                self.toocan_index_of_label[label] = i_toocan
        
        # Old method, to fix
        
#         # 0. initiate with nans
#         inds_labels_all = np.array(np.arange(0,np.nanmax(self.labels_toocan)+1),dtype=float)
#         inds_labels_all[0] = np.nan

#         # 1. Nans
#         # mask of valid label indices
#         mask_labels_valid = np.in1d(inds_labels_all,self.labels_valid)
#         # replace invalid with nans
#         inds_labels_all[~mask_labels_valid] = np.nan

#         # 2. Indices
#         # mask of valid toocan labels
#         mask_toocan_labels_valid = np.in1d(self.labels_toocan,self.labels_valid)
#         # find indices of valid toocan labels
#         ind_toocan_valid = np.where(mask_toocan_labels_valid)[0]
#         # assign 
#         inds_labels_all[mask_labels_valid] = ind_toocan_valid
#         # convert array to integer (nans become large negative numbers)
#         self.toocan_index_of_label = np.asarray(inds_labels_all,dtype=int)
        
        
    ###--- class Methods
    
    #--- i_t,label --> MCS age
        
    def indexOfLabel(self,label):
        """Returns index of MCS in TOOCAN list with corresponding label"""
    
        i_MCS = np.where(np.array(self.labels_toocan) == label)[0][0]
        
        return i_MCS
        
    def timeIndex2Timedelta(self,i_MCS,j_t_MCS):
    
        date_str = self.toocan[i_MCS].clusters.Utime[j_t_MCS]
        date_day = int(date_str)
        date_30mn = int(str(date_str).split('.')[-1])
        # compute time delta
        td = timedelta(days = int(date_str),seconds = date_30mn*30*60)
        
        # return
        return td
    
    def UtimeToRelTableIndex(self,Utime):
        """Convert Utime (in TOOCAN object for instance) to index in relation table."""
        
        i_t = np.where(self.relation_table.UTC == Utime)[0][0]
        
        return i_t
    
    def getAgeMCS(self,i_t,label,time_MCS='current'):
        """Compute the age of a MCS at a given time.

        i_t, label ---> 1 MCS age
        
        Arguments:
        - i_t: time index in simulation
        - label: MCS label
        - time_MCS: when to compute age in MCS lifetime (current --> i_t, end --> final MCS age)

        Returns:
        - age in hours"""

        # index in TOOCAN list
        i_MCS = self.indexOfLabel(label)
        # index of slice in TOOCAN lifecycle for current time
        if time_MCS == 'current':
            j_t_MCS = np.where(self.toocan[i_MCS].clusters.Utime == self.relation_table.loc[i_t].UTC)[0][0]
        elif time_MCS == 'end':
            j_t_MCS = -1
        # birth time
        time_birth = self.timeIndex2Timedelta(i_MCS,0)
        # current time
        time_now = self.timeIndex2Timedelta(i_MCS,j_t_MCS)
        # age delta
        age_delta = time_now - time_birth

        return age_delta.total_seconds()/3600
        
    
    def computeMCSAgeMetrics(self,i_t,label,segmask,metric='age'):
        """Compute age, normalized age, etc."""
        
        # index in list
        i_MCS = self.indexOfLabel(label)

        # check if current label occurs in current segmentation mask
        time_match = bool(np.any(segmask == self.toocan[i_MCS].label).data)

        if time_match:

            # turn it into system's age
            age_MCS_at_t = self.getAgeMCS(i_t,label) # in hours
            
            if metric == 'age':

                # to return
                out = age_MCS_at_t

            elif metric == 'norm_age':
                
                # total lifetime
                lifetime_MCS = self.getAgeMCS(i_MCS,label,'end') # in hours

                # normalized age
                norm_age_MCS_at_t = age_MCS_at_t/lifetime_MCS
                
                # to return 
                out = norm_age_MCS_at_t
                
            return out
                
    #--- i_t, segmentation mask --> mask of MCS ages
    
    def getLabelsInSegmentationMask(self,segmask):
        """Returns a list of labels appearing in segmentation mask"""

        seg_1D = segmask.values.flatten()
        seg_1D_nonans = seg_1D[~np.isnan(seg_1D)]

        return np.unique(np.array(seg_1D_nonans,dtype=int))
    
    def computeAgeMask(self,i_t,segmask,metric='age'):
        """Compute an array of MCS ages at MCS label locations in segmentation mask.
        
        i_t, segmentation mask ---> mask of MCS ages
        
        Arguments:
        - i_t: time index in relation table
        - seg_mask: TOOCAN segmentation mask at time i_t
        - metric: metric to compute (age, norm_age, etc.)
        """
    
        # all labels in segmask
        labels_in_segmask = self.getLabelsInSegmentationMask(segmask)
        
        #- compute sequence of ages in segmentation mask
        ages = {}
        for label in labels_in_segmask:
            
            if label in self.labels_valid:
            
                ages[label] = self.computeMCSAgeMetrics(i_t,label,segmask,metric=metric)
            
        #- create mask of ages
        age_mask = np.full(segmask.shape,np.nan)
        # fill
        for label in ages.keys():
            
            age_mask[segmask == label] = ages[label]
        
        return age_mask
    
    def sampleToMeanMCSAge(self,i_t,sample,normalize=False):
        """Compute an array of MCS ages from a sample for each percntile.
        
        i_t, precip values ---> MCS age composited onto precip distribution
        
        Arguments:
        - i_t: time index in relation table
        - sample: values at time i_t for the variable used to compute distribution
        """
        
        pass
    
    
        
        