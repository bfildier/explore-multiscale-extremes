"""@package lifecyclemapping
Documentation for module mapping.

Class LifeCycleMapping defines methods to establish correspondence between a statistical distribution from
class Distribution, and (Lagrangian) MCS objects.
"""
import numpy as np
from math import log10,ceil,floor,exp
import time
import sys
import re
import warnings

import datetime as dt

class LifeCycleMapping():
    
    ###--- Class constructor
    
    def __init__(self,relation_table,toocan_timetable,toocan,timestep=0.5):

        """Constructor for class LifeCycleMapping.
        
        Arguments:
        - timestep: model time step in hours
        - relation_table: data frame mapping simulation files with TOOCAN segmentation masks.
        - toocan_timetable: table for mapping TOOCAN objects and their birth and death time step and duration
        - toocan: list of TOOCAN objects.
        """
        
        # self.dist = dist
        self.timestep = timestep # in hours
        self.relation_table = relation_table
        self.toocan_timetable = toocan_timetable
        self.toocan = toocan
        self.N_MCS = len(self.toocan) if self.toocan is not None else 0
        
        # store list of toocan attributes to speed up later calculation
        self.labels_toocan = [self.toocan[i].label for i in range(self.N_MCS)] if self.toocan is not None else []
        
        # store list of toocan labels that appear for several MCSs
        self.findDuplicateLabels()
        self.findValidLabels()
        
        # # store list of indices accessible from labels --- TAKES TOO LONG
        # self.defineMappingLabelToIndex()
        
        # store indices i_t for birth and death and 
        
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
        
        # slow method, but that works
        
#         # init
#         self.toocan_index_of_label = np.full(np.nanmax(self.labels_toocan)+1,np.nan)
        
#         # fill 
#         for i_toocan in range(len(self.labels_toocan)):
            
#             label = self.labels_toocan[i_toocan]
            
#             if label in self.labels_valid:
                
#                 self.toocan_index_of_label[label] = i_toocan
        
#         # Old method (vector-based), to fix
        
        # 0. initiate with nans
        inds_labels_all = np.array(np.arange(0,np.nanmax(self.labels_toocan)+1),dtype=float)
        inds_labels_all[0] = np.nan

        # 1. Nans
        # a. mask of valid label indices
        mask_labels_valid = np.in1d(inds_labels_all,self.labels_valid)
        # b. replace invalid with nans
        inds_labels_all[~mask_labels_valid] = np.nan

        # 2. Indices
        # a. mask of valid toocan labels
        mask_toocan_labels_valid = np.in1d(self.labels_toocan,self.labels_valid)
        # b. extract order of valid toocan labels
        toocan_labels_valid_unordered = np.array(self.labels_toocan)[mask_toocan_labels_valid]
        order_toocan_valid = sorted(range(len(toocan_labels_valid_unordered)), key=lambda k: toocan_labels_valid_unordered[k])

        # find indices of valid toocan labels
        ind_toocan_valid = np.where(mask_toocan_labels_valid)[0]
        # assign them, in order as determined in 2.b
        inds_labels_all[mask_labels_valid] = np.take(np.array(ind_toocan_valid,dtype=int),order_toocan_valid)
        # convert array to integer (nans become large negative numbers)
        self.toocan_index_of_label = np.asarray(inds_labels_all,dtype=int)
        
        
# ## TO DEBUG method LifeCycleMapping.defineMappingLabelToIndex()

# ex_labels_toocan = [1,3,4,5,7,8,9,10]
# # max_lab_toocan = 
# ex_labels_valid = [3,5,7,8,9]

# # 0. initiate with nans
# # inds_labels_all = np.array(np.arange(0,np.nanmax(self.labels_toocan)+1),dtype=float)
# inds_labels_all = np.array(np.arange(0,np.nanmax(ex_labels_toocan)+1),dtype=float)
# print('inds_labels_all:',inds_labels_all)
# inds_labels_all[0] = np.nan
# print('inds_labels_all:',inds_labels_all)

# # 1. Nans
# # mask of valid label indices
# # mask_labels_valid = np.in1d(inds_labels_all,self.labels_valid)
# mask_labels_valid = np.in1d(inds_labels_all,ex_labels_valid)
# print('mask_labels_valid:',mask_labels_valid)

# # replace invalid with nans
# inds_labels_all[~mask_labels_valid] = np.nan
# print('inds_labels_all:',inds_labels_all)

# # 2. Indices
# # mask of valid toocan labels
# # mask_toocan_labels_valid = np.in1d(self.labels_toocan,self.labels_valid)
# mask_toocan_labels_valid = np.in1d(ex_labels_toocan,ex_labels_valid)
# print('mask_toocan_labels_valid:',mask_toocan_labels_valid)

# # find indices of valid toocan labels
# ind_toocan_valid = np.where(mask_toocan_labels_valid)[0]
# print('ind_toocan_valid:',ind_toocan_valid)
# # assign 
# inds_labels_all[mask_labels_valid] = np.array(ind_toocan_valid,dtype=int)
# print('inds_labels_all:',inds_labels_all)
# inds_labels_all = np.asarray(inds_labels_all,dtype=int)
# # inds_labels_all[inds_labels_all < 0] = np.nan
# print('inds_labels_all:',inds_labels_all)

# for lab in ex_labels_valid:
    
#     print(lab, inds_labels_all[lab])
    
    
        
    def buildTimeTable(self):
        """Construct a time table where rank i corresponds to MCS label i, and each column is a time metric:
        i_t_min, i_t_max (indices in DYAMOND relation table), duration. Stores in self.time_table.
        
        Ignores systems with identical labels.
        
        Takes 2mn46 for DYAMOND 1 - SAM, 347106 systems."""

        # initialize
        self.time_table = pd.DataFrame(columns=['label','i_t_min','i_t_max','duration'], index=np.arange(np.nanmax(labels_toocan)))

        # fill where label is found
        for i_MCS in range(len(toocan)):

            MCS = self.toocan[i_MCS]

            # label
            label = MCS.label

            if label in self.labels_valid:
                
                # birth
                i_t_min = np.where(self.relation_table.UTC == MCS.Utime_Init)[0][0]

                # death
                i_t_max = np.where(self.relation_table.UTC == MCS.Utime_End)[0][0]

                # duration (equal (i_t_max-i_t_min+1)*0.5 hrs)
                duration = MCS.duration

                # save
                self.time_table.loc[label] = pd.Series({'label':label,'i_t_min':i_t_min,'i_t_max':i_t_max,'duration':duration})
            
        
    ###--- class Methods
    
    #--- i_t,label --> MCS age
        
    def indexOfLabel(self,label):
        """Returns index of MCS in TOOCAN list with corresponding label"""
    
        i_MCS = np.where(np.array(self.labels_toocan) == label)[0][0]
        
        return i_MCS
        
    def timeIndex2Timedelta(self,i_MCS,j_t_MCS):
    
        date_str = self.toocan[i_MCS].clusters.Utime[j_t_MCS]
        date_day = int(date_str)
        date_30mn = int((str(date_str).split('.')[-1]).ljust(2,'0')) #ljust
        # compute time delta
        td = dt.timedelta(days = int(date_str),seconds = date_30mn*30*60)
        
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
    
    def getLabelsInSegMask(self,segmask):
        """Returns a list of labels appearing in segmentation mask"""

        seg_1D = segmask.flatten()
        seg_1D_nonans = seg_1D[~np.isnan(seg_1D)]

        return np.unique(np.array(seg_1D_nonans,dtype=int))
    
    
    # NOT EFFICIENT VERSION: 
    # execution time: 
    # 2016 labels in current segmentation mask
    # CPU times: user 9min 26s, sys: 3.57 s, total: 9min 30s
    # Wall time: 9min 30s
    
    def computeAgesFromSegMask_old(self,i_t,segmask,metric='age'):
        """Compute an array of MCS ages at MCS label locations in segmentation mask.
        
        i_t, segmentation mask ---> mask of MCS ages
        
        Arguments:
        - i_t: time index in relation table
        - seg_mask: TOOCAN segmentation mask at time i_t
        - metric: metric to compute (age, norm_age, etc.)
        """
        
        # all labels in segmask
        labels_in_segmask = self.getLabelsInSegMask(segmask)
        print('%d labels in current segmentation mask'%len(labels_in_segmask))
        # flatten segmentation mask to fill
        segmask_flat = segmask.values.flatten()
        
        #- compute sequence of ages in segmentation mask
        ages = {}
        age_mask_1D = np.full(segmask_flat.shape,np.nan)
        
        # fill
        for label in labels_in_segmask:
            
            if label in self.labels_valid:
                
                # print(label,end='..')
            
                ages[label] = self.computeMCSAgeMetrics(i_t,label,segmask,metric=metric)
                
                # fill age mask            
                age_mask_1D[segmask_flat == label] = ages[label]

        print()
            
        # reshape to 2D 
        age_mask = np.reshape(age_mask_1D,segmask.shape)
        
        return age_mask
    
    
    def computeAgesFromSegMask(self,i_t,segmask,metric='age'):
        """Compute an array of MCS ages at MCS label locations in segmentation mask.
        
        i_t, segmentation mask ---> mask of MCS ages
        
        Arguments:
        - i_t: time index in relation table
        - segmask: TOOCAN segmentation mask at time i_t
        - metric: metric to compute (age, norm_age, etc.)
        """
    
        #- step 1. Mask of nans and flatten segmentation mask without nans

        # flatten segmask
        segmask_1D = segmask.flatten()
        # mask of nans
        mask_nans_segmask_1D = np.isnan(segmask_1D)
        # remove nans
        segmask_1D_valid = segmask_1D[~mask_nans_segmask_1D]
        
        #- step 2. Take time metrics at labels ## CAREFUL, indices must be equal to labels in toocan_timetable

        # births
        i_t_min_valid = np.array(self.toocan_timetable['i_t_min'].take(segmask_1D_valid))
        # deaths
        i_t_max_valid = np.array(self.toocan_timetable['i_t_max'].take(segmask_1D_valid))
        # durations
        durations_valid = np.array(self.toocan_timetable['duration'].take(segmask_1D_valid))

        if metric == 'age':
            # age
            metric_valid = (i_t-i_t_min_valid)*self.timestep
        elif metric == 'norm_age':
            # normalized age
            metric_valid = (i_t-i_t_min_valid)*self.timestep/durations_valid
        elif metric == 'duration':
            # duration
            metric_valid = durations_valid
        
        #- step 3. Insert in 1D arrays with nans

        metric_1D = np.full(segmask_1D.shape,np.nan)
        metric_1D[~mask_nans_segmask_1D] = metric_valid
        
        #- step 4. reshape

        metric = np.reshape(metric_1D, segmask.shape)
        
        #- return result
        
        return metric
    
    
    
    def compositeMcsAgeOnDist(self,i_t,segmask,sample,dist_var,xymask=True,mask_T='allT',diag='mean',metric='age'):
        """Compute an array of MCS ages from a sample for each percentile.
        
        i_t, precip values ---> MCS age composited onto precip distribution
        
        Arguments:
        - i_t: time index in relation table
        - segmask: TOOCAN segmentation mask at time i_t
        - sample: values at time i_t for the variable used to compute distribution
        - xymask: spatial mask on which to restrict calculation of diagnostic (e.g. 'land','ocean', etc.) fed as a flattened 1D array
        - mask_T: name of mask on TOOCAN MCS lifetime on which to restrict calculation of diagnostic (e.g. allT, min5hr, etc.)
        - diagnostic: how to combine the metric quantified in a given bin (mean, var, etc.)
        - metric: metric to compute (age, norm_age, etc.)
        """

        # 1. compute array of metric
        metric = self.computeAgesFromSegMask(i_t,segmask,metric).flatten()
        
        # 1'. compute mask of durations (from mask_T)
        duration = self.computeAgesFromSegMask(i_t,segmask,'duration').flatten()
        
        if mask_T == 'allT':
            
            Tmask = np.full(duration.size,True)
            
        elif re.compile('min.*hr').match(mask_T):
            
            # get Tmin from argument
            Tmin = int(mask_T[3:][:-2])
            # extract mask
            Tmask = np.greater(duration,int(Tmin/self.timestep))
            
        elif re.compile('btw.*and.*hr').match(mask_T):
            
            # get Tbounds from argument
            Tmin = int(mask_T.split('and')[0][3:])
            Tmax = int(mask_T.split('and')[1][:-2])
            # extract mask
            Tmask_min = np.greater(duration,int(Tmin/self.timestep))
            Tmask_max = np.logical_not(np.greater(duration,int(Tmax/self.timestep)))
            Tmask = np.logical_and(Tmask_min,Tmask_max)
        
        # 2. digitize sample values in distribution bins
        digits = dist_var.getBinIndices(sample) # (flattened)
        N_dig = dist_var.ranks.size
        
        # 3. for each bin i, apply diagnostic to the metric and save sample size
        # init
        diag_all_bins = np.full(N_dig,np.nan)
        N_all_bins = np.full(N_dig,np.nan)
        N_valid_all_bins = np.full(N_dig,np.nan)
        # loop
        for i_bin in range(0,N_dig):

            # x-y mask in bin i
            mask_bin = (digits == i_bin)
            # take intersection with spatial mask
            mask_bin = np.logical_and(mask_bin,xymask)
            # take intersection with MCS duration mask
            mask_bin = np.logical_and(mask_bin,Tmask)
            
            # extract age in bin
            age_in_bin = metric[mask_bin]

            # get bin sample size
            N_all_bins[i_bin] = age_in_bin.size
            N_valid_all_bins[i_bin] = np.sum(~np.isnan(age_in_bin))
            # if N_valid_all_bins[i_bin] == 0:
            #     print('bin #%d = all nans'%i_bin)
            
            # apply diagnostic
            if diag in ['mean','max','var','min']:
                
                if N_valid_all_bins[i_bin] > 0:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        diag_all_bins[i_bin] = getattr(np,'nan%s'%diag)(age_in_bin)
                else:
                    diag_all_bins[i_bin] = 0

        return diag_all_bins, N_all_bins, N_valid_all_bins
    
        
        