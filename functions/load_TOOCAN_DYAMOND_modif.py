# -*- coding: utf-8 -*-

##########################################################################
# Histogrammes permettant de choisir les variables les plus pertinentes  #
# pour ensuite les utiliser dans les methodes de clustering              #
# modified for Python 3 and ascii decoding
##########################################################################

### Importation des bibliotheques
import sys
import matplotlib as mpl
mpl.use('agg')
import time
from math import *
import numpy as np
import gzip
import subprocess

from struct import unpack
from struct import *
#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.basemap import addcyclic

import matplotlib.pyplot as plt
import matplotlib.colors as mc
import matplotlib.gridspec as gridspec
import gzip

from datetime import *
from math import *
from mpl_toolkits.mplot3d import Axes3D
from xml.dom import minidom
from random import randint
from datetime import datetime

##################################################################################
#
#  Class convective system_IntParameters   : 
#		
#
#
##################################################################################

class MCS_IntParameters(object): 

    def __init__(self):
        self.label       		= 0
        self.qc_MCS	 			= 0
        self.duration			= 0
        self.classif            = 0
        self.Tmax               = 0.
        self.Utime_Init			= 0
        self.lonInit			    = 0
        self.latInit			    = 0
        self.Utime_End			= 0
        self.lonEnd				= 0
        self.latEnd				= 0
        self.lonmin				= 0
        self.latmin				= 0		
        self.lonmax				= 0
        self.latmax				= 0
        self.vavg				= 0
        self.dist				= 0
        self.olrmin				= 0
        self.surfmaxPix_172Wm2	= 0
        self.surfmaxkm2_172Wm2	= 0
        self.surfmaxkm2_132Wm2 	= 0
        self.surfmaxkm2_110Wm2 	= 0
        self.surfmaxkm2_90Wm2  	= 0
        self.surfcumkm2_172Wm2 	= 0
        self.surfcumkm2_132Wm2 	= 0
        self.surfcumkm2_110Wm2 	= 0
        self.surfcumkm2_90Wm2 	= 0
        self.precip_total   	= 0
        self.precip_max      	= 0
        self.maxSurf00mmh_km2 	= 0
        self.maxSurf02mmh_km2	= 0
        self.maxSurf05mmh_km2	= 0
        self.maxSurf10mmh_km2	= 0
        self.classif_JIRAK      = 0    

class MCS_Lifecycle(object):

    def __init__(self):
        self.qc_im				= []
        self.olrmin				= []
        self.olravg_172Wm2 		= []
        self.olravg_110Wm2      = []
        self.olravg_90Wm2       = []
        self.olr_90th           = []
        self.surfPix_172Wm2		= []
        self.surfPix_110Wm2		= []
        self.surfKm2			= []
        self.Utime				= []
        self.Localtime			= []
        self.lon				= []
        self.lat				= []
        self.x					= []
        self.y					= []
        self.velocity			= []
        self.semiminor_132Wm2	= []
        self.semimajor_132Wm2	= []
        self.orientation_132Wm2	= []
        self.excentricity_132Wm2= []
        self.semiminor_172Wm2		= []
        self.semimajor_172Wm2		= []
        self.orientation_172Wm2	= []
        self.excentricity_172Wm2	= []
        self.surfkm2_172Wm2		= []
        self.surfkm2_132Wm2		= []
        self.surfkm2_110Wm2		= []
        self.surfkm2_90Wm2		= []     

def load_TOOCAN_DYAMOND(FileTOOCAN):

    lunit=gzip.open(FileTOOCAN,'rb')
    print(FileTOOCAN)
    #
    # Read the Header
    ##########################
    header1 = lunit.readline()
    header2 = lunit.readline()
    header3 = lunit.readline()
    header4 = lunit.readline()
    header5 = lunit.readline()
    header6 = lunit.readline()
    header7 = lunit.readline()
    header8 = lunit.readline()
    header9 = lunit.readline()
    header10 = lunit.readline()
    header11 = lunit.readline()
    header12 = lunit.readline()
    header13 = lunit.readline()
    header14 = lunit.readline()
    header15 = lunit.readline()

    header16 = lunit.readline()
    header17 = lunit.readline()
    header18 = lunit.readline()
    header19 = lunit.readline()
    header20 = lunit.readline()
    header21 = lunit.readline()

    #print(header20)
    #print(header21)

    data = []
    iMCS = -1
    lines = lunit.readlines()

    #toggle = True

    for iline in lines:
        Values = iline.split()

    for i in range(len(Values)): 
        Values[i] = Values[i].decode('ascii') #modified!
        
    print(iline,Values)

    # if toggle : 
    #     print(Values[0])
    #     print(Values[1])
    # #if(Values[0] == '==>'): print("it works")
    # toggle = False

    if(Values[0] == '==>'):       
        #
        # Read the integrated parameters of the convective systems
        ###########################################################
        data.append(MCS_IntParameters())
        iMCS = iMCS+1
        data[iMCS].label 			= int(Values[1])	    # Label of the convective system in the segmented images
        data[iMCS].qc_MCS			= int(Values[2])	    # Quality control of the convective system 
        data[iMCS].classif			= int(Values[3])	    # Quality control of the convective system        
        data[iMCS].duration			= float(Values[4])/2	# duration of the convective system (slot)
        data[iMCS].Utime_Init		= float(Values[5])		# time TU of initiation of the convective system
        data[iMCS].localtime_Init	= float(Values[6])		# local time of inititiation
        data[iMCS].lonInit			= float(Values[7])		# longitude of the center of mass at inititiation
        data[iMCS].latInit			= float(Values[8])		# latitude of the center of mass at inititiation
        data[iMCS].Utime_End		= float(Values[9])		# time TU of dissipation of the convective system
        data[iMCS].localtime_End	= float(Values[10])		# local hour of dissipation
        data[iMCS].lonEnd			= float(Values[11])		# longitude of the center of mass at dissipation
        data[iMCS].latEnd			= float(Values[12])		# latitude of the center of mass at dissipation
        data[iMCS].vavg				= float(Values[13])		# average velocity during its life cycle(m/s)
        data[iMCS].dist				= float(Values[14])		# distance covered by the convective system during its life cycle(km)

        data[iMCS].lonmin			= float(Values[15])		# longitude min of the center of mass during its life cycle
        data[iMCS].lonmax			= float(Values[16])		# longitude max of the center of mass during its life cycle
        data[iMCS].latmin			= float(Values[17])		# latitude min of the center of mass during its life cycle
        data[iMCS].latmax			= float(Values[18])		# latitude max of the center of mass during its life cycle

        data[iMCS].olrmin			= float(Values[19])		#minimum Brigthness temperature (K)

        data[iMCS].surfmaxPix_172Wm2	= int(Values[20])		# maximum surface for a 235K threshold of the convective system during its life cycle (pixel)
        data[iMCS].surfmaxkm2_172Wm2	= float(Values[21])		#maximum surfacefor a 235K threshold of the convective system during its life cycle (km2)
        data[iMCS].surfmaxkm2_132Wm2	= float(Values[22])		#maximum surfacefor a 235K threshold of the convective system during its life cycle (km2)
        data[iMCS].surfmaxkm2_110Wm2	= float(Values[23])		#maximum surfacefor a 235K threshold of the convective system during its life cycle (km2)
        data[iMCS].surfmaxkm2_90Wm2	= float(Values[24])		#maximum surfacefor a 235K threshold of the convective system during its life cycle (km2        
        data[iMCS].surfcumkm2_172Wm2	= float(Values[25]) 	# integrated cumulated surface for a 235K threshold of the convective system during its life cycle (km2)
        data[iMCS].classif_JIRAK		= float(Values[26])        
        data[iMCS].clusters = MCS_Lifecycle()        
        inc = 0
    else:
        #
        # Read the parameters of the convective systems 
        #along their life cycles
        ##################################################
        data[iMCS].clusters.qc_im.append(int(Values[0]))	    			# quality control on the Infrared image
        data[iMCS].clusters.olrmin.append(float(Values[1]))	    			# min brightness temperature of the convective system at day TU (K)
        data[iMCS].clusters.olravg_172Wm2.append(float(Values[2]))	    			#  average brightness temperature of the convective system at day TU (K) 
        data[iMCS].clusters.olravg_110Wm2.append(float(Values[3]))            #  min brightness temperature of the convective system at day TU (K)
        data[iMCS].clusters.olravg_90Wm2.append(float(Values[4]))	    	#  min brightness temperature of the convective system at day TU (K)
        data[iMCS].clusters.olr_90th.append(float(Values[5]))	    		#  min brightness temperature of the convective system at day TU (K        
        data[iMCS].clusters.Utime.append(float(Values[6]))	    			# day TU 
        data[iMCS].clusters.Localtime.append(float(Values[7]))	    		# local hour (h)
        data[iMCS].clusters.lon.append(float(Values[8]))	    			# longitude of the center of mass (°)
        data[iMCS].clusters.lat.append(float(Values[9]))	    			# latitude of the center of mass (°)
        data[iMCS].clusters.x.append(int(Values[10]))		    			# column of the center of mass (pixel)
        data[iMCS].clusters.y.append(int(Values[11]))		    			# line of the center of mass(pixel)
        data[iMCS].clusters.velocity.append(float(Values[12]))	    		# instantaneous velocity of the center of mass (m/s        
        data[iMCS].clusters.semiminor_172Wm2.append(float(Values[13]))	    #
        data[iMCS].clusters.semimajor_172Wm2.append(float(Values[14]))	    #
        data[iMCS].clusters.excentricity_172Wm2.append(float(Values[15]))  	#
        data[iMCS].clusters.orientation_172Wm2.append(float(Values[16]))   	        
        data[iMCS].clusters.semiminor_132Wm2.append(float(Values[17]))	    #
        data[iMCS].clusters.semimajor_132Wm2.append(float(Values[18]))	    #
        data[iMCS].clusters.excentricity_132Wm2.append(float(Values[19]))  	#
        data[iMCS].clusters.orientation_132Wm2.append(float(Values[20]))   	        
        data[iMCS].clusters.surfPix_172Wm2.append(int(Values[21]))	    	# surface of the convective system at time day TU (pixel)
        data[iMCS].clusters.surfPix_110Wm2.append(int(Values[22]))	    	# surface of the convective system at time day TU (pixel)

        data[iMCS].clusters.surfkm2_172Wm2.append(float(Values[23]))  		# surface of the convective system for a 235K threshold
        data[iMCS].clusters.surfkm2_132Wm2.append(float(Values[24]))  		# surface of the convective system for a 200K threshold
        data[iMCS].clusters.surfkm2_110Wm2.append(float(Values[25]))  		# surface of the convective system for a 210K threshold
        data[iMCS].clusters.surfkm2_90Wm2.append(float(Values[26]))  		# surface of the convective system for a 220K threshol        
        #print float(Values[19]),float(Values[20]),float(Values[21]),float(Values[22]),float(Values[23])

    #data = MCS_Classif(data)
    #data = Compute_Tmax(data)

    return data    

# fileTracking="/data/fiolleau/DYAMOND/TOOCAN/TOOCAN_v2.07_olr/GLOBAL/2016/FileTracking/TOOCAN-GLOBAL-20160801-20160831.dat.gz"
#data           = load_TOOCAN_DYAMOND(fileTracking)
