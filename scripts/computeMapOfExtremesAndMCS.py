"""script computeMapOfExtremesAndMcs.py

load prec at one time slice from DYAMOND and a TOOCAN segmentation mask
to compute a map of frequency or rain amount above a extreme percentile
of the rain distribution, and the likelihood of MCS at each location.

B. Fildier November 2023
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

