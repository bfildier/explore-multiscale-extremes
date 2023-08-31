import argparse
import sys
import os

# own settings
from settings import *
workdir = os.getcwd()
moduledir, fcndir = defineDir(workdir)

# own functions
from casestudy import *
from fcns_load_DYAMOND_SAM import *
from PrecipGrid import *

#-- to redirect print output to standard output

class Unbuffered(object):

    def __init__(self, stream):
        self.stream = stream
    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
    def writelines(self, datas):
        self.stream.writelines(datas)
        self.stream.flush()
    def __getattr__(self, attr):
        return getattr(self.stream, attr)

sys.stdout = Unbuffered(sys.stdout)


def defineCaseStudy(region,varid,mask):
    
    # init
    cs = CaseStudy(name='DYAMOND-SAM',
                  region=region,
                  rel_tab_dyam_seg=loadRelTable('DYAMOND_SEG'))

    # define case study
    cs.setSimulationSpecs(i_t_min = 832+29,
                          i_t_max = 1917,
                          lat_slice=slice(-30,30),
                          lon_slice=slice(None))

    # where sliced distributions are stored
    cs.setDirectories(varid,mask)

    # load all distributions
    cs.loadDistSliced(varid,mask)

    # compute mean precip time series
    cs.computeMean(varid='Prec',mask='all')

    # find "missing" time indices
    cs.findTimeIndToIgnore()
    
    return cs


if __name__ == "__main__":
    
    ## Parse arguments
    parser = argparse.ArgumentParser(description='Compute the maximum precipitation for each MCS and the distance between the maximum and the center of the MCS')
    parser.add_argument('--day_i', type=int, default=1, help='First day to be simulated, included')
    parser.add_argument('--day_f', type=int, default=1, help='Last day to be simulated, excluded')
    parser.add_argument('--varid','-v',type=str, help='variable to regrid')
    parser.add_argument('--func',type=str,default='mean',help='What to compute in coarser grid (mean,max)')
    args = parser.parse_args()
    
    day_i = args.day_i
    day_f = args.day_f
    varid = args.varid
    func = args.func
    
    region = 'tropics'
    mask = 'all'

    # define case study
    cs = defineCaseStudy(region,varid='Prec',mask=mask)
    
    # Initialize
    grid = PrecipGrid(casestudy=cs,verbose_steps=True,overwrite=False)
    
    # prepare data ## -- to optimize
    grid.prepare_data()

    # compute surface
    grid.compute()
    
    # remove wrong variables
    for key in 'mean_prec','max_prec':
        if key in list(grid.ds.keys()):
            grid.ds = grid.ds.drop_vars(key)
    
    # regrid and save ## -- to optimize
    if func != 'all':
        
        grid.saveVarRegridded(day0=day_i,dayf=day_f,varid=varid,func=func)

    else:
        
        # mean
        grid.saveVarRegridded(day0=day_i,dayf=day_f,varid=varid,func='mean')
        # max
        grid.saveVarRegridded(day0=day_i,dayf=day_f,varid=varid,func='max')
        
    # exit
    print('Regridding completed.')
    sys.exit(0)
