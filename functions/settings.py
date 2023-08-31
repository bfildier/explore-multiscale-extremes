import os,glob,sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


#---- Paths ----#

DIR_DYAMOND = '/bdd/DYAMOND/SAM-4km/OUT_2D'
DIR_DYAMOND_DIAG2D = '/data/bfildier/DYAMOND_REGIONS/tropics/SAM/diagnostics_2D'
DIR_DYAMOND_PROCESSED = '/data/bfildier/DYAMOND_REGIONS'
# DIR_RCEMIP = '/bdd/MT_WORKSPACE/REMY/RCEMIP/SAM'
DIR_RCEMIP = '/scratchx/bfildier/RCEMIP/SAM'
DIR_TOOCANSEG_RCEMIP = '/bdd/MT_WORKSPACE/MCS/RCE/SAM/TOOCAN/TOOCAN_v2022_04/irtb'
DIR_DATA = '../input'
DIR_FIG = '../figures'
DIR_OUT = '../results'
DIR_TEMPDATA = '../temp'

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


#---- colors anb bounds ----#

clim_specs = {'prec':(1e-2,1e2), # mm (in 30mn)
              'PW':(10,70)}      # mm

cmap_specs = {'prec':plt.cm.ocean_r,   # alternative plt.cm.bone_r
              'PW':plt.cm.RdBu,
              'mcs':plt.cm.get_cmap('Accent', 10)}

norm_specs = {'prec':LogNorm(vmin=clim_specs['prec'][0], vmax=clim_specs['prec'][1]),
              'PW':None}

#---- regions of analysis ----#

# define boxes (xmin,xmax,ymin,ymax)
box_1 = [310,340,0,20] # Atlantic ITCZ
box_2 = [205,250,0,20] # Eastern Pacific ITCZ
box_3 = [130,165,0,20] # Pacific Warm Pool
box_4 = [-20,35,0,20] # Central Africa
