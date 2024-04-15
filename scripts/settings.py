import os,glob,sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import socket

#---- Paths ----#

if socket.gethostname() == 'clarity':
    
    DIR_DYAMOND = '/Users/bfildier/Data/vrac_DYAMOND_etc'
    DIR_DYAMOND_DIAG2D = '/Users/bfildier/Data/simulations/DYAMOND/DYAMOND_REGIONS/tropics/SAM/diagnostics_2D'
    DIR_DYAMOND_PROCESSED = '/Users/bfildier/Data/simulations/DYAMOND/DYAMOND_REGIONS'
    # DIR_TOOCANSEG_DYAMOND = '/Users/bfildier/Data/simulations/DYAMOND/TOOCAN/TOOCAN_v2.07/GLOBAL/2016'
    # DIR_TOOCAN_DYAMOND = '/Users/bfildier/Data/simulations/DYAMOND/TOOCAN/TOOCAN_v2.07/GLOBAL/2016/FileTracking'
    DIR_TOOCANSEG_DYAMOND = '/bdd/MT_WORKSPACE/lgouttes/MODELS/DYAMOND/Summer/SAM/235K/TOOCAN/TOOCAN_v2.08/GLOBAL/2016'
    DIR_TOOCAN_DYAMOND = '/bdd/MT_WORKSPACE/lgouttes/MODELS/DYAMOND/Summer/SAM/235K/TOOCAN/TOOCAN_v2.08/GLOBAL/2016/FileTracking'
    DIR_RCEMIP = None
    DIR_TOOCANSEG_RCEMIP = None
    
else:

    DIR_DYAMOND = '/bdd/DYAMOND/SAM-4km/OUT_2D'
    DIR_DYAMOND_DIAG2D = '/data/bfildier/DYAMOND_REGIONS/tropics/SAM/diagnostics_2D'
    DIR_DYAMOND_PROCESSED = '/data/bfildier/DYAMOND_REGIONS'
    # DIR_TOOCAN_DYAMOND = '/data/fiolleau/DYAMOND/TOOCAN/TOOCAN_v2.07/GLOBAL/2016/FileTracking'
    # DIR_TOOCANSEG_DYAMOND = '/data/fiolleau/DYAMOND/TOOCAN/TOOCAN_v2.07/GLOBAL/2016'
    DIR_TOOCANSEG_DYAMOND = '/bdd/MT_WORKSPACE/lgouttes/MODELS/DYAMOND/Summer/SAM/235K/TOOCAN/TOOCAN_v2.08/GLOBAL/2016'
    DIR_TOOCAN_DYAMOND = '/bdd/MT_WORKSPACE/lgouttes/MODELS/DYAMOND/Summer/SAM/235K/TOOCAN/TOOCAN_v2.08/GLOBAL/2016/FileTracking'
    # DIR_RCEMIP = '/bdd/MT_WORKSPACE/REMY/RCEMIP/SAM'
    DIR_RCEMIP = '/scratchx/bfildier/RCEMIP/SAM'
    DIR_TOOCANSEG_RCEMIP = '/bdd/MT_WORKSPACE/MCS/RCE/SAM/TOOCAN/TOOCAN_v2022_04/irtb'
    
    # DIRECTORIES = {'TOOCANSEG':
    #                    {'v2.08':
    #                         {'DYAMOND':
    #                              {'Summer':
    #                                   {'SAM':'/bdd/MT_WORKSPACE/lgouttes/MODELS/DYAMOND/Summer/SAM/235K/TOOCAN/TOOCAN_v2.08/GLOBAL/2016',}}},
    #                     'v2.07':
    #                         {'DYAMOND':
    #                              {'Summer':
    #                                   {'SAM':}}}},
    #                     'v2022_04':
    #                         {'RCEMIP':
    #                              {'SAM':}}}},
    #                'model_output':
    #                     {'DYAMOND':
    #                         {'Summer':
    #                             {'SAM':}},}
                
    DIR_ROOTS = {'MCSMIP':{'TOOCANSEG':'/bdd/MT_WORKSPACE/lgouttes/MODELS/MCSMIP/DYAMOND/Summer/%s'+\
                           '/TOOCAN/TOOCAN_%s/GLOBAL/2016', # expect %(<model or obs>,<toocan_version>)
                           'model_output':'/bdd/MT_WORKSPACE/MCS/MODELS/MCSMIP/DYAMOND/Summer/%s/olr_pcp_instantaneous'}, # expect %(<model or obs>)
                 '4km-30mn':{'TOOCANSEG':'/bdd/MT_WORKSPACE/lgouttes/MODELS/DYAMOND/Summer/%s'+\
                             '/235K/TOOCAN/TOOCAN_%s/GLOBAL/2016', # expect %(<model or obs>,<toocan_version>)
                             'model_output':'/bdd/MT_WORKSPACE/lgouttes/MODELS/DYAMOND/Summer/%s/235K/TOOCAN/TOOCAN_%s/GLOBAL/2016'}, # expect %(<model or obs>)
                }
    
    def getDirDYAMOND(resolution,data_type,model_name,toocan_version='v2.07'):
        """Return the full path for model or observation data, for OLR or precipitation data or TOOCAN segmentation.
        
        Arguments:
        - resolution: [MCSMIP, 4km-30mn]
        - data_type: [TOOCANSEG, TOOCAN, model_output]
        - model_name: [ARPEGE FV3 IFS MPAS NICAM OBS SAM UM]
        - toocan_version: [v2.07, v2.08, v2022_04]
        """
        
        # get root
        if data_type == 'TOOCAN':
            root = DIR_ROOTS[resolution]['TOOCANSEG']
        else:
            root = DIR_ROOTS[resolution][data_type]
        if model_name == 'OBS':
            root.replace('olr_pcp_instantaneous','olr_pcp')
                 
        # insert expected strings
        if data_type == 'TOOCANSEG':
            full_path = root%(model_name,toocan_version)
        elif data_type == 'TOOCAN':
            full_path = os.path.join(root%(model_name,toocan_version),'FileTracking')
        elif data_type == 'model_output':
            full_path = root%(model_name)
            
        return full_path
    
    def getDirRCEMIP(data_type,model_name):
        """Return the full path for model or observation data, for OLR or precipitation data or TOOCAN segmentation.
        
        Arguments:
        - data_type: [TOOCANSEG, TOOCAN, model_output]
        - model_name: [ICON  MESONH  SAM]
        """
        
        model_names_4_data = {'ICON':'ICON2','MESONH':'MNH','SAM':'SAM'}
        
        if data_type == 'model_output':
            full_path = '/bdd/MT_WORKSPACE/REMY/RCEMIP/%s/300K/'%(model_names_4_data[model_name])
        elif data_type == 'TOOCAN':
            full_path = '/bdd/MT_WORKSPACE/MCS/RCE/%s/TOOCAN/TOOCAN_v2022_04/irtb'%model_name
            
        
                 

DIR_DATA = '../input'
DIR_FIG = '../figures'
# DIR_OUT = '../results'
DIR_OUT = '/data/bfildier/multiscale-extremes'
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
box_0 = [0,360,-30,30] # whole tropics
box_1 = [310,340,0,20] # Atlantic ITCZ
box_2 = [205,250,0,20] # Eastern Pacific ITCZ
box_3 = [130,165,0,20] # Pacific Warm Pool
box_4 = [-20,35,0,20] # Central Africa

# coord slices
coord_slices = {'lat':{'tropics':slice(-30,30)},
                'lon':{'tropics':None}}
