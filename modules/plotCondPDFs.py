
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np


def subplotConditionalProfiles(ax,z,profiles,refbins=[0],refmin=None,refmax=None,ls='-',labels=None,normfunction='Normalize'):
    
    # Number of profiles
    Nprof = 1
    if len(profiles.shape) > 1:
        Nprof = profiles.shape[1]
    
    # min-max for refvariable
    if refmin is None: refmin = refbins[0] 
    if refmax is None: refmax = refbins[-1]
    
    # Color scale
    colorInd = range(Nprof)
    cm = plt.get_cmap('Spectral')
    cNorm = getattr(colors,normfunction)(vmin=refmin, vmax=refmax)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

    # Plot
    for i in np.flipud(range(Nprof)):
        
        # define label
        if labels is not None: # suppose it's a list of the right size
            lab = labels[i]
        else:
            lab = None
            
        colorVal = scalarMap.to_rgba(refbins[i])
        if Nprof > 1:
            ax.plot(profiles[:,i],z,c=colorVal,linestyle=ls,label=lab)
        else:
            ax.plot(profiles,z,c=colorVal,linestyle=ls,label=lab)

# def subplotDistributions(ax,distributions,perc,SSTs,refmin=None,refmax=None,ls='-',labels=None,normfunction='Normalize'):
    
#     # Number of profiles
#     Nd = 1
#     if len(distributions.shape) > 1:
#         Nd = distributions.shape[1]
    
#     # min-max for refvariable
#     if refmin is None: refmin = refbins[0] 
#     if refmax is None: refmax = refbins[-1]
    
#     # Color scale
#     colorInd = range(Nd)
#     cm = plt.get_cmap('Spectral') 
#     cNorm = getattr(colors,normfunction)(vmin=refmin, vmax=refmax)
#     scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

#     # Plot
#     for i in np.flipud(range(Nd)):
        
#         # define label
#         if labels is not None: # suppose it's a list of the right size
#             lab = labels[i]
#         else:
#             lab = None
            
#         colorVal = scalarMap.to_rgba(SSTs[i])
#         if len(refbins.shape) > 1:
#             ax.plot(perc,distributions[:,i],c=colorVal,linestyle=ls,label=lab)
#         else:
#             ax.plot(perc,distributions,c=colorVal,linestyle=ls,label=lab)

            
def showColorBar(fig,axs,values,vmin=None,vmax=None,cbar_factor=11,label='',normfunction='Normalize'):
    
    if axs.__class__ != np.ndarray:
        axs_all = axs,
    else:
        axs_all = axs.flatten()

    # min-max for refvariable
    if vmin is None: vmin = values[0] 
    if vmax is None: vmax = values[-1]
    
    x_left = np.nan
    x_right = np.nan
    y_bot = np.nan
    y_top = np.nan

    for ax in axs_all:

        # Save boundaries for legend
        x,y,w,h = ax.get_position().bounds
        x_left = np.nanmin(np.array([x,x_left]))
        x_right = np.nanmax(np.array([x+w,x_right]))
        y_bot = np.nanmin(np.array([y,y_bot]))
        y_top = np.nanmax(np.array([y+h,y_top]))


    # Color bar    
    cm = plt.get_cmap('Spectral') 
    cNorm = getattr(colors,normfunction)(vmin=vmin, vmax=vmax)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

    dy = (y_top-y_bot)/60
    cax = plt.axes([x_left,y_bot-cbar_factor*dy,x_right-x_left,dy])
    cbar = fig.colorbar(scalarMap, cax=cax, orientation='horizontal')

    return cbar.ax.set_xlabel(label)
