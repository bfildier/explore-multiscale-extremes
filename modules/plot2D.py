#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 17:05:02 2023

Module plot2D

Includes functions to 

- plot vertical profiles on transformed x axis or profiles conditioned on x variables
- plot joint histograms on inverse-logarithmic axes (extremes)

@author: bfildier
"""

#---- Modules ----#

from scipy.interpolate import griddata,interp2d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LogNorm

# ## Plot vertical data (transect, or vertical profiles over time)
# def subplotVerticalData(ax,x,y,Z,cmap=plt.cm.seismic,vmin=None,vmax=None,cbar=True):
    
#     """Arguments:
#         - x and y are coordinate values
#         - Z are the data values flattened"""
    
#     xs0,ys0 = np.meshgrid(x,y)
#     xmin = min(x[0],x[-1])
#     xmax = max(x[0],x[-1])
#     ymin = min(y[0],y[-1])
#     ymax = max(y[0],y[-1])
#     X0 = np.vstack([xs0.flatten(),ys0.flatten()]).T
  
#     extent = (xmax,xmin,ymin,ymax)

#     # New coordinates
#     xs,ys = np.meshgrid(np.linspace(xmin,xmax,num=len(x)),
#                         np.linspace(ymin,ymax,num=len(y)))
#     X = np.vstack([xs.flatten(),ys.flatten()]).T

#     # New values
#     resampled = griddata(X0,Z,X, method='cubic')
#     resampled_2D = np.reshape(resampled,xs.shape)

#     im = ax.imshow(np.flipud(resampled_2D),extent=extent,interpolation='bilinear',
#                    cmap=cmap,vmin=vmin,vmax=vmax,aspect='auto',origin='upper')
#     if cbar:
#         plt.colorbar(im,ax=ax)
    
#     if y[0] > y[-1]:
#         plt.gca().invert_yaxis()



#---- Functions ----#


#---- Transect or transformed x coordinate

# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work their way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def subplotSmooth2D(ax,x,y,Z,fplot='contourf',xmin=None,xmax=None,nx=50,nlev=50,vmin=None,vmax=None,**kwargs):
    """
    Plot 2D contours (exact method is defined by fplot) with user-defined Z-range and x range.
    """
    
    # set levels
    levels = nlev
    if vmin is not None and vmax is not None:
        levels = np.linspace(vmin,vmax,nlev+1)
    
    # only keep x_values where Z is not a nan
    notanan_x = np.any(~np.isnan(Z)&np.isfinite(Z),axis=0)
    notanan2d = np.vstack([notanan_x[None,:]]*Z.shape[0])
    x_valid = x[notanan_x]
    y_valid = y
    Nx,Ny = len(x_valid), len(y_valid)
    Z_valid=Z[notanan2d].reshape(Ny,Nx)

    if xmin is not None and xmax is not None:
        x_new = np.linspace(xmin,xmax,nx+1)
        f_interp = interp2d(x_valid,y_valid,Z_valid,kind='cubic')
        Z_new = f_interp(x_new,y_valid)
    # else:
    #     x_new = x_valid
    X,Y = np.meshgrid(x_new,y_valid)


    # remove values outside xrange
    X_out = np.logical_or(X < np.nanmin(x_valid),X > np.nanmax(x_valid))
    Z_new[X_out] = np.nan

    # plot
    return getattr(ax,fplot)(X,Y,Z_new,levels=levels,**kwargs)


#---- 2D joint PDFs

def computeTickLabels(xranks):
    
    # nonlinear transformation of ranks
    x = np.flipud(1./(1-xranks/100.))
    # linear scale
    k_all = np.log10(np.round(x,5))
    # tick positions
    xtick_pos = np.mod(k_all,1) == 0
    # positions of ticks on transformed axis
    xticks = x[xtick_pos]
    # linear ranks of ticks
    ks = k_all[xtick_pos]
    # n_digits in final tick labels
    ndigits = [int(max(k-2,1)) for k in ks]
    # labels
    xlab_floats = np.array(xranks[xtick_pos])
    # convert labels to strings
    xticklabels = np.array([('%'+('2.%d'%ndig)+'f')%s for ndig,s in zip(np.flipud(ndigits),xlab_floats)])
    
    return xticks, xticklabels

def setFrameIL(ax,xranks,yranks,aspect='1'):
    """Set inverse-logarithmic axes on x and y axes"""
    
    ##-- create inverse-log frame
    # define axes values
    x = np.flipud(1./(1-xranks/100.))
    y = np.flipud(1./(1-yranks/100.))
    Nx = len(x)
    Ny = len(y)
    # h = subplotSmooth2D(ax,x,y,Z,fplot='contourf',xmin=ymin,xmax=ymin,nx=50,nlev=50,vmin=None,vmax=None)
    ax.matshow(np.full((Nx,Ny),np.nan),origin='lower',extent=[x[0],x[-1],y[0],y[-1]],aspect=aspect)
    ax.set_xscale('log')
    ax.set_yscale('log')
    # axes labels and positions
    ax.xaxis.set_ticks_position('bottom')

    # set xticks
    xticks, xticklabels = computeTickLabels(xranks)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    # set yticks
    yticks, yticklabels = computeTickLabels(yranks)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

    return ax


def showJointHistogram(ax,values,scale='linear',vmin=1e-3,vmax=1,cmap=None,**kwargs):
    """Show matrix data as it is, regardless of preset frame and ticks"""

    if scale == 'linear':
        h = ax.matshow(values,vmin=vmin,vmax=vmax,origin='lower',cmap=cmap,**kwargs)
    elif scale == 'log':
        h = ax.matshow(values,norm=LogNorm(vmin=vmin,vmax=vmax),origin='lower',cmap=cmap,**kwargs)

    ax.set_xticks([])
    ax.set_yticks([])

    return h