#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 10:28:06 2023

Module plot1D

Includes new functions to 

- plot curves on inverse-logarithmic x-axis (extremes). Improved from plot1DInvLog
- add shadings

@author: bfildier
"""

#---- Modules ----#

import matplotlib.pyplot as plt
import numpy as np
from math import log10,ceil
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm

#---- Functions ----#

#-- Inverse logarithmic axes to display extremes with correct ticks
#
# The trick here is to duplicate the axis and rescale the background frame
# so that the data, already correctly sampled on the transformed axis, can be displayed
# independently from the figure frame
#
# Example : start with the setting the figure and duplicate y axis
#
# fig,ax = plt.subplots(figsize=(6,4.5))
# ax_show = ax.twiny()


#-- create inverse-log frame

def setXaxisIL(ax,ranks,xtickrotation=0):
    
    # define axes values
    x = np.flipud(1./(1-ranks/100.))
    Nx = len(x)
    # display empty data
    h = ax.plot(x,np.full((Nx,),np.nan))
    ax.set_xscale('log')
    # invert direction
    ax.invert_xaxis()
    
    # axis position
    ax.xaxis.set_ticks_position('bottom')

    # set xticks
    xtick_pos = np.mod(np.round(np.log10(x),7),1) == 0
    xticks = x[xtick_pos]
    ks = np.round(np.log10(xticks),7)
    xticklabels = np.array([("%%2.%1df"%max(0,k-2))%((1-10**(-k))*100) for k in np.flipud(ks)])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels,rotation=xtickrotation)
    # set labelpad
    # ax.set_xlabel(labelpad=4.0)
    
def setYaxisIL(ax,ranks):
    
    # define axes values
    y = np.flipud(1./(1-ranks/100.))
    Ny = len(y)
    # display empty data
    h = ax.plot(np.full((Ny,),np.nan),y)
    ax.set_yscale('log')
    # invert direction
    ax.invert_yaxis()
    
    # axis position
    ax.yaxis.set_ticks_position('left')

    # set xticks
    ytick_pos = np.mod(np.log10(np.round(y,7)),1) == 0
    yticks = y[ytick_pos]
    ks = np.round(np.log10(yticks),7)    
    yticklabels = np.array([("%%2.%1df"%max(0,k-2))%((1-10**(-k))*100) for k in np.flipud(ks)])
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

def setFrame(ax,rankmin=0,rankmax=99.99,axisIL='x',xtickrotation=0):
    
    #- duplicate axes for rescaling frame
    ax_frame = ax.twiny()
    
    # set ranks
    k_min = -np.round(log10(1-rankmin/100))
    k_max = -np.round(log10(1-rankmax/100))
    dk = 0.1
    scale_invlog = np.arange(k_min,k_max+dk,dk)
    ranks_frame = (1 - np.power(10,-scale_invlog))*100
    
    #- set axis
    if axisIL == 'x':
        setXaxisIL(ax_frame,ranks_frame,xtickrotation=xtickrotation)
    elif axisIL == 'y':
        setXaxisIL(ax_frame,ranks_frame)
        
    return ax_frame
    
def showData(ax,ranks,values,axisIL='x',rankmin=0,rankmax=99.99,**kwargs):
    """Show data as it is, regardless of preset frame and ticks"""

    if axisIL == 'x':
        
        x = 1/(1-ranks/100.)
        xmin = xmax = None
        if rankmin is not None:
            xmin = 1/(1-rankmin/100.)
        if rankmax is not None:
            xmax = 1/(1-rankmax/100.)
        # show
        h = ax.plot(x,values,**kwargs)
        # be careful that the x bounds are precisely the same as the background frame
        ax.margins(x=0)
        # bounds
        ax.set_xlim(xmin,xmax)
        # log
        ax.set_xscale ('log')
        # remove ticks
        ax.set_xticks([])
        ax.set_xticks([], minor=True)
        
        return h

    elif axisIL == 'y':
        
        y = 1/(1-ranks/100.)
        if rankmin is not None:
            ymin = 1/(1-rankmin/100.)
        if rankmax is not None:
            ymax = 1/(1-rankmax/100.)
        # show
        h = ax.plot(values,y,**kwargs)
        # be careful that the y bounds are precisely the same as the background frame
        ax.margins(y=0)
        # bounds
        ax.set_ylim(ymin,ymax)
        # log
        ax.set_yscale ('log')
        # remove ticks
        ax.set_yticks([])
        ax.set_yticks([], minor=True)

        return h


def subplotRanksILog(ax,ranks,y,sl=slice(None,None),rankmin=0,rankmax=99.999,setframe=True,\
                     col=None,ltype=None,linewidth=None,alpha=None,labels=None,offset=0,xtickrotation=0):
    """Display one or several curves on inverse=logarithmic axis.
    To allow successive use of this method, set frame based on rankmin and rankmax, duplicate axis for frame, and use
    axes given in argument to display curve"""
    
    #- set frame
    ax_frame = None
    if setframe:
        ax_frame = setFrame(ax,rankmin=rankmin,rankmax=rankmax,xtickrotation=xtickrotation)
        
    #- show data
    # init handle list
    h_all = []
    # show
    if isinstance(y,list):
        for i in range(len(y)):
            lab = None
            if labels is not None:
                lab = labels[i]
            lt = ltype[i] if ltype is not None else '-'
            a = alpha[i] if alpha is not None else 1
            c = col[i] if col is not None else 'k'
            lw = linewidth[i] if linewidth is not None else 1.5
            h = showData(ax,ranks[i],y[i],axisIL='x',rankmin=rankmin,rankmax=rankmax,c=c,alpha=a,linestyle=lt,linewidth=lw,label=lab)
            h_all.append(h)
    else:
        h_all = showData(ax,ranks,y,axisIL='x',rankmin=rankmin,rankmax=rankmax,c=col,alpha=alpha,linestyle=ltype,linewidth=linewidth,label=labels)

    return ax_frame, h_all


def addXHatch(ax,ranks,i_xlim,color='gray',hatch='//',
    alpha=1,fill=False,closed=True,**kwargs):
    """Add vertical shading"""

    x = 1/(1-ranks/100.)
    
    ax.add_patch(Polygon([[x[i_xlim[0]], ax.get_ylim()[0]],\
                          [x[i_xlim[1]], ax.get_ylim()[0]],\
                          [x[i_xlim[1]], ax.get_ylim()[1]],\
                          [x[i_xlim[0]], ax.get_ylim()[1]]],\
                          closed=closed, fill=fill, hatch=hatch,linewidth=0,
                          color=color,alpha=alpha,**kwargs))
    
def addYHatch(ax,ranks,i_ylim,color='gray',hatch='//',
    alpha=1,fill=False,**kwargs):
    """Add horizontal shading"""

    y = 1/(1-ranks/100.)
    
    ax.add_patch(Polygon([[ax.get_xlim()[0],y[i_ylim[0]]],\
                          [ax.get_xlim()[0],y[i_ylim[1]]],\
                          [ax.get_xlim()[1],y[i_ylim[1]]],\
                          [ax.get_xlim()[1],y[i_ylim[0]]]],\
                          closed=True, fill=fill, hatch=hatch,linewidth=0,
                          color=color,alpha=alpha,**kwargs))
