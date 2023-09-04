import numpy as np

from settings import *
from scipy.optimize import curve_fit

##------ MULTISCALE EXTREMES -----##

def fitBranches(cont,N):

    def func(x, a, b, c):
        return a * np.exp(-b * x) + c
    
    if cont.__class__ is list:
        seg_1 = np.flip(cont[0],axis=1)
    else:
        seg_1 = cont.allsegs[0][0]
        
    # Branch 1 -- end of contour (upward branch)
    xdata_1 = seg_1[-N:,0]
    y_1 = ydata_1 = seg_1[-N:,1]

    # fit
    popt_1, pcov_1 = curve_fit(func, ydata_1, xdata_1,p0=(-10,1,0))
    x_1 = func(ydata_1, *popt_1)
    
    # Branch 2 -- start of contour
    x_2 = xdata_2 = seg_1[:N,0]
    ydata_2 = seg_1[:N,1]

    # fit
    popt_2, pcov_2 = curve_fit(func, xdata_2, ydata_2,p0=(-10,1,0))
    y_2 = func(xdata_2, *popt_2)
    
    return popt_1, x_1, y_1, popt_2, x_2, y_2, func