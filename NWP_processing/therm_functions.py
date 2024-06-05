#####
## Collection of functions for converting temperatures, 
##    calculating mixing heights etc...
##
#####
import numpy as np

def thl_to_ta(thl,ql,p,p0=100000):

    lv = 2.5E6
    T0 = 273.
    rcp = 0.286
    cp = 1004.
    th = thl
    n_it = 10
    for i in range(n_it):
        T = th/ ((p0/p)**(rcp))
        th = thl * (1. + ((lv * T0 * ql) / (cp * T)))
        # print(T,th)
    return T

def MixingHeight(turb,z):
    #
    # Adjusted: Christos' script to determine the mixing height for sigmaw2 > 0.069 and sigmaw2 < 0.129
    # outputs the mixing height 50, 75 and 25 percentile of the mixing height over the 21 intervals
    # !! for 1 timestep !!
    #
    mh = []
    for k in np.arange(0,20):
        alpha_v= 0.069+0.003*k
        tmp2 = []
        tmp = []
        tmp1 = []
        for i in np.arange(len(z[:])-2):
            if ((turb[i]<alpha_v) & (turb[i+1]<alpha_v) & (turb[i+2]<alpha_v)):
                tmp = i
            elif ((turb[i]<alpha_v) & (turb[i+1]<alpha_v)):
                tmp =  i
            elif (turb[i]<alpha_v):
                tmp =  i
            tmp1.append(tmp)
        if len(np.array(tmp1))==0:
            tmp1 = np.nan
        try:
            tmp2.append(z[np.nanmin(tmp1)])
        except ValueError:
            tmp2.append(np.nan)
        mh.append(tmp2)
    mh_50 = np.nanmedian(mh)
#    mh_75 = np.percentile(mh,75)
#    mh_25 = np.percentile(mh,25)
    return mh_50