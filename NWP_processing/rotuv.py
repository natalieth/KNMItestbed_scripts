# -*- coding: utf-8 -*-
"""
@author: Hylke de Vries. NB. This script does not come with any warranty/guaranty that it
@ is doing the right thing. Please check whether it behaves as you expect by
@ testing on synthetic data (for example by setting one component to 0)
#
# script to unrotate Harmonie wind field vectors to become oriented along
# north-south and east-west direction. The native grib output is directed
# along the grid directions. The idea is to locally compute the angle between 
# the grid-lines and the geo-coordinates and rotate the vectors accordingly.
# The assumption here is that the grid forms a locally orthogonal set (not sure
# completely whether this is true for all gridchoices, but over Europe it 
# seems to be the case). The new version (09/2020) implements this using two 
# angles, which together add up to 180 deg..
#  
# Created on Tue 20 Mar 2018. modification log: 
# 2020-09 / implemented a new version with two angles rather than one.
# both give almost identical results, except at the edge-rows..
# these edge-numbers shouldn't be trusted anyway.
# 
# 
"""

import os.path
import sys,getopt
import numpy as np
import xarray as xr

USE_NEW_VERSION=True

def rotuvb(u,v,lon,lat):
    # determine angle of rotation; assume grid is locally orthogonal...
    # lon and latitudes are assumed to be in radians...
    # cos(lat) weighting to account for reduction of distance between longitudes?
    print('... computing rotation factor')
    DLON = np.nan*lon
    DLAT = np.nan*lon
    ny = lon.shape[0]
    nx = lon.shape[1]
    for jy in range(0,ny):
        jyp1 = np.min([jy+1,ny-1])
        jym1 = np.max([jy-1,0])
        DLAT[jy,:] = (lat[jyp1,:]-lat[jym1,:])
        DLON[jy,:] = (lon[jyp1,:]-lon[jym1,:])*np.cos(lat[jy,:])
    alpha = np.arctan2(DLON,DLAT)
    print('... rotating actual vectors')
    ugeo =  u*np.cos(alpha) + v*np.sin(alpha) 
    vgeo =  v*np.cos(alpha) - u*np.sin(alpha) 
    return ugeo,vgeo


def rotuvbNEW(u,v,lon,lat):
    # determine angle of rotation; assume grid is locally orthogonal...
    # lon and latitudes are assumed to be in radians...
    # cos(lat) weighting to account for reduction of distance between longitudes?
    print('... computing rotation factor')
    DLON = np.nan*lon
    DLAT = np.nan*lon
    ny = lon.shape[0]
    nx = lon.shape[1]
    for jy in range(0,ny):
        jyp1 = np.min([jy+1,ny-1])
        jym1 = np.max([jy-1,0])
        DLAT[jy,:] = (lat[jyp1,:]-lat[jym1,:])
        DLON[jy,:] = (lon[jyp1,:]-lon[jym1,:])*np.cos(lat[jy,:])
    beta = np.arctan2(DLAT,DLON)
    for jx in range(0,nx):
        jxp1 = np.min([jx+1,nx-1])
        jxm1 = np.max([jx-1,0])
        DLAT[:,jx] = (lat[:,jxp1]-lat[:,jxm1])
        DLON[:,jx] = (lon[:,jxp1]-lon[:,jxm1])*np.cos(lat[:,jx])
    alpha = np.arctan2(DLAT,DLON)    
    #print((-alpha+beta)/np.pi) # almost a half... grid almost orthogonal...
        
    print('... rotating actual vectors')
    ugeo =  u*np.cos(alpha) + v*np.cos(beta) 
    vgeo =  v*np.sin(beta)  + u*np.sin(alpha) 
    return ugeo,vgeo


  

if __name__ == '__main__':
    argv = (sys.argv[1:])    
    try:
        opts, args = getopt.getopt(argv,"i:o:k:l:u:v",["uifile=","vifile=","uofile=","vofile=","upar=","vpar="])
    except getopt.GetoptError:
        print('use the correct input formattings please')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('sorry, no help at this stage....')
            sys.exit()
        elif opt in ("--uifile"):
            ufile = str(arg)
        elif opt in ("--vifile"):
            vfile = str(arg)
        elif opt in ("--uofile"):
            tufile = str(arg)
        elif opt in ("--vofile"):
            tvfile = str(arg)
        elif opt in ("--upar"):
            upar = str(arg)
        elif opt in ("--vpar"):
            vpar = str(arg)

    #%%
    print('read')
    unc  = xr.open_dataset(ufile)
    vnc  = xr.open_dataset(vfile)
    uu   = unc[upar]
    vv   = vnc[vpar]
    lons = unc['lon'].values * np.pi/180.
    lats = unc['lat'].values * np.pi/180.
    
    #%% rotate
    print('rotate')
    if USE_NEW_VERSION:
        ugeo,vgeo = rotuvbNEW(uu,vv,lons,lats)
    else:
        ugeo,vgeo = rotuvb(uu,vv,lons,lats)
      
    #%% copy and output
    print('save2nc')
    tunc = unc.copy(deep=True) 
    tvnc = vnc.copy(deep=True)
    tunc[upar] = ugeo[:]
    tvnc[vpar] = vgeo[:]
    tunc.to_netcdf(tufile)
    tvnc.to_netcdf(tvfile)
    
    # sometimes the script produces strange results. 
    # somehow then the lines below might work.
    #tunc.to_netcdf(tufile,unlimited_dims='time')
    #tvnc.to_netcdf(tvfile,unlimited_dims='time')
    
    # Finished
    print('done')
    
    

## LEGACY CODE FROM HERE...
    
    
#def OLD_rotuvb_OLD(u,v,lon,lat):
    ## determine angle of rotation; assume grid is locally orthogonal...
    ## NOTE the cos-lat weighting factor seems crucial...
    ## to get the end-results similar to output from gl...
    ## first dim=y,second dim=x...
    #print('... computing rotation factor')
    #DLON = np.nan*lon
    #DLAT = np.nan*lon
    #ny = lon.shape[0]
    #nx = lon.shape[1]
    #for jx in range(1,(nx-1)):
        #for jy in range(1,(ny-1)):
            #DLAT[jy,jx] = (lat[jy+1,jx]-lat[jy-1,jx]) #*(np.pi/180.)
            #DLON[jy,jx] = np.cos(lat[jy,jx])*(lon[jy+1,jx]-lon[jy-1,jx]) #*(np.pi/180.)
    #alpha = np.arctan2(DLON,DLAT)
    
    ## rotate the vectors to lat-lon coordinates
    #print('... rotating actual vectors')
    #urot = u*np.cos(alpha) + v*np.sin(alpha)
    #vrot = v*np.cos(alpha) - u*np.sin(alpha)
    #return urot,vrot,alpha
  
#def rotuvb_SIM2MPI(u,v,lon,lat):
    ## THIS VERSION GIVES ALMOST IDENTICAL RESULTS TO THE MPI_OM rotate-vecto
    ## BUT IS IT CORREC? I DOUBT IT... SIMPLE TEST WITH U=1 V=0 gives strange
    ## ugeo and vgeo...
    ## determine angle of rotation; assume grid is locally orthogonal...
    ## lon and latitudes are assumed to be in radians...
    ## cos(lat) weighting to account for reduction of distance between longitudes?
    #print('... computing rotation factor (1)')
    #DLON = np.nan*lon
    #DLAT = np.nan*lon
    #ny = lon.shape[0]
    #nx = lon.shape[1]
    ## seems to be working...
    #for jx in range(0,nx):
        #jxp1 = np.min([jx+1,nx-1])
        #jxm1 = np.max([jx-1,0])
        #DLAT[:,jx] = (lat[:,jxp1]-lat[:,jxm1])
        #DLON[:,jx] = (lon[:,jxp1]-lon[:,jxm1])*np.cos(lat[:,jx])
    #alpha1 = np.arctan2(DLON,DLAT) # working
    ## rotate the vectors to lat-lon coordinates
    #print('... rotating actual vectors')
    #urot1 = -u*np.cos(alpha1) + v*np.sin(alpha1) # working
    #vrot1 =  v*np.cos(alpha1) + u*np.sin(alpha1) # working 
    
    ## a similar solution but with the other index...
    #print('... computing rotation factor (2)')
    #for jy in range(0,ny):
        #jyp1 = np.min([jy+1,ny-1])
        #jym1 = np.max([jy-1,0])
        #DLAT[jy,:] = (lat[jyp1,:]-lat[jym1,:])
        #DLON[jy,:] = (lon[jyp1,:]-lon[jym1,:])*np.cos(lat[jy,:])
    #alpha2 = np.arctan2(DLAT,DLON) # seems to be pi-alpha1
    ## rotate the vectors to lat-lon coordinates
    #print('... rotating actual vectors')
    #urot2 =   u*np.cos(alpha2) + v*np.sin(alpha2)
    #vrot2 =  -v*np.cos(alpha2) + u*np.sin(alpha2)
    ## taking the mean of both ...
    #urot = ( urot1 + urot2 )/2.0
    #vrot = ( vrot1 + vrot2 )/2.0
    #return urot,vrot


### FROM A fortran SCRIPT MPI_OM Rotate vectors.
## but is this one correct? It seems to be doing funny things with unit winds.
#def rotuvb2(u,v,lon,lat):
   #PI=np.pi
   #u_lon=np.nan*u
   #v_lat=np.nan*v
   #IX = lon.shape[0]
   #IY = lon.shape[1]
   #for i in range(0,IX):
       #print(i)
       #ip1 = np.min([i+1,IX-1])
       #im1 = np.max([0,i-1])
       #for j in range(0,IY):
           #jp1 = np.min([j+1,IY-1])
           #jm1 = np.max([0,j-1])
           
           ##! difference in latitudes
           #dlat_i = lat[ip1,j] - lat[im1,j]
           #dlat_j = lat[i,jp1] - lat[i,jm1]
           
           ##! difference in longitudes
           #dlon_i = lon[ip1,j] - lon[im1,j]
           #dlon_j = lon[i,jp1] - lon[i,jm1]

           #lat_factor = np.cos(lat[i,j])
           #dlon_i = dlon_i * lat_factor
           #dlon_j = dlon_j * lat_factor

           ##! projection by scalar product
           ##! ----------------------------
           #u_lon[:,i,j] = u[:,i,j]*dlon_i + v[:,i,j]*dlat_i
           #v_lat[:,i,j] = u[:,i,j]*dlon_j + v[:,i,j]*dlat_j

           #dist_i = np.sqrt(dlon_i**2+dlat_i**2)
           #dist_j = np.sqrt(dlon_j**2+dlat_j**2)
           #if (dist_i != 0. and dist_j != 0.):
               #u_lon[:,i,j] = u_lon[:,i,j]/dist_i
               #v_lat[:,i,j] = v_lat[:,i,j]/dist_j
           #else:
               #u_lon[:,i,j] = 0.0
               #v_lat[:,i,j] = 0.0
   #return u_lon,v_lat




