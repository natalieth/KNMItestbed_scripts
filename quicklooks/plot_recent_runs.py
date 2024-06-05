import os
import sys
import matplotlib.pyplot as plt
import xarray as xr 
import datetime 
import pandas as pd
import numpy as np 

import quicklooks as ql 

from therm_functions import *
from check_avail import * 

def create_time_df(date,what,D_exp):
    DALES_path = "/nobackup/users/theeuwes/testbed/DALES/{0:04d}/{1:02d}/{2:02d}/".format(date.year,date.month,date.day)
    OBS_path = "/nobackup/users/theeuwes/tmp_data/"
    HARM_path = "/nobackup/users/theeuwes/testbed/HARMONIE/nc/"
    
    dummy_times = [datetime.datetime(date.year,date.month,date.day,0,0)+datetime.timedelta(seconds=3600*i) for i in range(24)]
    dummy_var = [np.nan for i in range(len(dummy_times))]

    if what == "EB":
        # Get DALES EB componants
        try:
            DALES = xr.open_dataset(DALES_path + "tmser.{0:03d}.{1:04d}{2:02d}{3:02d}.nc".format(D_exp,date.year,date.month,date.day))
            DALES['time'] = [datetime.datetime(date.year,date.month,date.day,0,0) 
                         + datetime.timedelta(seconds=int(x)) for x in DALES.time.values]
            DALES = DALES.resample(time='1H').nearest()
            DALES_ar = {'H_D': DALES.H.values, 'LE_D': DALES.LE.values, 'G_D': DALES.G.values*-1.}
            df1 = pd.DataFrame(DALES_ar, index = DALES.time.values)

        except FileNotFoundError:
            DALES_ar = {'H_D': dummy_var, 'LE_D': dummy_var, 'G_D': dummy_var}
            df1 = pd.DataFrame(DALES_ar, index = dummy_times)

        # Get OBS EB componants
        try:
            OBS = xr.open_dataset(OBS_path + 
                    "cesar_surface_flux_la1_t10_v1.0_{0:04d}{1:02d}{2:02d}.nc".format(date.year,date.month,date.day))
            OBS = OBS.resample(time='1H').nearest()
            OBS_ar = {'H_O': OBS.HSON.values,'LE_O': OBS.LEED.values, 'G_O': OBS.FG0.values }
            df2 = pd.DataFrame(OBS_ar, index = OBS.time.values)
        except FileNotFoundError:
            OBS_ar = {'H_O': dummy_var,'LE_O': dummy_var, 'G_O': dummy_var }
            df2 = pd.DataFrame(OBS_ar, index = dummy_times)


        # Get HARMONIE EB componants
        try:
            HAP2 = xr.open_dataset(HARM_path + "CABAUW_HAP2_{0:04d}{1:02d}{2:02d}00.nc".format(date.year,date.month,date.day))
            # HAP2 = HAP2.resample(time='10min').mean()
            HAP2_ar = {'H_HAP2': HAP2.hfssde_surf.values}
            df3 = pd.DataFrame(HAP2_ar, index = HAP2.time.values)
        except FileNotFoundError:
            HAP2_ar = {'H_HAP2': dummy_var}
            df3 = pd.DataFrame(HAP2_ar, index = dummy_times)
        try:
            HA43 = xr.open_dataset(HARM_path + "CABAUW_HA43_{0:04d}{1:02d}{2:02d}00.nc".format(date.year,date.month,date.day))
            # HA43 = HAP2.resample(time='10min').mean()
            HA43_ar = {'H_HA43': HA43.hfssde_surf.values}
            df4 = pd.DataFrame(HA43_ar, index = HA43.time.values)
        except FileNotFoundError:
            HA43_ar = {'H_HA43': dummy_var}
            df4 = pd.DataFrame(HA43_ar, index = dummy_times)

        df  = pd.concat([df1, df2, df3, df4], axis=1)

    elif what=="RB":
        # Get DALES RB componants
        try:
            DALES = xr.open_dataset(DALES_path + "profiles.{0:03d}.{1:04d}{2:02d}{3:02d}.nc".format(D_exp,date.year,date.month,date.day))
            DALES['time'] = [datetime.datetime(date.year,date.month,date.day,0,0)
                         + datetime.timedelta(seconds=int(x)) for x in DALES.time.values]
#            DALES = DALES.resample(time='10min').mean()
            DALES_ar = {'SWD_D': DALES.swd.values[:,0]*-1., 'SWU_D': DALES.swu.values[:,0], 'LWD_D': DALES.lwd.values[:,0]*-1.,
                        'LWU_D': DALES.lwu.values[:,0]}

            df1 = pd.DataFrame(DALES_ar, index = DALES.time.values)

        except FileNotFoundError:
            DALES_ar = {'SWD_D': dummy_var, 'SWU_D': dummy_var, 'LWD_D': dummy_var,
                        'LWU_D': dummy_var}

            df1 = pd.DataFrame(DALES_ar, index = dummy_times)

        # Get OBS RB componants
        try: 
            OBS = xr.open_dataset(OBS_path +
                        "cesar_surface_radiation_la1_t10_v1.0_{0:04d}{1:02d}{2:02d}.nc".format(date.year,date.month,date.day))
#            OBS = OBS.resample(time='10min').mean()
            OBS_ar = {'SWD_O': OBS.SWD.values,'SWU_O': OBS.SWU.values, 'LWD_O': OBS.LWD.values, 'LWU_O':OBS.LWU.values }
            df2 = pd.DataFrame(OBS_ar, index = OBS.time.values)
        except FileNotFoundError:
            OBS_ar = {'SWD_O': dummy_var,'SWU_O': dummy_var, 'LWD_O': dummy_var, 'LWU_O':dummy_var }
            df2 = pd.DataFrame(OBS_ar, index = dummy_times)
        try:
            print(OBS_path + "MSGCPP_cabauw{0:04d}-{1:02d}-{2:02d}.nc".format(date.year,date.month,date.day))
            OBS_MSG = xr.open_dataset(OBS_path + "MSGCPP_cabauw{0:04d}-{1:02d}-{2:02d}.nc".format(date.year,date.month,date.day))
            OBS2_ar = {'SWD_MSG': OBS_MSG.sds.values, 'SWD_max_MSG': OBS_MSG.sds_max.values, 'SWD_min_MSG': OBS_MSG.sds_min.values}
            df3 = pd.DataFrame(OBS2_ar, index = OBS_MSG.time.values)
        except FileNotFoundError:
            OBS2_ar = {'SWD_MSG': dummy_var, 'SWD_max_MSG': dummy_var,'SWD_min_MSG': dummy_var,}
            df3 = pd.DataFrame(OBS2_ar, index = dummy_times)

        df  = pd.concat([df1, df2, df3], axis=1)
    elif what=="cwp":
        try: 
            DALES = xr.open_dataset(DALES_path + "tmser.{0:03d}.{1:04d}{2:02d}{3:02d}.nc".format(D_exp,date.year,date.month,date.day))
            DALES['time'] = [datetime.datetime(date.year,date.month,date.day,0,0)
                         + datetime.timedelta(seconds=int(x)) for x in DALES.time.values]
#            DALES = DALES.resample(time='10min').mean()
            DALES = DALES.resample(time='1H').nearest()
            DALES_ar = {'CWP_D': DALES.lwp_bar.values*1000,'CWP_max_D': DALES.lwp_max.values*1000.}
            df1 = pd.DataFrame(DALES_ar, index = DALES.time.values)
        except FileNotFoundError:
            DALES_ar = {'CWP_D': dummy_var,'CWP_max_D': dummy_var}
            df1 = pd.DataFrame(DALES_ar, index = dummy_times)

        try: 
            OBS_MSG = xr.open_dataset(OBS_path +
                        "MSGCPP_cabauw{0:04d}-{1:02d}-{2:02d}.nc".format(date.year,date.month,date.day))
            OBS_MSG = OBS_MSG.resample(time='1H').nearest()
            OBS_ar = {'CWP_MSG': OBS_MSG.cwp.values,'CWP_max_MSG': OBS_MSG.cwp_max.values,'CWP_min_MSG': OBS_MSG.cwp_min.values}
            df2 = pd.DataFrame(OBS_ar, index = OBS_MSG.time.values)
        except FileNotFoundError:
            OBS_ar = {'CWP_MSG': dummy_var,'CWP_max_MSG': dummy_var,'CWP_min_MSG': dummy_var }
            df2 = pd.DataFrame(OBS_ar, index = dummy_times)
        try:
            OBS_CN = xr.open_dataset(OBS_path +
                        "{0:04d}{1:02d}{2:02d}_cabauw_hatpro.nc".format(date.year,date.month,date.day))
            OBS_CN = OBS_CN.resample(time='1min').mean()
            print(OBS_CN)
            OBS_CNA = OBS_CN.resample(time='1H').mean()
            OBS_CNB = OBS_CN.resample(time='1H').min()
            OBS_CNC = OBS_CN.resample(time='1H').max()
            OBS_ar = {'CWP_CN': OBS_CNA.lwp.values,'CWP_max_CN': OBS_CNC.lwp.values,'CWP_min_CN': OBS_CNB.lwp.values}
            df2b = pd.DataFrame(OBS_ar, index = OBS_CNA.time.values)
        except FileNotFoundError:
            OBS_ar = {'CWP_CN': dummy_var,'CWP_max_CN': dummy_var,'CWP_min_CN': dummy_var }
            df2b = pd.DataFrame(OBS_ar, index = dummy_times)

        try: 
            HAP2 = xr.open_dataset(HARM_path + "CABAUW_HAP2_{0:04d}{1:02d}{2:02d}00.nc".format(date.year,date.month,date.day))
            HAP2_ar = {'CWP_HAP2': HAP2.cwp.values*1000.} 
            df3 = pd.DataFrame(HAP2_ar, index = HAP2.time.values)
        except FileNotFoundError:
            HAP2_ar = {'CWP_HAP2':dummy_var}
            df3 = pd.DataFrame(HAP2_ar, index = dummy_times)
        try:
            HA43 = xr.open_dataset(HARM_path + "CABAUW_HA43_{0:04d}{1:02d}{2:02d}00.nc".format(date.year,date.month,date.day))
            HA43_ar = {'CWP_HA43': HA43.cwp.values*1000.}
            df4 = pd.DataFrame(HA43_ar, index = HA43.time.values)
        except FileNotFoundError:
            HA43_ar = {'CWP_HA43':dummy_var}
            df4 = pd.DataFrame(HA43_ar, index = dummy_times)

        df  = pd.concat([df1, df2, df2b, df3, df4], axis=1)
#        df  = pd.concat([df1, df2, df3], axis=1)
        print(df)


    elif what=="2_10mDIAG":
        # Get DALES 10 and 2 m diagnostics
        try:
            DALES = xr.open_dataset(DALES_path + "profiles.{0:03d}.{1:04d}{2:02d}{3:02d}.nc".format(D_exp,date.year,date.month,date.day))
            DALES['time'] = [datetime.datetime(date.year,date.month,date.day,0,0)
                         + datetime.timedelta(seconds=int(x)) for x in DALES.time.values]
            DALES = DALES[['u','v']]
            DALES = DALES.resample(time='1H').mean()
            DALES['ws'] = np.sqrt(DALES.u**2 + DALES.v**2)
            DALES['wd'] = 180.+np.arctan2(DALES.u,DALES.v)*180./np.pi
            DALES_ar = {'WS_D': DALES.ws.values[:,0], 'WD_D': DALES.wd.values[:,0]}
            df1 = pd.DataFrame(DALES_ar, index = DALES.time.values)
        except FileNotFoundError:
            DALES_ar = {'WS_D': dummy_var, 'WD_D': dummy_var}
            df1 = pd.DataFrame(DALES_ar, index = dummy_times)


        # Get OBS 10 and 2 m diagnostics
        try: 
            OBS = xr.open_dataset(OBS_path +
                        "cesar_tower_meteo_la1_t10_v1.2_{0:04d}{1:02d}{2:02d}.nc".format(date.year,date.month,date.day))
            OBS = OBS[['D','F','TA','TD']]
            OBS= OBS.astype(float)
            OBS['U'] = -1.*OBS.F*np.sin(OBS.D*np.pi/180.)
            OBS['V'] = -1.*OBS.F*np.cos(OBS.D*np.pi/180.)

            OBS = OBS.resample(time='1H').mean()
    
            OBS['ws'] = np.sqrt(OBS.U**2 + OBS.V**2)
            OBS['wd'] = 180.+np.arctan2(OBS.U,OBS.V)*180./np.pi
       
            OBS_ar = {'T2_O': OBS.TA.values[:,-1],'TD2_O': OBS.TD.values[:,-1], 'WS_O': OBS.ws.values[:,-2], 'WD_O':OBS.wd.values[:,-2]}
            df2 = pd.DataFrame(OBS_ar, index = OBS.time.values)
        except FileNotFoundError:
            OBS_ar = {'T2_O': dummy_var,'TD2_O': dummy_var, 'WS_O': dummy_var, 'WD_O':dummy_var}
            df2 = pd.DataFrame(OBS_ar, index = dummy_times)

        # Get HARMONIE EB componants
        try:
            HAP2 = xr.open_dataset(HARM_path + "CABAUW_HAP2_{0:04d}{1:02d}{2:02d}00.nc".format(date.year,date.month,date.day))
            HAP2 = HAP2.resample(time='1H').mean()
            HAP2['ws'] = np.sqrt(HAP2.ua**2 + HAP2.va**2)
            HAP2['wd'] = 180.+np.arctan2(HAP2.ua,HAP2.va)*180./np.pi
            
            HAP2_ar = {'TD2_HAP2': HAP2.td_2m.values, 'T2_HAP2': HAP2.ta_2m.values, 
                       'WS_HAP2': HAP2.ws.values[:,-1], 'WD_HAP2': HAP2.wd.values[:,-1]} ## wind speed and direction now at 12 m not rotated!!
            df3 = pd.DataFrame(HAP2_ar, index = HAP2.time.values)
        except FileNotFoundError:
            HAP2_ar = {'TD2_HAP2': dummy_var, 'T2_HAP2': dummy_var,
                       'WS_HAP2': dummy_var, 'WD_HAP2': dummy_var} 
            df3 = pd.DataFrame(HAP2_ar, index = dummy_times)

        try:
            HA43 = xr.open_dataset(HARM_path + "CABAUW_HA43_{0:04d}{1:02d}{2:02d}00.nc".format(date.year,date.month,date.day))
            HA43 = HA43.resample(time='1H').mean()
            HA43['ws'] = np.sqrt(HA43.ua**2 + HA43.va**2)
            HA43['wd'] = 180.+np.arctan2(HA43.ua,HA43.va)*180./np.pi

            HA43_ar = {'TD2_HA43': HA43.td_2m.values, 'T2_HA43': HA43.ta_2m.values,  
                       'WS_HA43': HA43.ws.values[:,-1], 'WD_HA43': HA43.wd.values[:,-1]} ## wind speed and direction now at 12 m not rotated!!
            df4 = pd.DataFrame(HA43_ar, index = HA43.time.values)
        except FileNotFoundError:
            HA43_ar = {'TD2_HA43': dummy_var, 'T2_HA43': dummy_var,
                       'WS_HA43': dummy_var, 'WD_HA43': dummy_var}
            df4 = pd.DataFrame(HA43_ar, index = dummy_times)

        df  = pd.concat([df1, df2, df3, df4], axis=1)

    return df

def create_prof_df(date,D_exp,incl_turb=False):
    DALES_path = "/nobackup/users/theeuwes/testbed/DALES/{0:04d}/{1:02d}/{2:02d}/".format(date.year,date.month,date.day)
    OBS_path = "/nobackup/users/theeuwes/tmp_data/"
    LID_path =  "/usr/people/theeuwes/sensordata/WLIDARNET/Windcube200S/testbed/"
    HARM_path = "/nobackup/users/theeuwes/testbed/HARMONIE/nc/"

    df_ar = []
    df_str = []
    try:
        df_ar.append(xr.open_dataset(OBS_path +
                            "cesar_tower_meteo_la1_t10_v1.2_{0:04d}{1:02d}{2:02d}.nc".format(date.year,date.month,date.day)))
        df_str.append("OBS")
    except FileNotFoundError:
        dummy = 0
    try:
        df_ar.append(xr.open_dataset(OBS_path +
                            "ceilonet_chm15k_backsct_la1_t12s_v1.0_06348_A{0:04d}{1:02d}{2:02d}.nc".format(date.year,date.month,date.day)))
        df_str.append("OBS_ceil")
    except FileNotFoundError:
        dummy = 0
    try: 
        f = pd.read_csv(LID_path + 
                     "{0:04d}{1:02d}{2:02d}_windspeed.csv".format(date.year,date.month,date.day), index_col=0,parse_dates=True)
        f2 = pd.read_csv(LID_path +
                     "{0:04d}{1:02d}{2:02d}_winddirection.csv".format(date.year,date.month,date.day), index_col=0,parse_dates=True)
        z = [np.float(f.columns[i][1:]) for i in range(len(f.columns))]
        arr = np.array(f.values)
        arr2 = np.array(f2.values)
        df_ar.append(xr.Dataset(data_vars=dict(WS =(["time","z"],arr),WD = (["time","z"],arr2)),coords = dict(time=f.index,z=np.array(z))))

        df_str.append("OBS_lid")
    except FileNotFoundError: 
        dummy = 0
    if incl_turb :
        try: 
            if date < datetime.date(2022,6,21): 
                f3 = pd.read_csv(LID_path +
                     "{0:04d}{1:02d}{2:02d}_vertical_variance_30T.csv".format(date.year,date.month,date.day), index_col=0,parse_dates=True)
            else:
                f3 = pd.read_csv(LID_path +
                     "{0:04d}{1:02d}{2:02d}_vertical_variance_60T.csv".format(date.year,date.month,date.day), index_col=0,parse_dates=True)
            arr3 = np.array(f3.values)
            z = [np.float(f3.columns[i][3:]) for i in range(len(f3.columns))]
            mh = np.zeros([len(f3.index)]) * np.nan

            for t in range(len(f3.index)):
                mh[t] = MixingHeight(arr3[t,:],z)

            df_ar.append(xr.Dataset(data_vars=dict(WVAR = (["time","z"],arr3), MH = (["time"],mh)),
                  coords = dict(time=f3.index,z=np.array(z))))

            df_str.append("OBS_lidT")
        except FileNotFoundError:
            dummy = 0
    try: 
        DALES = xr.open_dataset(DALES_path + "profiles.{0:03d}.{1:04d}{2:02d}{3:02d}.nc".format(D_exp,date.year,date.month,date.day), engine='netcdf4')
        DALES['time'] = [datetime.datetime(date.year,date.month,date.day,0,0)
                         + datetime.timedelta(seconds=int(x)) for x in DALES.time.values]
        DALES['TA'] = xr.DataArray(np.array([thl_to_ta(DALES.thl.values[:,k],DALES.ql.values[:,k]*1E-3,DALES.presh.values[:,k])
                                   for k in np.arange(len(DALES.zt))]).T, coords=[("time", DALES.time.values), ("zt", DALES.zt.values)])
        DALES['ws'] = np.sqrt(DALES.u**2 + DALES.v**2) 
        DALES['wd'] = 180.+np.arctan2(DALES.u,DALES.v)*180./np.pi
        if incl_turb :
            DALES['MH'] = xr.DataArray(np.array([MixingHeight(DALES.w2r.values[t,:]+ DALES.w2s.values[t,:], DALES.zm.values)
                                   for t in np.arange(len(DALES.time))]), coords=[("time", DALES.time.values)])
        df_ar.append(DALES)
        df_str.append("DALES")
    except FileNotFoundError:
        dummy = 0
    try:
        HAP2 = xr.open_dataset(HARM_path + "CABAUW_HAP2_{0:04d}{1:02d}{2:02d}00.nc".format(date.year,date.month,date.day))
        HAP2['ws'] = np.sqrt(HAP2.ua**2 + HAP2.va**2)
        HAP2['wd'] = 180.+np.arctan2(HAP2.ua,HAP2.va)*180./np.pi
        if incl_turb :
            HAP2['sigma2w'] = (2./3.)*(HAP2.tke)   # assuming isotropic turbulence
            HAP2['MH'] = xr.DataArray(np.array([MixingHeight(HAP2.sigma2w.values[t,::-1], HAP2.z.values[t,::-1])
                                   for t in np.arange(len(HAP2.time))]), coords=[("time", HAP2.time.values)])
        df_ar.append(HAP2)
        df_str.append("HAP2")
    except FileNotFoundError:
        dummy = 0
    try:
        HA43 = xr.open_dataset(HARM_path + "CABAUW_HA43_{0:04d}{1:02d}{2:02d}00.nc".format(date.year,date.month,date.day))
        HA43['ws'] = np.sqrt(HA43.ua**2 + HA43.va**2)
        HA43['wd'] = 180.+np.arctan2(HA43.ua,HA43.va)*180./np.pi
        if incl_turb :
            HA43['sigma2w'] = (2./3.)*(HA43.tke)   # assuming isotropic turbulence
            HA43['MH'] = xr.DataArray(np.array([MixingHeight(HA43.sigma2w.values[t,::-1], HA43.z.values[t,::-1])
                                   for t in np.arange(len(HA43.time))]), coords=[("time", HA43.time.values)])
        df_ar.append(HA43)
        df_str.append("HA43")
    except FileNotFoundError:
        dummy = 0
#    try:  # discontinued 
#        HAP3 = xr.open_dataset(HARM_path + "CABAUW_HAP3_{0:04d}{1:02d}{2:02d}00.nc".format(date.year,date.month,date.day))
#        HAP3['ws'] = np.sqrt(HAP3.ua**2 + HAP3.va**2)
#        HAP3['wd'] = 180.+np.arctan2(HAP3.ua,HAP3.va)*180./np.pi
#        df_ar.append(HAP3)
#        df_str.append("HAP3")
#    except FileNotFoundError:
#        dummy = 0
    try:
        HAP1 = xr.open_dataset(HARM_path + "CABAUW_HAP1_{0:04d}{1:02d}{2:02d}00.nc".format(date.year,date.month,date.day))
        HAP1['ws'] = np.sqrt(HAP1.ua**2 + HAP1.va**2)
        HAP1['wd'] = 180.+np.arctan2(HAP1.ua,HAP1.va)*180./np.pi
        if incl_turb :
            HAP1['sigma2w'] = (2./3.)*(HAP1.tke)   # assuming isotropic turbulence
            HAP1['MH'] = xr.DataArray(np.array([MixingHeight(HAP1.sigma2w.values[t,::-1], HAP1.z.values[t,::-1])
                                   for t in np.arange(len(HAP1.time))]), coords=[("time", HAP1.time.values)])
        df_ar.append(HAP1)
        df_str.append("HAP1")
    except FileNotFoundError:
        dummy = 0

    return (df_ar, df_str)



if __name__ == '__main__':
    # use 3 days ago?? 
#    date = datetime.date.today() - datetime.timedelta(days=4)
    start = datetime.date(2022,6,19)
#    end = datetime.date(2022,5,3)
    start =  datetime.date.today() - datetime.timedelta(days=30) #datetime.datetime(2021,8,1)
    end =    datetime.date.today() - datetime.timedelta(days=1)  #datetime.datetime(2021,8,26)
    dates = [start + datetime.timedelta(days=i) for i in range((end-start).days + 1)]
#    date = datetime.date(2021,8,6)
    for date in dates:
        exp = 1
        hh = 0
        print(f"working on {date} and experiment {exp}")

    ### check if DALES data is there if not get from BULL
        check_DALES_exists(date,exp)

    ### check if CABAUW data is there, if not get from KDP
        check_cabauw_exists(date)

    ### add HARMONIE when available 
        check_HARMONIE_exists(date,hh)
    
    ### plot the energy balance partitioning

        df_EB = create_time_df(date,"EB",exp)
        ql.plot_EB(date,df_EB)

        df_RB = create_time_df(date,"RB",exp)
        ql.plot_RB(date,df_RB)

        df_DIA = create_time_df(date,"2_10mDIAG",exp)
        ql.plot_DIA(date,df_DIA)

        df_CWP = create_time_df(date,"cwp",exp)
        ql.plot_CWP(date,df_CWP)
    
        df_ar, df_str = create_prof_df(date,exp,incl_turb=True)
        ql.plot_TPROF(date,df_ar,df_str)
        ql.plot_QPROF(date,df_ar,df_str)
        ql.plot_QLPROF(date,df_ar,df_str)
        ql.plot_WDPROF(date,df_ar,df_str)
        ql.plot_WSPROF(date,df_ar,df_str)
        ql.plot_WSBLPROF(date,df_ar,df_str)
        ql.plot_WDBLPROF(date,df_ar,df_str)
        ql.plot_TURBPROF(date,df_ar,df_str)

 #   date = datetime.date(2021,9,5)
 #   df_ar,df_str = create_prof_df(date,exp,incl_turb=True)
 #   ql.plot_TURBPROF(date,df_ar,df_str)
