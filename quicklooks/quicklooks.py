import matplotlib.pyplot as plt
import xarray as xr 
import datetime
import numpy as np

def make_colormap(colors):
    from matplotlib.colors import LinearSegmentedColormap, ColorConverter
    from numpy import sort
    
    z  = np.array(sorted(colors.keys()))
    n  = len(z)
    z1 = min(z)
    zn = max(z)
    x0 = (z - z1) / (zn - z1)
    
    CC = ColorConverter()
    R = []
    G = []
    B = []
    for i in range(n):
        Ci = colors[z[i]]
        if type(Ci) == str:
            RGB = CC.to_rgb(Ci)
        else:
            RGB = Ci
        R.append(RGB[0])
        G.append(RGB[1])
        B.append(RGB[2])
    
    cmap_dict = {}
    cmap_dict['red']   = [(x0[i],R[i],R[i]) for i in range(len(R))]
    cmap_dict['green'] = [(x0[i],G[i],G[i]) for i in range(len(G))]
    cmap_dict['blue']  = [(x0[i],B[i],B[i]) for i in range(len(B))]
    mymap = LinearSegmentedColormap('mymap',cmap_dict)
    return mymap

test1 = make_colormap({0.:'#ffffff',0.2:'cornflowerblue',0.4:'limegreen',0.6:'yellow',0.8:'darkorange',1.:'crimson'})

def plot_EB(date,df):
    import matplotlib.dates as md

    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/EB/"
    figs, ax = plt.subplots(1,1)
    ax.set_title(date)
    ax.scatter(df.index,df.H_O,marker='o',color='k', label ='H obs')
    ax.plot(df.index,df.H_D,color='k',label = 'H DALES')
    ax.plot(df.index,df.H_HAP2,color='k',label = 'H HAP2',ls=':')
    ax.plot(df.index,df.H_HA43,color='k',label = 'H HA43',ls='-.')

    ax.scatter(df.index,df.LE_O,marker='o',color='forestgreen', label ='LE obs')
    ax.plot(df.index,df.LE_D,color='forestgreen',label = 'LE DALES')
    plt.legend()
    ax.set_xlabel('Time')
    ax.set_ylabel('Flux [W m$^{-2}$]')
    date_form = md.DateFormatter("%H:%M")
    ax.xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax.set_xlim(datet,datet+datetime.timedelta(days=1))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    print(f"Saving EB plot for {date} to {ql_dir}")

    plt.savefig("{0}EB_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))

def plot_RB(date,df):
    import matplotlib.dates as md

    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/RB/"
    figs, ax = plt.subplots(2,1,figsize=(6,8),sharex=True)
    ax[0].set_title(date)
    ax[0].scatter(df.index,df.SWD_O,marker='o',color='k', label ='SWD obs')
    ax[0].plot(df.index,df.SWD_D,color='k',label = 'SWD DALES')
    ax[0].scatter(df.index,df.SWU_O,marker='*',color='forestgreen', label ='SWU obs')
    ax[0].plot(df.index,df.SWU_D,color='forestgreen',label = 'SWU DALES')
    #plt.legend()
    ax[0].set_ylabel('Shortwave radiation [W m$^{-2}$]')
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)

    ax[1].scatter(df.index,df.LWD_O,marker='o',color='k', label ='LWD obs')
    ax[1].plot(df.index,df.LWD_D,color='k',label = 'LWD DALES')
    ax[1].scatter(df.index,df.LWU_O,marker='*',color='forestgreen', label ='LWU obs')
    ax[1].plot(df.index,df.LWU_D,color='forestgreen',label = 'LWU DALES')
    plt.legend()
    ax[1].set_xlabel('Time')
    ax[1].set_ylabel('Longwave radiation [W m$^{-2}$]')
    date_form = md.DateFormatter("%H:%M")
    ax[1].xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax[1].set_xlim(datet,datet+datetime.timedelta(days=1))
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    print(f"Saving RB plot for {date} to {ql_dir}")

    plt.savefig("{0}RB_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))

def plot_DIA(date,df):
    # plot near-surface diagnostics
    import matplotlib.dates as md

    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/DIAG/"
    figs, ax = plt.subplots(4,1,sharex=True,figsize=(6,18))
    ax[0].set_title(date)
    ax[0].scatter(df.index,df.T2_O,marker='o',color='k', label ='T2 obs')
    #ax[0].plot(df.index,df.T2_D,color='k',label = 'H DALES')  # should still be diagnosed
    ax[0].plot(df.index,df.T2_HAP2,color='k',label = 'T2 HAP2',ls=':')
    ax[0].plot(df.index,df.T2_HA43,color='k',label = 'T2 HA43',ls='-.')
    plt.legend()
    ax[0].set_ylabel('2m Temperature [K]')
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)

    ax[1].scatter(df.index,df.TD2_O,marker='o',color='k', label ='TD2 obs')
    #ax[1].plot(df.index,df.T2_D,color='k',label = 'H DALES')  # should still be diagnosed
    ax[1].plot(df.index,df.TD2_HAP2,color='k',label = 'TD2 HAP2',ls=':')
    ax[1].plot(df.index,df.TD2_HA43,color='k',label = 'TD2 HA43',ls='-.')
    plt.legend()
    ax[1].set_ylabel('2m Dew Point Temperature [K]')
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)

    ax[2].scatter(df.index,df.WS_O,marker ='o',color='k',label = 'U10 DALES')  # should still be diagnosed
    ax[2].plot(df.index,df.WS_D,color='k',label = 'U10 DALES')  # should still be diagnosed
    ax[2].plot(df.index,df.WS_HAP2,color='k',label = 'U10 HAP2',ls=':')
    ax[2].plot(df.index,df.WS_HA43,color='k',label = 'U10 HA43',ls='-.')
    plt.legend()
    ax[2].set_ylabel('10m wind speed [m s-1]')
    ax[2].spines['right'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
    
    ax[3].scatter(df.index,df.WD_O,marker='o',color='k', label = 'obs')
    ax[3].plot(df.index,df.WD_D,color='k',label = 'DALES')  # should still be diagnosed
    ax[3].plot(df.index,df.WD_HAP2,color='k',label = 'HAP2',ls=':')
    ax[3].plot(df.index,df.WD_HA43,color='k',label = 'HA43',ls='-.')
    plt.legend()
    ax[3].set_xlabel('Time')
    ax[3].set_ylabel('10m wind direction [$^{\circ}$]')
    date_form = md.DateFormatter("%H:%M")
    ax[3].xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax[3].set_xlim(datet,datet+datetime.timedelta(days=1))
    ax[3].spines['right'].set_visible(False)
    ax[3].spines['top'].set_visible(False)
    print(f"Saving DIAG plot for {date} to {ql_dir}")

    plt.savefig("{0}DIAG_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))

def plot_CWP(date,df):
    import matplotlib.dates as md

    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/CWP/"
    figs, ax = plt.subplots(1,1)
    print(df.CWP_MSG.values[:])
    ax.set_title(date)
    ax.scatter(df.index,df.CWP_MSG,marker='o',color='crimson',edgecolor='k',label ='MSG CPP')
    ax.fill_between(df.index,df.CWP_min_MSG,df.CWP_max_MSG,facecolor='crimson',alpha=0.1,edgecolor=None)
    ax.scatter(df.index,df.CWP_CN,marker='v',color='violet',edgecolor='k',label ='Microwave Radiometer')
    ax.fill_between(df.index,df.CWP_min_CN,df.CWP_max_CN,facecolor='violet',alpha=0.1,edgecolor=None)
    ax.plot(df.index,df.CWP_D,color='royalblue',label = 'DALES')
    ax.fill_between(df.index,df.CWP_D,df.CWP_max_D,facecolor='royalblue',alpha=0.1,edgecolor=None)
    ax.plot(df.index,df.CWP_HAP2,color='k',label = 'HAP2',ls=':')
    ax.plot(df.index,df.CWP_HA43,color='k',label = 'HA43',ls='-.')

    plt.legend()
    ax.set_xlabel('Time')
    ax.set_ylabel('Cloud water/ice path [g m$^{-3}$]')
    date_form = md.DateFormatter("%H:%M")
    ax.xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax.set_xlim(datet,datet+datetime.timedelta(days=1))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    print(f"Saving CWP plot for {date} to {ql_dir}")

    plt.savefig("{0}CWP_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))



def plot_TPROF(date,df_arr,df_str):
    import matplotlib.dates as md
    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/TPROF/"
    pltmax = 5
    obs_there = False

    fig, ax = plt.subplots(pltmax,1,sharey=True,sharex=True, figsize=(7,15))
    tmin = 270 
    tmax = 300
    fig.suptitle(date)
    for df,mo in zip(df_arr,df_str) :
        if mo == 'OBS':
            tmin = np.nanmin(df.TA.values)-1.
            tmax = np.nanmax(df.TA.values)+1.
            obs_there = True
            im = ax[0].pcolormesh(df.time,df.z,df.TA.values.T,cmap='Spectral_r',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'DALES':
            #print(df)
            if (~obs_there==False):
                tmin = np.nanmin(df.TA.values[:,0:20])-1
                tmax = np.nanmax(df.TA.values[:,0:20])+1.
            im = ax[1].pcolormesh(df.time.values,df.zt.values,df.TA.values.T,cmap='Spectral_r',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'HAP2':
            #df = df.where(df.z<300.,drop=True)
            im = ax[2].pcolormesh(df.time.values,df.z.values[0,:],df.ta.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HA43':
            #df = df.where(df.z<300.,drop=True)
            im = ax[3].pcolormesh(df.time.values,df.z.values[0,:],df.ta.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')
        #elif mo == 'HAP3':
        #    #df = df.where(df.z<300.,drop=True)
        #    im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.ta.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HAP1':
            #df = df.where(df.z<300.,drop=True)
            im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.ta.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')
    
    for i in range(pltmax):
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].set_ylabel(r'$z$ [m]')

    ax[0].set_ylim(0.,300.)
    
    date_form = md.DateFormatter("%H:%M")
    ax[3].xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax[3].set_xlim(datet,datet+datetime.timedelta(days=1))
    ax[4].set_xlabel('Time')
    plt.tight_layout()
    fig.subplots_adjust(top=0.925, right=0.8,hspace=0.2)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax,label=r'$T_a$ [K]')
    ax[0].set_title("Tower observations")
    ax[1].set_title("DALES")
    ax[2].set_title("HARMONIE-AROME HAP2")
    ax[3].set_title("HARMONIE-AROME HA43")
   # ax[4].set_title("HARMONIE-AROME HAP3")
    ax[4].set_title("HARMONIE-AROME HAP1")
    print(f"Saving TPROF plot for {date} to {ql_dir}")

    plt.savefig("{0}TPROF_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))


def plot_QPROF(date,df_arr,df_str):
    import matplotlib.dates as md
    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/QPROF/"
    pltmax = 5
    obs_there = False
    fig, ax = plt.subplots(pltmax,1,sharey=True,sharex=True,figsize=(7,15))
    tmin = 0
    tmax = 15
    fig.suptitle(date)
    for df,mo in zip(df_arr,df_str) :
        if mo == 'OBS':
            tmin = np.nanmin(df.Q.values)-0.5
            tmax = np.nanmax(df.Q.values)+0.5
            obs_there = True
            im = ax[0].pcolormesh(df.time,df.z,df.Q.values.T,cmap='GnBu',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'DALES':
            if (obs_there==False):
                tmin = np.nanmin(df.qt.values * 1E3)-0.5
                tmax = np.nanmax(df.qt.values * 1E3)+0.5
            im = ax[1].pcolormesh(df.time.values,df.zt.values,df.qt.values.T*1E3,cmap='GnBu',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'HAP2':
            #df = df.where(df.z<300.,drop=True)
            im = ax[2].pcolormesh(df.time.values,df.z.values[0,:],df.hus.values.T*1E3,cmap='GnBu',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HA43':
            #df = df.where(df.z<300.,drop=True)
            im = ax[3].pcolormesh(df.time.values,df.z.values[0,:],df.hus.values.T*1E3,cmap='GnBu',vmin=tmin,vmax=tmax,shading='auto')
        #elif mo == 'HAP3':
        #    #df = df.where(df.z<300.,drop=True)
        #    im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.hus.values.T*1E3,cmap='GnBu',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HAP1':
            #df = df.where(df.z<300.,drop=True)
            im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.hus.values.T*1E3,cmap='GnBu',vmin=tmin,vmax=tmax,shading='auto')

    for i in range(pltmax):
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].set_ylabel(r'$z$ [m]')

    ax[0].set_ylim(0.,300.)

    date_form = md.DateFormatter("%H:%M")
    ax[3].xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax[3].set_xlim(datet,datet+datetime.timedelta(days=1))
    ax[4].set_xlabel('Time')
    plt.tight_layout()
    fig.subplots_adjust(top= 0.925,right=0.8,hspace=0.2)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax,label=r'$q_t$ [g kg$^{-1}$]')
    ax[0].set_title("Tower observations")
    ax[1].set_title("DALES")
    ax[2].set_title("HARMONIE-AROME HAP2")
    ax[3].set_title("HARMONIE-AROME HA43")    
   # ax[4].set_title("HARMONIE-AROME HAP3")
    ax[4].set_title("HARMONIE-AROME HAP1")

    print(f"Saving qPROF plot for {date} to {ql_dir}")

    plt.savefig("{0}QPROF_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))

def plot_WSPROF(date,df_arr,df_str):
    import matplotlib.dates as md
    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/WSPROF/"
    pltmax = 5
    obs_there = False
    fig, ax = plt.subplots(pltmax,1,sharey=True,sharex=True,figsize=(7,15))
    tmin = 0
    tmax = 15
    fig.suptitle(date)
    for df,mo in zip(df_arr,df_str) :
        if mo == 'OBS':
            tmin = np.nanmin(df.F.values)-1.0
            tmax = np.nanmax(df.F.values)+1.0
            obs_there = True
            im = ax[0].pcolormesh(df.time,df.z,df.F.values.T,cmap='Spectral_r',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'DALES':
            #tmin = np.nanmin(df.ws.values[:,0:20])-1.0
            #tmax = np.nanmax(df.ws.values[:,0:20])+1.0
            im = ax[1].pcolormesh(df.time.values,df.zt.values,df.ws.values.T,cmap='Spectral_r',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'HAP2':
            #df = df.where(df.z<300.,drop=True)
            im = ax[2].pcolormesh(df.time.values,df.z.values[0,:],df.ws.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HA43':
            #df = df.where(df.z<300.,drop=True)
            im = ax[3].pcolormesh(df.time.values,df.z.values[0,:],df.ws.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')
        #elif mo == 'HAP3':
        #    #df = df.where(df.z<300.,drop=True)
        #    im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.ws.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HAP1':
            #df = df.where(df.z<300.,drop=True)
            im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.ws.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')

    for i in range(pltmax):
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].set_ylabel(r'$z$ [m]')

    ax[0].set_ylim(0.,300.)

    date_form = md.DateFormatter("%H:%M")
    ax[3].xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax[3].set_xlim(datet,datet+datetime.timedelta(days=1))
    ax[4].set_xlabel('Time')
    plt.tight_layout()
    fig.subplots_adjust(top= 0.925,right=0.8,hspace=0.2)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax,label=r'$U$ [m s$^{-1}$]')
    ax[0].set_title("Tower observations")
    ax[1].set_title("DALES")
    ax[2].set_title("HARMONIE-AROME HAP2")
    ax[3].set_title("HARMONIE-AROME HA43")
   # ax[4].set_title("HARMONIE-AROME HAP3")
    ax[4].set_title("HARMONIE-AROME HAP1")

    print(f"Saving WS-PROF plot for {date} to {ql_dir}")

    plt.savefig("{0}WSPROF_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))

def plot_WSBLPROF(date,df_arr,df_str):
    import matplotlib.dates as md
    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/WSBLPROF/"
    pltmax = 5
    obs_there = False
    fig, ax = plt.subplots(pltmax,1,sharey=True,sharex=True,figsize=(7,15))
    tmin = 0
    tmax = 30
    fig.suptitle(date)
    for df,mo in zip(df_arr,df_str) :
        if mo == 'OBS_lid':
            obs_there = True
            im = ax[0].pcolormesh(df.time.values,df.z.values,df.WS.values.T,cmap='Spectral_r',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'DALES':
            #tmin = np.nanmin(df.ws.values[:,0:60])-1.0
            #tmax = np.nanmax(df.ws.values[:,0:60])+1.0
            im = ax[1].pcolormesh(df.time.values,df.zt.values,df.ws.values.T,cmap='Spectral_r',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'HAP2':
            #df = df.where(df.z<300.,drop=True)
            im = ax[2].pcolormesh(df.time.values,df.z.values[0,:],df.ws.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HA43':
            #df = df.where(df.z<300.,drop=True)
            im = ax[3].pcolormesh(df.time.values,df.z.values[0,:],df.ws.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')
        #elif mo == 'HAP3':
        #    #df = df.where(df.z<300.,drop=True)
        #    im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.ws.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HAP1':
            #df = df.where(df.z<300.,drop=True)
            im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.ws.values.T,cmap='Spectral_r',vmin=tmin,vmax=tmax,shading='auto')

    for i in range(pltmax):
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].set_ylabel(r'$z$ [m]')

    ax[0].set_ylim(0.,3000.)

    date_form = md.DateFormatter("%H:%M")
    ax[3].xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax[3].set_xlim(datet,datet+datetime.timedelta(days=1))
    ax[4].set_xlabel('Time')
    plt.tight_layout()
    fig.subplots_adjust(top= 0.925,right=0.8,hspace=0.2)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax,label=r'$U$ [m s$^{-1}$]')
    ax[0].set_title("Doppler Lidar observations")
    ax[1].set_title("DALES")
    ax[2].set_title("HARMONIE-AROME HAP2")
    ax[3].set_title("HARMONIE-AROME HA43")
    ax[4].set_title("HARMONIE-AROME HAP1")

    print(f"Saving WS_BL-PROF plot for {date} to {ql_dir}")

    plt.savefig("{0}WSBLPROF_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))

def plot_WDPROF(date,df_arr,df_str):
    import matplotlib.dates as md
    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/WDPROF/"
    pltmax = 5
    fig, ax = plt.subplots(pltmax,1,sharey=True,sharex=True,figsize=(7,15))
    tmin = 0
    tmax = 360
    fig.suptitle(date)
    for df,mo in zip(df_arr,df_str) :
        if mo == 'OBS':
           # tmin = np.nanmin(df.Q.values)-0.5
           # tmax = np.nanmax(df.Q.values)+0.5
            im = ax[0].pcolormesh(df.time,df.z,df.D.values.T,cmap='twilight',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'DALES':
           # dfn = df.where(df.zt<300.)
            im = ax[1].pcolormesh(df.time.values,df.zt.values,df.wd.values.T,cmap='twilight',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'HAP2':
            #df = df.where(df.z<300.,drop=True)
            im = ax[2].pcolormesh(df.time.values,df.z.values[0,:],df.wd.values.T,cmap='twilight',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HA43':
            #df = df.where(df.z<300.,drop=True)
            im = ax[3].pcolormesh(df.time.values,df.z.values[0,:],df.wd.values.T,cmap='twilight',vmin=tmin,vmax=tmax,shading='auto')
#        elif mo == 'HAP3':
#            #df = df.where(df.z<300.,drop=True)
#            im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.wd.values.T,cmap='twilight',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HAP1':
            #df = df.where(df.z<300.,drop=True)
            im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.wd.values.T,cmap='twilight',vmin=tmin,vmax=tmax,shading='auto')

    for i in range(pltmax):
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].set_ylabel(r'$z$ [m]')

    ax[0].set_ylim(0.,300.)

    date_form = md.DateFormatter("%H:%M")
    ax[3].xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax[3].set_xlim(datet,datet+datetime.timedelta(days=1))
    ax[4].set_xlabel('Time')
    plt.tight_layout()
    fig.subplots_adjust(top= 0.925,right=0.8,hspace=0.2)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax,label=r'Wind direction [$^{\circ}$]')
    ax[0].set_title("Tower observations")
    ax[1].set_title("DALES")
    ax[2].set_title("HARMONIE-AROME HAP2")
    ax[3].set_title("HARMONIE-AROME HA43")
  #  ax[4].set_title("HARMONIE-AROME HAP3")
    ax[4].set_title("HARMONIE-AROME HAP1")

    print(f"Saving WD-PROF plot for {date} to {ql_dir}")

    plt.savefig("{0}WDPROF_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))
def plot_WDBLPROF(date,df_arr,df_str):
    import matplotlib.dates as md
    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/WDBLPROF/"
    pltmax = 5
    fig, ax = plt.subplots(pltmax,1,sharey=True,sharex=True,figsize=(7,15))
    tmin = 0
    tmax = 360
    fig.suptitle(date)
    for df,mo in zip(df_arr,df_str) :
        if mo == 'OBS_lid':
            im = ax[0].pcolormesh(df.time.values,df.z.values,df.WD.values.T,cmap='twilight',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'DALES':
           # dfn = df.where(df.zt<300.)
            im = ax[1].pcolormesh(df.time.values,df.zt.values,df.wd.values.T,cmap='twilight',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'HAP2':
            #df = df.where(df.z<300.,drop=True)
            im = ax[2].pcolormesh(df.time.values,df.z.values[0,:],df.wd.values.T,cmap='twilight',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HA43':
            #df = df.where(df.z<300.,drop=True)
            im = ax[3].pcolormesh(df.time.values,df.z.values[0,:],df.wd.values.T,cmap='twilight',vmin=tmin,vmax=tmax,shading='auto')
#        elif mo == 'HAP3':
#            #df = df.where(df.z<300.,drop=True)
#            im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.wd.values.T,cmap='twilight',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HAP1':
            #df = df.where(df.z<300.,drop=True)
            im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.wd.values.T,cmap='twilight',vmin=tmin,vmax=tmax,shading='auto')

    for i in range(pltmax):
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].set_ylabel(r'$z$ [m]')

    ax[0].set_ylim(0.,3000.)

    date_form = md.DateFormatter("%H:%M")
    ax[3].xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax[3].set_xlim(datet,datet+datetime.timedelta(days=1))
    ax[4].set_xlabel('Time')
    plt.tight_layout()
    fig.subplots_adjust(top= 0.925,right=0.8,hspace=0.2)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax,label=r'Wind direction [$^{\circ}$]')
    ax[0].set_title("Doppler lidar observations")
    ax[1].set_title("DALES")
    ax[2].set_title("HARMONIE-AROME HAP2")
    ax[3].set_title("HARMONIE-AROME HA43")
  #  ax[4].set_title("HARMONIE-AROME HAP3")
    ax[4].set_title("HARMONIE-AROME HAP1")

    print(f"Saving WD_BL-PROF plot for {date} to {ql_dir}")
    plt.savefig("{0}WDBLPROF_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))


def plot_QLPROF(date,df_arr,df_str):
    import matplotlib.dates as md
    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/QLPROF/"
    pltmax = 5
    fig, ax = plt.subplots(pltmax,1,sharey=True,sharex=True,figsize=(7,15))
    tmin = 0
    tmax = 0.5
    im = None
    im2 = None
    fig.suptitle(date)
    for df,mo in zip(df_arr,df_str) :
        if mo == 'OBS_ceil':
            print(df,mo)
            bmin = 10#np.nanmin(np.log(df.beta_raw.values))
            bmax = 20#np.nanmax(np.log(df.beta_raw.values))
            im2 = ax[0].pcolormesh(df.time,df.range,np.log(df.beta_raw.values.T),cmap='cubehelix_r',vmin=bmin, vmax=bmax,shading='auto')
            cbar_ax = fig.add_axes([0.85, 0.78, 0.02, 0.15])
            fig.colorbar(im2, cax=cbar_ax,label=r'log($\beta$) [-]')
        elif mo == 'DALES':
            tmin = np.nanmin(df.ql.values)*1E3
            tmax = np.nanmax(df.ql.values)*1E3
            im = ax[1].pcolormesh(df.time.values,df.zt.values,df.ql.values.T*1E3,cmap='cubehelix_r',vmin=tmin, vmax=tmax,shading='auto')
        elif mo == 'HAP2':
            #df = df.where(df.z<300.,drop=True)
            im = ax[2].pcolormesh(df.time.values,df.z.values[0,:],(df.clw.values.T+df.cli.values.T)*1E3,cmap='cubehelix_r',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HA43':
            #df = df.where(df.z<300.,drop=True)
            im = ax[3].pcolormesh(df.time.values,df.z.values[0,:],(df.clw.values.T+df.cli.values.T)*1E3,cmap='cubehelix_r',vmin=tmin,vmax=tmax,shading='auto')
  #      elif mo == 'HAP3':
  #          #df = df.where(df.z<300.,drop=True)
  #          im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],(df.clw.values.T+df.cli.values.T)*1E3,cmap='cubehelix_r',vmin=tmin,vmax=tmax,shading='auto')
        elif mo == 'HAP1':
            #df = df.where(df.z<300.,drop=True)
            im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],(df.clw.values.T+df.cli.values.T)*1E3,cmap='cubehelix_r',vmin=tmin,vmax=tmax,shading='auto')


    for i in range(pltmax):
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].set_ylabel(r'$z$ [m]')

    ax[0].set_ylim(0.,9000.)

    date_form = md.DateFormatter("%H:%M")
    ax[3].xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax[3].set_xlim(datet,datet+datetime.timedelta(days=1))
    ax[4].set_xlabel('Time')
    plt.tight_layout()
    fig.subplots_adjust(top= 0.925,right=0.8,hspace=0.2)
    cbar_ax = fig.add_axes([0.85, 0.10, 0.02, 0.55])
    fig.colorbar(im, cax=cbar_ax,label=r'$q_l$ [g kg$^{-1}$]')
    ax[0].set_title("Ceilometer observations")
    ax[1].set_title("DALES")
    ax[2].set_title("HARMONIE-AROME HAP2")
    ax[3].set_title("HARMONIE-AROME HA43")
#    ax[4].set_title("HARMONIE-AROME HAP3")
    ax[4].set_title("HARMONIE-AROME HAP1")
    print(f"Saving qlPROF plot for {date} to {ql_dir}")

    plt.savefig("{0}QLPROF_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))

def plot_TURBPROF(date,df_arr,df_str):
    import matplotlib.dates as md
    cmap = make_colormap({0.:'#ffffff',0.2:'cornflowerblue',0.4:'limegreen',0.6:'yellow',0.8:'darkorange',1.:'crimson'})
    ql_dir = "/usr/people/theeuwes/public_html/quicklooks/recent_plots/plots/TURBPROF/"
    pltmax = 5
    obs_there = False
    fig, ax = plt.subplots(pltmax,1,sharey=True,sharex=True,figsize=(7,15))
    tmin = 0
    tmax = 2
    fig.suptitle(date)
    for df,mo in zip(df_arr,df_str) :
        if mo == 'OBS_lidT':
            obs_there = True
            im = ax[0].pcolormesh(df.time.values,df.z.values,df.WVAR.values.T,cmap=cmap,vmin=tmin, vmax=tmax,shading='auto')
            ax[0].scatter(df.time.values,df.MH.values,c='k',s=5)
        elif mo == 'DALES':
            #tmin = np.nanmin(df.ws.values[:,0:60])-1.0
            #tmax = np.nanmax(df.ws.values[:,0:60])+1.0
            im = ax[1].pcolormesh(df.time.values,df.zt.values,df.w2s.values.T+df.w2r.values.T,cmap=cmap,vmin=tmin, vmax=tmax,shading='auto')
            ax[1].scatter(df.time.values,df.MH.values,c='k',s=5)
        elif mo == 'HAP2':
            #df = df.where(df.z<300.,drop=True)
            im = ax[2].pcolormesh(df.time.values,df.z.values[0,:],df.sigma2w.values.T,cmap=cmap,vmin=tmin,vmax=tmax,shading='auto')
            ax[2].scatter(df.time.values,df.MH.values,c='k',s=5)
            ax[2].scatter(df.time.values,df.pblh_surf.values,c='crimson',s=5)
        elif mo == 'HA43':
            #df = df.where(df.z<300.,drop=True)
            im = ax[3].pcolormesh(df.time.values,df.z.values[0,:],df.sigma2w.values.T,cmap=cmap,vmin=tmin,vmax=tmax,shading='auto')
            ax[3].scatter(df.time.values,df.MH.values,c='k',s=5)
            ax[3].scatter(df.time.values,df.pblh_surf.values,c='crimson',s=5)
        elif mo == 'HAP1':
            #df = df.where(df.z<300.,drop=True)
            im = ax[4].pcolormesh(df.time.values,df.z.values[0,:],df.sigma2w.values.T,cmap=cmap,vmin=tmin,vmax=tmax,shading='auto')
            ax[4].scatter(df.time.values,df.MH.values,c='k',s=5)
            ax[4].scatter(df.time.values,df.pblh_surf.values,c='crimson',s=5)
    for i in range(pltmax):
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
        ax[i].set_ylabel(r'$z$ [m]')

    ax[0].set_ylim(0.,3000.)

    date_form = md.DateFormatter("%H:%M")
    ax[3].xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax[3].set_xlim(datet,datet+datetime.timedelta(days=1))
    ax[4].set_xlabel('Time')
    plt.tight_layout()
    fig.subplots_adjust(top= 0.925,right=0.8,hspace=0.2)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax,label=r'$\sigma^2_W$ / $TKE$ [m$^2$ s$^{-2}$]')
    ax[0].set_title("Doppler Lidar observations")
    ax[1].set_title("DALES")
    ax[2].set_title("HARMONIE-AROME HAP2")
    ax[3].set_title("HARMONIE-AROME HA43")
    ax[4].set_title("HARMONIE-AROME HAP1")
    print(f"Saving TURB_BL-PROF plot for {date} to {ql_dir}")
    plt.savefig("{0}TURBPROF_recent_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))

def plot_WSBLPROF_test(date,df_arr,df_str):
    import matplotlib.dates as md
    ql_dir = ""
    pltmax = 1
    obs_there = False
    fig, ax = plt.subplots(pltmax,1,sharey=True,sharex=True,figsize=(7,3))
    tmin = 0
    tmax = 35
    #fig.suptitle(date)
    for df,mo in zip(df_arr,df_str) :
        if mo == 'HAP2':
            time = df.time +  np.timedelta64(1, 'h')# datetime.timedelta(seconds=3600)
            im = ax.pcolormesh(time,df.z.values[0,:],df.ws.values.T,cmap='jet',vmin=tmin,vmax=tmax,shading='auto')

    for i in range(pltmax):
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_ylabel(r'Hoogte [m]')

    ax.set_ylim(0.,2500.)

    date_form = md.DateFormatter("%H:%M")
    ax.xaxis.set_major_formatter(date_form)
    datet = datetime.datetime(date.year,date.month,date.day,0,0)
    ax.set_xlim(datet,datet+datetime.timedelta(days=1))
    ax.set_xlabel('Lokale tijd')
    plt.tight_layout()
    fig.subplots_adjust(top= 0.925,right=0.8,hspace=0.2)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax,label=r'Wind snelheid [m/s]')
    ax.set_title("HARMONIE-AROME")

    print(f"Saving WS_BL-PROF plot for {date} to {ql_dir}")

    plt.savefig("{0}WSBLPROF_HAP2_{1:04d}{2:02d}{3:02d}".format(ql_dir,date.year,date.month,date.day))