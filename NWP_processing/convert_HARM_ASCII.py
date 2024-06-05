import os
import numpy as np 
import datetime 
import pandas as pd
import xarray as xr

def unzip_HARM(file,location):
   
    os.system("tar -xvf {} --directory {}".format(file,location))

def add_pres_height(df):
    ##### calculate pressure #####

    ahalf= (pd.read_csv('/usr/people/theeuwes/testbed/recent_runs/H43_65lev.txt',
                         header=None,index_col=[0],delim_whitespace=True))[1].values[:]
    bhalf= (pd.read_csv('/usr/people/theeuwes/testbed/recent_runs/H43_65lev.txt',
                         header=None,index_col=[0],delim_whitespace=True))[2].values[:]
    ph = np.array([ahalf + (p * bhalf) for p in df['psurf_surf'].values])
    dp = np.diff(ph,axis=1)
    print(dp)
    p = np.zeros((df.ta.values).shape)
    z = np.zeros((df.ta.values).shape)
    for k in range(0,len(df.level)):
        p[:,k] = 0.5 * (ph[:,k] + ph[:,k+1])
        p0p = (df['psurf_surf'].values/p[:,k])**(1/5.257)
        z[:,k] = ((p0p -1)*df.ta.values[:,k]) / 0.0065
    
    df['p'] =  xr.DataArray(data=p,dims = dict(time = df.time, level = df.level))
    df['z'] =  xr.DataArray(data=z,dims = dict(time = df.time, level = df.level))
    df['dp']=  xr.DataArray(data=dp,dims= dict(time = df.time, level = df.level))

    return df

def add_cwp(df):
    #### calculate the cloud water path #####
    g = 9.81 # gravitational constant
    print(df)
    dp = df.dp.values
    qc = df.clw.values + df.cli.values

    cwpc = 0
    for i in np.arange(dp.shape[1]):
        cwpc += (dp[:,i]*qc[:,i])
    cwp=cwpc/g
    print(cwp)
    
    df['cwp'] =  xr.DataArray(data=cwp,dims = dict(time = df.time))

    return df

def deaccumulate(df):
    dt = 3600
    h_deaccu = np.array(df.hfss_surf.values).copy()
    rain_deaccu = np.array(df.rain_surf.values).copy()
    snow_deaccu = np.array(df.snow_surf.values).copy()
    grau_deaccu = np.array(df.grau_surf.values).copy()

    for t in range(1,len(df.time.values)):
        h_deaccu[t] = (df.hfss_surf[t] - (df.hfss_surf[t-1]))/dt
        rain_deaccu[t] = (df.rain_surf[t] - df.rain_surf[t-1])
        snow_deaccu[t] = (df.snow_surf[t] - df.snow_surf[t-1])
        grau_deaccu[t] = (df.grau_surf[t] - df.grau_surf[t-1])
    h_deaccu[0] = df.hfss_surf[0]/dt
  
    df['hfssde_surf'] =  xr.DataArray(data=-1.*h_deaccu,dims = dict(time = df.time))
    df['rainde_surf'] =  xr.DataArray(data=rain_deaccu,dims = dict(time = df.time))
    df['snowde_surf'] =  xr.DataArray(data=snow_deaccu,dims = dict(time = df.time))
    df['graude_surf'] =  xr.DataArray(data=grau_deaccu,dims = dict(time = df.time))

    return df

def read_ASCII(date,suite,hh):
    zip_dir = "/nobackup/users/theeuwes/testbed/HARMONIE/"
    tmp_dir = "/nobackup/users/theeuwes/tmp_data/"
    LT_HARM = pd.read_csv('/usr/people/theeuwes/testbed/recent_runs/lookuptable_HARM.csv',index_col=0)

    if os.path.isfile(tmp_dir+"fc{0:04d}{1:02d}{2:02d}{3:02d}+001ascii".format(date.year,date.month,date.day,hh)):
        print("file already unzipped")
    else: 
        unzip_HARM(zip_dir+"Ruisdael_{0}_{1:04d}{2:02d}{3:02d}{4:02d}.tar.gz"
               .format(suite,date.year,date.month,date.day,hh), tmp_dir)
    df = xr.Dataset(None)
    ntime = 49
    nlev = 65
    ndims = [ntime,nlev]
    for t in range(1,ntime):  ##    t = 1  # loop over time at some point
        name = []
        long = []
        units = []
        level = []
 
        filename = tmp_dir+"fc{0:04d}{1:02d}{2:02d}{3:02d}+{4:03d}ascii".format(date.year,date.month,date.day,hh,t)
        try: 
            with open(filename, 'r') as input:
                lines = input.readlines()[1:] # skip first line
    
                for line in lines:
                
              #  if '1    51.971     4.927  134  130' in line: # when gets to cabauw data
              #      break
                
                    line = line.strip()
                    row = line.split()
                    try: 
                        variable = int(row[1])
                        lev = int(row[2])
                        vtype = int(row[3])
     
                        if vtype == 109: 
                            # variables at model levels 
                            tbl = LT_HARM[LT_HARM['var'] == variable]
                            name.append(tbl['name'].values[0])
                            long.append(tbl['long'].values[0])
                            units.append(tbl['units'].values[0])
                            level.append(lev)
                        elif vtype == 105: 
                            # variabples at specific height
                            tbl = LT_HARM[LT_HARM['var'] == variable]
                            if lev>1:
                                name.append(tbl['name'].values[0]+"_{0:d}m".format(lev))
                                long.append(tbl['long'].values[0]+" at {0:d} meter".format(lev))
                            else:
                                name.append(tbl['name'].values[0]+"_surf")
                                long.append(tbl['long'].values[0]+" at the surface")
                            units.append(tbl['units'].values[0])
                            level.append(-1)
    
                    except (IndexError,ValueError):
                        dummy = 0
                        #print("Reached end of file header")
            with open(filename, 'r') as f:
                last_line = f.readlines()[-1]
                values = last_line.split()
                values = [float(i) for i in values]
    
        
            name = np.array(name)
            long = np.array(long)
            units = np.array(units)
            values = np.array(values)
            level = np.array(level)
            
            time = np.array([ datetime.datetime(date.year,date.month,date.day,hh,0) + datetime.timedelta(seconds=3600*t)])
    
            variables = np.unique(name)
            for vrs in variables:
                ll  = level[np.where(name==vrs)]
                if any(ll==-1):
                    bla = values[np.where(name==vrs)]
                    ds =  xr.Dataset(data_vars=dict(
                        var =(["time"], bla)),
                        coords=dict(time=time))
                    ds = ds.rename({'var' : str(vrs)})
                else: 
                    bla = np.expand_dims(values[np.where(name==vrs)],axis=0)
                    ds =  xr.Dataset(data_vars=dict(
                        var =(["time","level"], bla)),
                        coords=dict(time=time,level=ll))
                    ds = ds.rename({'var' : str(vrs)})
    
                df = df.merge(ds,compat='no_conflicts')

        except FileNotFoundError:
            break    
    return df



def make_HARM_nc(date,hh,suite):
    # makes a netcdf file
    testbed ="/nobackup/users/theeuwes/testbed/HARMONIE/nc/"
    tmpdir = "/nobackup/users/theeuwes/tmp_data/"

    df = read_ASCII(date,suite,hh)
    try:
        print(df.time) # will cause a crash if dataset is empty
        df = add_pres_height(df)
        df = add_cwp(df)
        df = deaccumulate(df)
        print(df)
        df.to_netcdf(testbed+"CABAUW_{0}_{1:04d}{2:02d}{3:02d}{4:02d}.nc".format(suite,date.year,date.month,date.day,hh),
                     mode='w',format ="NETCDF4")
        os.system("rm {0}fc{1:04d}{2:02d}{3:02d}{4:02d}*ascii".format(tmpdir,date.year,date.month,date.day,hh))

    except AttributeError:
        print("Missing forecast {} HH: {} suite {}".format(date,hh,suite))


if __name__ == '__main__':
    date = datetime.date(2022,3,6)

    testbed ="/nobackup/users/theeuwes/testbed/HARMONIE/nc/"
    tmpdir = "/nobackup/users/theeuwes/tmp_data/"
    
    su = ["HAP1","HAP2","HAQ2","HAP3","HA43"]
    forecasts = [0,3,6,9,12,15,18,21]
    for suite in su:
        for hh in forecasts:
            df = read_ASCII(date,suite,hh)
            try: 
                print(df.time) # will cause a crash if dataset is empty
                df = add_pres_height(df)
                df = add_cwp(df)
                print(df.cwp)
                df = deaccumulate(df)
                print(df.hfss_surf)
                df.to_netcdf(testbed+"CABAUW_{0}_{1:04d}{2:02d}{3:02d}{4:02d}.nc".format(suite,date.year,date.month,date.day,hh),
                         mode='w',format ="NETCDF4")
                os.system("rm {0}fc{1:04d}{2:02d}{3:02d}{4:02d}*ascii".format(tmpdir,date.year,date.month,date.day,hh))

            except AttributeError:
                print("Missing forecast {} HH: {} suite {}".format(date,hh,suite)) 




