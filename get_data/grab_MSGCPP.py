#!/usr/bin/env python3
import os
import datetime
import wget
import numpy as np
import xarray as xr

def grab_MSG_CPP(date):
    # create times of the day array

    ws_path= "/nobackup/users/theeuwes/tmp_data/"
 
    dt = np.array([datetime.datetime.combine(date,datetime.time()) + datetime.timedelta(minutes=15*x) for x in range(0,96)])
    CWP = []
    CWP_min = []
    CWP_max = []
    SWD = []
    SWD_min = []
    SWD_max = []
    times = []

    yr = str(date.year).zfill(4)
    mo = str(date.month).zfill(2)
    da = str(date.day).zfill(2)
    for t in dt:
        hr = str(t.hour).zfill(2)
        mi = str(t.minute).zfill(2)
        url1 = f"https://msgcpp-adaguc.knmi.nl/adaguc-server?dataset=msgrt&service=wcs&request=getcoverage&coverage=surface_downwelling_shortwave_flux_in_air&FORMAT=NetCDF4&CRS=EPSG%3A4326&BBOX=4.8509,51.8944,4.9891,52.0656&RESX=0.046&RESY=0.057&time={yr}-{mo}-{da}T{hr}:{mi}:00Z"
        url2 = f"https://msgcpp-adaguc.knmi.nl/adaguc-server?dataset=msgrt&service=wcs&request=getcoverage&coverage=atmosphere_cloud_condensed_water_content&FORMAT=NetCDF4&CRS=EPSG%3A4326&BBOX=4.8509,51.8944,4.9891,52.0656&RESX=0.046&RESY=0.057&time={yr}-{mo}-{da}T{hr}:{mi}:00Z"
        
        file1 = ws_path+f"MSG_CPP_SWD_{yr}-{mo}-{da}-T{hr}{mi}00Z.nc"
        file2 = ws_path+f"MSG_CPP_CWP_{yr}-{mo}-{da}-T{hr}{mi}00Z.nc"
        wget.download(url1,file1)
        wget.download(url2,file2)
        
        # check file size 
        s1 = os.path.getsize(file1)
        s2 = os.path.getsize(file2)

        if (s1<1000): 
            SWD.append(np.nan)
            SWD_min.append(np.nan)
            SWD_max.append(np.nan)
        else:
            df = xr.open_dataset(file1)
            swd_tmp = np.nanmean(df.sds.values)
            swdm_tmp = np.max(df.sds.values)
            swdi_tmp = np.min(df.sds.values)
            swd_tmp = np.nan if (~np.isfinite(swd_tmp)) else swd_tmp
            swdm_tmp = np.nan if (~np.isfinite(swdm_tmp)) else swdm_tmp
            swdi_tmp = np.nan if (~np.isfinite(swdi_tmp)) else swdi_tmp
            SWD.append(swd_tmp)
            SWD_max.append(swdm_tmp)
            SWD_min.append(swdi_tmp)

        if (s2<1000):
            CWP.append(np.nan)
            CWP_min.append(np.nan)
            CWP_max.append(np.nan)
        else:
            df = xr.open_dataset(file2)
            cwp_tmp = np.nanmean(df.cwp.values)
            cwpm_tmp = np.max(df.cwp.values)
            cwpi_tmp = np.min(df.cwp.values)
            cwp_tmp = np.nan if (~np.isfinite(cwp_tmp)) else cwp_tmp
            cwpm_tmp = np.nan if (~np.isfinite(cwpm_tmp)) else cwpm_tmp
            cwpi_tmp = np.nan if (~np.isfinite(cwpi_tmp)) else cwpi_tmp
            CWP.append(cwp_tmp)
            CWP_max.append(cwpm_tmp)
            CWP_min.append(cwpi_tmp)
            times.append(df.time.values)

    times = np.reshape(np.array(times),-1)   
    ds = xr.Dataset(data_vars= dict(sds=(["time"],np.array(SWD),dict(long_name=("Surface downwelling solar radiation"),
                   units=("W m-2"))), sds_max=(["time"],np.array(SWD_max),dict(long_name=("Maximum Surface downwelling solar radiation"),
                   units=("W m-2"))),sds_min=(["time"],np.array(SWD_min),dict(long_name=("Minimum Surface downwelling solar radiation"),
                   units=("W m-2"))),cwp=(["time"],np.array(CWP),dict(long_name=("Cloud Water Path"),
                   units=("g m-2"))),cwp_max=(["time"],np.array(CWP_max),dict(long_name=("Maximum Cloud Water Path"),
                   units=("g m-2"))),cwp_min=(["time"],np.array(CWP_min),dict(long_name=("Minimum Cloud Water Path"),
                   units=("g m-2")))), coords = dict(time=dt), attrs = df.attrs)
    #units = 'days since 1900-01-01'
    #ds.time.encoding['units'] = units
    #ds.time.attrs.update(units = units)
    #ds.time.encode['units'] = 'days since 1900-01-01'
    ds.to_netcdf(ws_path+f"MSGCPP_cabauw{yr}-{mo}-{da}.nc",encoding={'time':{'units':'days since 1900-01-01'}})
    
    # clean up
    os.system(f"rm {ws_path}MSG_CPP*")


if __name__ == '__main__':
   # argv = sys.argv[1]  # date yyyymmdd
   # date = datetime.datetime.strptime(argv,'%Y%m%d')
    date = datetime.date.today() - datetime.timedelta(days=2)
    print(date)

    grab_MSG_CPP(date)
