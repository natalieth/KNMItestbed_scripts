#!/usr/bin/env python3

import os
import sys
import glob
import shutil
import datetime

def grab_DALES(date,exp):
    ws_path= "/nobackup/users/theeuwes/testbed/DALES/{0:04d}/{1:02d}/{2:02d}/".format(date.year,date.month,date.day)
    bull_path = "theeuwes@bxshnr02:/nfs/home/users/theeuwes/work/testbed/DALES/{0:04d}/{1:02d}/{2:02d}/".format(date.year,date.month,date.day)

    # make directories

    if not os.path.exists(ws_path):
        os.makedirs(ws_path)

    FilesA_from = [bull_path+'profiles.{0:03d}.{1:04d}{2:02d}{3:02d}.nc'.format(exp,date.year,date.month,date.day),
                   bull_path+'tmser.{0:03d}.{1:04d}{2:02d}{3:02d}.nc'.format(exp,date.year,date.month,date.day),
                   bull_path+'namoptions.{0:03d}.{1:04d}{2:02d}{3:02d}'.format(exp,date.year,date.month,date.day)]
    FilesA_to   = [ws_path+'profiles.{0:03d}.{1:04d}{2:02d}{3:02d}.nc'.format(exp,date.year,date.month,date.day),
                   ws_path+'tmser.{0:03d}.{1:04d}{2:02d}{3:02d}.nc'.format(exp,date.year,date.month,date.day),
                   ws_path+'namoptions.{0:03d}.{1:04d}{2:02d}{3:02d}'.format(exp,date.year,date.month,date.day)]

    for i,j in zip(FilesA_from,FilesA_to):
        print("Copying :", i)
        os.system("rsync -vau {} {}".format(i,j))


if __name__ == '__main__':
   # argv = sys.argv[1]  # date yyyymmdd
   # date = datetime.datetime.strptime(argv,'%Y%m%d')
    date = datetime.date.today() - datetime.timedelta(days=4)
    print(date)

  #  exp = int(sys.argv[2])  # experiment number:: 1
    exp = int(sys.argv[1])  # experiment number:: 1


    grab_DALES(date,exp)
