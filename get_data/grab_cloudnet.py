#!/usr/bin/env python3
import datetime

def get_cloudnet_cabauw(date):
    import requests
    
    datadir = '/nobackup/users/theeuwes/tmp_data/'
    url = 'https://cloudnet.fmi.fi/api/files'
    payload = {
        'date': '{0:04d}-{1:02d}-{2:02d}'.format(date.year,date.month,date.day),
        'site': 'cabauw',
        'product': 'mwr'
        }
    metadata = requests.get(url, payload).json()
    print(metadata)
    for row in metadata:
        res = requests.get(row['downloadUrl'])
        filename = datadir + row['filename']
        with open(filename, 'wb') as f:
            f.write(res.content)


if __name__ == '__main__':
    date = datetime.date.today() - datetime.timedelta(days=4)

    
    get_cloudnet_cabauw(date)
