# -*- coding: utf-8 -*-
"""
Created on Wed Mar 07 15:34:40 2012

@author: geophysics
"""

import numpy as np
import h5py
import os
import MTtools as mt

#dirpath=r"d:\InjectionJuly2011"
#slst=['PB01','PB04','PB13','PB15','PB24','PB25','PB27','PB49','PB54','PBRT1']
dirpath=r"g:\Paralana\BaseMarch2011"
slst=['pb12','pb34','pb35','pb45','pb51','pb54','pb04','pb05','pb06',
      'pb07','pb08','pb09']

day='095'
start='050000'
stop='220000'
crate='001000'
d=100
df=500.
complst=['EX','EY','BX','BY']

nc=len(complst)
ns=len(slst)*nc
na=int(df/d)*3600*(int(stop[0:2])-int(start[0:2]))


cstr='Comb'+start[0:2]+'to'+stop[0:2]+'d'+str(d)
sstr=start[0:2]+'to'+stop[0:2]+'d'+str(d)

#make database
db=h5py.File(os.path.join(dirpath,'TimeSeriesDB.hdf5'))

#make a subgroup for the day
sgdb=db.create_group(day)

#create data set in the subgroup
dset=sgdb.create_dataset('TS',(ns,na))

#load the timeseries into the database
for ii,station in enumerate(slst):
    print 'Loading '+station
    for jj,comp in enumerate(complst):
        try:
            dset[ii*nc+jj,:]=np.loadtxt(os.path.join(dirpath,station,day,cstr,
                                        station+sstr+'.'+comp))
            print '\t Loaded '+comp
        except IOError:
            cfile,clst=mt.combineFewFiles(os.path.join(dirpath,station,day),
                                          station,start,stop,crate,
                                          complst=[comp],d=d)
            dset[ii*nc+jj,:]=np.loadtxt(cfile[0])
        
        
db.close()
                                    