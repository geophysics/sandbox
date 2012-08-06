# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Jared Peacock)s
"""

import os
import fnmatch
import Z
import LatLongUTMconversion as utm2ll
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np

#dirpath=r"C:\Peacock\PHD\Geothermal\Paralana\jared\edi17Aug\dr"
#dirpath=r"C:\Peacock\My Dropbox\Paralana\EDIFilesBaseSurvey"
#dirpath=r"C:\Peacock\PHD\Geothermal\Paralana\EDIFilesBaseSurvey\CFA"
#dirpath=r'C:\Peacock\PHD\Geothermal\Paralana\ForwardModels\ParalanaBase\Base'
#dirpath=r'C:\Peacock\PHD\Geothermal\Paralana\Quantec'
#dirpath=r'C:\Peacock\My Dropbox\Dongara_GRE\EDIFiles\CF'
#dirpath=r"d:\Data\EDIfiles"
#dirpath=r"C:\Peacock\My Dropbox\Paralana\InjectionEDIfiles\CFA\SS\DR"
#dirpath=r"C:\Peacock\PHD\Geothermal\Paralana\EDIFilesBaseSurvey\InjectionPlot"
#dirpath=r"c:\Users\Own er\Documents\PHD\Geothermal\Paralana\EDIFilesBaseSurvey\CFA"
#dirpath=r"c:\Users\Own er\Documents\PHD\Geothermal\Paralana\ForwardModels\ParalanaBase\Base"
dirpath=r"c:\Peacock\PHD\OCCAM\Antarctica\Sline\DR"

spa='\t'
llfmt='%2.6f'
nefmt='%.2f'

#edifid=open(r"C:\Peacock\PHD\Geothermal\QuantecLine.txt",'w')
#edifid=open(r"C:\Peacock\My Dropbox\Dongara_GRE\StationLocations2.txt",'w')
#edifid=open(r"C:\Peacock\PHD\Geothermal\StationLocationsInj.txt",'w')
##edifid=open(r"c:\Users\Own er\StationLocations.txt",'w')
#edifid.write('Station'+spa+'Latitude'+spa+'Longitude'+spa+'Easting'+spa+'Northing'
#            +spa+'Zone'+'\n')

latlst=[]
lonlst=[]
stationlst=[]
eastlst=[]
northlst=[]
for filename in os.listdir(dirpath):
    if fnmatch.fnmatch(filename,'*.edi'):
        imp=Z.Z(os.path.join(dirpath,filename))
        lat=imp.lat
        lon=imp.lon
        station=imp.station[0:6]
        zone,east,north=utm2ll.LLtoUTM(23,lat,lon)
        latlst.append(imp.lat)
        lonlst.append(imp.lon)
        stationlst.append(imp.station)
        zone,east,north=utm2ll.LLtoUTM(23,imp.lat,imp.lon)
        eastlst.append(east)
        northlst.append(north)
#        edifid.write(station+spa+llfmt % lat+spa+llfmt % lon+spa+nefmt % east+spa+
#                    nefmt % north+spa+str(zone)+'\n')
#        edifid.write(station+spa+llfmt % lat+spa+llfmt % lon+'\n')
#edifid.close()

plt.rcParams['font.size']=14

#plot stations
fig=plt.figure(1)
plt.plot(lonlst,latlst,'vk',ms=16)
for ii,stat in enumerate(stationlst):
    plt.text(lonlst[ii],latlst[ii]+.009,stat[:-2],horizontalalignment='center',
             verticalalignment='bottom',fontdict={'size':18,'weight':'bold'})
plt.xlabel('Longitude',fontsize=16,fontweight='bold')
plt.ylabel('Latititude',fontsize=16,fontweight='bold')
plt.grid('on')

#plot a bestfitting line
#p=sp.polyfit(lonlst,latlst,1)
#plt.plot(lonlst,sp.polyval(p,lonlst),'-b')
#theta=np.arctan(p[0])
#
#for ii in range(len(stationlst)):
#    d=(latlst[ii]-sp.polyval(p,lonlst[ii]))*np.cos(theta)
#    x0=lonlst[ii]+d*np.sin(theta)
#    y0=latlst[ii]-d*np.cos(theta)
#    plt.plot(x0,y0,'x',color='r',ms=5,mew=3)
#    plt.text(x0,y0,stationlst[ii][0:7],horizontalalignment='center',
#             verticalalignment='baseline')