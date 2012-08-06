# -*- coding: utf-8 -*-
"""
Created on Fri May 25 12:01:15 2012

@author: a1185872
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sps
#import MTPlotTools as mtplot
#import matplotlib.gridspec as gridspec
#import matplotlib.pylab as pylab
#import MTtools as mt
from matplotlib.ticker import FormatStrFormatter,MultipleLocator,FixedLocator
import Z
#from matplotlib.colors import LinearSegmentedColormap,Normalize
#from matplotlib.colorbar import *
#from matplotlib.patches import Ellipse,Rectangle
import pickle
#
#ptcmapdict={'red':((0.0,1.0,1.0),(1.0,1.0,1.0)),
#            'green':((0.0,0.0,1.0),(1.0,0.0,1.0)),
#            'blue':((0.0,0.0,0.0),(1.0,0.0,0.0))}
#ptcmap=LinearSegmentedColormap('ptcmap',ptcmapdict,256)

#drillhole lat and long
dhll=(139.72851,-30.2128)
#rotation of phase tensors so y points north
rot=180
#median filter size
mks=(3,3)
esized=3

ctype='data'
ttype='rt'


llst=[':k',':k',':k',':k']*2
mlst=['*','o','s','v']*2
clst=['blue','maroon','purple','green']*2

#===============================================================================
# Initialize parameters and locate directory path
#===============================================================================
if ctype=='data':
#    edipath=r"c:\Users\Own er\Documents\PHD\Geothermal\Paralana\EDIFilesInjection\EDIfiles\AdvPro24Hrs"
#    edipath=r"c:\Peacock\PHD\Geothermal\Paralana\EDIFilesInjection\EDIfiles\AdvPro24Hrs"
    edipath=r"c:\Peacock\PHD\Geothermal\Paralana\EDIFilesInjection\EDIfiles\AdvPro4Hr"
    stationlst=['pb01','pb04','pb24','pb27']
    
#    daylst=['','192','193','194','195','196','197']
    daylst=['191','192','193','194','195','196','197']
    
#    hrlst=['0'+str(hh) for hh in range(10)]+[str(hh) for hh in range(10,25)]
    hrlst=['0'+str(hh) for hh in range(0,10,4)]+[str(hh) for hh in range(12,25,4)]
    
    dmlst=[day+hr for day in daylst for hr in hrlst][12:]
    nf=36
#    nd=len(daylst)
    nd=len(dmlst)
    ns=(len(stationlst))
    nrow=3
    ncol=3
#    prange=[8,12,15,17,18,21,22,25,28]
    prange=[12,18,20,22,24,26,28,30,31]
    fignum=3

pkfn=r"c:\Peacock\PHD\ParalanaTimeLapseZErr"+ctype+ttype+'.pkl'

#==============================================================================
# get estimates of error in Z 
#==============================================================================

if os.path.isfile(pkfn)==False:
    zerr=np.zeros((nf,nd,ns,2,2))
    
    for ss,station in enumerate(stationlst):
        for ii,dm in enumerate(dmlst):
            try:
                z1=Z.Z(os.path.join(edipath,station+dm+'.edi'))
                zerr[:len(z1.period),ii,ss,:,:]=z1.zvar
            except IOError:
                pass
            
#        edilst=[os.path.join(edipath,edi) for edi in os.listdir(edipath) 
#                if edi.find(station)==0]
##        edilst=[os.path.join(edipath,station+'.edi')]+\
##                [os.path.join(edipath,station+day+'00.edi') for day in daylst[1:]]
#        for ii,edi in enumerate(edilst):
#            z1=Z.Z(edi)
#            zerr[:len(z1.period),ii,ss,:,:]=z1.zvar
    period=z1.period
    
    pkfid=file(pkfn,'w')
    pickle.dump((zerr,period),pkfid)
    pkfid.close()
    
else:
    pkfid=file(pkfn,'r')
    zerr,period=pickle.load(pkfid)
    pkfid.close()
    
#==============================================================================
# Plot Errors 
#==============================================================================

fig1=plt.figure(fignum,dpi=200)
llst=[]

zerr=np.nan_to_num(zerr)
zerr[np.where(zerr>100)]=100
for ii,ff in enumerate(prange):
    if ii==0:
        ax=fig1.add_subplot(nrow,ncol,ii+1)
    else:
        ax=fig1.add_subplot(nrow,ncol,ii+1,sharex=ax)

    for ss in range(ns):
        l1=ax.plot(np.arange(nd),zerr[ff,:,ss,1,1])
        if ii==0:        
            llst.append(l1)
    ax.plot(np.arange(nd),np.median(zerr[ff,:,:,1,1],axis=1),'orange')
    
    ax.set_ylim(0,np.median(zerr)*10)
    ax.text(.25,np.median(zerr)*10-.05,'T={0:.2g}(s)'.format(period[ff]),
            fontdict={'size':12,'weight':'bold'},
            bbox={'facecolor':'white'},verticalalignment='top',
            horizontalalignment='left')
    ax.xaxis.set_major_locator(FixedLocator([12,24,48,72,96,120,144,168]))
    ax.xaxis.set_ticklabels(daylst)
    ax.xaxis.set_minor_locator(MultipleLocator(1))

        


        






