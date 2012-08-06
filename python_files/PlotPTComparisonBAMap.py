# -*- coding: utf-8 -*-
"""
Created on Wed May 25 16:39:30 2011

@author: a1185872
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sps
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LinearSegmentedColormap
import Z
from matplotlib.colorbar import *
from matplotlib.patches import Ellipse,Rectangle,Arrow
import pickle
import LatLongUTMconversion as utm2ll

ctype='data'
refe=23

ptcmapdict={'red':((0.0,1.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,0.0,1.0),(1.0,0.0,1.0)),
            'blue':((0.0,0.0,0.0),(1.0,0.0,0.0))}
ptcmap=LinearSegmentedColormap('ptcmap',ptcmapdict,256)

ptcmapdict3={'red':[(0.0,1.0,0.0),
                    (.49,1.0,1.0),
                    (0.5,0.0,1.0),
                    (1.0,1.0,1.0)],
            'green':[(0.0,0.0,1.0),
                     (.49,1.0,1.0),
                     (0.5,0.0,1.0),
                     (1.0,0.0,1.0)],
            'blue':[(0.0,0.0,1.0),
                    (.49,0.0,1.0),
                    (0.5,0.0,0.0),
                    (1.0,0.0,0.0)]}
#ptcmapdict3={'red':((0.0,0.0,0.0),(1.0,0.0,0.0)),
#            'green':((0.0,1.0,1.0),(1.0,0.0,0.0)),
#            'blue':((0.0,1.0,1.0),(1.0,1.0,1.0))}
ptcmap3=LinearSegmentedColormap('ptcmap3',ptcmapdict3,256)

bhll=(139.72851,-30.2128)
bhz,bhe,bhn=utm2ll.LLtoUTM(refe,bhll[1],bhll[0])
#===============================================================================
#Initialize parameters
#===============================================================================
if ctype=='data':    
    #set edipaths
    edipathi=r"C:\Peacock\My Dropbox\Paralana\InjectionEDIfiles\CFA"
    edipathb=r"C:\Peacock\My Dropbox\Paralana\EDIFilesBaseSurvey\CFA"
    
    pkfn=r'C:\Peacock\Python\PTBAComparison.pkl'
#    pkfn=r'c:\Users\Own er\Documents\PHD\Geothermal\Paralana\PTBAComparison.pkl'
#    
#    edipathi=r"c:\Users\Own er\Documents\PHD\Geothermal\Paralana\EDIFilesInjection\EDIfiles\CFA\DR"
#    edipathb=r"c:\Users\Own er\Documents\PHD\Geothermal\Paralana\EDIFilesBaseSurvey\CFA\SS\DR"

    
    #make list of existing edifiles
    edilst=[[os.path.join(edipathb,edib),os.path.join(edipathi,edii)] 
            for edib in os.listdir(edipathb) 
            for edii in os.listdir(edipathi)
            if edib.find('.')>0  if edib==edii]
    
    pstart=16
    pstop=36
    pstep=1
    #prange=[7,9,14,17,18,22,24,26,28,30,31,32,33,34,35]
    prange=[17,24,28,30,32,35]
    ncols=3
    mfs=(3,3)
    #set parameters for plotting
    esize=1.2
    emax=1
    cmax=.5
    a=1
    #number of frequencies
    nf=43
    ns=len(edilst)
    xlimits=(-3.2,3.2)
    ylimits=(-2.8,2.8)
    noise=None
    fignum=3

elif ctype=='fm':
    #set edipaths
    edipathb=r"C:\Peacock\PHD\Geothermal\Paralana\ForwardModels\ParalanaBase\base"
    edipathi=r"C:\Peacock\PHD\Geothermal\Paralana\ForwardModels\ParalanaBase\Step6"
    
#    edipathb=r"c:\Users\Own er\Documents\PHD\Geothermal\Paralana\ForwardModels\ParalanaBase\Base"
#    edipathi=r"c:\Users\Own er\Documents\PHD\Geothermal\Paralana\ForwardModels\ParalanaBase\Step6"
    
    pkfn=r'C:\Peacock\Python\PTBAComparisonFM.pkl'
#    pkfn=r'c:\Users\Own er\Documents\PHD\Geothermal\Paralana\PTBAComparisonFM.pkl'
    
    #make list of existing edifiles
    edilst=[[os.path.join(edipathb,edib),os.path.join(edipathi,edii)] 
            for edib in os.listdir(edipathb) 
            for edii in os.listdir(edipathi)
            if edib.find('.')>0  if edib[0:-4]==edii[0:-4]]
    
    pstart=4
    pstop=10
    pstep=1
    #prange=[7,9,14,17,18,22,24,26,28,30,31,32,33,34,35]
    prange=[5,6,7,8,9,10]
    ncols=3
    mfs=(1,1)
    #set parameters for plotting
    cmax=.0022
    #number of frequencies
    nf=14
    ns=len(edilst)
    esize=.6
    a=1000
#    xlimits=(139.19,139.75)
#    ylimits=(-30.23,-30.18)
    xlimits=(-3.4,4)
    ylimits=(-3.4,3.4)
#    xlimits=(-6,6)
#    ylimits=(-6,6)
    bhe=bhe
    bhn=bhn
    noise=None
    fignum=1
    

if not os.path.isfile(pkfn):
    azimutharr=np.zeros((nf,ns))
    phimaxarr=np.zeros((nf,ns))
    phiminarr=np.zeros((nf,ns))
    betaarr=np.zeros((nf,ns))
    colorarr=np.zeros((nf,ns))
    
    latlst=np.zeros(ns)
    lonlst=np.zeros(ns)

    stationlst=[]
    
    for ss,station in enumerate(edilst):
        #make a data type Z      
        imp1=Z.Z(station[0])
        imp2=Z.Z(station[1])
        
        stationlst.append(imp1.station)        
        
        sz,se,sn=utm2ll.LLtoUTM(refe,imp1.lat,imp1.lon)
        latlst[ss]=(sn-bhn)/1000.
        lonlst[ss]=(se-bhe)/1000.
        
        #get the phase tensor information
        pt1=imp1.getPhaseTensor(rotate=180)
        pt2=imp2.getPhaseTensor(rotate=180)
    
        
        #loop over period plotting the difference between phase tensors
        period=imp1.period
        n=len(period)
        for ii in range(nf):
            
            if noise!=None:
                sigman=np.sqrt(abs(pt1.phi[ii,0,1]*pt1.phi[ii,1,0]))*noise
                pt1.phi[ii]=pt1.phi[ii]+sigman*np.random.normal(size=(2,2))
                pt2.phi[ii]=pt2.phi[ii]+sigman*np.random.normal(size=(2,2))
            #calculate the difference between the two phase tensor ellipses
            phi=np.eye(2)-\
                    (np.dot(np.linalg.inv(pt2.phi[ii]),pt1.phi[ii])+
                    np.dot(pt1.phi[ii],np.linalg.inv(pt2.phi[ii])))/2       
            
            #compute the trace        
            tr=phi[0,0]+phi[1,1]
            #Calculate skew of phi and the cooresponding error
            skew=phi[0,1]-phi[1,0]
            #calculate the determinate and determinate error of phi
            phidet=abs(np.linalg.det(phi))
            
            #calculate reverse trace and error
            revtr=phi[0,0]-phi[1,1]
            
            #calculate reverse skew and error
            revskew=phi[1,0]+phi[0,1]
            
            beta=.5*np.arctan2(skew,tr)*(180/np.pi)
            alpha=.5*np.arctan2(revskew,revtr)*(180/np.pi)
            
            #need to figure out why azimuth is off by 90 deg
            azimuth=(alpha-beta)                   
           
            #calculate phimax
            phimax=np.sqrt(abs((.5*tr)**2+(.5*skew)**2))+\
                    np.sqrt(abs((.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2))
                
            #calculate minimum value for phi
            if phidet>=0:
                phimin=np.sqrt(abs((.5*tr)**2+(.5*skew)**2))-\
                np.sqrt(abs((.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2))
            elif phidet<0:
                phimin=-1*np.sqrt(abs((.5*tr)**2+(.5*skew)**2))-np.sqrt(abs(
                            (.5*tr)**2+(.5*skew)**2-(np.sqrt(phidet))**2))
#            ecolor=(abs(phi.min())+abs(phi.max()))/2
            ecolor=np.sign(pt1.phimax[ii]-pt2.phimin[ii])*\
                    (abs(phi.min())+abs(phi.max()))/2
            
            #put things into arrays
            phimaxarr[ii,ss]=phimax
            phiminarr[ii,ss]=phimin
            azimutharr[ii,ss]=azimuth
            betaarr[ii,ss]=beta
            colorarr[ii,ss]=ecolor
            
    #===============================================================================
    # Filter the arrays if desired
    #===============================================================================
    phimaxarr=sps.medfilt2d(phimaxarr,kernel_size=mfs)
    phiminarr=sps.medfilt2d(phiminarr,kernel_size=mfs)
    azimutharr=sps.medfilt2d(azimutharr,kernel_size=mfs)
    betaarr=sps.medfilt2d(betaarr,kernel_size=mfs)
    colorarr=sps.medfilt2d(colorarr,kernel_size=mfs)
    
    #normalize ecolor
    #cmax=colorarr.max()

#    ecolorarr=colorarr/cmax
#    cpass=np.where(abs(ecolorarr)>1)
#    ecolorarr[cpass]=1

    #===============================================================================
    # Pickle results so don't have to reload them everytime
    #===============================================================================
    fid=file(pkfn,'w')
    pickle.dump((phimaxarr,phiminarr,azimutharr,betaarr,colorarr,latlst,lonlst,
                 stationlst,imp1.period),fid)
    fid.close()


#===============================================================================
# Plot ellipses in map view
#===============================================================================

#load pickled file
pkfid=file(pkfn,'r')
phimaxarr,phiminarr,azimutharr,betarr,ecolorarr,latlst,lonlst,stationlst,period=\
                                                            pickle.load(pkfid)
pkfid.close()

#ecolorarr=ecolorarr/cmax
#cpass=np.where(abs(ecolorarr)>1)
#ecolorarr[cpass]=1

ecolorarr=np.nan_to_num(ecolorarr)
ecmax=abs(ecolorarr).max()
#ecmax=90.
#ecmax=180.
nrows=len(prange)/ncols

plt.rcParams['font.size']=6
plt.rcParams['figure.subplot.left']=.1
plt.rcParams['figure.subplot.right']=.91
plt.rcParams['figure.subplot.bottom']=.08
plt.rcParams['figure.subplot.top']=.95
plt.rcParams['figure.subplot.hspace']=.005
plt.rcParams['figure.subplot.wspace']=.005


emax=2
fig=plt.figure(fignum,[14,14],dpi=300)
for ii,ff in enumerate(prange,1):
    ax1=fig.add_subplot(nrows,ncols,ii,aspect='equal')
#    esized=a*np.mean([phimaxarr[ff,:],phiminarr[ff,:]])
#    emax=5*np.median([phimaxarr[ff,:],phiminarr[ff,:]])*esized
    for ss in range(ns):
#        eheightd=phimaxarr[ff,ss]*esized
#        ewidthd=phiminarr[ff,ss]*esized
        eheightd=phimaxarr[ff,ss]/phimaxarr[ff,:].max()*esize
        ewidthd=phiminarr[ff,ss]/phimaxarr[ff,:].max()*esize
#        eheightd=phiminarr[ff,ss]/phimaxarr[ff,:].max()*esize
#        ewidthd=phimaxarr[ff,ss]/phimaxarr[ff,:].max()*esize
        
        if eheightd>emax or ewidthd>emax:
            print 'face'
            pass
            
        else:
            ellipd=Ellipse((lonlst[ss],latlst[ss]),
                           width=ewidthd,
                           height=eheightd,
                           angle=azimutharr[ff,ss])
#            cvar=betarr[ff,ss]/ecmax
            cvar=ecolorarr[ff,ss]/ecmax
#            cvar=azimutharr[ff,ss]/ecmax
#            cvar=phiminarr[ff,ss]/abs(phiminarr).max()
#            if ecolorarr[ff,ss]<0:
#                ellipd.set_facecolor((1-cvar,1,cvar))
#            else:
#                ellipd.set_facecolor((1,1-cvar,.1))
            if cvar<0:
                ellipd.set_facecolor((1-abs(cvar),1,abs(cvar)))
            else:
                ellipd.set_facecolor((1,1-abs(cvar),.1))
#            if phiminarr[ff,ss]<0:
#                ellipd.set_facecolor((1-abs(cvar),1,.5+abs(cvar)/2))
#            else:
#                ellipd.set_facecolor((1,1-abs(cvar),.1))
#            if betarr[ff,ss]<0:
#                ellipd.set_facecolor((1-cvar,1,cvar))
#            else:
#                ellipd.set_facecolor((1,1-cvar,.1))
            ax1.add_artist(ellipd)
            arrow=Arrow(lonlst[ss],latlst[ss],
                        np.sin(-azimutharr[ff,ss]*np.pi/180)*eheightd,
                        np.cos(-azimutharr[ff,ss]*np.pi/180)*eheightd,
                        width=.05)
            ax1.add_patch(arrow)
            
    
    ax1.set_xlim(xlimits)
    ax1.set_ylim(ylimits)

    ax1.text(xlimits[0]+.05,ylimits[1]-.1,'T={0:.2g}'.format(period[ff]),
             verticalalignment='top',horizontalalignment='left',
             fontdict={'size':8,'weight':'bold'})
    ax1.text(0,0,'X',
             verticalalignment='center',
             horizontalalignment='center',
             fontdict={'size':9,'weight':'bold'})
    ellips=Ellipse((xlimits[1]-esize*1.2,ylimits[0]+esize/2),width=esize,height=esize,
                   angle=0)
    ellips.set_facecolor((.1,.1,1.))
    ax1.add_artist(ellips)
    ax1.grid(alpha=.2)
    
    ax1.text(xlimits[1]-esize*1.3,ylimits[0]+esize*1.4,
             '$\Delta$={0:.2g}'.format(esize*phimaxarr[ff,:].max()*a),
             horizontalalignment='center',
             verticalalignment='baseline')

    if ii==len(prange)-1:
        ax1.set_xlabel('easting (km)',fontdict={'size':9,'weight':'bold'})
    if ii<4:
        ax1.xaxis.set_ticklabels(['' for hh in 
                                        range(len(ax1.xaxis.get_ticklabels()))])
    if ii==1 or ii==ncols+1:
        pass
    else:
        ax1.yaxis.set_ticklabels(['' for hh in 
                                        range(len(ax1.yaxis.get_ticklabels()))])
    if ii==1:
        ax1.set_ylabel('northing (km)',fontdict={'size':9,'weight':'bold'})

#add colorbar
ax2=fig.add_subplot(1,1,1)
ax2.set_visible(False)
cbax=make_axes(ax2,shrink=.90,fraction=.01,pad=.2)
#cbx=ColorbarBase(cbax[0],cmap=ptcmap,norm=Normalize(vmin=0,vmax=cmax*a),
#                orientation='vertical',format='%.2g')
cbx=ColorbarBase(cbax[0],cmap=ptcmap3,norm=Normalize(vmin=-ecmax,vmax=ecmax),
                orientation='vertical',format='%.2g')
#cbx.set_label('Beta',fontdict={'size':7,'weight':'bold'})
cbx.set_label('sgn($\Phi_{max}^1-\Phi_{max}^2$)(|$\Delta_{max}$|+|$\Delta_{min}$|)/2 ',
                  fontdict={'size':7,'weight':'bold'})
plt.show()


#===============================================================================
# Plot data for one period
#===============================================================================

#plt.rcParams['font.size']=7
#plt.rcParams['figure.subplot.left']=.11
#plt.rcParams['figure.subplot.right']=.98
#plt.rcParams['figure.subplot.bottom']=.08
#plt.rcParams['figure.subplot.top']=.98
#plt.rcParams['figure.subplot.hspace']=.005
#plt.rcParams['figure.subplot.wspace']=.005
#
#if ctype=='fm':
#    ff=8
#else:
#    ff=28
#xlimits=(-3,3)
#ylimits=(-2.4,1.75)
#
#
#
#
#emax=1
#fig=plt.figure(1,[14,14],dpi=300)
#ax1=fig.add_subplot(1,1,1,aspect='equal')
#for ss in range(ns):
#    eheightd=phimaxarr[ff,ss]/phimaxarr[ff,:].max()*esize
#    ewidthd=phiminarr[ff,ss]/phimaxarr[ff,:].max()*esize
#    
#    if eheightd>emax or ewidthd>emax:
#        print 'face'
#        
#    else:
#        ellipd=Ellipse((lonlst[ss],latlst[ss]),
#                       width=ewidthd,
#                       height=eheightd,
#                       angle=azimutharr[ff,ss])
##        ellipd.set_facecolor((1,1-ecolorarr[ff,ss],.1))
#        ellipd.set_facecolor((1,1-abs(betarr[ff,ss])/90.,.1))
##       if cvars<0:
##           ellip.set_facecolor((.1,1+cvars,abs(cvars)))
##       else:
##           ellip.set_facecolor((1,1-cvars,.1))
#           
##        if ecolorarr[ff,ss]<0:
##            ellipd.set_facecolor((.1,1-abs(ecolorarr[ff,ss]),1-abs(ecolorarr[ff,ss])))
##        else:
##            
#        ax1.add_artist(ellipd) 
#        
#
#ax1.set_xlim(xlimits)
#ax1.set_ylim(ylimits)
#
#ax1.text(xlimits[0]+.05,ylimits[1]-.1,'T={0:.2g}'.format(period[ff]),
#         verticalalignment='top',horizontalalignment='left',
#         fontdict={'size':8,'weight':'bold'})
#ax1.text(0,0,'X',
#         verticalalignment='center',
#         horizontalalignment='center',
#         fontdict={'size':9,'weight':'bold'})
#ellips=Ellipse((2.45,-2.),width=esize/2,
#               height=esize/2,angle=0)
#ellips.set_facecolor((.1,.1,1.))
#ax1.add_artist(ellips)
#ax1.grid(alpha=.2)
#
#ax1.text(2.45,-2.+esize/3.3,
#         '$\Delta$={0:.2g}'.format(esize/2*phimaxarr[ff].max()),
#         horizontalalignment='center',
#         verticalalignment='baseline')
#
#ax1.set_xlabel('easting (km)',fontdict={'size':9,'weight':'bold'})
#ax1.set_ylabel('northing (km)',fontdict={'size':9,'weight':'bold'})
#ax1.xaxis.set_minor_locator(MultipleLocator(.25))
#ax1.yaxis.set_minor_locator(MultipleLocator(.25))
#
##add colorbar
#cbax=make_axes(ax1,shrink=.40,orientation='horizontal')
#cbx=ColorbarBase(cbax[0],cmap=ptcmap,norm=Normalize(vmin=0,vmax=cmax*a),
#                orientation='horizontal',format='%.2g')
#cbx.set_ticks(MultipleLocator(.1))
#cbx.set_label('(|$\Delta_{max}$|+|$\Delta_{min}$|)/2 ',
#                  fontdict={'size':7,'weight':'bold'})
#plt.show()

