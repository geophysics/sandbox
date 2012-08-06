# -*- coding: utf-8 -*-
"""
Created on Wed May 25 16:39:30 2011

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
#from matplotlib.ticker import FormatStrFormatter,MultipleLocator
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
ttype='pt'


llst=[':k',':k',':k',':k']*2
mlst=['*','o','s','v','h','x']*2
clst=['blue','maroon','purple','green','orange','red','gray','brown']*2

#===============================================================================
# Initialize parameters and locate directory path
#===============================================================================
if ctype=='data':
#    edipath=r"c:\Users\Own er\Documents\PHD\Geothermal\Paralana\EDIFilesInjection\EDIfiles\AdvPro24Hrs"
    edipath=r"c:\Peacock\PHD\Geothermal\Paralana\EDIFilesInjection\EDIfiles\AdvPro24Hrs"
    stationlst=['pb13','pb15','pb01','pb04','pb24','pb25','pb27','pb49']
    daylst=['','192','193','194','195','196','197']
    
    nf=36
    nd=len(daylst)
    ns=(len(stationlst))
    
    pstart=15
    pstop=24
    if ttype=='pt':
#        prange=[15,16,17,18,19,20,21,22,23]
#        prange=[8,17,22,28]
        prange=[8,17,21,23]
#        prange=np.arange(start=6,stop=33,step=3)
#        prange=[8,12,15,17,18,21,22,25,28]
        fignum=4
        ptm='phimax'
        nrow=2
        ncol=2
    elif ttype=='rt':
#        prange=[20,21,22,23,24,25,26,27,28]
        prange=[22,24,26,29]
        fignum=3
        ptm='rhodet'

elif ctype=='fm':
    edipath=r"c:\Peacock\PHD\Geothermal\Paralana\ForwardModels\ParalanaBase\EDIFiles"
    stationlst=['par17ew','par19ew','par21ew','par76','par87','par98','par67','par68',
                'par69','par64','par62','par47']#,'par05ew']
    daylst=['0{0}'.format(ii) for ii in range(7)]
    
    nf=14
    nd=len(daylst)
    ns=(len(stationlst))
    
    pstart=4
    pstop=10

pkfn=r"c:\Peacock\PHD\ParalanaTimeLapsePT"+ctype+ttype+'.pkl'
#===============================================================================
# Compute difference and put into arrays
#===============================================================================
if not os.path.isfile(pkfn):

    latlst=np.zeros(ns)
    lonlst=np.zeros(ns)
    phiminarr=np.zeros((nd+2,nf,ns,2))
    phimaxarr=np.zeros((nd+2,nf,ns,2))
    azimutharr=np.zeros((nd+2,nf,ns,2))
    betaarr=np.zeros((nd+2,nf,ns,2))
    colorarr=np.zeros((nd+2,nf,ns))
    zerr=np.zeros((nd+2,nf,ns))
    perr=np.zeros((nd+2,nf,ns))
    
    for ss,station in enumerate(stationlst):
        for dd,day in enumerate(daylst):
            
            edi1=os.path.join(edipath,station+daylst[0]+'.edi')
            if os.path.exists(edi1)==True:
                print edi1
            else:
                edi1=os.path.join(edipath,station+daylst[1]+'00.edi')
                if os.path.exists(edi1)==True:
                    print edi1
                else:
                    raise IOError('File '+edi1+'does not exist')
            
            edi2=os.path.join(edipath,station+day+'.edi')
            if os.path.exists(edi2)==True:
                print edi2
            else:
                edi2=os.path.join(edipath,station+day+'00.edi')
                if os.path.exists(edi2)==True:
                    print edi2
                else:
                    print 'Skipping File '+edi2
                    edi2=edi1
            print '-'*10
                    
            z1=Z.Z(edi1)
            z2=Z.Z(edi2)
            
            p1=z1.period
            p2=z2.period
            
            if ttype=='pt':
                pt1=z1.getPhaseTensor(rotate=rot)
                pt2=z2.getPhaseTensor(rotate=rot)
            elif ttype=='rt':
                pt1=z1.getResTensor(rotate=rot)
                pt2=z2.getResTensor(rotate=rot)
            
#            zerr[dd+1,:len(z2.period),ss]=np.array([abs(np.median(
#                                               np.dot(z2.zvar[kk],
#                                               np.linalg.inv(z2.z[kk]))))
#                                               for kk in range(len(z2.period))])            
            zerr[dd+1,:len(z2.period),ss]=np.array([abs(np.linalg.det(
                                               np.dot(z2.zvar[kk],
                                               np.linalg.inv(z2.z[kk]))))**.5
                                               for kk in range(len(z2.period))])+\
                                           np.array([abs(np.linalg.det(
                                               np.dot(z1.zvar[kk],
                                               np.linalg.inv(z1.z[kk]))))**.5
                                               for kk in range(len(z1.period))])            
            
            for nn,f1 in enumerate(p1):
                for mm,f2 in enumerate(p2):
                    if f1==f2:
                        if ttype=='pt':
#                            phi=np.eye(2)-\
#                                (np.dot(pt2.phi[mm],
#                                        np.linalg.inv(pt1.phi[nn]))+
#                                np.dot(pt1.phi[nn],
#                                       np.linalg.inv(pt2.phi[mm])))/2
                            phi=np.eye(2)-np.dot(np.linalg.inv(pt1.phi[mm]),
                                                   pt2.phi[nn])
                                       
                            dphi=pt1.phivar[mm]+pt2.phivar[mm]
                        elif ttype=='rt':
                            phi=np.eye(2)-\
                                (np.dot(pt2.rho[mm],
                                        np.linalg.inv(pt1.rho[nn]))+
                                np.dot(pt1.rho[nn],
                                       np.linalg.inv(pt2.rho[mm])))/2
                            dphi=np.repeat(zerr[dd,nn,ss],4)
                            dphi=dphi.reshape(2,2)
#                        phi=np.eye(2)-\
#                            (np.dot(pt1.phi[nn],
#                                    np.linalg.inv(pt2.phi[mm])))
                        
#                        dphi=np.eye(2)-np.dot(pt1.phivar[nn],
#                                            np.linalg.inv(pt2.phivar[mm]))
                        
                         
                        
#                        #put variance into standard deviation
#                        zvar=z1.zvar[nn]**2-z2.zvar[mm]**2
#                        #create a matrix for errors to be calculated
#                        dphi=np.zeros(np.shape(z[ii]))
#                        #compute determinate of X
#                        detX=np.linalg.det(X)
#                        #calculate the deteriminate of the error matrix 
#                        ddet=np.sqrt((X[0,0]*X[1,1])**2*((zvar[0,0]/X[0,0])**2+
#                            (zvar[1,1]/X[1,1])**2)+(X[1,0]*X[0,1])**2*(
#                            (zvar[0,1]/X[0,1])**2+(zvar[1,0]/X[1,0])**2))
#                        #calculate errors for each component of the matrix ala Caldwell 
#                        #2004
#                        dphi[0,0]=np.sqrt((X[1,1]*Y[0,0])**2*((zvar[1,1]/X[1,1])**2+
#                            (zvar[0,0]/Y[0,0])**2)+(X[0,1]*Y[1,0])**2*(
#                            (zvar[0,1]/X[0,1])**2+(zvar[1,0]/X[1,0])**2))
#                        dphi[0,1]=np.sqrt((X[1,1]*Y[0,1])**2*((zvar[1,1]/X[1,1])**2+
#                            (zvar[0,1]/Y[0,1])**2)+(X[0,1]*Y[1,1])**2*(
#                            (zvar[0,1]/X[0,1])**2+(zvar[1,1]/X[1,1])**2))
#                        dphi[1,0]=np.sqrt((X[0,1]*Y[1,0])**2*((zvar[0,0]/X[0,0])**2+
#                            (zvar[1,0]/Y[1,0])**2)+(X[1,0]*Y[0,0])**2*(
#                            (zvar[1,0]/X[1,0])**2+(zvar[0,0]/X[0,0])**2))
#                        dphi[1,1]=np.sqrt((X[0,0]*Y[1,1])**2*((zvar[0,0]/X[0,0])**2+
#                            (zvar[1,1]/Y[1,1])**2)+(X[1,0]*Y[0,1])**2*(
#                            (zvar[1,0]/X[1,0])**2+(zvar[0,1]/X[0,1])**2))
#                        #rotate the error matrix
#                        dphi=np.rot90(dphi,2)
#                        #finish calculating the errors
#                        dphi[0,0]=(phi[0,0]/detX)**2*np.sqrt((dphi[0,0]/phi[0,0])**2+
#                                                                        (ddet/detX)**2)
#                        dphi[0,1]=(phi[0,1]/detX)**2*np.sqrt((dphi[0,1]/phi[0,1])**2+
#                                                                        (ddet/detX)**2)
#                        dphi[1,0]=(phi[1,0]/detX)**2*np.sqrt((dphi[1,0]/phi[1,0])**2+
#                                                                        (ddet/detX)**2)
#                        dphi[1,1]=(phi[1,1]/detX)**2*np.sqrt((dphi[1,1]/phi[1,1])**2+
#                                                                        (ddet/detX)**2)
#                        
                        #Calculate Trace of Phi and error of trace of phi
                        tr=phi[0,0]+phi[1,1]
                        trvar=np.sqrt(dphi[0,0]**2+dphi[1,1]**2)
                        
                        #Calculate skew of phi and the cooresponding error
                        skew=phi[0,1]-phi[1,0]
                        skewvar=np.sqrt(dphi[0,1]**2+dphi[1,1]**2)
                        
                        #calculate the determinate and determinate error of phi
                        phidet=abs(np.linalg.det(phi))
                        phidetvar=np.sqrt((np.sqrt((dphi[0,0]/phi[0,0])**2+(
                            dphi[1,1]/phi[1,1])**2)*phi[0,0]*phi[1,1])**2+(
                            np.sqrt((dphi[0,1]/phi[0,1])**2+(
                            dphi[1,0]/phi[1,0])**2)*phi[0,1]*phi[1,0])**2)
                        #calculate reverse trace and error
                        revtr=phi[0,0]-phi[1,1]
                        revtrvar=np.sqrt(dphi[0,0]**2+dphi[1,1]**2)
                        
                        #calculate reverse skew and error
                        revskew=phi[1,0]+phi[0,1]
                        revskewvar=np.sqrt(phi[0,1]**2+dphi[1,0]**2)
                        
                        #calculate skew angle beta and error
                        beta=.5*(np.arctan2(skew,tr)*180/np.pi)
                        betavar=abs(np.arctan(skew*tr*np.sqrt((skewvar/skew)**2+(
                                                            trvar/tr)**2))*180/np.pi)
                        
                        #calculate angle alpha corresponding to phase tensor's 
                        #dependence on coordinate system
                        alpha=.5*(np.arctan2(revskew,revtr)*180/np.pi)
                        alphavar=abs(.5*np.arctan(revskew*revtr*np.sqrt(
                                (revskewvar/revskew)**2+(revtrvar/revtr)**2))*180/np.pi)
                        
                        #calculate azimuth as angle between major axis and x-axis
                        azimuth=alpha-beta
                        azimuthvar=np.sqrt(alphavar**2+betavar**2)
                        
                        #calulate maximum value for phi
                        phimax=np.sqrt((.5*tr)**2+(.5*skew)**2)+\
                            np.sqrt(abs((.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2))
                        phimaxvar=.5*np.sqrt(2*trvar**2+2*skewvar**2)+.5*np.sqrt(
                                                    2*trvar**2+2*skewvar**2+phidetvar)
                        
                        #calculate minimum value for phi
                        if np.linalg.det(phi)>=0:
                            phimin=np.sqrt((.5*tr)**2+(.5*skew)**2)-\
                            np.sqrt((.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2)
                        elif np.linalg.det(phi)<0:
                            phimin=-1*np.sqrt((.5*tr)**2+(.5*skew)**2)-np.sqrt(
                                        (.5*tr)**2+(.5*skew)**2-
                                        (np.sqrt(abs(phidet)))**2)
                        phiminvar=phimaxvar
                        ecolor=abs((phi.min()+phi.max())/2)
                        
                        #put things into arrays
                        phimaxarr[dd+1,nn,ss,0]=phimax
                        phiminarr[dd+1,nn,ss,0]=phimin
                        azimutharr[dd+1,nn,ss,0]=azimuth
                        betaarr[dd+1,nn,ss,0]=beta
                        phimaxarr[dd+1,nn,ss,1]=phimaxvar
                        phiminarr[dd+1,nn,ss,1]=phiminvar
                        azimutharr[dd+1,nn,ss,1]=azimuthvar
                        betaarr[dd+1,nn,ss,1]=betavar
                        if ttype=='rt':
                            colorarr[dd+1,nn,ss]=abs(1-pt2.rhodet[mm]/pt1.rhodet[nn])
                        else:
                            colorarr[dd+1,nn,ss]=ecolor
    
    for ss,station in enumerate(stationlst):
        #get base edifile        
        try:
            z1=Z.Z(os.path.join(edipath,station+daylst[0]+'.edi'))
            
        except IOError:
            try:
                z1=Z.Z(os.path.join(edipath,station+daylst[1]+'00.edi'))
            except IOError:
                raise IOError('No base file found')                    
        latlst[ss]=z1.lat
        lonlst[ss]=z1.lon
    period=z1.period
        
        
    pkfid=file(pkfn,'w')
    pickle.dump((latlst,lonlst,stationlst,phiminarr,phimaxarr,betaarr,azimutharr,
                colorarr,period,zerr),pkfid)
    pkfid.close()
#===========================================================================
# Plot ellipses
#===========================================================================

else:
    pkfid=file(pkfn,'r')
    latlst,lonlst,stationlst,phiminarr,phimaxarr,betaarr,azimutharr,colorarr,period,zerr=pickle.load(pkfid)
    pkfid.close()    

#pad the arrays to remove edge effects from median filter
#phiminarr[0,:,:,:]=phiminarr[1,:,:,:]
#phimaxarr[0,:,:,:]=phimaxarr[1,:,:,:]
#betaarr[0,:,:,:]=betaarr[1,:,:,:]
#azimutharr[0,:,:,:]=azimutharr[1,:,:,:]
#colorarr[0,:,:]=colorarr[1,:,:,]
#zerr[0,:,:]=zerr[1,:,:,]
#
phiminarr[-1,:,:,:]=phiminarr[-2,:,:,:]
phimaxarr[-1,:,:,:]=phimaxarr[-2,:,:,:]
betaarr[-1,:,:,:]=betaarr[-2,:,:,:]
azimutharr[-1,:,:,:]=azimutharr[-2,:,:,:]
colorarr[-1,:,:]=colorarr[-2,:,:,]
zerr[-1,:,:]=zerr[-2,:,:,]

#filter data if desires
for ss in range(ns):
    phiminarr[:,:,ss,0]=sps.medfilt2d(phiminarr[:,:,ss,0],kernel_size=mks)
    phimaxarr[:,:,ss,0]=sps.medfilt2d(phimaxarr[:,:,ss,0],kernel_size=mks)
    betaarr[:,:,ss,0]=sps.medfilt2d(betaarr[:,:,ss,0],kernel_size=mks)
    azimutharr[:,:,ss,0]=sps.medfilt2d(azimutharr[:,:,ss,0],kernel_size=mks)
    #error
    phiminarr[:,:,ss,1]=sps.medfilt2d(phiminarr[:,:,ss,1],kernel_size=mks)
    phimaxarr[:,:,ss,1]=sps.medfilt2d(phimaxarr[:,:,ss,1],kernel_size=mks)
    betaarr[:,:,ss,1]=sps.medfilt2d(betaarr[:,:,ss,1],kernel_size=mks)
    azimutharr[:,:,ss,1]=sps.medfilt2d(azimutharr[:,:,ss,1],kernel_size=mks)
    colorarr[:,:,ss]=sps.medfilt2d(colorarr[:,:,ss],kernel_size=mks)
    zerr[:,:,ss]=sps.medfilt2d(zerr[:,:,ss],kernel_size=mks)
#for ff in range(nf):
#    phiminarr[:,ff,:,0]=sps.medfilt2d(phiminarr[:,ff,:,0],kernel_size=mks)
#    phimaxarr[:,ff,:,0]=sps.medfilt2d(phimaxarr[:,ff,:,0],kernel_size=mks)
#    betaarr[:,ff,:,0]=sps.medfilt2d(betaarr[:,ff,:,0],kernel_size=mks)
#    azimutharr[:,ff,:,0]=sps.medfilt2d(azimutharr[:,ff,:,0],kernel_size=mks)
#    colorarr[:,ff,:]=sps.medfilt2d(colorarr[:,ff,:],kernel_size=mks)
#    zerr[:,ff,:]=sps.medfilt2d(zerr[:,ff,:],kernel_size=mks)

cmax=3*colorarr.mean()
colorarr=colorarr/cmax
   
   
#===============================================================================
# Plot errors for pt parameters   
#===============================================================================

pslst=list(stationlst)
try:
    pslst.remove('pb49')
except ValueError:
    pass


plt.rcParams['font.size']=12
plt.rcParams['figure.subplot.left']=.08
plt.rcParams['figure.subplot.right']=.97
plt.rcParams['figure.subplot.bottom']=.12
plt.rcParams['figure.subplot.top']=.90
plt.rcParams['figure.subplot.wspace']=.09
plt.rcParams['figure.subplot.hspace']=.09


pumping=np.zeros(nd)
pumping[3:-2]=20

fig2=plt.figure(fignum,dpi=250)

for jj,ff in enumerate(prange,1):
    ax=fig2.add_subplot(nrow,ncol,jj)
    
    if ttype=='pt':
        ax.axvline(194.5,ymax=20,lw=120,color='orange',alpha=.25)
        mederrd=ax.plot(np.arange(nd)+191,
                        (np.max(zerr,axis=2)[1:nd+1,ff]*4*100)**2,
                        '-',color='k',lw=1)
        ax.fill_between(np.arange(nd)+191,0,
                        (np.max(zerr,axis=2)[1:nd+1,ff]*4*100)**2+.16,
                         facecolor='gray', alpha=0.5)
#        ax.fill_between(np.arange(nd)+191,0,pumping,
#                        facecolor='orange',alpha=.25)
        
        ax.set_ylim(0,19)
        ax.text(.2+191,18.,'T={0:.2g}(s)'.format(period[ff]),
            fontdict={'size':16,'weight':'bold'},
            horizontalalignment='left',
            verticalalignment='top',
            bbox={'facecolor':'white'})
    elif ttype=='rt':
        mederrd=ax.plot(np.arange(nd)+191,
                        (np.median(zerr,axis=2)[1:nd+1,ff]**2)*100,
                        '-',color='k',lw=1)
        ax.fill_between(np.arange(nd)+191,0,
                        (np.median(zerr,axis=2)[1:nd+1,ff]**2)*100,
                         facecolor='gray', alpha=0.5)
#        mederrd=ax.plot(np.arange(nd)+191,
#                        (.2*period[ff]*np.median(zerr,axis=2)[1:nd+1,ff]**2)*5000,
#                        '-',color='k',lw=1)
#        ax.fill_between(np.arange(nd)+191,0,
#                        (.2*period[ff]*np.median(zerr,axis=2)[1:nd+1,ff]**2)*5000,
#                         facecolor='gray', alpha=0.5)
        ax.set_ylim(0,99)
        ax.text(.2+191,3.8,'T={0:.2g}(s)'.format(period[ff]),
            fontdict={'size':16,'weight':'bold'},
            horizontalalignment='left',
            verticalalignment='top',
            bbox={'facecolor':'white'})
#    e1=ax.errorbar(np.arange(nd)+191,
#                   (180/np.pi)*np.median(phimaxarr[1:nd+1,ff,:,0],axis=1),
#                    yerr=1./(np.median(phimaxarr[1:nd+1,ff,:,1],axis=1)*180/np.pi),
#                    marker='h',ms=4,mfc='None',mec='orange',
#                    mew=1,ls='--',ecolor='orange',color='orange')

    ax.set_xlim(191,197)

    lines=[]
    for ss in range(len(stationlst)):
        if stationlst[ss]=='pb49':
            pass
        else:
            if ptm=='phimax':
#                e1=ax.errorbar(np.arange(nd)+191,(180/np.pi)*phimaxarr[1:nd+1,ff,ss,0],
#                               yerr=1./(phimaxarr[1:nd+1,ff,ss,1]*180/np.pi),
#                                marker=mlst[ss],ms=4,mfc='None',mec='k',
#                                mew=1,ls=':',ecolor='k',color='k')
#                e1=ax.errorbar(np.arange(nd)+191,
#                               (180/np.pi)*phimaxarr[1:nd+1,ff,ss,0],
#                               yerr=phimaxarr[1:nd+1,ff,ss,1]*180/np.pi,
#                                marker=mlst[ss],ms=4,mfc='None',mec=clst[ss],
#                                mew=1,ls=':',ecolor=clst[ss],color=clst[ss])
#                e1=ax.errorbar(np.arange(nd)+191,
#                               np.sqrt(abs(phimaxarr[1:nd+1,ff,ss,0]*
#                                       phiminarr[1:nd+1,ff,ss,0]))*100,
#                               yerr=phimaxarr[1:nd+1,ff,ss,1]*100,
#                                marker=mlst[ss],ms=4,mfc='None',mec=clst[ss],
#                                mew=1,ls=':',ecolor=clst[ss],color=clst[ss])
                e1=ax.errorbar(np.arange(nd)+191,
                               phimaxarr[1:nd+1,ff,ss,0]*100,
                               yerr=phimaxarr[1:nd+1,ff,ss,1]*100,
                                marker=mlst[ss],ms=8,mfc='None',mec=clst[ss],
                                mew=2,ls=':',ecolor=clst[ss],color=clst[ss],
                                lw=3)
            elif ptm=='rhodet':
                e1=ax.errorbar(np.arange(nd)+191,(colorarr[1:nd+1,ff,ss])*100,
                               yerr=.2*period[ff]*abs(phimaxarr[1:nd+1,ff,ss,1])**2,
                                marker=mlst[ss],ms=4,mfc='None',mec=clst[ss],
                                mew=1,ls=':',ecolor=clst[ss],color=clst[ss])
            elif ptm=='phimin':
                e1=ax.errorbar(np.arange(nd)+191,
                               (180/np.pi)*phiminarr[:,ff,ss,0],
                               yerr=1./(phiminarr[:,ff,ss,1]*180/np.pi))
            elif ptm=='beta':
                e1=ax.errorbar(np.arange(nd)+191,abs(betaarr[1:nd+1,ff,ss,0]),
                               yerr=1./(betaarr[1:nd+1,ff,ss,1]))
                
            lines.append(e1[0])
    
#    ax.set_ylim(0,90)
#    ax.set_xlim(191,197)
#    medferr=ax.plot(np.arange(nd)+191,
#            np.repeat(np.median(np.median(zerr,axis=2),axis=0)[ff],nd)
#            *2*180/np.pi,'-.k')
#    lines.append(medferr)
#    pslst.append('Median_err (f)')
#    
#    meanerr=ax.plot(np.arange(nd)+191,
#            np.repeat(np.mean(zerr),nd)
#            *2*180/np.pi,':k')
#    lines.append(meanerr)
#    pslst.append('Mean_err')
#    
#    mederr=ax.plot(np.arange(nd)+191,
#            np.repeat(np.median(zerr),nd)
#            *2*180/np.pi,'--k')
#    lines.append(mederr)
#    pslst.append('Median_err')
#    
    
#    lines.append(mederrd)
#    pslst.append('Max_err(day)')
    
    ax.grid(alpha=.25)
    if jj==1 or jj==(nrow-1)+ncol:
        if ttype=='pt':
            ax.set_ylabel('$\mathbf{\Delta \Phi_{max}}$ $\mathbf{(\%)}$',
                          fontdict={'size':18,'weight':'bold'})
        elif ttype=='rt':
            ax.set_ylabel('Det(rho) ($\Omega \cdot$m)',
                          fontdict={'size':9,'weight':'bold'})
    else:
        plt.setp(ax.yaxis.get_ticklabels(),visible=False)

    if jj<((nrow-1)*ncol)+1:
        plt.setp(ax.xaxis.get_ticklabels(),visible=False)
    else:
        ax.set_xlabel('day (UTC)',fontdict={'size':16,'weight':'bold'})
fig2.legend(lines,pslst,ncol=len(pslst),loc='upper center',columnspacing=.8)
plt.show()
    


#===============================================================================
# Plot maps for each period    
#===============================================================================

#plt.rcParams['font.size']=12
#plt.rcParams['figure.subplot.left']=.08
#plt.rcParams['figure.subplot.right']=.99
#plt.rcParams['figure.subplot.bottom']=.06
#plt.rcParams['figure.subplot.top']=.96
#plt.rcParams['figure.subplot.wspace']=.005
#plt.rcParams['figure.subplot.hspace']=.005
#
#nrows=pstop-pstart
#period=p1
#    
#fig=plt.figure(1,[16,16])
#for ii,ff in enumerate(range(pstart,pstop)):
#    emax=1000*np.median([phiminarr[:,ff,:],phimaxarr[:,ff,:]])*esized
#    for jj,dd in enumerate(range(1,nd-1)):
#        ax1=fig.add_subplot(nrows,nd-1,(jj+1)+((nd-1)*ii),aspect='equal')
#        
#        for ss in range(ns):
#            if phimaxarr[dd,ff,ss,0]*esized>emax or phiminarr[dd,ff,ss,0]*esized>emax:
#                pass
#            else:
#                ewidth=phimaxarr[dd,ff,ss,0]*esized
#                eheight=phiminarr[dd,ff,ss,0]*esized
#                ellipd=Ellipse((lonlst[ss],latlst[ss]),
#                               width=ewidth,
#                               height=eheight,
#                               angle=azimutharr[dd,ff,ss,0])
#                ax1.add_artist(ellipd)
#                if colorarr[dd,ff,ss]>1:
#                    colorarr[dd,ff,ss]=0
#                ellipd.set_facecolor((1,1-colorarr[dd,ff,ss],.1))
#        if ii==0:
#            ax1.set_title(daylst[dd],fontdict={'size':18,'weight':'bold'})
#        if ii!=nrows-1:
#            ax1.xaxis.set_ticklabels(['' for hh in 
#                                        range(len(ax1.xaxis.get_ticklabels()))])
#        else:
#            ax1.xaxis.set_major_formatter(FormatStrFormatter('%3.2f'))
#            ax1.xaxis.set_major_locator(MultipleLocator(.015))
#        if jj==0:
##            ax1.set_ylabel('%.2g' % period[ff],
##                           fontdict={'size':12,'weight':'bold'})
#            ax1.text(max(lonlst),min(latlst),'T=%.2g(s)' % period[ff],
#                           fontdict={'size':12,'weight':'bold'},
#                    verticalalignment='bottom',horizontalalignment='center')
#            ax1.yaxis.set_major_formatter(FormatStrFormatter('%3.2f'))
#            ax1.yaxis.set_major_locator(MultipleLocator(.015))              
#        else:
#            ax1.yaxis.set_ticklabels(['' for hh in 
#                                        range(len(ax1.yaxis.get_ticklabels()))])
#
#        if ii==(nrows/2) and jj==0:
#            ax1.set_ylabel('Latitude',fontdict={'size':18,'weight':'bold'})
#        
#            
#            
#        if ii==nrows-1 and jj==2:
#            ax1.set_xlabel('Longitude',fontdict={'size':18,'weight':'bold'})
#        
#        #mark borehole
#        ax1.text(dhll[0],dhll[1],'X',fontdict={'size':12,'weight':'bold'},
#                 verticalalignment='center',horizontalalignment='center')
#        ax1.set_xlim(lonlst.min()-.005,lonlst.max()+.005)
#        ax1.set_ylim(latlst.min()-.005,latlst.max()+.005)
#        ax1.grid()
#        
#
#ax2=fig.add_subplot(1,1,1)
#ax2.set_visible(False)
#cbax=make_axes(ax2,shrink=.99)
#cbx=ColorbarBase(cbax[0],cmap=ptcmap,norm=Normalize(vmin=0,vmax=cmax),
#                orientation='vertical')
#cbx.set_label('(|$\Delta_{max}$|+|$\Delta_{min}$|)/2',fontdict={'size':18})
#plt.show()
