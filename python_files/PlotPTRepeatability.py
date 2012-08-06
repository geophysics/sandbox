# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Jared Peacock)s
"""

import numpy as np
import Z
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colorbar import *
from matplotlib.patches import Ellipse

ptcmapdict={'red':((0.0,1.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,0.0,1.0),(1.0,0.0,1.0)),
            'blue':((0.0,0.0,0.0),(1.0,0.0,0.0))}
ptcmap=LinearSegmentedColormap('ptcmap',ptcmapdict,256)

#edipath1=r"c:\Users\Own er\Documents\PHD\Geothermal\Paralana\EDIFilesBaseSurvey\March2010"
#edipath2=r"c:\Users\Own er\Documents\PHD\Geothermal\Paralana\EDIFilesBaseSurvey\May2010"
#edipath3=r"c:\Users\Own er\Documents\PHD\Geothermal\Paralana\EDIFilesBaseSurvey\March2011"

edipath1=r"c:\Peacock\PHD\Geothermal\Paralana\EDIFilesBaseSurvey\March2010\DR"
edipath2=r"c:\Peacock\PHD\Geothermal\Paralana\EDIFilesBaseSurvey\May2010\DR"
edipath3=r"c:\Peacock\PHD\Geothermal\Paralana\EDIFilesBaseSurvey\March2011\DR"


ediext='cdr.edi'
#station='pb23'
#
#stationlst=['pb0'+str(ii) for ii in range(1,10)]+\
#            ['pb'+str(ii) for ii in range(10,55)]
#            
stationlst=['pb'+nn for nn in ['03','06','09','14','15','16','17','23','24',
                               '25','26','29','30','35','38','41','42','52',
                               '54']]

#stationlst=['pb'+nn for nn in ['03','05','09','15','17','24','25','30','36','38',
#                               '41','42','44','52']]


#===============================================================================
plt.rcParams['font.size']=6
plt.rcParams['figure.subplot.left']=.11
plt.rcParams['figure.subplot.right']=.98
plt.rcParams['figure.subplot.bottom']=.1
plt.rcParams['figure.subplot.top']=.98
plt.rcParams['figure.subplot.hspace']=.005
plt.rcParams['figure.subplot.wspace']=.005

fig=plt.figure(1,[14,14],dpi=300)
ax1=fig.add_subplot(1,1,1,aspect='equal')

esize=4
emax=8
xspacing=6
yspacing=3
cmax=.15
ns,nf=len(stationlst),29

phiminarr=np.zeros((ns,nf))
phimaxarr=np.zeros((ns,nf))
azimutharr=np.zeros((ns,nf))
ecolorarr=np.zeros((ns,nf))
betaarr=np.zeros((ns,nf))
phiarr=np.zeros((ns,nf,2,2))
reparr=np.zeros((ns,nf,4))
kk=-1
pslst=[]

for ss,station in enumerate(stationlst):
    implst=[]
    if os.path.isfile(os.path.join(edipath1,station+ediext))==True:
        imp1=os.path.join(edipath1,station+ediext)
        implst.append(imp1)
        
    if os.path.isfile(os.path.join(edipath2,station+ediext))==True:
        imp2=os.path.join(edipath2,station+ediext)
        implst.append(imp2)
    if os.path.isfile(os.path.join(edipath3,station+ediext))==True:
        imp3=os.path.join(edipath3,station+ediext)
        implst.append(imp3)
    
    if implst==[] or len(implst)==1: 
        pass
    else:
        kk+=1
        pslst.append(station[2:4])
        implst=[Z.Z(implst[0]),Z.Z(implst[1])]
        period=implst[0].period
        nf=len(period)

        pt1=implst[0].getPhaseTensor()
        pt2=implst[1].getPhaseTensor()
        
        
        for ii in range(nf):
#            reparr[kk,ii,0]=(1-pt2.phimin[ii]/pt1.phimin[ii])*100            
#            reparr[kk,ii,1]=(1-pt2.phimax[ii]/pt1.phimax[ii])*100 
            reparr[kk,ii,0]=(1-(pt2.phi[ii][0,0]+pt2.phi[ii][1,1])/ 
                            (pt1.phi[ii][0,0]+pt1.phi[ii][1,1]))*100
            reparr[kk,ii,1]=(1-(pt2.phi[ii][0,1]-pt2.phi[ii][1,0])/ 
                            (pt1.phi[ii][0,1]-pt1.phi[ii][1,0]))*100
            reparr[kk,ii,2]=(1-np.linalg.det(pt2.phi[ii])/\
                            np.linalg.det(pt1.phi[ii]))*100  
#            reparr[kk,ii,2]=(1-np.sqrt(abs(pt2.phimin[ii]*pt2.phimax[ii]))/
#                            np.sqrt(abs(pt1.phimin[ii]*pt1.phimax[ii])))*100    
#            reparr[kk,ii,3]=(1-pt2.alpha[ii]/pt1.alpha[ii])*100  
#            phi=np.eye(2)-\
#                    (np.dot(np.linalg.inv(pt2.phi[ii]),pt1.phi[ii])+
#                    np.dot(pt1.phi[ii],np.linalg.inv(pt2.phi[ii])))/2       
            phi=np.eye(2)-np.dot(np.linalg.inv(pt1.phi[ii]),pt2.phi[ii]) 
            
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
#            ecolor=abs((phi.min()+phi.max())/2)/cmax
            ecolor=np.sqrt(abs((phi.min()*phi.max())))*cmax
            if ecolor>1:
                ecolor=1
                
            phiminarr[kk,ii]=phimin
            phimaxarr[kk,ii]=phimax
            azimutharr[kk,ii]=azimuth
            ecolorarr[kk,ii]=ecolor
            betaarr[kk,ii]=beta
            phiarr[kk,ii,:,:]=phi
                
            
            #put things into arrays
            eheightd=phimax*esize
            ewidthd=phimin*esize
            
            if eheightd>emax or ewidthd>emax:
                pass
                
            else:
                ellipd=Ellipse((xspacing*kk,yspacing*(nf-ii)),
                               width=ewidthd,
                               height=eheightd,
                               angle=azimuth)
                ellipd.set_facecolor((1,1-ecolor,.1))
                ax1.add_artist(ellipd) 
        

ax1.set_xlim(-esize/2,kk*xspacing+esize/2)
ax1.set_ylim(-esize/2,nf*yspacing+esize/2)
#
#ax1.text(xlimits[0]+.05,ylimits[1]-.1,'T={0:.2g}'.format(period[ff]),
#         verticalalignment='top',horizontalalignment='left',
#         fontdict={'size':8,'weight':'bold'})
#ax1.text(0,0,'X',
#         verticalalignment='center',
#         horizontalalignment='center',
#         fontdict={'size':9,'weight':'bold'})
ellips=Ellipse((kk*xspacing,yspacing/2),width=esize,
               height=esize,angle=0)
ellips.set_facecolor((.1,.1,1.))
ax1.add_artist(ellips)
ax1.grid(alpha=.2)
#
#ax1.text(2.45,-2.+esize/3.3,
#         '$\Delta$={0:.2g}'.format(esize/2*phimaxarr[ff].max()),
#         horizontalalignment='center',
#         verticalalignment='baseline')
#
ax1.set_xlabel('station',fontdict={'size':9,'weight':'bold'})
ax1.set_ylabel('period (s)',fontdict={'size':9,'weight':'bold'})
ax1.xaxis.set_major_locator(MultipleLocator(xspacing))
ax1.yaxis.set_major_locator(MultipleLocator(yspacing*2))
ax1.xaxis.set_ticklabels(pslst)
ax1.yaxis.set_ticklabels(['{0:.2g}'.format(jj) for jj in period[range(-1,-nf-2,-2)]])

#add colorbar
cbax=make_axes(ax1,shrink=.40,orientation='vertical')
cbx=ColorbarBase(cbax[0],cmap=ptcmap,norm=Normalize(vmin=0,vmax=cmax),
                orientation='vertical',format='%.2g')
cbx.set_label('(|$\Delta_{max}$|+|$\Delta_{min}$|)/2 ',
                  fontdict={'size':7,'weight':'bold'})
plt.show()

    

#===============================================================================
# Plot data
#===============================================================================
#plt.rcParams['font.size']=6
#plt.rcParams['figure.subplot.left']=.11
#plt.rcParams['figure.subplot.right']=.98
#plt.rcParams['figure.subplot.bottom']=.08
#plt.rcParams['figure.subplot.top']=.98
#plt.rcParams['figure.subplot.hspace']=.005
#plt.rcParams['figure.subplot.wspace']=.005
#
#ff=28
#xlimits=(-3,3)
#ylimits=(-2.4,1.75)
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
#        pass
#        
#    else:
#        ellipd=Ellipse((lonlst[ss],latlst[ss]),
#                       width=ewidthd,
#                       height=eheightd,
#                       angle=azimutharr[ff,ss])
#        ellipd.set_facecolor((1,1-ecolorarr[ff,ss],.1))
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
#cbx.set_label('(|$\Delta_{max}$|+|$\Delta_{min}$|)/2 ',
#                  fontdict={'size':7,'weight':'bold'})
#plt.show()
#

