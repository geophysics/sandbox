# -*- coding: utf-8 -*-
"""
Created on Sat May 26 17:33:34 2012

@author: a1185872
"""

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import Z 
import os
from matplotlib.ticker import FormatStrFormatter,MultipleLocator,FixedLocator
#import matplotlib.gridspec as gridspec

#==============================================================================
# input parameters
#==============================================================================

edipath=r"c:\Peacock\PHD\Geothermal\Paralana\EDIFilesInjection\EDIfiles\AdvPro24Hrs"
#edipath=r"c:\Peacock\PHD\Geothermal\Paralana\EDIFilesInjection\EDIfiles\AdvPro4Hr"

stationlst=['pb01','pb02','pb04','pb24','pb25','pb27','pb49','pb54','pbrt1']
#stationlst=['pb01']

edis=[]
for station in stationlst:    
    edilst=[os.path.join(edipath,edi) for edi in os.listdir(edipath) 
            if edi.find(station)==0]
    edis.append(edilst)
daylst=[os.path.basename(edi)[4:7] for edi in edis[0]]
        
nedi=len(edilst)
#==============================================================================
# Plot data
#==============================================================================

plt.rcParams['font.size']=8
plt.rcParams['figure.subplot.left']=.08
plt.rcParams['figure.subplot.right']=.98
plt.rcParams['figure.subplot.bottom']=.1
plt.rcParams['figure.subplot.top']=.905
plt.rcParams['figure.subplot.hspace']=.00
plt.rcParams['figure.subplot.wspace']=.1

glst=[(ii,ii,ii) for ii in np.arange(start=0,stop=1,step=.1)]
clst=[(0,1-jj,0) for jj in np.arange(0,1,.1)]
clstb=[(0,0,1-jj) for jj in np.arange(0,1,.1)]

llst=[]

fig1=plt.figure(1,dpi=250)

#ylimlst=[(1,45),(2,9),(8,27),(26,57)]
ylimlst=[(10,45),(8,21),(26,57)]
#ylimlst=[(100,130),(20,60),(26,50),(26,61)]
#textlst=['$r_{xy}$','$r_{yx}$','$\phi_{xy}$','$\phi_{yx}$']
textlst=['$r_{xy}$','$\phi_{xy}$']
mc=(.5,.5,.5)
fd={'size':14,'weight':'bold'}

kk=0
#for edilst in edis[:]:
zb=Z.Z(edis[kk][0])
rb=zb.getResPhase()
pb=zb.getPhaseTensor()

zf=Z.Z(edis[kk][-1])
rf=zf.getResPhase()
pf=zf.getPhaseTensor()
kk+=1

pbi=np.where(np.round(zb.period)!=14)
pbf=np.where(np.round(zf.period)!=14)

#axr1=fig1.add_subplot(2,2,1)
#r1=axr1.errorbar(zb.period[pbi],rb.resxy[pbi],yerr=rb.resxyerr[pbi],mec=glst[kk],color=glst[kk],
#                marker='+',ms=3,ls='-')
#r2=axr1.errorbar(zf.period[pbf],rf.resxy[pbf],yerr=rf.resxyerr[pbf],mec=clstb[kk],color=clstb[kk],
#                marker='s',ms=1,ls=':',mew=.5)
#                
#axr2=fig1.add_subplot(2,2,2)
#r3=axr2.errorbar(zb.period[pbi],rb.resyx[pbi],yerr=rb.resyxerr[pbi],mec=glst[kk],color=glst[kk],
#                marker='+',ms=3,ls='-')
#r4=axr2.errorbar(zf.period[pbf],rf.resyx[pbf],yerr=rf.resyxerr[pbf],mec=clst[kk],color=clst[kk],
#                marker='o',ms=1,ls=':',mew=.5)
#                
#axp1=fig1.add_subplot(2,2,3)
#r5=axp1.errorbar(zb.period[pbi],rb.phasexy[pbi],yerr=rb.phasexyerr[pbi],mec=glst[kk],color=glst[kk],
#                marker='+',ms=3,ls='-')
#r6=axp1.errorbar(zf.period[pbf],rf.phasexy[pbf],yerr=rf.phasexyerr[pbf],mec=clstb[kk],color=clstb[kk],
#                marker='s',ms=3,ls=':')
#                
#axp2=fig1.add_subplot(2,2,4)
#r7=axp2.errorbar(zb.period[pbi],rb.phaseyx[pbi]+180,yerr=rb.phaseyxerr[pbi],mec=glst[kk],color=glst[kk],
#                marker='+',ms=3,ls='-')
#r8=axp2.errorbar(zf.period[pbf],rf.phaseyx[pbf]+180,yerr=rf.phaseyxerr[pbf],mec=clst[kk],color=clst[kk],
#                marker='o',ms=3,ls=':')
axr1=fig1.add_subplot(2,1,1)
r1=axr1.errorbar(zb.period[pbi],rb.resxy[pbi],yerr=rb.resxyerr[pbi],mec='k',
                 color='k',marker='+',ms=3,ls='-')
r2=axr1.errorbar(zf.period[pbf],rf.resxy[pbf],yerr=rf.resxyerr[pbf],
                 mec=mc,color=mc,marker='s',ms=1,ls=':',mew=.5)
                 
#axr1=fig1.add_subplot(2,1,1)
#r1=axr1.errorbar(zb.period[pbi],pb.phimaxang[pbi],yerr=pb.phimaxangvar[pbi],
#                 mec='k',color='k',marker='+',ms=3,ls='-')
#r2=axr1.errorbar(zf.period[pbf],pf.phimaxang[pbf],yerr=pb.phimaxangvar[pbf],
#                 mec=mc,color=mc,marker='s',ms=1,ls=':',mew=.5)
                
#axr2=fig1.add_subplot(2,1,2)
#r3=axr2.errorbar(zb.period[pbi],rb.resyx[pbi],yerr=rb.resyxerr[pbi],mec='k',
#                 color='k', marker='+',ms=3,ls='-')
#r4=axr2.errorbar(zf.period[pbf],rf.resyx[pbf],yerr=rf.resyxerr[pbf],
#                 mec=mc,color=mc,marker='o',ms=1,ls=':',mew=.5)
                
axp1=fig1.add_subplot(2,1,2,sharex=axr1)
r5=axp1.errorbar(zb.period[pbi],rb.phasexy[pbi],yerr=rb.phasexyerr[pbi],
                 mec='k',color='k',marker='+',ms=3,ls='-')
r6=axp1.errorbar(zf.period[pbf],rf.phasexy[pbf],yerr=rf.phasexyerr[pbf],
                 mec=mc,color=mc,marker='s',ms=.5,ls=':')
                 
#r5=axp1.errorbar(zb.period[pbi],pb.phiminang[pbi],yerr=pb.phiminangvar[pbi],
#                 mec='k',color='k',marker='+',ms=3,ls='-')
#r6=axp1.errorbar(zf.period[pbf],pf.phiminang[pbf],yerr=pb.phiminangvar[pbf],
#                 mec=mc,color=mc,marker='s',ms=1,ls=':',mew=.5)
#                
#axp2=fig1.add_subplot(2,2,4)
#r7=axp2.errorbar(zb.period[pbi],rb.phaseyx[pbi]+180,yerr=rb.phaseyxerr[pbi],
#                 mec='k',color='k',marker='+',ms=3,ls='-')
#r8=axp2.errorbar(zf.period[pbf],rf.phaseyx[pbf]+180,yerr=rf.phaseyxerr[pbf],
#                 mec=mc,color=mc,marker='o',ms=.5,ls=':')

#for jj,ax in enumerate([axr1,axr2,axp1,axp2]):
for jj,ax in enumerate([axr1,axp1]):
    ax.set_xscale('log')
    ax.set_xlim(1,30)
#    ax.set_ylim(ylimlst[jj])
    ax.tick_params(which='both')
    if jj>0:
        ax.xaxis.set_minor_locator(FixedLocator([2,3,4,5,6,7,8,9,20,30]))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        ax.xaxis.set_minor_formatter(FormatStrFormatter('%.2g'))
        ax.set_xlabel('period(s)',fontdict=fd)
    else:
        plt.setp(ax.xaxis.get_ticklabels(),visible=False)
    
    if jj==0:
        ax.set_ylabel('app. res. ($\Omega \cdot$m)',fontdict=fd)
    if jj==1:
        ax.set_ylabel('phase (deg)',fontdict=fd)
        
    ax.grid(alpha=.2,which='both')
    ax.text(1.1,ylimlst[jj][1]*.92,textlst[jj],fontdict=fd,
            bbox={'facecolor':'white'})

    llst.append(r2[0])
    
leg=fig1.legend([r1[0],r2[0]],['UT Day 192','UT Day 196'],
            loc='upper center',ncol=2,borderpad=.1,columnspacing=.2,
            prop=fd)

#for ii,edi in enumerate(edilst[1:],1):
#    zi=Z.Z(edi)
#    ri=zi.getResPhase()
#    if ii==1:    
#        axr=fig1.add_subplot(nedi-1,2,2*ii-1)
#        axp=fig1.add_subplot(nedi-1,2,2*ii,sharex=axr)
#    else:
#        axr=fig1.add_subplot(nedi-1,2,2*ii-1,sharex=axr)
#        axp=fig1.add_subplot(nedi-1,2,2*ii,sharex=axr)
#    
##    r1=axr.errorbar(zb.period,rb.resxy,yerr=rb.resxyerr,mec=glst[kk],color=glst[kk],
##                    marker='+',ms=3,ls='-')
##    r2=axr.errorbar(zi.period,ri.resxy,yerr=ri.resxyerr,mec=clstb[kk],color=clstb[kk],
##                    marker='s',ms=1,ls=':',mew=.5)
#                    
#    r1=axr.errorbar(zb.period,rb.resyx,yerr=rb.resyxerr,mec=glst[kk],color=glst[kk],
#                    marker='+',ms=3,ls='-')
#    r2=axr.errorbar(zi.period,ri.resyx,yerr=ri.resyxerr,mec=clst[kk],color=clst[kk],
#                    marker='o',ms=1,ls=':',mew=.5)
#                    
##    r3=axp.errorbar(zb.period,rb.phasexy,yerr=rb.phasexyerr,mec=glst[kk],color=glst[kk],
##                    marker='+',ms=3,ls='-')
##    r4=axp.errorbar(zi.period,ri.phasexy,yerr=ri.phasexyerr,mec=clstb[kk],color=clstb[kk],
##                    marker='s',ms=3,ls=':')
##                    
#    r3=axp.errorbar(zb.period,rb.phaseyx+180,yerr=rb.phaseyxerr,mec=glst[kk],color=glst[kk],
#                    marker='+',ms=3,ls='-')
#    r4=axp.errorbar(zi.period,ri.phaseyx+180,yerr=ri.phaseyxerr,mec=clst[kk],color=clst[kk],
#                    marker='o',ms=3,ls=':')
#    
##    axr.set_yscale('log')
#    axr.set_xscale('log')
##    axp.set_yscale('log')
#    
#    axp.set_xscale('log')
#    axr.set_xlim(2,30)
##    axr.set_ylim(8,60)
##    axp.set_ylim(2,50)
#    axp.set_ylim(3,11)
#    axp.set_ylim(26,40)
##    axp.set_ylim(5,18)
#    axr.grid(alpha=.25)
#    axp.grid(alpha=.25)
#    
#    
#    axr.text(.85,42,daylst[ii],fontdict={'size':9,'weight':'bold'},
#             horizontalalignment='left',verticalalignment='top',
#             bbox={'facecolor':'white'})
#             
#    if ii<nedi-2:
#        plt.setp(axr.get_xticklabels(),visible=False)
#        plt.setp(axp.get_xticklabels(),visible=False)
#    else:
#        axr.set_xlabel('period (s)',fontdict={'size':10,'weight':'bold'})
#        axp.set_xlabel('period (s)',fontdict={'size':10,'weight':'bold'})
#        
##    axr.yaxis.set_major_locator(MultipleLocator(3))
##    axp.yaxis.set_major_locator(MultipleLocator(3))
##    axr.xaxis.set_major_locator(MultipleLocator(1))
##    axp.xaxis.set_major_locator(MultipleLocator(3))
#    
##    axr.set_ylabel('phase (deg)',fontdict={'size':10,'weight':'bold'})
##    axp.set_ylabel('phase (deg)',fontdict={'size':10,'weight':'bold'})
#    
#
##fig1.legend([p1[0],p2[0],p3[0],p4[0]],['Pre_xy','ii_xy','Pre_yx','ii_yx'],
##            loc='upper center',ncol=4,borderpad=.01,columnspacing=.2)
#fig1.legend([r1[0],r2[0],r3[0],r4[0]],['Pre_xy','ii_xy','Pre_yx','ii_yx'],
#            loc='upper center',ncol=4,borderpad=.01,columnspacing=.2)
    
                    
    