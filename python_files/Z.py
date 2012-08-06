# -*- coding: utf-8 -*-
"""
Created on Mon May 03 13:44:51 2010

@author: a1185872
"""

import numpy as np
import MTtools as mt
import os

class Z:
    """Z is an impedance tensor class that can calculate invariants according to
    Weaver et al. 2003 and phase tensor components according to Caldwell et al.
    2004. Input is an .edi or .imp file or a dictionary with keys:
    impvar['station']
    impvar['lat']
    impvar['lon']
    impvar['periodlst']
    impvar['zlst']
    impvar['zvar']
    impvar['tip']
    impvar['tipvar']"""
    def __init__(self,impvar):
        if type(impvar) is str:
            if impvar[-3:].find('imp')>=0:
                self.filename=impvar
                self.station,self.periodlst,self.zlst=mt.readimp(impvar)
                self.lat='NaN'
                self.lon='NaN'
                self.zvar=np.zeros_like(self.zlst)
                self.tip=np.zeros((len(self.zlst),2))
                self.tipvar=np.zeros_like(self.tip)
            elif impvar[-3:].find('edi')>=0:
                self.filename=impvar
                self.station=os.path.basename(impvar)[:-4]
                impdict=mt.readedi(impvar)
                self.lat=impdict['lat']
                self.lon=impdict['lon']
                #print impdict['frequency']
                #flip things around to order from short period to long period
                if impdict['frequency'][0]<impdict['frequency'][-1]:
                    self.freqlst=impdict['frequency'][::-1]
                    self.periodlst=1/impdict['frequency'][::-1] 
                    self.zlst=impdict['z'][::-1,:,:]    
                    self.zvar=impdict['zvar'][::-1,:,:] 
                    self.tip=impdict['tipper'][::-1,:]      
                    self.tipvar=impdict['tippervar'][::-1,:] 
                    print 'Flipped arrays so frequency is decending'
                else:
                    self.freqlst=impdict['frequency']
                    self.periodlst=1/impdict['frequency']
                    self.zlst=impdict['z']
                    self.zvar=impdict['zvar']
                    self.tip=impdict['tipper']
                    self.tipvar=impdict['tippervar']
            print 'Opened impedance file: '+impvar
        else:
            try:
                impvar['station']
            except KeyError:
                impvar['station']='n/a'
            try:
                impvar['lat']
            except KeyError:
                impvar['lat']='NaN'
            try:
                impvar['lon']
            except KeyError:
                impvar['lon']='NaN'
            try:
                impvar['periodlst']
            except KeyError:
                impvar['periodlst']=range(len(impvar['zlst']))
            try:
                impvar['zvar']
            except KeyError:
                impvar['zvar']=np.zeros_like(impvar['zlst'])
            try:
                impvar['tip']
            except KeyError:
                impvar['tip']=np.zeros((len(impvar['zlst']),2))
            try:
                impvar['tipvar']
            except KeyError:
                impvar['tipvar']=np.zeros((len(impvar['zlst']),2))
  
            self.station=impvar['station']
            self.lat=impvar['lat']
            self.lon=impvar['lon']
            self.periodlst=impvar['periodlst']
            self.zlst=impvar['zlst']
            self.zvar=impvar['zvar']
            self.tip=impvar['tip']
            self.tipvar=impvar['tipvar']
            print 'Impedance tensor read.'
    
    def getInvariants(self):
        """Calculate the invariants according to Weaver et al. (2003) output is:
        a dictionary with keys:
            inv1,
            inv2,
            inv3,
            inv4,
            inv5,
            inv6,
            inv7,
            q,
            strike,
            strikeerr"""
        
        z=np.array(self.zlst)
        
        if z.size<4:
            raise ValueError('Size of impedance tensor is less than 4.')
        if z.size==4 and z.shape==(2,2):
            #print 'Computing invariants for single impedance tensor...'
            #compute invariant parameters
            x1=.5*(z[0,0].real+z[1,1].real)
            x2=.5*(z[0,1].real+z[1,0].real)
            x3=.5*(z[0,0].real-z[1,1].real)
            x4=.5*(z[0,1].real-z[1,0].real)
            e1=.5*(z[0,0].imag+z[1,1].imag)
            e2=.5*(z[0,1].imag+z[1,0].imag)
            e3=.5*(z[0,0].imag-z[1,1].imag)
            e4=.5*(z[0,1].imag-z[1,0].imag)
            ex=x1*e1-x2*e2-x3*e3+x4*e4
            d12=(x1*e2-x2*e1)/ex
            d34=(x3*e4-x4*e3)/ex
            d13=(x1*e3-x3*e1)/ex
            d24=(x2*e4-x4*e2)/ex
            d41=(x4*e1-x1*e4)/ex
            d23=(x2*e3-x3*e2)/ex
            s12=(x1*e2+x2*e1)/ex
            s34=(x3*e4+x4*e3)/ex
            s13=(x1*e3+x3*e1)/ex
            s24=(x2*e4+x4*e2)/ex
            s41=(x4*e1+x1*e4)/ex
            s23=(x2*e3+x3*e2)/ex
            inv1=np.sqrt(x4**2+x1**2)
            inv2=np.sqrt(e4**2+e1**2)
            inv3=np.sqrt(x2**2+x3**2)/inv1
            inv4=np.sqrt(e2**2+e3**2)/inv2
            inv5=s41*ex/(inv1*inv2)
            inv6=d41*ex/(inv1*inv2)
            q=np.sqrt((d12-d34)**2+(d13+d24)**2)
            inv7=(d41-d23)/q
            strikeang=.5*np.arctan2(d12-d34,d13+d24)*(180/np.pi)
            strikeangerr=abs(.5*np.arcsin(inv7))*(180/np.pi)
            invdict={}
            invdict['inv1']=inv1
            invdict['inv2']=inv2
            invdict['inv3']=inv3
            invdict['inv4']=inv4
            invdict['inv5']=inv5
            invdict['inv6']=inv6
            invdict['inv7']=inv7
            invdict['q']=q
            invdict['strike']=strikeang
            invdict['strikeerr']=strikeangerr
            
            return invdict
        elif z.size>4:
            invdict={}
            invdict['inv1']=[]
            invdict['inv2']=[]
            invdict['inv3']=[]
            invdict['inv4']=[]
            invdict['inv5']=[]
            invdict['inv6']=[]
            invdict['inv7']=[]
            invdict['q']=[]
            invdict['strike']=[]
            invdict['strikeerr']=[]
            #print 'Calculating invariants for '+str(len(z))+' impedance tensors'
            for ii in range(len(z)):
                #compute invariant parameters
                x1=.5*(z[ii,0,0].real+z[ii,1,1].real)
                x2=.5*(z[ii,0,1].real+z[ii,1,0].real)
                x3=.5*(z[ii,0,0].real-z[ii,1,1].real)
                x4=.5*(z[ii,0,1].real-z[ii,1,0].real)
                e1=.5*(z[ii,0,0].imag+z[ii,1,1].imag)
                e2=.5*(z[ii,0,1].imag+z[ii,1,0].imag)
                e3=.5*(z[ii,0,0].imag-z[ii,1,1].imag)
                e4=.5*(z[ii,0,1].imag-z[ii,1,0].imag)
                ex=x1*e1-x2*e2-x3*e3+x4*e4
                d12=(x1*e2-x2*e1)/ex
                d34=(x3*e4-x4*e3)/ex
                d13=(x1*e3-x3*e1)/ex
                d24=(x2*e4-x4*e2)/ex
                d41=(x4*e1-x1*e4)/ex
                d23=(x2*e3-x3*e2)/ex
                s12=(x1*e2+x2*e1)/ex
                s34=(x3*e4+x4*e3)/ex
                s13=(x1*e3+x3*e1)/ex
                s24=(x2*e4+x4*e2)/ex
                s41=(x4*e1+x1*e4)/ex
                s23=(x2*e3+x3*e2)/ex
                inv1=np.sqrt(x4**2+x1**2)
                inv2=np.sqrt(e4**2+e1**2)
                inv3=np.sqrt(x2**2+x3**2)/inv1
                inv4=np.sqrt(e2**2+e3**2)/inv2
                inv5=s41*ex/(inv1*inv2)
                inv6=d41*ex/(inv1*inv2)
                q=np.sqrt((d12-d34)**2+(d13+d24)**2)
                inv7=(d41-d23)/q
                strikeang=.5*np.arctan2(d12-d34,d13+d24)*(180/np.pi)
                strikeangerr=abs(.5*np.arcsin(inv7))*(180/np.pi)
                invdict['inv1'].append(inv1)
                invdict['inv2'].append(inv2)
                invdict['inv3'].append(inv3)
                invdict['inv4'].append(inv4)
                invdict['inv5'].append(inv5)
                invdict['inv6'].append(inv6)
                invdict['inv7'].append(inv7)
                invdict['q'].append(q)
                invdict['strike'].append(strikeang)
                invdict['strikeerr'].append(strikeangerr)

            return invdict
    
    def getPhaseTensor(self,rotate=180):
        """Calculate phase tensor elements following Caldwell et al. 2004. 
        Returns a dictionary of lists with keys:
            ptdict['phi']
            ptdict['phivar']
            ptdict['phimin']
            ptdict['phiminvar']
            ptdict['phimax']
            ptdict['phimaxvar']
            ptdict['alpha']
            ptdict['alphavar']
            ptdict['beta']
            ptdict['betavar']
            ptdict['azimuth']
            ptdict['azimuthvar']
            ptdict['phiminang']
            ptdict['phiminangvar']
            ptdict['phimaxang']
            ptdict['phimaxangvar']
            ptdict['ellipticity']
            ptdict['ellipticityvar']"""
        
        z=np.array(self.zlst)
                
        if z.size<4:
            raise ValueError('Shape of impedance tensor is less than 4.')
        if z.size==4 and z.shape==(2,2):
            X=z.real
            Y=z.imag
            #calculate phase tensor to align x pointing 
            #east and y pointing north Caldwell (2004)
            if rotate==90:
                phi=np.rot90(np.dot(np.linalg.inv(X),Y),1)
            elif rotate==180:
                phi=np.rot90(np.dot(np.linalg.inv(X),Y),2)
            elif rotate==270:
                phi=np.rot90(np.dot(np.linalg.inv(X),Y),3)
            else:
                phi=np.dot(np.linalg.inv(X),Y)
            #put variance into standard deviation
            zvar=self.zvar**2
            #create a matrix for errors to be calculated
            dphi=np.zeros_like(z)
            #compute determinate of X
            detX=np.linalg.det(X)
            #calculate the deteriminate of the error matrix 
            ddet=np.sqrt((X[0,0]*X[1,1])**2*((zvar[0,0]/X[0,0])**2+
                        (zvar[1,1]/X[1,1])**2)+(X[1,0]*X[0,1])**2*
                        ((zvar[0,1]/X[0,1])**2+(zvar[1,0]/X[1,0])**2))
            #calculate errors for each component of the matrix ala Caldwell 2004
            dphi[0,0]=np.sqrt((X[1,1]*Y[0,0])**2*((zvar[1,1]/X[1,1])**2+
                (zvar[0,0]/Y[0,0])**2)+(X[0,1]*Y[1,0])**2*((zvar[0,1]/X[0,1])**2
                +(zvar[1,0]/X[1,0])**2))
            dphi[0,1]=np.sqrt((X[1,1]*Y[0,1])**2*((zvar[1,1]/X[1,1])**2+
                (zvar[0,1]/Y[0,1])**2)+(X[0,1]*Y[1,1])**2*((zvar[0,1]/X[0,1])**2
                +(zvar[1,1]/X[1,1])**2))
            dphi[1,0]=np.sqrt((X[0,1]*Y[1,0])**2*((zvar[0,0]/X[0,0])**2+
                (zvar[1,0]/Y[1,0])**2)+(X[1,0]*Y[0,0])**2*((zvar[1,0]/X[1,0])**2
                +(zvar[0,0]/X[0,0])**2))
            dphi[1,1]=np.sqrt((X[0,0]*Y[1,1])**2*((zvar[0,0]/X[0,0])**2+
                (zvar[1,1]/Y[1,1])**2)+(X[1,0]*Y[0,1])**2*((zvar[1,0]/X[1,0])**2
                +(zvar[0,1]/X[0,1])**2))
            #rotate the error matrix
            dphi=np.rot90(dphi)
            #finish calculating the errors
            dphi[0,0]=(phi[0,0]/detX)**2*np.sqrt((dphi[0,0]/phi[0,0])**2+
                                                                (ddet/detX)**2)
            dphi[0,1]=(phi[0,1]/detX)**2*np.sqrt((dphi[0,1]/phi[0,1])**2+
                                                                (ddet/detX)**2)
            dphi[1,0]=(phi[1,0]/detX)**2*np.sqrt((dphi[1,0]/phi[1,0])**2+
                                                                (ddet/detX)**2)
            dphi[1,1]=(phi[1,1]/detX)**2*np.sqrt((dphi[1,1]/phi[1,1])**2+
                                                                (ddet/detX)**2)
            
            #Calculate Trace of Phi and error of trace of phi
            tr=phi[0,0]+phi[1,1]
            trvar=np.sqrt(dphi[0,0]**2+dphi[1,1]**2)
            
            #Calculate skew of phi and the cooresponding error
            skew=phi[0,1]-phi[1,0]
            skewvar=np.sqrt(dphi[0,1]**2+dphi[1,1]**2)
            
            #calculate the determinate and determinate error of phi
            phidet=abs(np.linalg.det(phi))
            pd1=np.sqrt((dphi[0,0]/phi[0,0])**2+
                                (dphi[1,1]/phi[1,1])**2)*phi[0,0]*phi[1,1]
            pd2=np.sqrt((dphi[0,1]/phi[0,1])**2+
                                (dphi[1,0]/phi[1,0])**2)*phi[0,1]*phi[1,0]
            phidetvar=np.sqrt(pd1**2+pd2**2)
            
            #calculate reverse trace and error
            revtr=phi[0,0]-phi[1,1]
            revtrvar=np.sqrt(dphi[0,0]**2+dphi[1,1]**2)
            
            #calculate reverse skew and error
            revskew=phi[1,0]+phi[0,1]
            revskewvar=np.sqrt(phi[0,1]**2+dphi[1,0]**2)
            
            #calculate skew angle beta and error
            beta=.5*np.arctan2(skew,tr)*180/np.pi
            betavar=abs(np.arctan(skew*tr*np.sqrt((skewvar/skew)**2+
                                                    (trvar/tr)**2))*180/np.pi)
            
            #calculate angle alpha corresponding to phase tensor's dependence on 
            #coordinate system
            alpha=.5*np.arctan2(revskew,revtr)*180/np.pi
            alphavar=abs(.5*np.arctan(revskew*revtr*np.sqrt(
                        (revskewvar/revskew)**2+(revtrvar/revtr)**2))*180/np.pi)
            
            #calculate azimuth as angle between major axis and x-axis
            azimuth=alpha-beta
            azimuthvar=np.sqrt(alphavar**2+betavar**2)
            
            #calulate maximum value for phi
            phimax=np.sqrt((.5*tr)**2+(.5*skew)**2)+np.sqrt((.5*tr)**2+
                                            (.5*skew)**2-(np.sqrt(phidet))**2)
            phimaxvar=.5*np.sqrt(2*trvar**2+2*skewvar**2)+.5*np.sqrt(2*trvar**2+
                                                        2*skewvar**2+phidetvar)
            
            #calculate minimum value for phi
            if np.linalg.det(phi)>=0:
                phimin=np.sqrt((.5*tr)**2+(.5*skew)**2)-np.sqrt((.5*tr)**2+
                                            (.5*skew)**2-(np.sqrt(phidet))**2)
            elif np.linalg.det(phi)<0:
                phimin=-1*np.sqrt((.5*tr)**2+(.5*skew)**2)-np.sqrt((.5*tr)**2+
                                            (.5*skew)**2-(np.sqrt(phidet))**2)
            
            phiminvar=phimaxvar
            
            #calculate ellipticity
            phiminang=(180/np.pi)*np.arctan(np.array(phimin))
            phiminangvar=(180/np.pi)*np.arctan(np.array(phiminvar))
            phimaxang=(180/np.pi)*np.arctan(np.array(phimax))
            phimaxangvar=(180/np.pi)*np.arctan(np.array(phimaxvar))
            
            ellipticity=(phimaxang-phiminang)/(phimaxang+phiminang)
            ellipticityvar=ellipticity*np.sqrt(phimaxangvar**2+
                        phiminangvar**2)*np.sqrt((1/(phimaxang-phiminang))**2+
                        (1/(phimaxang+phiminang))**2)
                
            
            #put things into a dictionary for easy mining
            ptdict={}
            ptdict['phi']=np.array(phi)
            ptdict['phivar']=np.array(dphi)
            ptdict['phimin']=np.array(phimin)
            ptdict['phiminvar']=np.array(phiminvar)
            ptdict['phimax']=np.array(phimax)
            ptdict['phimaxvar']=np.array(phimaxvar)
            ptdict['alpha']=np.array(alpha)
            ptdict['alphavar']=np.array(alphavar)
            ptdict['beta']=np.array(beta)
            ptdict['betavar']=np.array(betavar)
            ptdict['azimuth']=np.array(azimuth)
            ptdict['azimuthvar']=np.array(azimuthvar)
            ptdict['phiminang']=np.array(phiminang)
            ptdict['phiminangvar']=np.array(phiminangvar)
            ptdict['phimaxang']=np.array(phimaxang)
            ptdict['phimaxangvar']=np.array(phimaxangvar)
            ptdict['ellipticity']=np.array(ellipticity)
            ptdict['ellipticityvar']=np.array(ellipticityvar)
            
            return ptdict
        
        elif z.size>4:
            #create a dictionary of lists to put things into and keep easier 
            #track of things
            
#            ptdict={}
#            self.phi=[]
#            self.phivar=[]
#            self.phimin=[]
#            self.phiminvar=[]
#            self.phimax=[]
#            self.phimaxvar=[]
#            self.alpha=[]
#            self.alphavar=[]
#            self.beta=[]
#            self.betavar=[]
#            self.azimuth=[]
#            self.azimuthvar=[]
#            self.phiminang=[]
#            self.phiminangvar=[]
#            self.phimaxang=[]
#            self.phimaxangvar=[]
#            self.ellipticity=[]
#            self.ellipticityvar=[]
            ptdict={}
            ptdict['phi']=[]
            ptdict['phivar']=[]
            ptdict['phimin']=[]
            ptdict['phiminvar']=[]
            ptdict['phimax']=[]
            ptdict['phimaxvar']=[]
            ptdict['alpha']=[]
            ptdict['alphavar']=[]
            ptdict['beta']=[]
            ptdict['betavar']=[]
            ptdict['azimuth']=[]
            ptdict['azimuthvar']=[]
            ptdict['phiminang']=[]
            ptdict['phiminangvar']=[]
            ptdict['phimaxang']=[]
            ptdict['phimaxangvar']=[]
            ptdict['ellipticity']=[]
            ptdict['ellipticityvar']=[]
            
            for ii in range(len(z)):
                X=z[ii].real
                Y=z[ii].imag
                #calculate phase tensor and rotate  to align x pointing 
                #east and y pointing north Caldwell (2004)
                if rotate==90:
                    phi=np.rot90(np.dot(np.linalg.inv(X),Y),1)
                elif rotate==180:
                    phi=np.rot90(np.dot(np.linalg.inv(X),Y),2)
                elif rotate==270:
                    phi=np.rot90(np.dot(np.linalg.inv(X),Y),3)
                else:
                    phi=np.dot(np.linalg.inv(X),Y)
                #put variance into standard deviation
                zvar=self.zvar[ii]**2
                #create a matrix for errors to be calculated
                dphi=np.zeros(np.shape(z[ii]))
                #compute determinate of X
                detX=np.linalg.det(X)
                #calculate the deteriminate of the error matrix 
                ddet=np.sqrt((X[0,0]*X[1,1])**2*((zvar[0,0]/X[0,0])**2+
                    (zvar[1,1]/X[1,1])**2)+(X[1,0]*X[0,1])**2*(
                    (zvar[0,1]/X[0,1])**2+(zvar[1,0]/X[1,0])**2))
                #calculate errors for each component of the matrix ala Caldwell 
                #2004
                dphi[0,0]=np.sqrt((X[1,1]*Y[0,0])**2*((zvar[1,1]/X[1,1])**2+
                    (zvar[0,0]/Y[0,0])**2)+(X[0,1]*Y[1,0])**2*(
                    (zvar[0,1]/X[0,1])**2+(zvar[1,0]/X[1,0])**2))
                dphi[0,1]=np.sqrt((X[1,1]*Y[0,1])**2*((zvar[1,1]/X[1,1])**2+
                    (zvar[0,1]/Y[0,1])**2)+(X[0,1]*Y[1,1])**2*(
                    (zvar[0,1]/X[0,1])**2+(zvar[1,1]/X[1,1])**2))
                dphi[1,0]=np.sqrt((X[0,1]*Y[1,0])**2*((zvar[0,0]/X[0,0])**2+
                    (zvar[1,0]/Y[1,0])**2)+(X[1,0]*Y[0,0])**2*(
                    (zvar[1,0]/X[1,0])**2+(zvar[0,0]/X[0,0])**2))
                dphi[1,1]=np.sqrt((X[0,0]*Y[1,1])**2*((zvar[0,0]/X[0,0])**2+
                    (zvar[1,1]/Y[1,1])**2)+(X[1,0]*Y[0,1])**2*(
                    (zvar[1,0]/X[1,0])**2+(zvar[0,1]/X[0,1])**2))
                #rotate the error matrix
                dphi=np.rot90(dphi,2)
                #finish calculating the errors
                dphi[0,0]=(phi[0,0]/detX)**2*np.sqrt((dphi[0,0]/phi[0,0])**2+
                                                                (ddet/detX)**2)
                dphi[0,1]=(phi[0,1]/detX)**2*np.sqrt((dphi[0,1]/phi[0,1])**2+
                                                                (ddet/detX)**2)
                dphi[1,0]=(phi[1,0]/detX)**2*np.sqrt((dphi[1,0]/phi[1,0])**2+
                                                                (ddet/detX)**2)
                dphi[1,1]=(phi[1,1]/detX)**2*np.sqrt((dphi[1,1]/phi[1,1])**2+
                                                                (ddet/detX)**2)
                
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
                    np.sqrt((.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2)
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
                
                #calculate ellipticity
                phiminang=(180/np.pi)*np.arctan(np.array(phimin))
                phiminangvar=(180/np.pi)*np.arctan(np.array(phiminvar))
                phimaxang=(180/np.pi)*np.arctan(np.array(phimax))
                phimaxangvar=(180/np.pi)*np.arctan(np.array(phimaxvar))
                
                ellipticity=(phimaxang-phiminang)/(phimaxang+phiminang)
                ellipticityvar=ellipticity*np.sqrt(phimaxangvar**2+
                        phiminangvar**2)*np.sqrt((1/(phimaxang-phiminang))**2+
                        (1/(phimaxang+phiminang))**2)
                
                
                #append things to a dictionary of lists
                
#                self.phi.append(phi)
#                self.phivar.append(np.array(dphi))
#                self.phimin.append(float(phimin))
#                self.phiminvar.append(float(phiminvar))
#                self.phimax.append(float(phimax))
#                self.phimaxvar.append(float(phimaxvar))
#                self.alpha.append(float(alpha))
#                self.alphavar.append(float(alphavar))
#                self.beta.append(float(beta))
#                self.betavar.append(float(betavar))
#                self.azimuth.append(float(azimuth))
#                self.azimuthvar.append(float(azimuthvar))
#                self.phiminang.append(float(phiminang))
#                self.phiminangvar.append(float(phiminangvar))
#                self.phimaxang.append(float(phimaxang))
#                self.phimaxangvar.append(float(phimaxangvar))
#                self.ellipticity.append(float(ellipticity))
#                self.ellipticityvar.append(float(ellipticityvar))
#            
#            self.phi=np.array(self.phi)
#            self.phivar=np.array(self.phivar)
                ptdict['phi'].append(phi)
                ptdict['phivar'].append(np.array(dphi))
                ptdict['phimin'].append(float(phimin))
                ptdict['phiminvar'].append(float(phiminvar))
                ptdict['phimax'].append(float(phimax))
                ptdict['phimaxvar'].append(float(phimaxvar))
                ptdict['alpha'].append(float(alpha))
                ptdict['alphavar'].append(float(alphavar))
                ptdict['beta'].append(float(beta))
                ptdict['betavar'].append(float(betavar))
                ptdict['azimuth'].append(float(azimuth))
                ptdict['azimuthvar'].append(float(azimuthvar))
                ptdict['phiminang'].append(float(phiminang))
                ptdict['phiminangvar'].append(float(phiminangvar))
                ptdict['phimaxang'].append(float(phimaxang))
                ptdict['phimaxangvar'].append(float(phimaxangvar))
                ptdict['ellipticity'].append(float(ellipticity))
                ptdict['ellipticityvar'].append(float(ellipticityvar))
            
            ptdict['phi']=np.array(ptdict['phi'])
            ptdict['phivar']=np.array(ptdict['phivar'])
            return ptdict

    def getTipper(self):
        """Get tipper information and return a dictionary of lists with keys:
        tipdict['magreal']
        tipdict['magimag']
        tipdict['anglereal']
        tipdict['angleimag']"""
        
        tip=np.array(self.tip)
        tdict={}
        if tip.size==2:
            tdict['magreal']=np.sqrt(tip[0].real**2+tip[1].real**2)
            tdict['magimag']=np.sqrt(tip[0].imag**2+tip[1].imag**2)
            tdict['anglereal']=np.arctan2(-tip[1].real,-tip[0].real)*180/np.pi
            tdict['angleimag']=np.arctan2(-tip[1].imag,-tip[0].imag)*180/np.pi
            
            return tdict
        
        else:
            tdict['magreal']=[]
            tdict['magimag']=[]
            tdict['anglereal']=[]
            tdict['angleimag']=[]
            for ii in range(len(tip)):
                tdict['magreal'].append(np.sqrt(tip[ii,0].real**2+
                                                            tip[ii,1].real**2))
                tdict['magimag'].append(np.sqrt(tip[ii,0].imag**2+
                                                            tip[ii,1].imag**2))
                tdict['anglereal'].append(np.arctan2(-tip[ii,1].real,
                                                    -tip[ii,0].real)*180/np.pi)
                tdict['angleimag'].append(np.arctan2(-tip[ii,1].imag,
                                                    -tip[ii,0].imag)*180/np.pi)
            
            return tdict
   
    def removeDistortion(self):
        """removeDistortion(self) will remove the distortion from the impedance
        tensor as prescribed by Bibby et al. [2005].  Saves the edifile as 
        station+dr.edi. Returns saved filename with full path."""
        
        zlst=np.array(self.zlst)
        zvar=np.array(self.zvar)
        #calculate ellipticity to figure out the distortion from 1d structure
        ptdict=Z.getPhaseTensor(self)
        zd=[]
        #rotated identity matrix for static shift removal
        im=np.array([[0,-1],[1,0]])
        for ii,ellip in enumerate(ptdict['ellipticity']):
            if ellip<.1:
                X=zlst[ii].real
                Y=zlst[ii].imag
                #calculate static shift factor
                gx=np.sqrt(abs(np.linalg.det(X)))
                gy=np.sqrt(abs(np.linalg.det(Y)))
                #remove static shift from real and imaginary parts
                #need to use the dot for matric multiplication
                Xs=np.dot(X,im)/gx
                Ys=np.dot(Y,im)/gy
                #append real and imaginary parts sequentially
                zd.append(Xs)
                zd.append(Ys)
            else:
                pass
        if len(zd)==0:
            print 'There is no 1D structure for this station'
        else:
            #estimate distortion matrix 
            zd=np.array(zd)
            D=np.array([[zd[:,0,0].mean(),zd[:,0,1].mean()],
                         [zd[:,1,0].mean(),zd[:,1,1].mean()]])
            Dvar=np.array([[zd[:,0,0].std(),zd[:,0,1].std()],
                            [zd[:,1,0].std(),zd[:,1,1].std()]])
            Dinv=np.linalg.inv(D)
            #remove distortion as (D^-1)Z[ii]
            zdr=np.array([np.dot(np.linalg.inv(D),zlst[ii]) for ii in range(len(zlst))])
            
            #estimate errors
            zvardr=np.zeros_like(zlst)
            for jj in range(len(zlst)):
                X=zlst[jj].real
                Xvar=zvar[jj].real
                zvardr[jj,0,0]=np.sqrt((Dinv[0,0]*X[0,0]*
                            np.sqrt((Dvar[1,1]/Dinv[0,0])**2+
                            (Xvar[0,0]/X[0,0])**2))**2+(Dinv[0,1]*X[1,0]
                            *np.sqrt((Dvar[0,1]/Dinv[0,1])**2+
                            (Xvar[0,1]/X[0,1])**2))**2)
                zvardr[jj,0,1]=np.sqrt((Dinv[0,0]*X[0,1]*
                            np.sqrt((Dvar[1,1]/Dinv[0,0])**2+
                            (Xvar[0,1]/X[0,1])**2))**2+(Dinv[0,1]*X[1,1]
                            *np.sqrt((Dvar[0,1]/Dinv[0,1])**2+
                            (Xvar[1,1]/X[1,1])**2))**2)
                zvardr[jj,1,0]=np.sqrt((Dinv[1,0]*X[0,0]*
                            np.sqrt((Dvar[1,0]/Dinv[1,0])**2+
                            (Xvar[0,0]/X[0,0])**2))**2+(Dinv[1,1]*X[1,0]
                            *np.sqrt((Dvar[0,0]/Dinv[1,1])**2+
                            (Xvar[1,0]/X[1,0])**2))**2)
                zvardr[jj,1,1]=np.sqrt((Dinv[1,0]*X[1,0]*
                            np.sqrt((Dvar[1,0]/Dinv[1,0])**2+
                            (Xvar[0,1]/X[0,1])**2))**2+(Dinv[1,1]*X[1,1]
                            *np.sqrt((Dvar[1,1]/Dinv[1,1])**2+
                            (Xvar[1,1]/X[1,1])**2))**2)
            #make new edi file. need to flip zdr and zvardr back to normal order 
            #for edi file.
            newedifn=mt.rewriteedi(self.filename,zdr,
                                   zvardr.real)
            
            return newedifn
  
 
    def getResPhase(self,ffactor=1):
        """imp2resphase(z,zvar,freq,df=100.0) will convert impedances 
        to resistivities (ohm-m) and phases (deg) as well as the 
        errors of each.  The output is dictioinary: """
    
        z=self.zlst
        zvar=self.zvar
        freq=1/self.periodlst
        
        resxy=[]
        resxyerr=[]
        resyx=[]
        resyxerr=[]
        resxx=[]
        resxxerr=[]
        resyy=[]
        resyyerr=[]
        
        phasexy=[]
        phasexyerr=[]
        phaseyx=[]
        phaseyxerr=[]
        phasexx=[]
        phasexxerr=[]
        phaseyy=[]
        phaseyyerr=[]
        for jj in range(len(z)):
            wt=.2/(freq[jj])*ffactor
            resxx.append(wt*abs(z[jj,0,0])**2)
            resxy.append(wt*abs(z[jj,0,1])**2)
            resyx.append(wt*abs(z[jj,1,0])**2)
            resyy.append(wt*abs(z[jj,1,1])**2)
            
            resxxerr.append(wt*(abs(z[jj,0,0])+zvar[jj,0,0])**2
                                                                -resxx[jj])
            resxyerr.append(wt*(abs(z[jj,0,1])+zvar[jj,0,1])**2
                                                                -resxy[jj])
            resyxerr.append(wt*(abs(z[jj,1,0])+zvar[jj,1,0])**2
                                                                -resyx[jj])
            resyyerr.append(wt*(abs(z[jj,1,1])+zvar[jj,1,1])**2
                                                                -resyy[jj])
            
            phasexx.append(np.arctan(z[jj,0,0].imag/z[jj,0,0].real)\
                                                                *(180/np.pi))
            phasexy.append(np.arctan(z[jj,0,1].imag/z[jj,0,1].real)\
                                                                *(180/np.pi))
            phaseyx.append(np.arctan(z[jj,1,0].imag/z[jj,1,0].real)\
                                                                *(180/np.pi))
            phaseyy.append(np.arctan(z[jj,1,1].imag/z[jj,1,1].real)\
                                                                *(180/np.pi))
            
            phasexxerr.append(np.arcsin(zvar[jj,0,0]/abs(z[jj,0,0]))\
                                                                *(180/np.pi))
            phasexyerr.append(np.arcsin(zvar[jj,0,1]/abs(z[jj,0,1]))\
                                                                *(180/np.pi))
            phaseyxerr.append(np.arcsin(zvar[jj,1,0]/abs(z[jj,1,0]))\
                                                                *(180/np.pi))
            phaseyyerr.append(np.arcsin(zvar[jj,1,1]/abs(z[jj,1,1]))\
                                                                *(180/np.pi))
        rpdict={}
        rpdict['resxx']=resxx
        rpdict['resxxerr']=resxxerr
        rpdict['resxy']=resxy
        rpdict['resxyerr']=resxyerr
        rpdict['resyx']=resyx
        rpdict['resyxerr']=resyxerr
        rpdict['resyy']=resyy
        rpdict['resyyerr']=resyyerr
        rpdict['phasexx']=phasexx
        rpdict['phasexxerr']=phasexxerr
        rpdict['phasexy']=phasexy
        rpdict['phasexyerr']=phasexyerr
        rpdict['phaseyx']=phaseyx
        rpdict['phaseyxerr']=phaseyxerr
        rpdict['phaseyy']=phaseyy
        rpdict['phaseyyerr']=phaseyyerr
        
        return rpdict



            