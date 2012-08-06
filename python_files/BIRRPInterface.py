"""
An interface to help with BIRPP (Chave and Thompson (2004)).  Need to have 
Notepad++ installed into program files.

Author: Jared Peacock 2010

Updated Feb 15,2011

    1)can input the default path at the beginning
    2)can change survey name at the beginning

Updated Feb 28, 2011
    1)input a scaling factor when computing apparent resistivity and phases
    2)plots all components of impedance tensor
    
"""


#===============================================================================
# Import necessary packages
#===============================================================================

import os

os.chdir(r'c:\Peacock\Python\MTPro\src')

import easygui
import subprocess
import shutil
import BIRRPTools as brp
import MTPlotTools as mtplot
import MTtools as mt
import datetime
import time
import numpy as np

#set default path
defaultpath=r'H:\Paralana\BaseMarch2011'

#survey Name to go into header file
surveyname='Paralana'

#notepad++ path
notepadpath=r"c:\Program Files\Notepad++\notepad++.exe"

#if you find your resistivities are off by a factor change the ffactor to that
#else just leave it at 1
ffactor=1

#===============================================================================
# Find out what station and load .ini file if one exists if not load default
#===============================================================================
cwd=os.getcwd()

stationdir=easygui.diropenbox(msg='Locate the folder you are about to process', 
                            title='Station Find',default=defaultpath)
stationdir=str(stationdir)
stationi=os.path.basename(stationdir)
#===============================================================================
# Get Station, Header and Orientation information
#===============================================================================
#if there is a file with header information pass if not make one from info file
#the .ini,.hdr,.ori files need to reside in the station folder
#check to see if there is already an .ini file. If there is read it in
todaysdate=datetime.date.today()
if os.path.isfile(os.path.join(stationdir,stationi+'.ini'))==True:
    inidict=brp.readini(os.path.join(stationdir,stationi+'.ini'))
    magtype=inidict['magtype']  #magnetic type
    mcomps=inidict['mcomps']    #measured componets: 4-bx,by, 5-bx,by,bz,ex,ey
    magdec=inidict['magdec']    #magnetic declintation
    df=inidict['df']            #sampling frequency (Hz)
    lpyn=inidict['lpyn']        #long period conversion to mV/nTy or n
    eyn=inidict['eyn']          #electric conversion to mV/m y or n
    dlgain=inidict['dlgain']    #data logger gain
    egain=inidict['egain']      #interface box gain for electrics
    dlength=inidict['dlength']  #dipole lengths ex,ey (m)
    lpbzcor=inidict['lpbzcor']  #long period Bz conversion
    bbfile=inidict['bbfile']    #broadband calibration file
    magori=inidict['magori'].upper()    #magnetic channel orientation
    #check to see if there is a .hdr and .ori file. If there is do nothing, if
    #there isn't get info to make one.
    if os.path.isfile(os.path.join(stationdir,stationi+'.hdr'))==True and \
            os.path.isfile(os.path.join(stationdir,stationi+'.ori'))==True:
        stationinfofile=inidict['stationinfofile']
    else:
        #ask if there is a file with the station information in it, if there is
        #read in the info, if there isn't manually enter it
        stationinfofile=easygui.fileopenbox(msg='Choose file',
                                            title='Station Info File',
                                            default=stationdir,
                                            filetypes=['*.txt','*.csv'])
        if stationinfofile!=None:
            statinfodict=mt.getStationInfo(stationinfofile,stationi)
            try:
				df=statinfodict['df']
            except KeyError:
				statinfodict['df']=df
            try:
				magtype=statinfodict['magtype']
            except KeyError:
				statinfodict['magtype']=magtype
            try:
				magdec=statinfodict['magdec']
            except KeyError:
				statinfodict['magdec']=magdec
			
            hdrfields=['Data ID',
                       'Acquired By',
                        'Acquired Date (DD/MM/YY)',
                        'File Date (DD/MM/YY)',
                        'Prospect Location',
                        'Station Name',
                        'Latitude',
                        'Longitude',
                        'Elevation (m)',
                        'HX sensor position (x,y,angle from geographic North)',
                        'HX sensor position (x,y,angle from geographic North)',
                        'EX electrode location (x1,y1,x2,y2)',
                        'EY electrode location (x1,y1,x2,y2)',
                        'Remote Reference HX (x,y,angle from geographic North)',
                        'Remote Reference HX (x,y,angle from geographic North)',
                        'Magnetic Measurement type lp=long period, bb=broadband',
                        'Measured components',
                        'Magnetic declination',
                        'Sampling Frequency (Hz)',
                        'Data logger gain (very low=2.5,low=1.0,high=0.1)',
                        'Interface Box gain',
                        'Dipole lengths (m) Ex,Ey',
                        'Longperiod Bz correction',
                        'Broadband calibration file location (full path)',
                        'Measured magnetic orientation (BX,BY,BZ)',
                        'Electric fields converted to microV/m (y or n)',
                        'Longperiod magnetics converted to microV/m (y or n)']
            hdrvalues=[statinfodict['station'],
                        'Adelaide University',
                        statinfodict['date'],
                        todaysdate.strftime("%d/%m/%y"),
                        surveyname,
                        statinfodict['station'],
                        statinfodict['lat'],
                        statinfodict['long'],
                        statinfodict['elev'],
                        '0 0 0',
                        '0 0 90',
                        '0 0 '+statinfodict['ex']+' 0',
                        '0 0 0 '+statinfodict['ey'],
                        '0 0 0',
                        '0 0 90',
                        statinfodict['magtype'],
                        'EX,EY,'+statinfodict['magori'].upper(),
                        statinfodict['magdec'],
                        statinfodict['df'],
                        statinfodict['dlgain'],
                        statinfodict['egain'],
                        statinfodict['ex']+','+statinfodict['ey'],
                        statinfodict['lpbzcor'],
                        bbfile,
                        statinfodict['magori'].upper(),
                        eyn,
                        lpyn]
            
            dataid,acqby,acqd,filedate,prospectloc,stationname,lat,long,elev,\
            hxloc,hyloc,exloc,eyloc,rxloc,ryloc,magtype,mcomps,magdec,df,\
            dlgain,egain,dlength,lpbzcor,bbfile,magori,eyn,\
            lpyn=easygui.multenterbox(msg='Enter station information for the header file and orientation file',
                                      title='Station, Header and Orientation File Information',
                                      fields=hdrfields,
                                      values=hdrvalues)
            #write header file to station folder as stationname.hdr 
            hdrfid=open(os.path.join(stationdir,stationname+'.hdr'),'w')
            hdrfid.write('DATAID="'+dataid+'"'+'\n')
            hdrfid.write('ACQBY="'+acqby+'"'+'\n')
            hdrfid.write('ACQDATE='+acqd+'\n')
            hdrfid.write('FILEDATE='+filedate+'\n')
            hdrfid.write('PROSPECT="'+prospectloc+'"'+'\n')
            hdrfid.write('LOC="'+stationname+'"'+'\n')
            hdrfid.write('LAT='+lat+'\n')
            hdrfid.write('LONG='+long+'\n')
            hdrfid.write('ELEV='+elev)
            hdrfid.close()
            #write orientation file to station folder as stationname.ori
            orifid=open(os.path.join(stationdir,stationname+'.ori'),'w')
            orifid.write(hxloc+'\n')
            orifid.write(hyloc+'\n')
            orifid.write(exloc+'\n')
            orifid.write(eyloc+'\n')
            orifid.write(rxloc+'\n')
            orifid.write(ryloc)
            orifid.close()
            
        #if there is no station info file   
        else:
            station=stationi
            stationhdr=stationi
            latstr='-00:00:00.0'
            lonstr='00:00:00.0'
            elev='000'
            startdate=todaysdate.strftime("%d/%m/%y")
                
            hdrfields=['Data ID',
                       'Acquired By',
                        'Acquired Date (DD/MM/YY)',
                        'File Date (DD/MM/YY)',
                        'Prospect Location',
                        'Station Name',
                        'Latitude (HH:MM:SS.AS)',
                        'Longitude ((HH:MM:SS.AS)',
                        'Elevation (m)',
                        'HX sensor position (x,y,angle from geographic North)',
                        'HX sensor position (x,y,angle from geographic North)',
                        'EX electrode location (x1,y1,x2,y2)',
                        'EY electrode location (x1,y1,x2,y2)',
                        'Remote Reference HX (x,y,angle from geographic North)',
                        'Remote Reference HX (x,y,angle from geographic North)',
                        'Magnetic Measurement type lp=long period, bb=broadband',
                        'Measured components',
                        'Magnetic declination',
                        'Sampling Frequency (Hz)',
                        'Data logger gain (very low=2.5,low=1.0,high=0.1)',
                        'Interface Box gain',
                        'Dipole lengths (m) Ex,Ey',
                        'Longperiod Bz correction',
                        'Broadband calibration file location (full path)',
                        'Measured magnetic orientation (BX,BY,BZ)',
                        'Electric fields converted to microV/m (y or n)',
                        'Longperiod magnetics converted to microV/m (y or n)']
            hdrvalues=[station,
                        'Adelaide University',
                        startdate,
                        startdate,
                        surveyname,
                        station,
                        latstr,
                        lonstr,
                        elev,
                        '0 0 0',
                        '0 0 90',
                        '0 0 100 0',
                        '0 0 0 100',
                        '0 0 0',
                        '0 0 90',
                        magtype,
                        'EX,EY,BX,BY',
                        magdec,
                        df,
                        dlgain,
                        egain,
                        dlength,
                        lpbzcor,
                        bbfile,
                        magori,
                        'y',
                        'y']
            
            dataid,acqby,acqd,filedate,prospectloc,stationname,lat,long,elev,\
            hxloc,hyloc,exloc,eyloc,rxloc,ryloc,magtype,mcomps,magdec,df,\
            dlgain,egain,dlength,lpbzcor,bbfile,magori,eyn,\
            lpyn=easygui.multenterbox(msg='Enter information for the header file and orientation file',
                                      title='Station, Header and Orientation File Information',
                                      fields=hdrfields,
                                      values=hdrvalues)
            #write header file to station folder as stationname.hdr
            hdrfid=open(os.path.join(stationdir,stationname+'.hdr'),'w')
            hdrfid.write('DATAID="'+dataid+'"'+'\n')
            hdrfid.write('ACQBY="'+acqby+'"'+'\n')
            hdrfid.write('ACQDATE='+acqd+'\n')
            hdrfid.write('FILEDATE='+filedate+'\n')
            hdrfid.write('PROSPECT="'+prospectloc+'"'+'\n')
            hdrfid.write('LOC="'+stationname+'"'+'\n')
            hdrfid.write('LAT='+lat+'\n')
            hdrfid.write('LONG='+long+'\n')
            hdrfid.write('ELEV='+elev)
            hdrfid.close()
            #write orientation file to station folder as stationname.ori
            orifid=open(os.path.join(stationdir,stationname+'.ori'),'w')
            orifid.write(hxloc+'\n')
            orifid.write(hyloc+'\n')
            orifid.write(exloc+'\n')
            orifid.write(eyloc+'\n')
            orifid.write(rxloc+'\n')
            orifid.write(ryloc)
            orifid.close()

#if there is no .ini file found read in a default one and enter appropriate values    
else:
    inidict=brp.readini(os.path.join(cwd,'BIRRPinterface.ini'))
    magtype=inidict['magtype']
    mcomps=inidict['mcomps']
    magdec=inidict['magdec']
    df=inidict['df']
    lpyn=inidict['lpyn']
    eyn=inidict['eyn']
    dlgain=inidict['dlgain']
    egain=inidict['egain']
    dlength=inidict['dlength']
    lpbzcor=inidict['lpbzcor']
    bbfile=inidict['bbfile']
    magori=inidict['magori']
    #if there is .hdr an .ori file pass, else make each
    if os.path.isfile(os.path.join(stationdir,stationi+'.hdr'))==True and \
            os.path.isfile(os.path.join(stationdir,stationi+'.ori'))==True:
        pass
    else:
        #ask if there is a station file with info in it
        stationinfofile=easygui.fileopenbox(msg='Choose file',
                                            title='Station Info File',
                                            default=stationdir,
                                            filetypes=['*.txt','*.csv'])
        if stationinfofile!=None:
            statinfodict=mt.getStationInfo(stationinfofile,stationi)
            try:
				df=statinfodict['df']
            except KeyError:
				statinfodict['df']=df
            try:
				magtype=statinfodict['magtype']
            except KeyError:
				statinfodict['magtype']=magtype
            try:
				magdec=statinfodict['magdec']
            except KeyError:
				statinfodict['magdec']=magdec
            
            hdrfields=['Data ID',
                       'Acquired By',
                        'Acquired Date (DD/MM/YY)',
                        'File Date (DD/MM/YY)',
                        'Prospect Location',
                        'Station Name',
                        'Latitude (HH:MM:SS.AS)',
                        'Longitude ((HH:MM:SS.AS)',
                        'Elevation (m)',
                        'HX sensor position (x,y,angle from geographic North)',
                        'HX sensor position (x,y,angle from geographic North)',
                        'EX electrode location (x1,y1,x2,y2)',
                        'EY electrode location (x1,y1,x2,y2)',
                        'Remote Reference HX (x,y,angle from geographic North)',
                        'Remote Reference HX (x,y,angle from geographic North)',
                        'Magnetic Measurement type lp=long period, bb=broadband',
                        'Measured components',
                        'Magnetic declination',
                        'Sampling Frequency (Hz)',
                        'Data logger gain (very low=2.5,low=1.0,high=0.1)',
                        'Interface Box gain',
                        'Dipole lengths (m) Ex,Ey',
                        'Longperiod Bz correction',
                        'Broadband calibration file location (full path)',
                        'Measured magnetic orientation (BX,BY,BZ)',
                        'Electric fields converted to microV/m (y or n)',
                        'Longperiod magnetics converted to microV/m (y or n)']
            hdrvalues=[statinfodict['station'],
                        'Adelaide University',
                        statinfodict['date'],
                        todaysdate.strftime("%d/%m/%y"),
                        surveyname,
                        statinfodict['station'],
                        statinfodict['lat'],
                        statinfodict['long'],
                        statinfodict['elev'],
                        '0 0 0',
                        '0 0 90',
                        '0 0 '+statinfodict['ex']+' 0',
                        '0 0 0 '+statinfodict['ey'],
                        '0 0 0',
                        '0 0 90',
                        statinfodict['magtype'],
                        'EX,EY,'+statinfodict['magori'].upper(),
                        magdec,
                        df,
                        statinfodict['dlgain'],
                        statinfodict['egain'],
                        statinfodict['ex']+','+statinfodict['ey'],
                        statinfodict['lpbzcor'],
                        bbfile,
                        statinfodict['magori'].upper(),
                        eyn,
                        lpyn]
            
            dataid,acqby,acqd,filedate,prospectloc,stationname,lat,long,elev,\
            hxloc,hyloc,exloc,eyloc,rxloc,ryloc,magtype,mcomps,magdec,df,\
            dlgain,egain,dlength,lpbzcor,bbfile,magori,eyn,\
            lpyn=easygui.multenterbox(msg='Enter station information for the header file and orientation file',
                                      title='Station, Header and Orientation File Information',
                                      fields=hdrfields,
                                      values=hdrvalues)
            #write header file 
            hdrfid=open(os.path.join(stationdir,stationname+'.hdr'),'w')
            hdrfid.write('DATAID="'+dataid+'"'+'\n')
            hdrfid.write('ACQBY="'+acqby+'"'+'\n')
            hdrfid.write('ACQDATE='+acqd+'\n')
            hdrfid.write('FILEDATE='+filedate+'\n')
            hdrfid.write('PROSPECT="'+prospectloc+'"'+'\n')
            hdrfid.write('LOC="'+stationname+'"'+'\n')
            hdrfid.write('LAT='+lat+'\n')
            hdrfid.write('LONG='+long+'\n')
            hdrfid.write('ELEV='+elev)
            hdrfid.close()
            #write orientation file
            orifid=open(os.path.join(stationdir,stationname+'.ori'),'w')
            orifid.write(hxloc+'\n')
            orifid.write(hyloc+'\n')
            orifid.write(exloc+'\n')
            orifid.write(eyloc+'\n')
            orifid.write(rxloc+'\n')
            orifid.write(ryloc)
            orifid.close()
            
        #if there is no station info file   
        else:
            station=stationi
            stationhdr=stationi
            latstr='-00:00:00.0'
            lonstr='00:00:00.0'
            elev='000'
            startdate=todaysdate.strftime("%d/%m/%y")
                
            hdrfields=['Data ID',
                       'Acquired By',
                        'Acquired Date (DD/MM/YY)',
                        'File Date (DD/MM/YY)',
                        'Prospect Location',
                        'Station Name',
                        'Latitude (HH:MM:SS.AS)',
                        'Longitude ((HH:MM:SS.AS)',
                        'Elevation (m)',
                        'HX sensor position (x,y,angle from geographic North)',
                        'HX sensor position (x,y,angle from geographic North)',
                        'EX electrode location (x1,y1,x2,y2)',
                        'EY electrode location (x1,y1,x2,y2)',
                        'Remote Reference HX (x,y,angle from geographic North)',
                        'Remote Reference HX (x,y,angle from geographic North)',
                        'Magnetic Measurement type lp=long period, bb=broadband',
                        'Measured components',
                        'Magnetic declination',
                        'Sampling Frequency (Hz)',
                        'Data logger gain (very low=2.5,low=1.0,high=0.1)',
                        'Interface Box gain',
                        'Dipole lengths (m) Ex,Ey',
                        'Longperiod Bz correction',
                        'Broadband calibration file location (full path)',
                        'Measured magnetic orientation (BX,BY,BZ)',
                        'Electric fields converted to microV/m (y or n)',
                        'Longperiod magnetics converted to microV/m (y or n)']
            hdrvalues=[station,
                        'Adelaide University',
                        startdate,
                        startdate,
                        surveyname,
                        station,
                        latstr,
                        lonstr,
                        elev,
                        '0 0 0',
                        '0 0 90',
                        '0 0 100 0',
                        '0 0 0 100',
                        '0 0 0',
                        '0 0 90',
                        magtype,
                        'EX,EY,BX,BY',
                        magdec,
                        df,
                        dlgain,
                        egain,
                        dlength,
                        lpbzcor,
                        bbfile,
                        magori,
                        eyn,
                        lpyn]
            
            dataid,acqby,acqd,filedate,prospectloc,stationname,lat,long,elev,\
            hxloc,hyloc,exloc,eyloc,rxloc,ryloc,magtype,mcomps,magdec,df,\
            dlgain,egain,dlength,lpbzcor,bbfile,magori,eyn,\
            lpyn=easygui.multenterbox(msg='Enter information for the header file and orientation file',
                                      title='Station, Header and Orientation File Information',
                                      fields=hdrfields,
                                      values=hdrvalues)
            #write header file 
            hdrfid=open(os.path.join(stationdir,stationname+'.hdr'),'w')
            hdrfid.write('DATAID="'+dataid+'"'+'\n')
            hdrfid.write('ACQBY="'+acqby+'"'+'\n')
            hdrfid.write('ACQDATE='+acqd+'\n')
            hdrfid.write('FILEDATE='+filedate+'\n')
            hdrfid.write('PROSPECT="'+prospectloc+'"'+'\n')
            hdrfid.write('LOC="'+stationname+'"'+'\n')
            hdrfid.write('LAT='+lat+'\n')
            hdrfid.write('LONG='+long+'\n')
            hdrfid.write('ELEV='+elev)
            hdrfid.close()
            #write orientation file
            orifid=open(os.path.join(stationdir,stationname+'.ori'),'w')
            orifid.write(hxloc+'\n')
            orifid.write(hyloc+'\n')
            orifid.write(exloc+'\n')
            orifid.write(eyloc+'\n')
            orifid.write(rxloc+'\n')
            orifid.write(ryloc)
            orifid.close()
    
#dipole lengths
dlen=dlength.split(',')
exlen=float(dlen[0])
eylen=float(dlen[1])
#magnetic orientation list
magorilst=magori.split(',')
cacherate=inidict['cacherate']

#===============================================================================
# START GETTING BIRRP INFO: Ask for script file first
#===============================================================================
scriptfile=easygui.fileopenbox(msg='Open script file, if one does not exist hit cancel',
                               title='Get .script file',
                               filetypes=['*.script'],
                               default=os.path.join(stationdir,stationi))
if scriptfile!=None:
    #generate a directory to put birrp output files from script file
    bfpathi=os.path.dirname(scriptfile)
    fileslst=['files already combined']
    fileslstr=['files already combined']
    
    #get a few parameters from the user
    paramfields=['Location of birrp.exe',
                 'Location of broadband calibration file',
                 '3 letter outfile','Directory for output files',
                 'Are Electric fields converted from counts to units (y or n)?',
                 'Are Magnetic fields converted from counts to units(y or n)?']
    paramvalues=[inidict['birrploc'],
                 inidict['bbfile'],
                 stationi,
                 bfpathi,
                 eyn,
                 lpyn]
    birrploc,bbfile,ofil,bfpath,eyn,lpyn=easygui.multenterbox(msg='Few input parameters',
                                                              title='Few Parameters',
                                                              fields=paramfields,
                                                              values=paramvalues)
    station=ofil
    
    #create a directory for output files if not already there
    if not os.path.exists(bfpath):
        os.mkdir(bfpath)
        print 'Made directory: ', bfpath
    
    #open the script file to make sure everything is peachy
    subprocess.call([notepadpath,scriptfile])
    fid=file(scriptfile,'r')
    fidlines=fid.readlines()
    #convert files if need be
    #convert electric channels to microV/m
    
    cfilenlst=[]
    for ff in range(len(fidlines)):
        if fidlines[ff].find('EX')>=0:
            cfilenlst.append(fidlines[ff].rstrip())
        elif fidlines[ff].find('EY')>=0:
            cfilenlst.append(fidlines[ff].rstrip())
        elif fidlines[ff].find('BX')>=0:
            cfilenlst.append(fidlines[ff].rstrip())
        elif fidlines[ff].find('BY')>=0:
            cfilenlst.append(fidlines[ff].rstrip())
        elif fidlines[ff].find('BZ')>=0:
            cfilenlst.append(fidlines[ff].rstrip())
    if eyn=='n' or lpyn=='n':
        mt.convertCounts2Units(cfilenlst,eyn=eyn,lpyn=lpyn,magtype=magtype,
                               egain=egain,dlgain=dlgain,exlen=exlen,
                               eylen=eylen,zadj=lpbzcor)
        lpynf='y'
        eynf='y'
    else:
        eynf=eyn
        lpynf=lpyn
        
                

    #get some birrp parameters to put into .log and .ini file            
    ilev=fidlines[0].rstrip()
    nout=fidlines[1].rstrip()
    ninp=fidlines[2].rstrip()
    tbw=fidlines[3].rstrip()
    nfftline=fidlines[5].rstrip()
    nfftline=nfftline.split(' ')
    nfft=nfftline[0]
    nsctmax=nfftline[1]
    uinline=fidlines[7].rstrip()
    uinline=uinline.split(' ')
    uin=uinline[0]
    ainuin=uinline[1]
    c2threshe=fidlines[8].rstrip()
    #if there is a Bz component different locations for values
    if len(fidlines[9].rstrip())>=3:
        ofil=os.path.basename(fidlines[9].rstrip())
        nlev=fidlines[10].rstrip()
        nar=fidlines[12].rstrip()
        imode=fidlines[13].rstrip()
        jmode=fidlines[14].rstrip()
        nfil=fidlines[16].rstrip()
        nz='None'
        c2threshe1='None'
    else:
        nz=fidlines[9].rstrip()
        c2threshe1=fidlines[10].rstrip()
        ofil=os.path.basename(fidlines[11].rstrip())
        nlev=fidlines[12].rstrip()
        nar=fidlines[14].rstrip()
        imode=fidlines[15].rstrip()
        jmode=fidlines[16].rstrip()
        nfil=fidlines[17].rstrip()
    thetae=fidlines[len(fidlines)-3].rstrip()
    thetab=fidlines[len(fidlines)-2].rstrip()
    thetaf=fidlines[len(fidlines)-1].rstrip()
    
    rbasename=os.path.basename(fidlines[len(fidlines)-5].rstrip())
    if rbasename.find('d')>=0:
        rfind=rbasename.find('d')
        rstation=rbasename[:rfind-6]
    else:
        rstation=rbasename[:-9]
    fid.close()
    complstr=magori
    
    #===========================================================================
    #  Run BIRRP with Script file already written   
    #===========================================================================
        #Give the user time to  check the file
    time.sleep(15)
    #run BIRRP
    stimeburp=time.time()
    os.chdir(birrploc)
    os.system('birrp5<'+scriptfile)

    #check to see if BIRRP ran correctly
    etimeburp=time.time()
    dtburp=etimeburp-stimeburp
    if nfft<=2**14:
        print 'BIRRP Run time (m:s) for '+station \
                +': %.2f' % int(np.floor(dtburp/60)) \
                +':%.2f' % np.round(np.remainder(dtburp,60))
    elif nfft>2**14:
        if dtburp<10.:
            raise RuntimeError('BIRRP did not run properly, check script file')
        else:
            print 'BIRRP Run time (m:s) for ' \
                +station+': %.2f' % int(np.floor(dtburp/60)) \
                +':%.2f' % np.round(np.remainder(dtburp,60))


#===============================================================================
#if no script file make one

else:
    # Get values for appropriate parameters for site and remote reference  
    #ask how many time series segments to be processed

    nfilesfields=['Number of time series segments (1<= <=3)',
                  'Number of remote reference segments (1<= <=3)']
    nfilevalues=[1,1]
    nfiles=easygui.multenterbox(msg='How many segments of data (excluding components and number of cache files) to be processed and\or combined. \n ex.  for 2 segments of 4 hours of data each with remote reference and combine to get 4 hour blocks -> 2, 2, 2, 2',
                                title='Number of Time Series',
                                fields=nfilesfields, 
                                values=nfilevalues)
    ntseries=int(nfiles[0])
    nrseries=int(nfiles[1])
    
    fancynumstr=['1st','2nd','3rd']
    #Combine files time series
    tfiles=[]
    cdirpatht=[]
    fileslst=[]
    for ii in range(ntseries):
        cdirpath=easygui.diropenbox(msg='Locate the folder where the '+fancynumstr[ii]+' time series resides',
                                    title='Combine '+fancynumstr[ii]+' time series path',
                                    default=stationdir)
        cdirpatht.append(cdirpath)
        combfields=['Station',
                    'Start Time (HHMMSS)',
                    'End Time (HHMMSS)',
                    'Cache Length (HHMMSS)',
                    'Components',
                    'Decimation Factor',
                    'Electrics converted (y or n)',
                    'Magnetics converted (y or n)']
        
        combvalues=[stationi,
                    '000000',
                    '240000',
                    inidict['cacherate'],
                    mcomps,
                    '0',
                    eyn,
                    lpyn]
        station,stime,etime,cacherate,combcomp,dec,eyn,lpyn=\
                easygui.multenterbox(msg='Enter values',
                                     title='Combine Files Parameters',
                                     fields=combfields,
                                     values=combvalues)
        combcomplst=combcomp.split(',')
        dec=int(dec)
        #Combined files
        cfilenlst,filelst=mt.combineFewFiles(cdirpath,station,stime,
                                                  etime,cacherate,
                                                  complst=combcomplst,d=dec)
        #append list of new filenames of the combines files
        cfilenlst=[str(cfilenlst[ii]) for ii in range(len(cfilenlst))]
        tfiles.append(cfilenlst)
        #append list of combined filenames to list
        fileslst.append(filelst)
        
       #convert electric  and lp magnetic channels if need be and put back into file of same name
        if eyn=='n' or lpyn=='n':
            mt.convertCounts2Units(cfilenlst,eyn=eyn,lpyn=lpyn,
                                        magtype=magtype,egain=egain,
                                        dlgain=dlgain,exlen=exlen,
                                        eylen=eylen,zadj=lpbzcor)
            lpynf='y'
            eynf='y'
        else:
            eynf=eyn
            lpynf=lpyn
            
            
    #combine files remote reference
    rfiles=[]
    cdirpathr=[]
    fileslstr=[]
    for ii in range(nrseries):
        cdirpathrr=easygui.diropenbox(
                    msg='Locate the folder where the '+fancynumstr[ii]+\
                    ' remote reference time series resides. Hit cancel if same as time series',
                    title='Combine '+fancynumstr[ii]+' remote reference time series path',
                    default=inidict['defdirpath'])
        if cdirpathrr!=None:
            cdirpathr.append(cdirpathrr)
            combfields=['Station',
                        'Start Time (HHMMSS)',
                        'End Time (HHMMSS)',
                        'Cache Length (HHMMSS)',
                        'Components',
                        'Decimation Factor',
                        'Electrics converted (y or n)',
                        'Magnetics converted (y or n)']
            try:
                rrstattest=float(str(os.path.basename(
                                        os.path.dirname(cdirpathrr))))
                rrstation=str(os.path.basename(os.path.dirname(
                                        os.path.dirname(cdirpathrr))))
            except ValueError:
                rrstation=str(os.path.basename(os.path.dirname(
                                os.path.dirname(os.path.dirname(cdirpathrr)))))
            combvalues=[rrstation,
                        stime,   
                        etime,
                        inidict['cacherate'],
                        'BX,BY',
                        dec,
                        inidict['eyn'],
                        inidict['lpyn']]
            rstation,stime,etime,cacherate,combcompr,dec,eyn,lpyn=\
                        easygui.multenterbox(msg='Enter values',
                                             title='Combine Files Parameters',
                                             fields=combfields,
                                             values=combvalues)
            combcomplstr=combcompr.split(',')
            dec=int(dec)
            #combine files
            cfilenlstr,filelstr=mt.combineFewFiles(cdirpathrr,rstation,
                                                        stime,etime,cacherate,
                                                        complst=combcomplstr,
                                                        d=dec)
            #append combined filenames to list
            cfilenlstr=[str(cfilenlstr[ii]) for ii in range(len(cfilenlstr))]
            rfiles.append(cfilenlstr)
            #append combined filenames to list 
            fileslstr.append(filelstr)
            
            #convert electric  and lp magnetic channels if need be and put back into file of same name
            if eyn=='n' or lpyn=='n':
                mt.convertCounts2Units(cfilenlstr,eyn=eyn,lpyn=lpyn,
                                            magtype=magtype,egain=egain,
                                            dlgain=dlgain,exlen=exlen,
                                            eylen=eylen,zadj=lpbzcor)
                lpynf='y'
                eynf='y'
            else:
                eynf=eyn
                lpynf=lpyn
            
        else:
            rstation=station
            rfiles.append(tfiles[ii])
            cdirpathr.append(cdirpatht[ii])
            fileslstr.append('Files already combined')
            combcomplstr=['BX','BY']
    
#===============================================================================
#GUI to get BIRRP Parameters 
#===============================================================================
    #calculate number of points to read for time series
    if dec=='0' or dec==0:
        dec=1
    nreadini=int(3600*(float(df)/dec)*(float(etime[0:2])-float(stime[0:2])))
    nreadstr=str(nreadini)
    for yy in range(int(ntseries)-1):
        nreadstr+=','+nreadini
    #calculate number of points to skip
    npow=np.floor(np.log2(float(nreadini)))-16
    if npow>6:
        nfftpow=17
    elif npow>=2 and npow<=6:
        nfftpow=16
    elif npow>=-2 and npow<2:
        nfftpow=15
    elif npow>=-6 and npow<-2:
        nfftpow=14
    
    #calculate nfft and nsctmax
    nfft=2**nfftpow
    nsctmax=nfftpow-4
    compsr=combcomplstr[0]
    for cc in range(1,len(combcomplstr)):
        compsr+=','+combcomplstr[cc]
    #put filenames into lists
    #timeseries
    tfilesarray=np.array(tfiles)
    fileexlst=[]
    fileeylst=[]
    filebxlst=[]
    filebylst=[]
    filebzlst=[]
    for ss in range(len(tfilesarray)):
        for fs in range(len(tfilesarray[ss])):
            if tfilesarray[ss,fs].find('EX')>=0:
                fileexlst.append(tfilesarray[ss,fs])
            elif tfilesarray[ss,fs].find('EY')>=0:
                fileeylst.append(tfilesarray[ss,fs])
            elif tfilesarray[ss,fs].find('BX')>=0:
                filebxlst.append(tfilesarray[ss,fs])
            elif tfilesarray[ss,fs].find('BY')>=0:
                filebylst.append(tfilesarray[ss,fs])
            elif tfilesarray[ss,fs].find('BZ')>=0:
                filebzlst.append(tfilesarray[ss,fs])
                
    #remote reference
    rfilesarray=np.array(rfiles)
    rfileexlst=[]
    rfileeylst=[]
    rfilebxlst=[]
    rfilebylst=[]
    rfilebzlst=[]
    for rr in range(len(rfilesarray)):
        for fr in range(len(rfilesarray[rr])):
            if rfilesarray[rr,fr].find('EX')>=0:
                rfileexlst.append(rfilesarray[rr,fr])
            elif rfilesarray[rr,fr].find('EY')>=0:
                rfileeylst.append(rfilesarray[rr,fr])
            elif rfilesarray[rr,fr].find('BX')>=0:
                rfilebxlst.append(rfilesarray[rr,fr])
            elif rfilesarray[rr,fr].find('BY')>=0:
                rfilebylst.append(rfilesarray[rr,fr])
            elif tfilesarray[rr,fr].find('BZ')>=0:
                rfilebzlst.append(rfilesarray[rr,fr])
        
    #put names into a string
    try:
        exstr=fileexlst[0]
    except IndexError:
        exstr=''
    try:
        eystr=fileeylst[0]
    except IndexError:
        eystr=''
    try:
        bxstr=filebxlst[0]
    except IndexError:
        bxstr=''
    try:
        bystr=filebylst[0]
    except IndexError:
        bystr=''
    try:
        bzstr=filebzlst[0]
    except IndexError:
        bzstr=''
    try:
        rexstr=rfileexlst[0]
    except IndexError:
        rexstr=''
    try:
        reystr=rfileeylst[0]
    except IndexError:
        reystr=''
    try:
        rbxstr=rfilebxlst[0]
    except IndexError:
        rbxstr=''
    try:
        rbystr=rfilebylst[0]
    except IndexError:
        rbystr=''
    try:
        rbzstr=rfilebzlst[0]
    except IndexError:
        rbzstr=''
    
    if len(tfilesarray)>1:
        for ff in range(1,len(tfilesarray)):
            try:
                exstr+=','+fileexlst[ff]
            except IndexError:
                exstr+=', '
            try:
                eystr+=','+fileeylst[ff]
            except IndexError:
                eystr+=', '
            try:
                bxstr+=','+filebxlst[ff]
            except IndexError:
                bxstr+=', '
            try:
                bystr+=','+filebylst[ff]
            except IndexError:
                bystr+=', '
            try:
                bzstr+=','+filebzlst[ff]
            except IndexError:
                bzstr+=', '
                
    if len(rfilesarray)>1:
        for ff in range(1,len(rfilesarray)):
            try:
                rexstr+=','+rfileexlst[ff]
            except IndexError:
                rexstr+=', '
            try:
                reystr+=','+rfileeylst[ff]
            except IndexError:
                reystr+=', '
            try:
                rbxstr+=','+rfilebxlst[ff]
            except IndexError:
                rbxstr+=', '
            try:
                rbystr+=','+rfilebylst[ff]
            except IndexError:
                rbystr+=', '
            try:
                rbzstr+=','+rfilebzlst[ff]
            except IndexError:
                rbzstr+=', '
    #need a 180 rotation when using the coils for phases to be correct
    #declination can be subtracted so everything is in geographic coordinates
    if magtype=='bb':
        inidict['thetae']='0,90,'+str(180-float(magdec))
    
    #suggest an initial path
    bfpathi=os.path.join(os.path.dirname(fileexlst[0]),'BF')
    
    fields=['BIRRP location',
            'Broadband Calibration File Location',
            'Interaction Level (0-run or 1-interactive)',
            'Number of Output time series (2 or 3-> for BZ)',
            'Numper of input time series for E-field (1,2,3)', 
            'Time bandwidth for Sepian sequence',
            'Sampling rate (+) for (s), (-) for (Hz)',
            'Length of FFT (should be even)',
            'Number of windows used in FFT',
            'Frequencies OK (y or n)',
            'Residual rejection factor low end (usually 0)',
            'Residual rejection factor high end (.95-.99)',
            'Coherence threshold (0 if not desired)',
            'Threshold for Bz (0=separate from E, 1=E threshold, 2=E and B) Input if 3 B components else None',
            'Squared coherence for Bz, input if NZ=0, Nout=3', 
            'Output file root(usually three letters, can add full path)',
            'Output files (0=Z; 1=Z,qq; 2=Z,qq,w; 3=Z,qq,w,d)',
            'Number of independent data to be processed (1 for one segement)',
            'Prewhitening Filter (3< >15) or 0 if not desired',
            'Output file mode (0=ascii; 1=binary; 2=headerless ascii; 3=ascii in TS mode',
            'Input file mode (0=user defined; 1=start time YYYY-MM-DD HH:MM:SS)',
            'Number of points to be read for each data set (if segments>1 -> npts1,npts2...)',
            'Filter parameters (0=none; >0=input parameters; <0=filename',
            'Component list (in order for BIRRP usually EX,EY,BX,BY,BZ)',
            'Skip number of points in time series (0) if no skip, (if segements >1 -> input1,input2...)',
            'Remote Reference component list (in order for BIRRP usually BX,BY or BX,BY,BZ)',
            'Number of points to skip over (0) if none,(if segements >1 -> input1,input2...)',
            'Rotation angles for electrics (relative to geomagnetic North)(N,E,rot)',
            'Rotation angles for magnetics (relative to geomagnetic North)(N,E,rot)',
            'Rotation angles for calculation (relative to geomagnetic North)(N,E,rot)']
    
    values=[inidict['birrploc'],
            inidict['bbfile'],
            inidict['ilev'],
            inidict['nout'],
            inidict['ninp'],
            inidict['tbw'],
            '-'+str(float(df)/float(dec)),
            inidict['nfft'],
            inidict['nsctmax'],
            'y',
            inidict['uin'],
            inidict['ainuin'],
            inidict['c2threshe'],
            inidict['nz'],
            inidict['c2threshe1'],
            os.path.join(bfpathi,stationi),
            inidict['nlev'],
            ntseries,
            inidict['nar'],
            inidict['imode'],
            inidict['jmode'],
            nreadstr,
            inidict['nfil'],
            mcomps,
            0,
            compsr,
            0,
            inidict['thetae'],
            inidict['thetab'],
            inidict['thetaf']]
    
    #===========================================================================
    # Get values via gui
    #===========================================================================
    birrploc,bbfile,ilev,nout,ninp,tbw,deltat,nfft,nsctmax,yesno,uin,ainuin,\
        c2threshe,nz,c2threshe1,ofil,nlev,npcs,nar,imode,jmode,nread,nfil,\
        mcomps,nskip,complstr,nskipr,thetae,thetab,thetaf=\
        easygui.multenterbox(msg='Enter values for BIRRP script file',
                             title='BIRRP Script File for basic mode',
                             fields=fields,values=values,windowsize=[900,950])
    
    #split file names into a list
    tfileex=exstr.split(',')
    tfileey=eystr.split(',')
    tfilebx=bxstr.split(',')
    tfileby=bystr.split(',')
    tfilebz=bzstr.split(',')
    
    rfileex=rexstr.split(',')
    rfileey=reystr.split(',')
    rfilebx=rbxstr.split(',')
    rfileby=rbystr.split(',')
    rfilebz=rbzstr.split(',')

    
    #split values into lists
    if nread.find(',')==-1:
        nreadlst=[nread]
    else:
        nreadlst=nread.split(',')
        
    if nskip.find(',')==-1:
        nskiplst=[nskip]
    else:
        nskiplst=nskip.split(',')
    
    if nskipr.find(',')==-1:
        nskiprlst=[nskipr]
    else:
        nskiprlst=nskipr.split(',')
    
    #===========================================================================
    # Print values to file
    #===========================================================================
    #make a birrp folder to put data
    bfpath=os.path.dirname(ofil)
    if not os.path.exists(bfpath):
        os.mkdir(bfpath)
        print 'Made directory: ', bfpath
    
    complst=mcomps.split(',')
    complstr=complstr.split(',')
    #write to a file
    fid=open(os.path.join(bfpath,ofil+'.script'),'w')
    fid.write(ilev+'\n')
    fid.write(nout+'\n')
    fid.write(ninp+'\n')
    fid.write(tbw+'\n')
    fid.write(deltat+'\n')
    fid.write(nfft+' '+nsctmax+'\n')
    fid.write(yesno+'\n')
    fid.write(uin+' '+ainuin+'\n')
    fid.write(c2threshe+'\n')
    if nz.find('None')>=0:
        pass
    else:
        fid.write(nz+'\n')
    if c2threshe1.find('None')>=0:
        pass
    else:
        fid.write(c2threshe1+'\n')
    fid.write(ofil+'\n')
    fid.write(nlev+'\n')
    fid.write(npcs+'\n')
    fid.write(nar+'\n')
    fid.write(imode+'\n')
    fid.write(jmode+'\n')
    for ll in range(ntseries):
        fid.write(nreadlst[ll]+'\n')
        #write filenames
        for ii in range(len(complst)):
            if ll==0:
                if complst[ii]=='EX':
                    fid.write(nfil+'\n')
                    fid.write(tfileex[ll]+'\n')
                    fid.write(str(nskip)+'\n')
                elif complst[ii]=='EY':
                    fid.write(nfil+'\n')
                    fid.write(tfileey[ll]+'\n')
                    fid.write(str(nskip)+'\n')
                elif complst[ii]=='BX':
                    fid.write(nfil+'\n')
                    fid.write(tfilebx[ll]+'\n')
                    fid.write(str(nskip)+'\n')
                elif complst[ii]=='BY':
                    fid.write(nfil+'\n')
                    fid.write(tfileby[ll]+'\n')
                    fid.write(str(nskip)+'\n')
                elif complst[ii]=='BZ':
                    fid.write(nfil+'\n')
                    fid.write(tfilebz[ll]+'\n')
                    fid.write(str(nskip)+'\n')
            else:
                if complst[ii]=='EX':
                    fid.write(tfileex[ll]+'\n')
                    fid.write(str(nskip)+'\n')
                elif complst[ii]=='EY':
                    fid.write(tfileey[ll]+'\n')
                    fid.write(str(nskip)+'\n')
                elif complst[ii]=='BX':
                    fid.write(tfilebx[ll]+'\n')
                    fid.write(str(nskip)+'\n')
                elif complst[ii]=='BY':
                    fid.write(tfileby[ll]+'\n')
                    fid.write(str(nskip)+'\n')
                elif complst[ii]=='BZ':
                    fid.write(tfilebz[ll]+'\n')
                    fid.write(str(nskip)+'\n')
    for rr in range(nrseries):
        #write filenames for remote referencing
        for jj in range(len(complstr)):
            if complstr[jj]=='EX':
                fid.write(nfil+'\n')
                fid.write(rfileex[rr]+'\n')
                fid.write(str(nskipr)+'\n')
            elif complstr[jj]=='EY':
                fid.write(nfil+'\n')
                fid.write(rfileey[rr]+'\n')
                fid.write(str(nskipr)+'\n')
            elif complstr[jj]=='BX':
                fid.write(nfil+'\n')
                fid.write(rfilebx[rr]+'\n')
                fid.write(str(nskipr)+'\n')
            elif complstr[jj]=='BY':
                fid.write(nfil+'\n')
                fid.write(rfileby[rr]+'\n')
                fid.write(str(nskipr)+'\n')
            elif complstr[jj]=='BZ':
                fid.write(nfil+'\n')
                fid.write(rfilebz[rr]+'\n')
                fid.write(str(nskipr)+'\n')
    
    fid.write(thetae.replace(',',' ')+'\n')
    fid.write(thetab.replace(',',' ')+'\n')
    fid.write(thetaf.replace(',',' '))    
    fid.close()
    
    #check to see if the script file was written properly
    subprocess.call([notepadpath,os.path.join(bfpath,ofil+'.script')])
    time.sleep(15)
    #===============================================================================
    # Run BIRRP
    #===============================================================================
    stimeburp=time.time()
    os.chdir(birrploc)
    subprocess.os.system('birrp5<'+os.path.join(bfpath,ofil+'.script'))
    
    #check to see if BIRRP ran correctly
    etimeburp=time.time()
    dtburp=etimeburp-stimeburp
    #print 'BIRRP Run time (m:s): ', str(int(np.floor(dtburp/60)))+':'+str(np.round(np.remainder(dtburp,60)))
    if int(nfft)<=2**14:
        print 'BIRRP Run time (m:s) for '+station+': ', str(int(np.floor(dtburp/60)))+':'+str(np.round(np.remainder(dtburp,60)))
    elif int(nfft)>2**14:
        if dtburp<10.:
            raise RuntimeError('BIRRP did not run properly, check script file')
        else:
            print 'BIRRP Run time (m:s) for '+station+': ', str(int(np.floor(dtburp/60)))+':'+str(np.round(np.remainder(dtburp,60)))

#write .ini file by putting parameters into a dictionary
#make complstr a string
complstrstr=complstr[0]
for bb in range(len(complstr)):
    complstrstr+=','+complstr[bb]

argsdict={}
argsdict['defdirpath']=stationdir
argsdict['station']=station
argsdict['magtype']=magtype
argsdict['lpyn']=lpynf
argsdict['eyn']=eynf
argsdict['mcomps']=mcomps
argsdict['magdec']=magdec
argsdict['df']=str(df)
argsdict['cacherate']=cacherate
argsdict['dlength']=dlength
argsdict['dlgain']=str(dlgain)
argsdict['egain']=egain
argsdict['lpbzcor']=lpbzcor
argsdict['bbfile']=bbfile
argsdict['magori']=magori
argsdict['birrploc']=birrploc
argsdict['ilev']=ilev
argsdict['nout']=nout
argsdict['ninp']=ninp
argsdict['tbw']=tbw
argsdict['nfft']=nfft
argsdict['nsctmax']=nsctmax
argsdict['uin']=uin
argsdict['ainuin']=ainuin
argsdict['c2threshe']=c2threshe
argsdict['nz']=nz
argsdict['c2threshe1']=c2threshe1
argsdict['ofil']=ofil
argsdict['nlev']=nlev
argsdict['nar']=nar
argsdict['imode']=imode
argsdict['jmode']=jmode
argsdict['nfil']=nfil
argsdict['complstr']=complstrstr
argsdict['thetae']=thetae
argsdict['thetab']=thetab
argsdict['thetaf']=thetaf
if stationinfofile==None:
    argsdict['stationinfofile']=' '
else:
    argsdict['stationinfofile']=stationinfofile

brp.writeini(stationdir,argsdict)

#===============================================================================
# read BIRRP output files into .dat,.coh,.edi
#===============================================================================
#write a coherence file .coh
dlgain=float(dlgain)
df=float(df)
if stationinfofile==' ':
    stationinfofile=None

cohfile=brp.writecoh(bfpath)
edifile=brp.writeedi(bfpath,station,stationinfofile=stationinfofile,
                     rrstation=rstation,birrpdict=argsdict,bbfile=bbfile,
                     ffactor=ffactor)
#edifile=brp.writeedi(bfpath, magtype, dlgain,bbfile=bbfile)
datfile=brp.writedat(bfpath,df=df,egain=egain,dlgain=dlgain,bbfile=bbfile,
                     dlen=[exlen,eylen],magtype=magtype,
                     ffactor=ffactor)
impfile=brp.writeimp(bfpath,egain=egain,dlgain=dlgain,dlen=[exlen,eylen],
                     magtype=magtype,bbfile=bbfile,ffactor=ffactor)


#make a folder to put all .edi files into
if not os.path.exists(os.path.join(os.path.dirname(stationdir),'EDIfiles')):
    os.mkdir(os.path.join(os.path.dirname(stationdir),'EDIfiles'))
    print 'Made path: '+os.path.join(os.path.dirname(stationdir),'EDIfiles')

#copy .edi file to common folder
shutil.copy(edifile,os.path.join(os.path.dirname(stationdir),'EDIfiles',
                                 station+'.edi'))

#===============================================================================
#Write an .ini file that has all the recent parameters and a .log file
#===============================================================================
logfid=file(os.path.join(bfpath,ofil+'.log'),'a')
todaysdate=datetime.datetime.today()
logfid.write('==========================================================='+'\n')
logfid.write('Processing log for station: ' +os.path.join(bfpath,station) +'\n')
logfid.write('Processed on: ' 
            +todaysdate.strftime("%b %d, %Y at %I:%M:%S %p local time")+'\n')
logfid.write('-----Files Combined----- \n')
fileslststr=str(fileslst)
flstr=fileslststr.replace('[','')
flstr=flstr.replace(']','')
flstr=flstr.replace('"','')
flstr=flstr.split(',')
for jj in range(len(flstr)):
    logfid.write(str(flstr[jj]) +'\n')
logfid.write('-----Remote Reference Files Combined----- \n')
fileslststrr=str(fileslstr)
flstrr=fileslststrr.replace('[','')
flstrr=flstrr.replace(']','')
flstrr=flstrr.replace('"','')
flstrr=flstrr.split(',')
for jj in range(len(flstrr)):
    logfid.write(str(flstrr[jj]) +'\n')
logfid.write('-----BIRRP Parameters----- \n')
logfid.write('BIRRP.exe location: '+birrploc+'\n')
logfid.write('BIRRP run time(m:s.ds): '+ str(int(np.floor(dtburp/60)))+':'
                                +str(np.round(np.remainder(dtburp,60)))+'\n')
birrpfid=file(os.path.join(bfpath,ofil+'.script'),'r')
birrplines=birrpfid.readlines()
for bb in range(len(birrplines)):
    line=birrplines[bb].rstrip()
    logfid.write(line+'\n')
birrpfid.close()
logfid.write('-----Impedance Calculated from BIRRP-----'+'\n')
impfid=file(impfile,'r')
implines=impfid.readlines()
for mm in range(len(implines)):
    logfid.write(implines[mm].rstrip()+'\n')
impfid.close()
logfid.write('-----Coherence Calculated from BIRRP-----'+'\n')
cohfid=file(cohfile,'r')
cohlines=cohfid.readlines()
for nn in range(len(cohlines)):
    logfid.write(cohlines[nn].rstrip()+'\n')
cohfid.close()
logfid.write('-----Resistivity and Phase Calculated from BIRRP-----'+'\n')
datfid=file(datfile,'r')
datlines=datfid.readlines()
for pp in range(len(datlines)):
    logfid.write(datlines[pp].rstrip()+'\n')
datfid.close()

logfid.close()


#===============================================================================
# Check the quality of the data by opening files and plotting coh, res, phase
#===============================================================================
#open the .coh, .imp, and .edi files to quality check
subprocess.call([notepadpath,
                 bfpath+os.sep+ofil+'.log'])
subprocess.call([notepadpath,cohfile])
subprocess.call([notepadpath,impfile])
subprocess.call([notepadpath,datfile])

#plot the coherence and apparent resistivity and phase
if not os.path.exists(os.path.join(os.path.dirname(stationdir),'Plots')):
    os.mkdir(os.path.join(os.path.dirname(stationdir),'Plots'))
    print 'Made path: '+os.path.join(os.path.dirname(stationdir),'Plots')

plotpath=os.path.join(os.path.dirname(stationdir),'Plots')

fig1=mtplot.plotcoh(cohfile,1,savefigfilename=os.path.join(plotpath,
                                                           station+'coh.pdf'))
fig2=mtplot.plotResPhase(edifile,plotnum=2,fignum=2,ffactor=ffactor,
                         savefigfilenames=[os.path.join(plotpath,
                                                        station+'Resxy.pdf'),
                                                        os.path.join(plotpath,
                                                                     station+\
                                                                     'Resxx.pdf')])
#tsn=np.shape(cfilenlst)
#try:
#    tsn[1]
#    for tt in range(tsn[0]):
#        mtplot.plotTS(cfilenlst[tt],)
#fig3=mtplot.plotTS(cfilenlst,
print 'Save plots: '+os.path.join(plotpath,station+'coh.pdf')
print '\t '+os.path.join(plotpath,station+'Resxy.pdf')
print '\t '+os.path.join(plotpath,station+'Resxx.pdf')

mtplot.plotcoh(cohfile,1)
mtplot.plotResPhase(edifile,plotnum=2,fignum=2,ffactor=ffactor)
